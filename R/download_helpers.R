# Download helpers for vignette_manchester.Rmd.
# Taken from the AI-generated final_manchester_2008_2012.Rmd, which used them
# successfully; they retry with several download methods because the ONS and
# ArcGIS servers intermittently drop large downloads.

download_if_missing <- function(url, destfile) {
  needs_download <- !file.exists(destfile) || isTRUE(file.info(destfile)$size == 0)
  if (!needs_download) {
    return(destfile)
  }

  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(destfile, ".part")
  methods <- c("libcurl", "curl", "auto")
  last_err <- NULL

  for (attempt in 1:4) {
    for (m in methods) {
      ok <- tryCatch({
        if (file.exists(tmp)) file.remove(tmp)
        suppressWarnings(
          download.file(
            url = url,
            destfile = tmp,
            mode = "wb",
            quiet = FALSE,
            method = m,
            headers = c("User-Agent" = "R/disagg notebook")
          )
        )
        file.exists(tmp) && file.info(tmp)$size > 0
      }, error = function(e) {
        last_err <<- e$message
        FALSE
      })

      if (isTRUE(ok)) {
        if (file.exists(destfile)) file.remove(destfile)
        file.rename(tmp, destfile)
        return(destfile)
      }
    }
    Sys.sleep(2 ^ attempt)
  }

  stop("Failed to download: ", url, if (!is.null(last_err)) paste0(" (last error: ", last_err, ")"))
}

fetch_arcgis_attrs <- function(base_query_url, out_fields, where = "1=1", page_size = 100) {
  fetch_json <- function(url) {
    tmp <- tempfile(fileext = ".json")
    on.exit(unlink(tmp), add = TRUE)
    methods <- c(getOption("download.file.method", "auto"), "curl", "auto")
    methods <- unique(methods)
    last_err <- NULL

    for (attempt in 1:4) {
      for (m in methods) {
        ok <- tryCatch({
          if (file.exists(tmp)) file.remove(tmp)
          suppressWarnings(
            download.file(
              url = url,
              destfile = tmp,
              mode = "wb",
              quiet = TRUE,
              method = m,
              headers = c("User-Agent" = "R/disagg notebook")
            )
          )
          file.exists(tmp) && file.info(tmp)$size > 0
        }, error = function(e) {
          last_err <<- e$message
          FALSE
        })
        if (isTRUE(ok)) {
          return(jsonlite::fromJSON(tmp, simplifyVector = FALSE))
        }
      }
      Sys.sleep(2 ^ attempt)
    }
    stop("ArcGIS download failed for URL: ", substr(url, 1, 140), if (!is.null(last_err)) paste0(" ... (last error: ", last_err, ")"))
  }

  attrs_to_dt <- function(features) {
    attrs <- lapply(features, function(f) f$attributes)
    attrs <- attrs[vapply(attrs, function(x) !is.null(x) && length(x) > 0, logical(1))]
    if (length(attrs) == 0) return(NULL)
    data.table::rbindlist(lapply(attrs, as.list), fill = TRUE, use.names = TRUE)
  }

  ids_url <- paste0(
    base_query_url,
    "?where=", URLencode(where, reserved = TRUE),
    "&returnIdsOnly=true&f=pjson"
  )
  ids_res <- fetch_json(ids_url)

  if (!is.null(ids_res$error)) {
    err_msg <- if (!is.null(ids_res$error$message)) ids_res$error$message else "Unknown error"
    stop("ArcGIS ID query failed: ", err_msg)
  }

  rows <- list()
  object_ids <- ids_res$objectIds

  if (!is.null(object_ids) && length(object_ids) > 0) {
    chunks <- split(object_ids, ceiling(seq_along(object_ids) / page_size))
    for (i in seq_along(chunks)) {
      oid_str <- paste(chunks[[i]], collapse = ",")
      q <- paste0(
        base_query_url,
        "?objectIds=", URLencode(oid_str, reserved = TRUE),
        "&outFields=", URLencode(paste(out_fields, collapse = ","), reserved = TRUE),
        "&f=pjson&returnGeometry=false"
      )
      res <- fetch_json(q)
      if (!is.null(res$error)) {
        err_msg <- if (!is.null(res$error$message)) res$error$message else "Unknown error"
        stop("ArcGIS chunk query failed (chunk ", i, "): ", err_msg)
      }
      if (!is.null(res$features) && length(res$features) > 0) {
        dt <- attrs_to_dt(res$features)
        if (!is.null(dt)) rows[[length(rows) + 1L]] <- dt
      }
    }
  } else {
    q <- paste0(
      base_query_url,
      "?where=", URLencode(where, reserved = TRUE),
      "&outFields=", URLencode(paste(out_fields, collapse = ","), reserved = TRUE),
      "&f=pjson&returnGeometry=false"
    )
    res <- fetch_json(q)
    if (!is.null(res$error)) {
      err_msg <- if (!is.null(res$error$message)) res$error$message else "Unknown error"
      stop("ArcGIS query failed: ", err_msg)
    }
    if (!is.null(res$features) && length(res$features) > 0) {
      dt <- attrs_to_dt(res$features)
      if (!is.null(dt)) rows[[length(rows) + 1L]] <- dt
    }
  }

  if (length(rows) == 0) {
    return(data.table::as.data.table(setNames(replicate(length(out_fields), character(), simplify = FALSE), out_fields)))
  }
  data.table::rbindlist(rows, fill = TRUE, use.names = TRUE)
}
