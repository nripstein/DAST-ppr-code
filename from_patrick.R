library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(aghq)
# devtools::install_github("nripstein/DAST", ref = "development")
# devtools::install_github("https://github.com/nripstein/DAST/tree/rand-betas-optional")
library(DAST)


######################## PRIOR AND POSTERIOR PLOT FUNCTION ########################

plot_priors_posteriors_manual <- function(fit_object) {
  
  # 1. Setup
  opt  <- fit_object$aghq_model$optresults
  mode <- opt$mode
  hess <- opt$hessian
  cov_mat <- tryCatch(solve(hess), error = function(e) MASS::ginv(hess))
  par_names <- names(mode)
  
  # 2. Direct Data Extraction
  # FIX: User identified the TMB object is located at fit_object$obj
  data_list <- fit_object$obj$env$data
  
  
  # Verify all required keys exist first
  req_vars <- c("prior_rho_min", "prior_rho_prob", "prior_sigma_max",
                "prior_sigma_prob", "prior_iideffect_sd_max", "prior_iideffect_sd_prob")
  
  missing <- setdiff(req_vars, names(data_list))
  if (length(missing) > 0) {
    stop(paste("Error: Missing required prior data in fit_object:", paste(missing, collapse = ", ")))
  }
  
  # Extract values directly
  rho_min  <- as.numeric(data_list$prior_rho_min)
  rho_prob <- as.numeric(data_list$prior_rho_prob)
  sig_max  <- as.numeric(data_list$prior_sigma_max)
  sig_prob <- as.numeric(data_list$prior_sigma_prob)
  tau_max  <- as.numeric(data_list$prior_iideffect_sd_max)
  tau_prob <- as.numeric(data_list$prior_iideffect_sd_prob)
  
  
  # Helper: Posterior Density
  get_post_dens <- function(param_name, transform_fun, jacobian_fun, x_range_sd=4) {
    idx <- which(par_names == param_name)
    if (length(idx) == 0) return(NULL)
    mu <- mode[idx]; sd <- sqrt(cov_mat[idx, idx])
    x_est <- seq(mu - x_range_sd * sd, mu + x_range_sd * sd, length.out = 1000)
    d_est <- dnorm(x_est, mean = mu, sd = sd)
    x_nat <- transform_fun(x_est)
    jac   <- jacobian_fun(x_est)
    d_nat <- d_est * (1 / abs(jac))
    data.frame(x = x_nat, y = d_nat, type = "Posterior")
  }
  
  plot_list <- list()
  
  # --- 3. HYPERPARAMETERS ---
  
  # A. Range (Rho)
  if ("log_rho" %in% par_names) {
    df_post <- get_post_dens("log_rho", exp, exp)
    
    # PC Prior Formula
    lam_rho <- -log(rho_prob) * rho_min
    x_grid  <- df_post$x[df_post$x > 0]
    y_prior <- (lam_rho / x_grid^2) * exp(-lam_rho / x_grid)
    
    df_prior <- data.frame(x = x_grid, y = y_prior, type = "Prior")
    plot_list[["rho"]] <- rbind(df_post %>% filter(x>0), df_prior) %>%
      mutate(param = "Range (Rho) [m]")
  }
  
  # B. Spatial SD (Sigma)
  if ("log_sigma" %in% par_names) {
    df_post <- get_post_dens("log_sigma", exp, exp)
    
    # Exponential Formula
    lam_sig <- -log(sig_prob) / sig_max
    x_grid  <- df_post$x
    y_prior <- dexp(x_grid, rate = lam_sig)
    
    df_prior <- data.frame(x = x_grid, y = y_prior, type = "Prior")
    plot_list[["sigma"]] <- rbind(df_post, df_prior) %>%
      mutate(param = "Spatial SD (Sigma)")
  }
  
  # C. NB Dispersion (Tau)
  if ("iideffect_log_tau" %in% par_names) {
    df_post <- get_post_dens("iideffect_log_tau", exp, exp)
    
    # Exponential Formula
    lam_tau <- -log(tau_prob) / tau_max
    x_grid  <- df_post$x
    y_prior <- dexp(x_grid, rate = lam_tau)
    
    df_prior <- data.frame(x = x_grid, y = y_prior, type = "Prior")
    plot_list[["tau"]] <- rbind(df_post, df_prior) %>%
      mutate(param = "NB Dispersion (Tau)")
  }
  
  # --- 4. FIXED EFFECTS ---
  beta_idx <- grep("intercept|slope|beta", par_names, ignore.case = TRUE)
  exclude  <- which(par_names %in% c("log_rho", "log_sigma", "iideffect_log_tau"))
  beta_idx <- setdiff(beta_idx, exclude)
  
  for (i in beta_idx) {
    pname <- par_names[i]
    df_post <- get_post_dens(pname, identity, function(x) 1)
    
    x_grid  <- df_post$x
    y_prior <- dnorm(x_grid, mean = 0, sd = 100) # Assuming SD=100 default
    
    dname <- tools::toTitleCase(gsub("_", " ", pname))
    df_prior <- data.frame(x = x_grid, y = y_prior, type = "Prior")
    plot_list[[pname]] <- rbind(df_post, df_prior) %>% mutate(param = dname)
  }
  
  # --- 5. PLOTTING ---
  
  full_df <- bind_rows(plot_list)
  hyper_names <- c("Range (Rho) [m]", "Spatial SD (Sigma)", "NB Dispersion (Tau)")
  
  p1 <- NULL
  df_hyper <- full_df %>% filter(param %in% hyper_names)
  if(nrow(df_hyper) > 0) {
    p1 <- ggplot(df_hyper, aes(x = x, y = y)) +
      geom_line(aes(color = type, linetype = type), linewidth = 1.1) +
      scale_color_manual(values = c("Posterior" = "navy", "Prior" = "firebrick")) +
      scale_linetype_manual(values = c("Posterior" = "solid", "Prior" = "dashed")) +
      facet_wrap(~param, scales = "free") +
      theme_bw() +
      labs(title = "Hyperparameters", subtitle = "PC Priors (Red) vs Posterior (Blue)",
           x = "Value", y = "Density")
  }
  
  p2 <- NULL
  df_fixed <- full_df %>% filter(!param %in% hyper_names)
  if(nrow(df_fixed) > 0) {
    p2 <- ggplot(df_fixed, aes(x = x, y = y)) +
      geom_line(aes(color = type, linetype = type), linewidth = 1.1) +
      scale_color_manual(values = c("Posterior" = "navy", "Prior" = "firebrick")) +
      scale_linetype_manual(values = c("Posterior" = "solid", "Prior" = "dashed")) +
      facet_wrap(~param, scales = "free") +
      theme_bw() +
      labs(title = "Fixed Effects", subtitle = "Priors (Red) vs Posterior (Blue)",
           x = "Value", y = "Density")
  }
  
  return(list(hyper = p1, fixed = p2))
}

######################## LOAD DATA ########################

## Choose resolution branch: "100m" or "1km"
res_choice <- "1km"  # or "100m"


data_dir <- "/home/noah/nuts_data/"

options(timeout = 3600)

Syears <- c(2003, 2006, 2010, 2013, 2016, 2021, 2024)
Surl <- paste(
  "https://gisco-services.ec.europa.eu/distribution/v2/nuts/gpkg/NUTS_RG_10M_",
  Syears, "_3035.gpkg",
  sep = ""
)
names(Surl) <- Syears

cat("About to download NUTS files...", "\n")
Sfiles <- Pmisc::downloadIfOld(Surl, path = data_dir)
nutsAll <- mapply(terra::vect, x = Sfiles)
names(nutsAll) <- Syears

nuts <- list()
for (D in names(nutsAll)) {
  x <- nutsAll[[D]]
  nuts[[D]] <- x[x$CNTR_CODE %in% c("NL", "BE") & x$LEVL_CODE == 3, ]
}

nutsYear <- list(
  "2018" = nuts[["2016"]],
  "2019" = nuts[["2016"]],
  "2020" = nuts[["2021"]],
  "2021" = nuts[["2021"]],
  "2022" = nuts[["2021"]],
  "2023" = nuts[["2024"]]
)

mort <- eurostat::get_eurostat("demo_r_magec3", type = "code", cache_dir = data_dir)
mort$year <- as.numeric(gsub("-.*", "", mort$TIME_PERIOD))

popDf <- expand.grid(year = names(nutsYear), country = c("NLD", "BEL"), age = c("00", "01"))



wp_base <- "https://data.worldpop.org/GIS/AgeSex_structures/Global_2015_2030/R2025A"

if (res_choice == "1km") {
  popDf$url <- paste0(
    wp_base, "/",
    popDf$year, "/",
    popDf$country,
    "/v1/1km_ua/constrained/",
    tolower(popDf$country),
    "_f_",
    popDf$age,
    "_",
    popDf$year,
    "_CN_1km_R2025A_UA_v1.tif"
  )
} else if (res_choice == "100m") {
  popDf$url <- paste0(
    wp_base, "/",
    popDf$year, "/",
    popDf$country,
    "/v1/100m/constrained/",
    tolower(popDf$country),
    "_f_",
    popDf$age,
    "_",
    popDf$year,
    "_CN_100m_R2025A_v1.tif"
  )
} else {
  stop("res_choice must be '100m' or '1km'")
}

cat("About to download population files...", "\n")
popDf$files <- Pmisc::downloadIfOld(popDf$url, path = data_dir)

popDfSplit <- split(popDf, popDf$country)
popSplit <- mapply(function(xx) terra::rast(xx$files), popDfSplit)

pop <- terra::mosaic(
  terra::tapp(popSplit[[1]], popDfSplit[[1]]$year, sum),
  terra::tapp(popSplit[[2]], popDfSplit[[2]]$year, sum)
)

popT1 <- terra::project(pop, terra::crs(nuts[[1]]))
pop_file <- file.path(data_dir, paste0("pop_", res_choice, ".tif"))
popT <- terra::writeRaster(popT1, pop_file, overwrite = TRUE)
popT <- terra::rast(pop_file)

mort1 <- as.data.frame(mort[mort$unit == "NR" & mort$age == "Y_LT5" & mort$sex == "F", ])
matchYear <- mort2 <- list()
for (D in names(nutsYear)) {
  mortHere <- mort1[mort1$year == as.integer(D), ]
  nutsYear[[as.character(D)]]$deaths <-
    mortHere[
      match(
        nutsYear[[as.character(D)]]$NUTS_ID,
        mortHere$geo
      ),
      "values"
    ]
}

lapply(nutsYear, function(xx) quantile(xx$deaths))

library("terra")
pCol <- mapmisc::colourScale(breaks = c(0, 0.1, 0.5, 1, 2, 20), style = "fixed")

mapmisc::map.new(popT)
plot(popT[[1]], breaks = pCol$breaks, col = pCol$col, legend = FALSE, maxcell = 1e7, add = TRUE)
plot(nuts[[6]], add = TRUE)
mapmisc::legendBreaks("topleft", pCol, inset = 0.1, bty = "n")

data("worldMap", package = "mapmisc")
world <- terra::unwrap(worldMap)
worldT <- terra::project(world, terra::crs(nuts[[1]]))

nuts3 <- nutsYear[c(1, 3, 6)]
Scol <- c("lightblue", "green", "blue")
Slwd <- seq(2, 0.8, len = length(nuts3))^2
Slty <- rep(1, length(nuts3))
mapmisc::map.new(worldT[grep("Belgium|Nether", worldT$NAME), ], buffer = 10 * 1000)
plot(popT[[1]], breaks = pCol$breaks, col = pCol$col, legend = FALSE, maxcell = 1e7, add = TRUE)
for (D in 1:length(nuts3)) plot(nuts3[[D]], lwd = Slwd[D], add = TRUE, border = Scol[D], lty = Slty[D])
legend("topleft", lwd = 3, col = Scol, legend = names(nuts3), bty = "n", inset = 0.1)

#Noah: fit model to nutsYear[["2017"]]$deaths, popT[["2017"]] (and other years)!

######################## NOAH START ########################

plot_sf <- function(sf_df, field_name) {
  ggplot2::ggplot(sf_df) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[field_name]])) +
    ggplot2::scale_fill_viridis_c(limits = range(sf_df[[field_name]], na.rm = TRUE))
}
# plot(nutsYear[["2018"]], "deaths")
# plot(popT[["X2018"]])

poly_list <- lapply(nutsYear, sf::st_as_sf)
offset_list <- setNames(
  lapply(names(popT), function(nm) popT[[nm]]),
  names(popT)
)
names(offset_list) <- sub("^X", "", names(popT))





limits <- sf::st_bbox(poly_list[[1]])
hypotenuse <- sqrt((limits$xmax - limits$xmin)^2 + (limits$ymax - limits$ymin)^2)
maxedge <- hypotenuse / 10



prep_time <- system.time({
  dis_data_mmap <- DAST::prepare_data_mmap(
    polygon_shapefile = poly_list,
    aggregation_rasters_list = offset_list,
    # mesh_args = list(
    #   max.edge = c(maxedge/5, maxedge),
    #   cutoff = 0.05*maxedge,
    #   offset = c(hypotenuse/2, hypotenuse),
    #   max.n=20*1000,
    #   resolution=600
    # ),
    mesh_args = list(
      max.edge = c(maxedge / 2, maxedge / 3 * 2 * 3),
      cutoff = 0.01,
      offset = c(hypotenuse / 5 * 5, hypotenuse / 5 * 5)
    ),
    response_var = "deaths",
    id_var = "NUTS_ID",
    na_action = TRUE,
    verbose = TRUE
  )
})["elapsed"]

mapmisc::map.new(terra::vect(poly_list[[1]]), buffer=100*1000)
plot(terra::vect(poly_list[[1]]), add=TRUE)
points(dis_data_mmap$mesh$loc[,1:2], pch=15, cex=0.4)     
print(paste0(nrow(dis_data_mmap$mesh$graph$tv), " basis functions"))

plot(dis_data_mmap)
# print(DAST::get_priors(dis_data_mmap))
crop_region <- terra::aggregate(nutsYear[["2018"]], dissolve = TRUE)


# get the default priors and modify them
my_priors <- DAST::get_priors(dis_data_mmap) # gets default priors encoded in the data object (they're dependent on the data)
my_priors$prior_rho_min  <- 25000  # Units are meters!
my_priors$prior_rho_prob  <- 0.5  # Units are meters!

my_priors$prior_iideffect_sd_max  = 0.1 
my_priors$prior_iideffect_sd_prob = 0.5

# my_priors$prior_rho_min <- 125000 # DEBUG


### time-pooled betas FALSE K=1 START

fm_nb_k1_pooled_priros <- DAST::disag_model_mmap(
  data = dis_data_mmap,
  # aghq_k = 1,
  family = "negbinomial",
  link = "log", # default
  iid = TRUE, # default
  verbose = TRUE, # FALSE is default
  priors = my_priors,
  engine = "AGHQ",
  engine.args = list(optimizer = "nlminb",
                     aghq_k = 1),
  fixed_effect_betas = TRUE,
  # optimizer = "nlminb", # "BFGS" is default
  time_varying_betas = FALSE # TRUE is default
)
pred_k1_pooled_priors <- predict(fm_nb_k1_pooled_priros, verbose=TRUE)


pred_mask = terra::mask(pred_k1_pooled_priors$mean_prediction$prediction$time_1, crop_region)
pred_col = mapmisc::colourScale(pred_mask, dec=1, breaks=8, 
                                style='equal', col='Spectral', 
                                transform='log', rev=TRUE)
mapmisc::map.new(pred_mask)
plot(pred_mask, col=pred_col$col, breaks=pred_col$breaks, legend=FALSE, add=TRUE)
mapmisc::legendBreaks('topleft', pred_col, bty='n')


plots_priors <- plot_priors_posteriors_manual(fm_nb_k1_pooled_priros)
print(plots_priors$hyper)
print(plots_priors$fixed)

#
print(summary(fm_nb_k1_pooled_priros))

### time varying FALSE END

### time varying TRUE K=1 START
fm_nb_k1 <- DAST::disag_model_mmap(
  data = dis_data_mmap,
  # aghq_k = 1,
  family = "negbinomial",
  link = "log",
  iid = TRUE,
  verbose = TRUE,
  priors = my_priors,
  engine = "AGHQ",
  engine.args = list(optimizer = "nlminb", aghq_k = 1),
  # fixed_effect_betas = TRUE,
  # optimizer = "nlminb",
  time_varying_betas = TRUE
)
pred_k1 <- predict(fm_nb_k1, verbose=TRUE)
plot(terra::mask(pred_k1$mean_prediction$prediction$time_1, crop_region), main = "Predicted deaths time 1 - K=1 nonpooled")
plot(terra::mask(pred_k1$mean_prediction$prediction$time_2, crop_region), main = "Predicted deaths time 2 - K=1 nonpooled")
summary(fm_nb_k1)

plots_priors_fm_nb_k1 <- plot_priors_posteriors_manual(fm_nb_k1)
print(plots_priors_fm_nb_k1$hyper)
print(plots_priors_fm_nb_k1$fixed)






# ?disag_model_mmap

