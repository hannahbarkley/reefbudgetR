#' Calculate production rates for each benthic component measured on
#' a carbonate budget transect
#'
#'@author Hannah Barkley
#'
#'@param substrate_class Type of substrate observed ("CORAL", "CCA", "TURF", "MA", ...).
#'@param substrate_code Taxa code observed (e.g. "PLOB", "MCAP", "PMEA", "CCA", ...).
#'@param morphology_code Taxa morphology observed (corals only; "BR", "MD", "EM", ...). Non-corals default to NA.
#'@param substrate_cover_cm Measure surface distance of benthic component, in cm.
#'@param region_code Survey region ("MHI", "MARIAN", ...).
#'#'@param dbase_type Production database to use, either Indo-Pacific ReefBudget ("IPRB")
#'or U.S. Pacific Islands NCRMP-specific database ("NCRMP"). The Indo-Pacific ReefBudget
#'database is derived from "IP Calcification and bioerosion rates database v.1.3",
#'downloaded from https://geography.exeter.ac.uk/reefbudget/indopacific/.
#'@param prod_dbase Production database to reference, either Indo-Pacific ReefBudget ("IPRB")
#'or NCRMP-specific ("NCRMP"). Defaults to "NCRMP".
#'
#'@details See included carbonate production databases `prod_dbase_iprb` and `prod_dbase_ncrmp`
#'for acceptable `substrate_class`, `substrate_code`, and `morphology_code` values.
#'
#'@export calc_prod
#'
#'@examples
#' calc_prod(
#' substrate_class = "CORAL",
#' substrate_code = "PLOB",
#' morphology_code = "MD",
#' substrate_cover_cm = 10,
#' dbase_type = "NCRMP",
#' prod_dbase = "NCRMP")
#'


calc_prod <- function(substrate_class,
                      substrate_code,
                      morphology_code,
                      substrate_cover_cm,
                      location_code = NULL,
                      dbase_type,
                      prod_dbase) {
  
  # Initialize results
  cp_i <- 0; cp_l95_i <- 0; cp_u95_i <- 0
  
  # Set non-calcifiers to 0 and exit function
  if (!(substrate_class %in% c("CORAL", "CCA"))) {
    return(data.frame(cp_i, cp_l95_i, cp_u95_i))
  }
  
  # Function to pull values from prod_dbase
  get_val <- function(code, morph = NULL, loc = NULL, field) {
    filt <- prod_dbase$SUBSTRATE_CODE == code
    if (!is.null(morph) && !is.na(morph)) filt <- filt & prod_dbase$MORPHOLOGYCODE == morph
    if (!is.null(loc) && !is.na(loc))     filt <- filt & prod_dbase$LOCATIONCODE == loc
    
    matches <- prod_dbase[[field]][filt]
    
    if (length(matches) == 0) stop(paste0("MISSING DATA: [", code, "] Morph:[", morph, "]"))
    if (length(matches) > 1)  stop(paste0("DUPLICATE DATA: Found ", length(matches), " rows for [", code, "]"))
    
    return(matches)
  }
  
  # Prod calculation for CCA
  if (substrate_class == "CCA") {
    loc_param <- if (dbase_type == "Custom") location_code else NULL
    g    <- get_val(substrate_code, field = "EXTENSION_CM_YR", loc = loc_param)
    g_ci <- get_val(substrate_code, field = "EXTENSION_CM_YR_CI", loc = loc_param)
    d    <- get_val(substrate_code, field = "DENSITY_G_CM3", loc = loc_param)
    
    cp_i     <- substrate_cover_cm * g * d
    cp_l95_i <- substrate_cover_cm * (g - g_ci) * d
    cp_u95_i <- substrate_cover_cm * (g + g_ci) * d
  }
  
  # Prod calculation for CORAL
  if (substrate_class == "CORAL") {
    loc_param <- if (dbase_type == "Custom") location_code else NULL
    df <- data.frame(x = seq(0, 135, 1))
    
    # --- NCRMP/Custom Math for LC ---
    if (morphology_code == "LC" && dbase_type != "IPRB") {
      g_la    <- get_val(substrate_code, "LC-LA", loc_param, "EXTENSION_CM_YR")
      g_la_ci <- get_val(substrate_code, "LC-LA", loc_param, "EXTENSION_CM_YR_CI")
      g_co    <- get_val(substrate_code, "LC-CO", loc_param, "EXTENSION_CM_YR")
      g_co_ci <- get_val(substrate_code, "LC-CO", loc_param, "EXTENSION_CM_YR_CI")
      d       <- get_val(substrate_code, "LC-CO", loc_param, "DENSITY_G_CM3")
      d_ci    <- get_val(substrate_code, "LC-CO", loc_param, "DENSITY_G_CM3_CI")
      c       <- get_val(substrate_code, "LC-CO", loc_param, "CONVERSION_FACTOR")
      
      prop_la <- 0.55; prop_co <- 0.45; h <- 2
      
      df$y <- (h * g_la * d + 0.1 * g_la * ((prop_la * df$x) + (h * g_la)) * d) + 
        ((((1 - c) * (prop_co * df$x)) * g_co * 0.1) + (c * (prop_co * df$x) * g_co)) * d
      
      df$y_l95 <- (h * (g_la-g_la_ci) * (d-d_ci) + 0.1*(g_la-g_la_ci) * ((prop_la*df$x)+(h*(g_la-g_la_ci))) * (d-d_ci)) + 
        ((((1 - c) * (prop_co * df$x)) * (g_co-g_co_ci) * 0.1) + (c*(prop_co*df$x)*(g_co-g_co_ci))) * (d-d_ci)
      
      df$y_u95 <- (h * (g_la+g_la_ci) * (d+d_ci) + 0.1*(g_la+g_la_ci) * ((prop_la*df$x)+(h*(g_la+g_la_ci))) * (d+d_ci)) + 
        ((((1 - c) * (prop_co * df$x)) * (g_co+g_co_ci) * 0.1) + (c*(prop_co*df$x)*(g_co+g_co_ci))) * (d+d_ci)
      
    } else {
      # --- Calculations for all other corals ---
      g    <- get_val(substrate_code, morphology_code, loc_param, "EXTENSION_CM_YR")
      g_ci <- get_val(substrate_code, morphology_code, loc_param, "EXTENSION_CM_YR_CI")
      d    <- get_val(substrate_code, morphology_code, loc_param, "DENSITY_G_CM3")
      d_ci <- get_val(substrate_code, morphology_code, loc_param, "DENSITY_G_CM3_CI")
      c    <- get_val(substrate_code, morphology_code, loc_param, "CONVERSION_FACTOR")
      
      if (morphology_code %in% c("MD", "ML", "FR")) {
        df$y     <- d * (((((df$x*2/pi)/2)+g)^2 * pi/2) - (((df$x*2/pi)/2)^2 * pi/2))
        df$y_l95 <- (d-d_ci) * (((((df$x*2/pi)/2)+(g-g_ci))^2 * pi/2) - (((df$x*2/pi)/2)^2 * pi/2))
        df$y_u95 <- (d+d_ci) * (((((df$x*2/pi)/2)+(g+g_ci))^2 * pi/2) - (((df$x*2/pi)/2)^2 * pi/2))
      } else if (morphology_code %in% c("EM", "PL", "FO", "TB")) {
        h <- 2
        df$y     <- h*g*d + 0.1*g*(df$x+(h*g))*d
        df$y_l95 <- h*(g-g_ci)*(d-d_ci) + 0.1*(g-g_ci)*(df$x+(h*(g-g_ci)))*(d-d_ci)
        df$y_u95 <- h*(g+g_ci)*(d+d_ci) + 0.1*(g+g_ci)*(df$x+(h*(g+g_ci)))*(d+d_ci)
      } else {
        df$y     <- ((((1-c)*df$x)*g*0.1) + (c*df$x*g)) * d
        df$y_l95 <- ((((1-c)*df$x)*(g-g_ci)*0.1) + (c*df$x*(g-g_ci))) * (d-d_ci)
        df$y_u95 <- ((((1-c)*df$x)*(g+g_ci)*0.1) + (c*df$x*(g+g_ci))) * (d+d_ci)
      }
    }
    
    # Final Linear Modeling
    df[1, ] <- 0
    cp_i     <- coef(lm(y ~ x, data = df))[2] * substrate_cover_cm + coef(lm(y ~ x, data = df))[1]
    cp_l95_i <- coef(lm(y_l95 ~ x, data = df))[2] * substrate_cover_cm + coef(lm(y_l95 ~ x, data = df))[1]
    cp_u95_i <- coef(lm(y_u95 ~ x, data = df))[2] * substrate_cover_cm + coef(lm(y_u95 ~ x, data = df))[1]
  }
  
  return(data.frame(cp_i, cp_l95_i, cp_u95_i))
}
