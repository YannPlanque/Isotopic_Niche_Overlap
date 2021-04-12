##--------------------------------------------------------------------------------------------------------
## SCRIPT : Identification of isotopic niche overlap between two species considering inter- and intra-individual variability in isotopic composition
## Specific content : - Fitting a model to identify isotopic niches of two seal species with isotopic data measured in individual whisker (i.e. with intra-individual variability)
##                    - Identifying isotopic niches with ellipses (similarly as Jackson et al.'s (2011) study)
##                    - Quantifying and characterising interspecific overlap of isotopic niches
##
## As part of : 
##    Planque et al. "Trophic niche overlap between sympatric harbour seals (Phoca vitulina) and grey seals (Halichoerus grypus)
##     at their Southern European limit range (Eastern English Channel)". Submitted in Ecology & Behaviour journal
##        >> Preprint version (1st version submitted in 2020-11): https://doi.org/10.22541/au.160508195.50224560/v1
##        >> Peer review status (2021-03-31): minor revisions
##
## Authors : Yann Planque(1)*, Matthieu Authier(2)(3)
## Affiliations : 
##    (1) Centre d'Etudes Biologiques de Chizé (CEBC, UMR 7372 CNRS - La Rochelle Université), La Rochelle, France
##    (2) Observatoire PELAGIS (UMS 3462 CNRS - La Rochelle Université), La Rochelle, France
##    (3) ADERA, Pessac, France
##
## Contact* : yann.planque@univ-lr.fr ; yann.planque@hotmail.fr
##
## First publication on GitHub : 2020-11-05
## Last update : 2021-04-12 (Version 1.3)
##
##
## R version 4.0.2 (2020-06-22) -- "Taking Off Again"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##-------------------------------------------------.-------------------------------------------------------


### 0 // Install & set up rstan ###########################################################################
# If rstan is not already installed on your computer, please follow these instructions:
# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

# If you work with 4.0 version of R, it may be required to make little modifications prior to use rstan, as some problem from recent R versions are not fixed for now.
library(withr) # Package 2.2.0 to run stan model properly https://cran.r-project.org/web/packages/withr/index.html

### Download Rtools 4.0
### /!\ for windows, in 'C:\...\Documents\.R\Makevars.win' modify the file with :
    # CXX14FLAGS += -O3 -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2
    # CXX14 = C:/Rtools40/mingw64/bin/g++ 

###########################################################################################################


### 0 // Packages ##########################################################################################

lapply(c("coda", "mgcv", "mvtnorm", "dplyr", "tidyr",
         "reshape", "ggplot2", "ggthemes", "ggrepel",
         "StanHeaders", "rstan", "sp", "raster",
         "sf", "viridis", "ggpubr", "MASS"), library, character.only=TRUE)
library(withr)

###########################################################################################################


### 0 // Set up rstan ######################################################################################

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
rm(list = ls())

###########################################################################################################


### 0 // Useful functions #################################################################################
# Get confidence intervals (default: CI95%)
get_ci <- function(x, alpha = 0.95){
  c(coda::HPDinterval(coda::as.mcmc(x), prob = alpha)[1], mean(x), coda::HPDinterval(coda::as.mcmc(x), prob = alpha)[2])
}

# Lower and upper bounds for confidence interval (CI, default 95%)
lower <- function(x, alpha = 0.95) { as.numeric(coda::HPDinterval(coda::as.mcmc(x), prob = alpha)[1]) }
upper <- function(x, alpha = 0.95) { as.numeric(coda::HPDinterval(coda::as.mcmc(x), prob = alpha)[2]) }


# Sample in a normal distribution, with lower and upper limits
mysamp <- function(n, m, s, lwr, upr) {
  samp <- rnorm(n, m, s)
  samp[samp < lwr] <- lwr
  samp[samp > upr] <- upr
  samp
}

# Generate posterior isotopic ellipses from the model to use SIBER indices
make_posteriorEllipse4SIBER <- function(stanfit, n_draws = 1000) {
  mu <- rstan::extract(stanfit, 'mu')$mu
  Omega <- rstan::extract(stanfit, 'Omega')$Omega
  if(n_draws < nrow(mu)) {
    writeLines(paste("Returning a posterior object of", n_draws, "lines", sep = " "))
    samp <- sample.int(nrow(mu), size = n_draws)
  }
  else {
    samp <- 1:nrow(mu)
    writeLines(paste("Returning a posterior object of", nrow(mu), "lines", sep = " "))
  }
  return(list("1.1" = as.matrix(cbind(Omega[samp, 1, 1, 1], Omega[samp, 1, 1, 2], Omega[samp, 1, 2, 1], Omega[samp, 1, 2, 2], mu[samp, 1, 1], mu[samp, 1, 2])),
              "2.1" = as.matrix(cbind(Omega[samp, 2, 1, 1], Omega[samp, 2, 1, 2], Omega[samp, 2, 2, 1], Omega[samp, 2, 2, 2], mu[samp, 2, 1], mu[samp, 2, 2]))
  )
  )
}

# Generate predict ellipses to plot them
pred_ell <- function(n_sim = 500, n_pts = 100, stanfit, alpha = 0.95) {
  
  mu <- rstan::extract(stanfit, "mu")$mu
  Omega <- rstan::extract(stanfit, "Omega")$Omega
  
  y_pred <- NULL
  for (z in 1:2){
    spp <- z
    SIM <- seq(1, n_sim*2)
    
    x_pred <- NULL
    for (j in 1:n_sim){
      
      if (spp==1) {
        SIM_selec <- SIM[n_sim+j]
      }else{
        SIM_selec <- SIM[j]
      }
      
      x <- ellipse::ellipse(x = as.matrix(Omega[j, spp, , ]), 
                            centre = mu[j, spp, ], npoints = n_pts
      )
      
      x <- as.data.frame(x); names(x) <- c("d13C", "d15N")
      x$param <- paste("iter", SIM_selec, sep = "")
      x$species <- spp
      x_pred <- rbind(x_pred, x)
      
    }
    
    y_pred <- rbind(y_pred, x_pred)
  }
  return(y_pred)  
}

###########################################################################################################


### I // Data #############################################################################################
  ## 1 / Direction 
  Direction <- ".../Script_Planque_et_al_Isotopic_Niche_Overlap_V1-3"
  #Direction <- "C:/Users/yplanq01/Documents/CEBC/Article_Planque_et_al_Niche_Overlap/Scripts/Isotopic_Niche_Overlap/Script_Planque_et_al_Isotopic_Niche_Overlap_V1-3"


  ## 2 / Import data
  #  Data available on SEANOE: 
  #   Planque Yann, Vincent Cécile, Guillou Gaël, Lebreton Benoit, Caurant Florence (2020).
  #   δ13C and δ15N stable isotope compositions of the whisker of 8 harbour seals (Phoca vitulina)
  #   and 10 grey seals (Halichoerus grypus) captured in the baie de Somme, France, in 2008 and 2012,
  #   for telemetry tracking. SEANOE. https://doi.org/10.17882/76528 

  #   Pv : Phoca vitulina (harbour seals)  //  Hg : Halichoerus grypus (grey seals)
  iso <- read.csv(file = paste(Direction, "Input", "SI_Seal_whiskers_BDS_data.csv", sep = "/"), header = TRUE, dec = ".", sep = ";")
  head(iso)
  

  for(j in 1:ncol(iso)) {
    if(is.factor(iso[, j])) { iso[, j] <- as.character(iso[, j]) }
  }; rm(j)
  
  
  ## 3 / Prepare data
  # Remove lacking data
  iso <- iso %>% dplyr::filter(!is.na(d13C), !is.na(d15N))  
  
  # Select used columns & rename them
  iso0 <- iso %>% dplyr::select(Species, Seal_ID, Segment_mm, d13C, d15N) %>%
    dplyr::rename(spp = Species, IND = Seal_ID, SEG = Segment_mm)
  
  colnames(iso0)
  # Numeric values for input model
  iso_data <- iso0
  iso_data$spp <- as.numeric(as.factor(iso_data$spp))
  iso_data$IND <- as.numeric(as.factor(iso_data$IND))
  
  iso_data <- sort_df(iso_data, vars = c("IND", "SEG"))
  
  theme_set(theme_bw(base_size = 14))
  
  
  ## 4 / Visualise raw data
  # Add label for each individual (according to their capture date)
  iso0$Ind_label <- factor(
                    ifelse(iso0$spp == "Pv",
                           "Harbour seals\n(--- Oct. 2008)",
                           ifelse(iso0$IND %in% c("G01", "G02", "G03", "G04",
                                                  "G05", "G07", "G08", "G09"),
                           "Grey seals\n(--- May-Jun. 2012)",
                           "Grey seals\n(--- Oct. 2012)")),
                    levels=c("Harbour seals\n(--- Oct. 2008)", "Grey seals\n(--- May-Jun. 2012)", "Grey seals\n(--- Oct. 2012)"))
  
  # Create individual labels for raw SI plot
  Labels_inds_plot <- 
  as.data.frame(
  left_join(iso0,
  iso0 %>% gather(key = "Isotope", value = "SI_measure",d13C, d15N) %>%
    group_by(spp, IND, Isotope) %>%
    summarise(SEG_max=max(SEG))) %>%
    filter(SEG == SEG_max) %>%
    dplyr::select(spp, IND, SEG, d13C, d15N, Ind_label) %>%
    mutate(spp = factor(spp, levels=c("Pv", "Hg"))))
  Labels_inds_plot <- Labels_inds_plot[seq(2, nrow(Labels_inds_plot), 2),]
  

  # Plot of d13C values along individual whiskers
  d13C_raw_plot <- 
  iso0 %>% mutate(spp = factor(spp, levels=c("Pv", "Hg"))) %>%
  ggplot(aes(x = -SEG, y = d13C, color = spp, fill = spp, group = IND)) +
    geom_line()+
    geom_point(aes(shape=IND), size=2.5) +
    geom_vline(xintercept=0, linetype="dashed", color = "black") +
    geom_text_repel(data = Labels_inds_plot,
                     mapping = aes(x = -SEG, y = d13C,
                                   label = IND),
                     #fill = 'white',
                     segment.colour = "black",
                     color="black",
                     #min.segment.length = unit(5, 'lines'),
                     box.padding = unit(0.75, "lines"), #3   2.9
                     size=3.5,
                     #label.size = 0.01,
                     fontface = 'bold',
                     segment.size = 0.5, segment.alpha = NULL,
                     direction = "both",
                     force=0.75, nudge_x = -13, nudge_y = 0) + #, nudge_y = 3  min.segment.length = unit(1, 'lines')    geom_line() +
    facet_grid(~Ind_label) +
    ### axis & legend
    xlab(" ") + #Whisker length in mm (tip to base)
    ylab(quote(delta^13*C)) +
    scale_x_continuous(labels = abs,  limits=c(-145,0), breaks=seq(-140, 0, by=20))+
    scale_shape_manual(values = c(19, 15, 17, 18, 0, 5, 2, 6,
                                  19, 18, 
                                  19, 15, 17, 18, 0, 5, 2, 6)) +
    scale_color_manual(name = "Seal species", values = c("#4DE600",  "#0000DE"),
                      labels = c("Harbour seals", "Grey seals")) +
    scale_fill_manual(name = "Seal species", values = c("#4DE600",  "#0000DE"),
                       labels = c("Harbour seals", "Grey seals")) +
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(size = 10, color="black")) +
    guides(fill = "none",
          color = "none", #guide_legend(order = 1)
          shape = "none")
  
  d13C_raw_plot
  
  WID <- 15
  HEI <- 4.75
  ggsave(d13C_raw_plot,
         filename = paste(Direction, "Plot", "00_SI_whisker_d13C.png", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )
  
  # Plot of d15N values along individual whiskers
  d15N_raw_plot <- 
    iso0 %>% mutate(spp = factor(spp, levels=c("Pv", "Hg"))) %>%
    ggplot(aes(x = -SEG, y = d15N, color = spp, fill = spp, group = IND)) +
    geom_line()+
    geom_point(aes(shape=IND), size=2.5) +
    geom_vline(xintercept=0, linetype="dashed", color = "black") +
    geom_text_repel(data = Labels_inds_plot,
                    mapping = aes(x = -SEG, y = d15N,
                                  label = IND),
                    #fill = 'white',
                    segment.colour = "black",
                    color="black",
                    #min.segment.length = unit(5, 'lines'),
                    box.padding = unit(0.75, "lines"), #3   2.9
                    size=3.5,
                    #label.size = 0.01,
                    fontface = 'bold',
                    segment.size = 0.5, segment.alpha = NULL,
                    direction = "both",
                    force=0.75, nudge_x = -12, nudge_y = 0.15) + #, nudge_y = 3  min.segment.length = unit(1, 'lines')    geom_line() +
    facet_grid(~Ind_label) +
    ### axis & legend
    xlab("Whisker length in mm (tip to base)") +
    ylab(quote(delta^15*N)) +
    scale_x_continuous(labels = abs,  limits=c(-145,0), breaks=seq(-140, 0, by=20))+
    scale_shape_manual(values = c(19, 15, 17, 18, 0, 5, 2, 6,
                                  19, 18, 
                                  19, 15, 17, 18, 0, 5, 2, 6)) +
    scale_color_manual(name = "Seal species", values = c("#4DE600",  "#0000DE"),
                       labels = c("Harbour seals", "Grey seals")) +
    scale_fill_manual(name = "Seal species", values = c("#4DE600",  "#0000DE"),
                      labels = c("Harbour seals", "Grey seals")) +
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(size = 10, color="black")) +
    guides(fill = "none",
           color = "none", #guide_legend(order = 1)
           shape = "none")
  
  d15N_raw_plot
  
  WID <- 15
  HEI <- 4.75
  ggsave(d15N_raw_plot,
         filename = paste(Direction, "Plot", "00_SI_whisker_d15N.png", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )
  
  
  # Plot of d13C & d15N values along individual whiskers
  SI_raw_plot <- ggarrange(d13C_raw_plot, d15N_raw_plot, nrow=2, align = "v")
  SI_raw_plot
  
  WID <- 15
  HEI <- 8.5
  ggsave(SI_raw_plot,
         filename = paste(Direction, "Plot", "00_SI_whisker_SI.png", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )
  
  ggsave(SI_raw_plot,
         filename = paste(Direction, "Plot", "00_SI_whisker_SI_PDF.pdf", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )
  
  
  ## 5 / Summarise data
  ID <- as.data.frame(
  iso_data %>% group_by(IND, spp) %>%
    summarize(nb_seg = n(), d13C_mean = mean(d13C), d15N_mean = mean(d15N), rho = cor(d13C, d15N)) )
  
  ID$n_ind_spp <- sapply(ID$spp, function(x) {nrow(subset(ID, spp == x))} )
  
###########################################################################################################


### II // Isotopic niche identification ###################################################################
######### Bayesian model taking intra- and inter-individual variability in iso values #####################

  ## 1 / The model (with rstan)
  hmodel <- '
  /*  Variable naming: 
  n_ind     = number of individuals 
  n_obs     = number of obs (isotopic measurements)
  n_species = number of species
  n_isotopes= number of isotopes
  ISO       = matrix of isotopic measurements
  SPECIES1  = species id for each ind
  SPECIES2  = species id for each measurement
  n_ind_samp= nb of measurement for each ind.
  */ 
  data { 
  int<lower = 1> n_ind;  
  int<lower = 1> n_species; 
  int<lower = 1> n_isotopes;
  int<lower = 1> n_obs;  // number of isotopic measurements
  int<lower = 1> n_ind_samp[n_ind];
  int<lower = 1, upper = n_ind> IND[n_obs];
  int<lower = 1, upper = n_species> SPECIES1[n_ind];
  int<lower = 1, upper = n_species> SPECIES2[n_obs];
  vector[n_isotopes] ISO[n_obs];
  vector[n_isotopes] prior_scale_mu;
  vector<lower = 0.0>[n_isotopes] prior_scale_ind;
  vector<lower = 0.0>[n_isotopes] prior_scale_res;
  }
  
  transformed data {
  real nu;
  real df;
  vector[n_isotopes] prior_var_ind;
  vector[n_isotopes] prior_var_res;
  nu = 2.0;
  df = n_isotopes + nu - 1;
  prior_var_ind = square(prior_scale_ind);
  prior_var_res = square(prior_scale_res);
  }
  
  parameters { 
  vector[n_isotopes] unscaled_mu[n_species];  
  vector[n_isotopes] alpha[n_ind];
  cov_matrix[n_isotopes] Sigma[n_species];      // residual-level covariance matrix for each species
  cov_matrix[n_isotopes] Omega[n_species];      // individual-level covariance matrix for each species
  vector<lower = 0.0>[n_isotopes] a[n_species]; // related to cov. matrix
  vector<lower = 0.0>[n_isotopes] b[n_species]; // related to cov. matrix
  } 
  
  transformed parameters {
  vector[n_isotopes] mu[n_species];
  vector[n_isotopes] sigma_res[n_species]; 
  vector[n_isotopes] sigma_ind[n_species];
  vector[n_isotopes] iso_hat[n_ind];
  for (j in 1:n_isotopes) {
  for (k in 1:n_species) {
  mu[k, j] = unscaled_mu[k, j] * prior_scale_mu[j];
  sigma_ind[k, j] = sqrt(Omega[k, j, j]);
  sigma_res[k, j] = sqrt(Sigma[k, j, j]);
  }
  for (i in 1:n_ind) {
  iso_hat[i, j] = mu[SPECIES1[i], j] + alpha[i, j];
  }
  }
  } 
  
  model { 
  // Priors
  // effets individuels
  for (i in 1:n_ind) {
  alpha[i] ~ multi_student_t(n_ind_samp[i], rep_vector(0.0, 2), Omega[SPECIES1[i]]);
  }
  // moyennes isotopiques
  for (k in 1:n_species) {
  unscaled_mu[k] ~ normal(0.0, 1.0);
  }
  
  for (k in 1:n_species) {
  a[k] ~ inv_gamma(0.5, 1.0); // ind. level covariance: Huang & Wand parametrization described in Bayesian Analysis 2013
  b[k] ~ inv_gamma(0.5, 1.0); // res. level covariance
  // degrees of freedom chosen to ensure marginal uniform distributions on corr. parameters
  Omega[k] ~ inv_wishart(df, diag_matrix(2 * a[k] .* prior_var_ind)); 
  Sigma[k] ~ inv_wishart(df, diag_matrix(2 * b[k] .* prior_var_res)); 
  }
  
  // Likelihood
  for (l in 1:n_obs) {
  ISO[l] ~ multi_normal(iso_hat[IND[l]], Sigma[SPECIES2[l]]);
  }
  }
  '

  ## 2 / Compiling model ...
  mystanmodel <- stan_model(model_code = hmodel, model_name = "Isotopes")

  
  ## 3 / Fitting model ... 
  fit <- sampling(mystanmodel, 
                  data = list(n_ind = nrow(ID),
                              n_obs = nrow(iso_data),
                              n_species = length(unique(iso_data$spp)),
                              n_isotopes = 2,
                              ISO = as.matrix(iso_data[, c("d13C", "d15N")]),
                              IND = iso_data$IND,
                              SPECIES1 = ID$spp,
                              SPECIES2 = iso_data$spp,
                              n_ind_samp = ID$nb_seg,
                              prior_scale_mu = c(20, 20),
                              prior_scale_ind = c(1, 1),
                              prior_scale_res = c(1, 1)
                  ),
                  pars = c("sigma_res", "sigma_ind", "alpha", "mu", "Omega", "Sigma", "iso_hat"),
                  iter = 2000,
                  warmup = 1000,
                  thin = 1,
                  chains = 4
  )
  
  
  ## 4 / Convergence
  stan_rhat(fit)
  traceplot(fit, pars = "mu")
  get_elapsed_time(fit)

  print(fit, dig = 3)
  
  save.image(paste(Direction, "Output", "hmodel_iso_Overlap1.RData", sep = "/"), safe = TRUE)
  #load(paste(Direction, "Output", "hmodel_iso_Overlap1.RData", sep = "/"))


  ## 5 / Extract posterior values (for observed inds.)
  
  # Posterior mean for observed individuals
  # d13C
  ID$pred_d13C <- apply(rstan::extract(fit, "iso_hat")$iso_hat[, ID$IND, 1], 2, mean)
  
  # d15N
  ID$pred_d15N <- apply(rstan::extract(fit, "iso_hat")$iso_hat[, ID$IND, 2], 2, mean)

  
  # Posterior sd for observed individuals
  # d13C
  ID$pred_d13C_SD <- apply(rstan::extract(fit, "iso_hat")$iso_hat[, ID$IND, 1], 2, sd)
  # d15N
  ID$pred_d15N_SD <- apply(rstan::extract(fit, "iso_hat")$iso_hat[, ID$IND, 2], 2, sd)

  
  # Posterior CI95% for observed individuals
  # lower boundary
  ID$pred_d13C_low <- apply(rstan::extract(fit, "iso_hat")$iso_hat[, ID$IND, 1], 2, lower, alpha = 0.95)
  # lower boundary
  ID$pred_d15N_low <- apply(rstan::extract(fit, "iso_hat")$iso_hat[, ID$IND, 2], 2, lower, alpha = 0.95)
  
  # higher boundary
  ID$pred_d13C_up <- apply(rstan::extract(fit, "iso_hat")$iso_hat[, ID$IND, 1], 2, upper, alpha = 0.95)
  # higher boundary
  ID$pred_d15N_up <- apply(rstan::extract(fit, "iso_hat")$iso_hat[, ID$IND, 2], 2, upper, alpha = 0.95)

  

  ## 6 / Extract interindividual variance in species' isotopic niches (posterior omega parameter)
  # For d13C
  Omega_Hg_d13C <- sqrt(rstan::extract(fit, "Omega")$Omega[,, 1 , 1][,1])
  Omega_Pv_d13C <- sqrt(rstan::extract(fit, "Omega")$Omega[,, 1 , 1][,2])
  
  get_ci(Omega_Hg_d13C)
  get_ci(Omega_Pv_d13C)

  # For d15N
  Omega_Hg_d15N <- sqrt(rstan::extract(fit, "Omega")$Omega[,, 2 , 2][,1])
  Omega_Pv_d15N <- sqrt(rstan::extract(fit, "Omega")$Omega[,, 2 , 2][,2])
  
  get_ci(Omega_Hg_d15N)
  get_ci(Omega_Pv_d15N)
  
  # Probabilities
  # Probability of interindividual variability in d13C to be higher for Pv (harbour seals) than for Hg (grey seals)
  mean(ifelse(Omega_Hg_d13C < Omega_Pv_d13C, 1, 0))
  
  # Probability of interindividual variability in d15N to be higher for Hg than for Pv
  mean(ifelse(Omega_Hg_d15N > Omega_Pv_d15N, 1, 0))

###########################################################################################################
  
  
### III // Isotopic niches with ellipses - quantify the overlap and plot niches  ##########################
  
  ## 1 / Generate ellipses from posterior parameters from the model
  # 1000 ellipses by species
  post <- make_posteriorEllipse4SIBER(stanfit = fit, n_draws = 1000)
  

  ## 2 / Overlap index (in SIBER package) (/!\ can last several minutes, ~ 3 min for 1000 ellipses)
  system.time(
    overlap <- SIBER::bayesianOverlap(ellipse1 = "1.1", ellipse2 = "2.1", ellipses.posterior = post, draws = NULL)
  )
  
  save.image(paste(Direction, "Output", "hmodel_iso_Overlap1.RData", sep = "/"), safe = TRUE)
  #load(paste(Direction, "Output", "hmodel_iso_Overlap1.RData", sep = "/"))
  
  # Summarise overlap results...
    # A - Isotopic niche size [in ‰]
    get_ci(overlap$area1, alpha = 0.95) # area for species 1 (here Hg) 
    get_ci(overlap$area2, alpha = 0.95) # area for species 2 (here Pv) 
    get_ci(overlap$overlap, alpha = 0.95) # overlapping area
    
    ggplot(overlap) +
      geom_density(aes(area1), colour="#0000DE") +
      geom_density(aes(area2), colour="#4DE600") +
      labs(x="Area in ‰", y="Density") +
      geom_vline(aes(xintercept=median(area1)),
                 colour="#0000DE",linetype="dashed") +
      geom_vline(aes(xintercept=median(area2)),
                 colour="#4DE600",linetype="dashed")
    
    # Probability of the niche of species 1 (Hg) to be larger than those of species 2 (Pv)
    mean(ifelse(overlap$area1 > overlap$area2, 1, 0))
    
    
    # B - Percentage of ellipses area overlapping
    proportion_overlap <- with(overlap, overlap / (area1 + area2 - overlap))
    hist(proportion_overlap*100)
    get_ci(proportion_overlap*100, alpha = 0.95)
    
    
    # C - Percentage of the niche of species 1 nested in those of species 2, and vice-versa
    # SPP 1 in SPP 2
    get_ci(overlap$overlap / overlap$area1, alpha = 0.95)*100
    
    # SPP 2 in SPP 1
    get_ci(overlap$overlap / overlap$area2, alpha = 0.95)*100
    
    # Probability of (SPP 2 in SPP 1) > (SPP 1 in SPP 2)
    mean(ifelse((overlap$overlap / overlap$area2) > (overlap$overlap / overlap$area1), 1, 0))
    
    
  ## 3 / Plot isotopic niches
    
    # A - Generate all ellipses (100 points by ellipses)
    my_pred <- pred_ell(stanfit = fit, n_sim = 1000, n_pts = 100, alpha = 0.95) # all ellipses
    my_pred$species <- factor(ifelse(my_pred$species==1, "Hg", "Pv"), levels=c("Pv", "Hg")) # rename and organise spp

    # Prepare ID for observed individuals
    ID$spp <- factor(ifelse(ID$spp==1, "Hg", "Pv"), levels=c("Pv", "Hg"))

    
    # B - Plot all ellipses characterising isotopic niches
    COLOR_PV <- "#4DE600"
    COLOR_HG <-"#0000DE"
    
    FILL_PV <- "#4DE600"
    FILL_HG <- "#0071FE"
    
    Plot_01_ellipses <- 
    ggplot(data = ID, aes(x = pred_d13C, y = pred_d15N)) + 
      # ellipses
      geom_path(data = my_pred, aes(x = d13C, y = d15N, color = species, group=param), alpha = 0.215, size=0.215) +
      # points
      geom_errorbarh(aes(xmax = pred_d13C_up, xmin = pred_d13C_low, height = 0.0), color="black", size=0.7) +
      geom_errorbar(aes(ymax = pred_d15N_up, ymin = pred_d15N_low, width = 0.0), color="black", size=0.7) +
      geom_point(aes(fill=spp), shape=21, size=3) +
      
      ### axis and legend
      xlab(quote(delta^13*C)) + ylab(quote(delta^15*N)) + 

      scale_color_manual(name = "Species", values = c(COLOR_PV, COLOR_HG),
                         labels = c("Harbour seals", "Grey seals")) +
      scale_fill_manual(name = "Species", values = c(FILL_PV, FILL_HG), 
                        labels = c("Harbour seals", "Grey seals")) +
      guides(fill = FALSE, color = FALSE
      ) +
      theme_bw(base_size = 14) +
      theme(legend.position = "right", 
            plot.title = element_text(lineheight = 0.8, face = "bold"), 
            axis.text = element_text(size = 14, color="black"),
            axis.title = element_text(size = 18)) +
      scale_y_continuous(breaks=seq(12,22,2)) + 
      scale_x_continuous(breaks=seq(-20, -8 ,2))
    
    Plot_01_ellipses
    
    WID <- 7
    HEI <- 6
    ggsave(Plot_01_ellipses,
           filename = paste(Direction, "Plot", "01_SI_niches_ellipses.png", sep = "/"), dpi = 300,
           width = WID, height = HEI
    )
    
    ggsave(Plot_01_ellipses,
           filename = paste(Direction, "Plot", "01_SI_niches_ellipses_PDF.pdf", sep = "/"), dpi = 300,
           width = WID, height = HEI
    )

    
###########################################################################################################
    
    
### IV // Probability belonging to the isotopic niche  ####################################################

  ## 1 / Generate uniform points in the isotopic space, each 0.05‰ for both d13C and d15N
  SEQ_d13C <- seq(floor(min(my_pred$d13C)), ceiling(max(my_pred$d13C)), 0.05)
  SEQ_d15N <- seq(floor(min(my_pred$d15N)), ceiling(max(my_pred$d15N)), 0.05)
  Points_unif <- NULL
  for (m in 1: length(SEQ_d15N)){
    SEQ_d15N[m] 
    Points_unif <- rbind(Points_unif, data.frame(d13C = SEQ_d13C, d15N = rep(SEQ_d15N[m], length(SEQ_d13C))))
  }

  nrow(Points_unif) # Number of  points in the isotopic space
  
  
  ## 2 / Spatialise these points for further analyses
  # Objective : selecting points in an ellipse
  coordinates(Points_unif) <- ~ d13C + d15N
  length(Points_unif) == length(SEQ_d15N) * length(SEQ_d13C)
  
  
  # Select each points in the isotopic ellipses (i.e. isotopic niches)
  # To ultimately calculate the probability of each point to be in all ellipses of a species
  # // ! \\ Last several minutes
  species_unique <- unique(my_pred$species)
  Points_unif$Number_point <- rep(0, length(Points_unif)) # Number of ellipses in which there are
  
  # Loops at species + ellipses levels to identify points in each ellipse
  Points_ellipses_ALL <- NULL
  for (i in 1:length(species_unique)){
    SPP <- species_unique[i]
    Points_prov <- Points_unif
    param_unique <- unique(my_pred[my_pred$species == SPP,]$param)
    
    for (j in 1:length(param_unique)){ # Spatial selection of points in the ellipse
      ITER <- param_unique[j]
      Ellipse_data <- my_pred[my_pred$param == ITER & my_pred$species == SPP,]
      
      p <- Polygon(Ellipse_data[,c(1:2)])
      ps <- Polygons(list(p),1)
      Polyg_1 <- SpatialPolygons(list(ps))
      
      # Binary presence / absence of points in the ellipse
      Points_prov$Number_point <- Points_prov$Number_point + ifelse(is.na(over(Points_unif, Polyg_1)), 0, 1)
      
    }
    
    Points_final <- as.data.frame(coordinates(Points_prov))
    Points_final$Number_point <- Points_prov$Number_point
    Points_final$Prob <- Points_final$Number_point / length(param_unique)
    Points_final$Species <- SPP
    #Points_final <- Points_final[Points_final$Number_point > 0,]
    
    Points_ellipses_ALL <- rbind(Points_ellipses_ALL, Points_final)
  }

  
  ## 3 / Plot of contour probabilities (at 5%, 25%, 50%, 75% and 100%)
  Plot_02_probabilities <- 
  ggplot(Points_ellipses_ALL) +
    geom_contour(data= Points_ellipses_ALL, 
                 aes(d13C, d15N, z=Prob, colour = Species), breaks=c(0.05, 0.25,0.5, 0.75, 1))+
    scale_color_manual(name = "Species", values = c(COLOR_PV,  COLOR_HG), breaks = c( "Pv",  "Hg"),
                       labels = c("Harbour seals", "Grey seals"))+
    xlab(quote(delta^13*C)) + ylab(quote(delta^15*N)) + 
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(size = 14, color="black"),
          axis.title = element_text(size = 18)) +
    scale_y_continuous(breaks=seq(12,22,1)) + 
    scale_x_continuous(breaks=seq(-20, -8 ,1))
  
  Plot_02_probabilities
  # ==>  contours of probabilities (belonging to isotopic niches) look like ellipses
  
  
  ALPHA <- 0.15 # fill transparency
  # Define breaks for iso plots
  Plot4_ranges <- 
  ggplot() +
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.05 & Points_ellipses_ALL$Species == "Hg",],
                 mapping = aes(d13C, d15N), fill = COLOR_HG,  color = COLOR_HG, level=0.95, geom="polygon", alpha=ALPHA)+
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.05 & Points_ellipses_ALL$Species == "Pv",],
                 mapping = aes(d13C, d15N), fill = COLOR_PV,  color = COLOR_PV, level=0.95, geom="polygon", alpha=ALPHA)
  Plot4_ranges <- ggplot_build(Plot4_ranges)
  
  Plot4_ranges_out <- rbind(Plot4_ranges[[1]][[1]], Plot4_ranges[[1]][[2]])
  
  Ranges_X <- c(min(Plot4_ranges_out$x), max(Plot4_ranges_out$x))
  Ranges_Y <- c(min(Plot4_ranges_out$y), max(Plot4_ranges_out$y))
  By_breaks <- 1
  Breaks_X <- seq(floor(Ranges_X[1]) - floor(Ranges_X[1]) %% 2,
                  ceiling(Ranges_X[2]) + ceiling(Ranges_X[2]) %% 2,  By_breaks)
  Breaks_Y <- seq(floor(Ranges_Y[1]) - floor(Ranges_Y[1]) %% 2,
                  ceiling(Ranges_Y[2]) + ceiling(Ranges_Y[2]) %% 2,  By_breaks)
  ##

  
  ## 4 / Plot of ranges of probability belonging to isotopic niches
  ALPHA <- 0.15 # fill transparency
  Plot_02_proba_ranges <- 
  ggplot() +
    ### Probability ranges
    ## Grey seals
    # 1- ellipse at 95% around points with a probability >= 0.05
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.05 & Points_ellipses_ALL$Species == "Hg",],
                 mapping = aes(d13C, d15N), fill = COLOR_HG,  color = COLOR_HG, level=0.95, geom="polygon", alpha=ALPHA)+
    # 2- prob >= 0.25
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.25 & Points_ellipses_ALL$Species == "Hg",],
                 mapping = aes(d13C, d15N), fill = COLOR_HG,  color = COLOR_HG, level=0.95, geom="polygon", alpha=ALPHA)+
    # 3- prob >= 0.5
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.5 & Points_ellipses_ALL$Species == "Hg",],
                 mapping = aes(d13C, d15N), fill = COLOR_HG,  color = COLOR_HG, level=0.95, geom="polygon", alpha=ALPHA)+
    # 4- prob >= 0.75
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.75 & Points_ellipses_ALL$Species == "Hg",],
                 mapping = aes(d13C, d15N), fill = COLOR_HG,  color = COLOR_HG, level=0.95, geom="polygon", alpha=ALPHA)+
    # 5- prob >= 1 (i.e. in ALL ellipses)
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob==1 & Points_ellipses_ALL$Species == "Hg",],
                 mapping = aes(d13C, d15N), fill = COLOR_HG,  color = COLOR_HG, level=0.95, geom="polygon", alpha=ALPHA)+
    
    ## Harbour seals
    # 1- prob >= 0.05
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.05 & Points_ellipses_ALL$Species == "Pv",],
                 mapping = aes(d13C, d15N), fill = COLOR_PV,  color = COLOR_PV, level=0.95, geom="polygon", alpha=ALPHA)+
    # 2- prob >= 0.25
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.25 & Points_ellipses_ALL$Species == "Pv",],
                 mapping = aes(d13C, d15N), fill = COLOR_PV,  color = COLOR_PV, level=0.95, geom="polygon", alpha=ALPHA)+
    # 3- prob >= 0.5
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.5 & Points_ellipses_ALL$Species == "Pv",],
                 mapping = aes(d13C, d15N), fill = COLOR_PV,  color = COLOR_PV, level=0.95, geom="polygon", alpha=ALPHA)+
    # 4- prob >= 0.75
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.75 & Points_ellipses_ALL$Species == "Pv",],
                 mapping = aes(d13C, d15N), fill = COLOR_PV,  color = COLOR_PV, level=0.95, geom="polygon", alpha=ALPHA)+
    # 5- prob >= 1 (i.e. in ALL ellipses)
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob==1 & Points_ellipses_ALL$Species == "Pv",],
                 mapping = aes(d13C, d15N), fill = COLOR_PV,  color = COLOR_PV, level=0.95, geom="polygon", alpha=ALPHA)+
    
    ### Observed individuals with pred values
    geom_errorbarh(ID, mapping=aes(y = pred_d15N, xmax = pred_d13C_up, xmin = pred_d13C_low, height = 0.0), color="black", size=0.7) +
    geom_errorbar(ID, mapping=aes(x = pred_d13C, ymax = pred_d15N_up, ymin = pred_d15N_low, width = 0.0), color="black", size=0.7) +
    geom_point(data = ID, mapping=aes(x = pred_d13C, y = pred_d15N, fill = spp), shape=21, size=3) +

    scale_color_manual(name = "Species", values = c(COLOR_PV,  COLOR_HG), breaks = c( "Pv",  "Hg"),
                       labels = c("Harbour seals", "Grey seals")) +
    scale_fill_manual(name = "Species", values = c(FILL_PV, FILL_HG), breaks = c( "Pv",  "Hg"),
                      labels = c("Harbour seals", "Grey seals")) +
    xlab(quote(delta^13*C)) + ylab(quote(delta^15*N)) +
    guides(fill = FALSE, color = FALSE) +
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(size = 14, color="black"),
          axis.title = element_text(size = 18)) +
    scale_x_continuous(limits = Ranges_X, breaks = Breaks_X) + 
    scale_y_continuous(limits = Ranges_Y, breaks = Breaks_Y)
  
  Plot_02_proba_ranges
  
  WID <- 7
  HEI <- 6
  ggsave(Plot_02_proba_ranges,
         filename = paste(Direction, "Plot", "02_SI_niches_proba_ranges.png", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )
  
  ggsave(Plot_02_proba_ranges,
         filename = paste(Direction, "Plot", "02_SI_niches_proba_ranges_PDF.pdf", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )

###########################################################################################################
  
  
### V // Probability of interspecific niche overlap  ######################################################

  ## 1 / Calculate this probability
  #   = probability belonging in the niche of species 1 (Hg) X probability belonging in the niche of species 2 (Pv)

  # Probabilities belong in harbour seals' niche (Pv)
  Points_ellipses_Pv <- Points_ellipses_ALL %>% filter(Species == "Pv") %>% dplyr::select(!Species) %>%
    dplyr::rename(Number_point_Pv = Number_point,  Prob_Pv = Prob)
  
  # Probabilities belong in grey seals' niche (Hg)
  Points_ellipses_Hg <- Points_ellipses_ALL %>% filter(Species == "Hg") %>% dplyr::select(!Species) %>%
    dplyr::rename(Number_point_Hg = Number_point,  Prob_Hg = Prob)
  
  Points_ellipses_over <- left_join(Points_ellipses_Pv, Points_ellipses_Hg, by = c("d13C", "d15N"))
  
  # Proba of interspecific overlap for each point:
  Points_ellipses_over$Proba_overlap <- Points_ellipses_over$Prob_Pv * Points_ellipses_over$Prob_Hg

  
  ## 2 / Plot the probability of overlap
  # Locate the hotspot in isotopic space...
  ggplot(Points_ellipses_over) +
    geom_point(aes(d13C, d15N, color=Proba_overlap)) +
    scale_colour_gradientn(name="Probability\nof overlap",
                           colors=c("white", "lightgoldenrod", "gold1", "orange", "red2", "red4")) +
    xlab(quote(delta^13*C)) + ylab(quote(delta^15*N)) +
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(size = 14, color="black"),
          axis.title = element_text(size = 18)) +
    scale_y_continuous(breaks=seq(12,22,1)) + 
    scale_x_continuous(breaks=seq(-20, -8 ,1))
  
  save.image(paste(Direction, "Output", "hmodel_iso_Overlap1.RData", sep = "/"), safe = TRUE)
  #load(paste(Direction, "Output", "hmodel_iso_Overlap1.RData", sep = "/"))
  
  # Add individual names to the data
  Individuals <- data.frame(Individuals = unique(iso0$IND),
                            IND = unique(as.numeric(factor(iso0$IND))))
  
  ID <- left_join(ID, Individuals) 
  
  # Plot the ranges of probability of overlap:
  FILL_overlap <- "red2"
  COL_overlap <- "red2"
  ALPHA <- 0.15
  
  Plot_03_overlap_ranges <- 
  ggplot() +
    # 1- ellipse at 95% around points with a probability >= 0.05
    stat_ellipse(Points_ellipses_over[Points_ellipses_over$Proba_overlap>=0.05 ,],
                 mapping = aes(d13C, d15N), fill = FILL_overlap,  color = COL_overlap, level=0.95,geom="polygon", alpha=ALPHA) +
    # 2- proba >= 0.25
    stat_ellipse(Points_ellipses_over[Points_ellipses_over$Proba_overlap>=0.25 ,],
                 mapping = aes(d13C, d15N), fill = FILL_overlap,  color = COL_overlap, level=0.95, geom="polygon", alpha=ALPHA) +
    # 3- proba >= 0.5
    stat_ellipse(Points_ellipses_over[Points_ellipses_over$Proba_overlap>=0.5 ,],
                 mapping = aes(d13C, d15N), fill = FILL_overlap,  color = COL_overlap, level=0.95, geom="polygon", alpha=ALPHA) +
    # 4- proba >= 0.75
    stat_ellipse(Points_ellipses_over[Points_ellipses_over$Proba_overlap>=0.75 ,],
                 mapping = aes(d13C, d15N), fill = FILL_overlap,  color = COL_overlap, level=0.95, geom="polygon", alpha=ALPHA) +
    # 5- proba >= 0.95
    stat_ellipse(Points_ellipses_over[Points_ellipses_over$Proba_overlap>=0.95 ,],
                 mapping = aes(d13C, d15N), fill = FILL_overlap,  color = COL_overlap, level=0.95, geom="polygon", alpha=ALPHA) +
    # 6- proba >= 0.99
    stat_ellipse(Points_ellipses_over[Points_ellipses_over$Proba_overlap>=0.99 ,],
                 mapping = aes(d13C, d15N), fill = FILL_overlap,  color = COL_overlap, level=0.95, geom="polygon", alpha=ALPHA) +
    
    ### Observed individuals with pred values
    geom_errorbarh(ID, mapping=aes(y = pred_d15N, xmax = pred_d13C_up, xmin = pred_d13C_low, height = 0.0), color="black", size=0.7) +
    geom_errorbar(ID, mapping=aes(x = pred_d13C, ymax = pred_d15N_up, ymin = pred_d15N_low, width = 0.0), color="black", size=0.7) +
    geom_point(data = ID, mapping=aes(x = pred_d13C, y = pred_d15N, fill = spp), shape=21, size=3) +
    geom_label_repel(data = ID, mapping = aes(x = pred_d13C, y = pred_d15N,
                                                       label = Individuals,
                                              fill = ifelse(spp=="Pv", COLOR_PV, COLOR_HG)),
                              color = 'white', segment.colour = "black",
                              #min.segment.length = unit(5, 'lines'),
                              box.padding = unit(2.9, "lines"), #3   2.9
                              size=5.3,
                              label.size = 0.01,
                              fontface = 'bold',
                              segment.size = 0.5, segment.alpha = NULL,
                              direction = "both",
                              force=2) + #, nudge_y = 3  min.segment.length = unit(1, 'lines')
    scale_fill_manual(name = "Species", values = c(FILL_PV, FILL_HG, COLOR_PV, COLOR_HG), breaks = c( "Pv",  "Hg", COLOR_PV, COLOR_HG),
                      labels = c("Harbour seals", "Grey seals", "Harbour seals", "Grey seals")) +
    scale_color_manual(name = "Species", values = c(COLOR_PV,  COLOR_HG), breaks = c( "Pv",  "Hg"),
                       labels = c("Harbour seals", "Grey seals")) +

    scale_fill_discrete(aesthetics = "segment.colour") +
    xlab(quote(delta^13*C)) + ylab(quote(delta^15*N)) +
    guides(fill = FALSE, color = FALSE) +
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(size = 14, color="black"),
          axis.title = element_text(size = 18)) +
    scale_x_continuous(limits = Ranges_X, breaks = Breaks_X) + 
    scale_y_continuous(limits = Ranges_Y, breaks = Breaks_Y)
  
  Plot_03_overlap_ranges
  
  WID <- 7
  HEI <- 6
  ggsave(Plot_03_overlap_ranges,
         filename = paste(Direction, "Plot", "03_SI_overlap_proba_ranges.png", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )
  
  ggsave(Plot_03_overlap_ranges,
         filename = paste(Direction, "Plot", "03_SI_overlap_proba_ranges_PDF.pdf", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )
  
  
### Export final plot  ####################################################################################
  
  Plot_01_02_03 <- 
      ggarrange(Plot_01_ellipses, Plot_02_proba_ranges, Plot_03_overlap_ranges, align = "v", ncol=1)
  
  Plot_01_02_03
  
  WID <- 7
  HEI <- 6*3
  ggsave(Plot_01_02_03,
         filename = paste(Direction, "Plot", "03_02_01_All_plots.png", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )
  
  ggsave(Plot_01_02_03,
         filename = paste(Direction, "Plot", "03_02_01_All_plots_PDF.pdf", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )

  
  
# save
save.image(paste(Direction, "Output", "hmodel_iso_Overlap1.RData", sep = "/"), safe = TRUE)
#load(paste(Direction, "Output", "hmodel_iso_Overlap1.RData", sep = "/"))
    
###########################################################################################################

  
### VI // OPTIONAL USE OF SUCH RESULTS: comparison with isotopic values of potential prey   ###############
  
### HERE: Comparison with isotopic values of fish and cephalopod species, published by Kopp et al. (2015) [http://dx.doi.org/10.1016/j.pocean.2014.11.001]
### We compared the isotopic values of potential prey (i.e. preferential prey identified in seals' diet) with seal isotopic niches
  
  ## 1 / Import data
  # Data extracted from Kopp et al. (2015) http://dx.doi.org/10.1016/j.pocean.2014.11.001 
  # with the file
  PREY_iso <- read.csv(file = paste(Direction, "Input", "SI_data_Kopp_et_al_EEC.csv", sep = "/"), header = TRUE, dec = ".", sep = ";")
  
  # OR...
  # without the file
  PREY_iso <- data.frame(Prey_latin	= c("Pleuronectes platessa", "Platichthys flesus", "Solea solea", "Microchirus variegatus",
                                             "Buglossidium luteum", "Merlangius merlangus", "Trisopterus luscus",
                                             "Clupea harengus", "Clupea harengus (0-20m)", "Clupea harengus (20-38m)", "Callionymus lyra",
                                             "Loligo vulgaris"),
                         Code = c("PP",	"PF",	"SS",	"MV",	"BL",	"MM",	"TL",	"CH",	"CH1",	"CH2",	"CL",	"LV"),
                         Prey_functionnal = c("Benthic flatfish", "Benthic flatfish", "Benthic flatfish", "Benthic flatfish", "Benthic flatfish", "Demersal fish",
                                              "Demersal fish", "Pelagic fish", "Pelagic fish", "Pelagic fish", "Benthic non flatfish", 	"Pelagic squids"),

                         Depth = c("All", "All", "All", "All", "All", "All", "All", "All", "0-20", "20-38", "All", "All"),
                         d13C_mean = c(-16.61,	-17.38,	-16.75,	-15.65,	-16.79,	-16.57,	-17.34,	-18.05,	-16.74,	-19.02,	-17.31,	-16.82),
                         d13C_SD = c(0.82, 0.52, 0.8, 0.21, 0.56, 0.44, 1.2, 1.72, 0.29, 1.48, 1.22,	0.97),
                         d15N_mean = c(13.41, 13.88, 13.9, 14.28, 13.65, 16.05, 14.78, 13.53, 12.92, 14.47, 13.37,	17.05),
                         d15N_SD = c(1.02, 0.78, 1, 0.2, 0.81, 0.53, 0.86, 1.69, 0.37, 2.54, 0.78, 0.35),
                         N = c(46, 10, 54, 6, 7, 48, 24, 10, 5, 4, 18, 7))
  
  PREY_iso <- PREY_iso[PREY_iso$Code != "CH",] # Remove Clupea harengus for all depths
  
  PREY_iso
  

  ## 2 / Apply a TEF, assuming a potential consumption by seals
  # Apply a trophic enrichment factor (TEF) on data, assuming a theoric consumption
  # Here, we used TEF values determined by Lerner et al. (2018) [https://doi.org/10.1371/journal.pone.0192241] for grey seal whiskers 
      # + 2.4 ± 1.3 for d13C (mean ± sd)
      # + 2.6 ± 1.2 for d15N (mean ± sd)

  TEF_vals <- data.frame(Iso = c("d13C", "d15N"), mean = c(2.4, 2.6), sd = c(1.3, 1.2))
  
  
  # Reconstruct mean and sd of isotopic values + TEF
  # Method: generate normal distributions of isotopic prey data and TEF to reconstruct that
  PREY_iso_TEF <- NULL
  for (i in 1:nrow(PREY_iso)){
    prey_data <- PREY_iso[i,]
    prey_data2 <- prey_data
    
    # Apply TEF for d13C
    Prey_d13C_rnorm <- rnorm(n = 1000, mean = prey_data$d13C_mean, sd = prey_data$d13C_SD) # reconstruct raw d13C distribution
    TEF_d13C_rnorm <- mysamp(n = 1000, m = TEF_vals[TEF_vals$Iso == "d13C",]$mean, s = TEF_vals[TEF_vals$Iso == "d13C",]$sd, lwr = 0, upr = Inf) # reconstruct TEF distribution

    prey_data2$d13C_mean <- mean(Prey_d13C_rnorm + TEF_d13C_rnorm) # calculate mean of d13C+TEF
    prey_data2$d13C_SD <- sd(Prey_d13C_rnorm + TEF_d13C_rnorm) # sd
    
    
    # Apply TEF for d15N
    Prey_d15N_rnorm <- rnorm(n = 1000, mean = prey_data$d15N_mean, sd = prey_data$d15N_SD) # reconstruct raw d15N distribution
    TEF_d15N_rnorm <- mysamp(n = 1000, m = TEF_vals[TEF_vals$Iso == "d15N",]$mean, s = TEF_vals[TEF_vals$Iso == "d15N",]$sd, lwr = 0, upr = Inf) # reconstruct TEF distribution
    
    prey_data2$d15N_mean <- mean(Prey_d15N_rnorm + TEF_d15N_rnorm) # calculate mean of d15N+TEF
    prey_data2$d15N_SD <- sd(Prey_d15N_rnorm + TEF_d15N_rnorm) # sd

    PREY_iso_TEF <- rbind(PREY_iso_TEF, prey_data2)
    rm(Prey_d13C_rnorm, TEF_d13C_rnorm, prey_data, prey_data2, i)
  }
  

  
  ## 3 / Visualise SI values (+TEF) of potential prey on isotopic niches

  ALPHA <- 0.15 # fill transparency
  Plot_04_potential_Prey_Niche <- 
    ggplot() +
    ### Probability ranges
    ## Grey seals
    # 1- ellipse at 95% around points with a probability >= 0.05
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.05 & Points_ellipses_ALL$Species == "Hg",],
                 mapping = aes(d13C, d15N), fill = COLOR_HG,  color = COLOR_HG, level=0.95, geom="polygon", alpha=ALPHA)+
    # 2- prob >= 0.25
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.25 & Points_ellipses_ALL$Species == "Hg",],
                 mapping = aes(d13C, d15N), fill = COLOR_HG,  color = COLOR_HG, level=0.95, geom="polygon", alpha=ALPHA)+
    # 3- prob >= 0.5
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.5 & Points_ellipses_ALL$Species == "Hg",],
                 mapping = aes(d13C, d15N), fill = COLOR_HG,  color = COLOR_HG, level=0.95, geom="polygon", alpha=ALPHA)+
    # 4- prob >= 0.75
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.75 & Points_ellipses_ALL$Species == "Hg",],
                 mapping = aes(d13C, d15N), fill = COLOR_HG,  color = COLOR_HG, level=0.95, geom="polygon", alpha=ALPHA)+
    # 5- prob >= 1 (i.e. in ALL ellipses)
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob==1 & Points_ellipses_ALL$Species == "Hg",],
                 mapping = aes(d13C, d15N), fill = COLOR_HG,  color = COLOR_HG, level=0.95, geom="polygon", alpha=ALPHA)+
    
    ## Harbour seals
    # 1- prob >= 0.05
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.05 & Points_ellipses_ALL$Species == "Pv",],
                 mapping = aes(d13C, d15N), fill = COLOR_PV,  color = COLOR_PV, level=0.95, geom="polygon", alpha=ALPHA)+
    # 2- prob >= 0.25
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.25 & Points_ellipses_ALL$Species == "Pv",],
                 mapping = aes(d13C, d15N), fill = COLOR_PV,  color = COLOR_PV, level=0.95, geom="polygon", alpha=ALPHA)+
    # 3- prob >= 0.5
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.5 & Points_ellipses_ALL$Species == "Pv",],
                 mapping = aes(d13C, d15N), fill = COLOR_PV,  color = COLOR_PV, level=0.95, geom="polygon", alpha=ALPHA)+
    # 4- prob >= 0.75
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob>=0.75 & Points_ellipses_ALL$Species == "Pv",],
                 mapping = aes(d13C, d15N), fill = COLOR_PV,  color = COLOR_PV, level=0.95, geom="polygon", alpha=ALPHA)+
    # 5- prob >= 1 (i.e. in ALL ellipses)
    stat_ellipse(Points_ellipses_ALL[Points_ellipses_ALL$Prob==1 & Points_ellipses_ALL$Species == "Pv",],
                 mapping = aes(d13C, d15N), fill = COLOR_PV,  color = COLOR_PV, level=0.95, geom="polygon", alpha=ALPHA)+
    
    ### Observed individuals with pred values
    geom_errorbarh(data = PREY_iso_TEF, mapping=aes(y = d15N_mean, xmax = d13C_mean + d13C_SD, xmin = d13C_mean - d13C_SD, color = Prey_functionnal , height = 0.0), size=0.8) +
    geom_errorbar(data = PREY_iso_TEF, mapping=aes(x = d13C_mean, ymax = d15N_mean + d15N_SD, ymin = d15N_mean - d15N_SD, color = Prey_functionnal , width = 0.0), size=0.8) +
    geom_point(data = PREY_iso_TEF, mapping=aes(x = d13C_mean, y = d15N_mean, color = Prey_functionnal ), shape=15, size=7)+
    
    scale_color_manual(name = "Potential prey", values = c("goldenrod",  "darkorchid2", "green3", "royalblue4","deepskyblue"),
                       breaks = c( "Benthic flatfish",  "Benthic non flatfish", "Demersal fish", "Pelagic fish", "Pelagic squids"),
                       labels = c( "Benthic flatfish",  "Benthic non flatfish", "Demersal fish", "Pelagic fish", "Pelagic squids"), guide =F
    )+
    scale_fill_manual(name = "Potential prey", values = c( "goldenrod",  "darkorchid2", "green3", "royalblue4","deepskyblue"),
                      breaks = c( "Benthic flatfish",  "Benthic non flatfish", "Demersal fish", "Pelagic fish", "Pelagic squids"),
                      labels = c( "Benthic flatfish",  "Benthic non flatfish", "Demersal fish", "Pelagic fish", "Pelagic squids"), guide =F
    )+
    xlab(quote(delta^13*C)) + ylab(quote(delta^15*N))+
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(size = 14, color="black"),
          axis.title.x = element_text(size = 18)) + 
    scale_x_continuous(limits = c(-18.6, -11.9), breaks=Breaks_X) +
    scale_y_continuous(limits = c( Ranges_Y[1], 21), breaks=Breaks_Y) +
    geom_text(data = PREY_iso_TEF, mapping=aes(x = d13C_mean, y = d15N_mean, label= Code, fontface=2),
              size=3, color="white") +
    labs(caption = "Source of isotopic data for potential prey: Kopp et al., (2015)\nApplied TEF from prey to predator: d13C = + 2.4 ± 1.3‰, d15N = + 2.6 ± 1.4‰ (Lerner et al., 2018)") 
  
  Plot_04_potential_Prey_Niche
  
  WID <- 7
  HEI <- 6
  ggsave(Plot_04_potential_Prey_Niche,
         filename = paste(Direction, "Plot", "04_SI_niche_potential_prey.png", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )
  
  ggsave(Plot_04_potential_Prey_Niche,
         filename = paste(Direction, "Plot", "04_SI_niche_potential_prey_PDF.pdf", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )
  
  
  ## 4 / Compare SI values (+TEF) of potential prey with ranges of overlap probability
  
  ALPHA <- 0.15 # fill transparency
  Plot_05_potential_prey_overlap <- 
    ggplot() +
    ### Probability ranges
    # 1- ellipse at 95% around points with a probability >= 0.05
    stat_ellipse(Points_ellipses_over[Points_ellipses_over$Proba_overlap>=0.05 ,],
                 mapping = aes(d13C, d15N), fill = FILL_overlap,  color = COL_overlap, level=0.95,geom="polygon", alpha=ALPHA) +
    # 2- proba >= 0.25
    stat_ellipse(Points_ellipses_over[Points_ellipses_over$Proba_overlap>=0.25 ,],
                 mapping = aes(d13C, d15N), fill = FILL_overlap,  color = COL_overlap, level=0.95, geom="polygon", alpha=ALPHA) +
    # 3- proba >= 0.5
    stat_ellipse(Points_ellipses_over[Points_ellipses_over$Proba_overlap>=0.5 ,],
                 mapping = aes(d13C, d15N), fill = FILL_overlap,  color = COL_overlap, level=0.95, geom="polygon", alpha=ALPHA) +
    # 4- proba >= 0.75
    stat_ellipse(Points_ellipses_over[Points_ellipses_over$Proba_overlap>=0.75 ,],
                 mapping = aes(d13C, d15N), fill = FILL_overlap,  color = COL_overlap, level=0.95, geom="polygon", alpha=ALPHA) +
    # 5- proba >= 0.95
    stat_ellipse(Points_ellipses_over[Points_ellipses_over$Proba_overlap>=0.95 ,],
                 mapping = aes(d13C, d15N), fill = FILL_overlap,  color = COL_overlap, level=0.95, geom="polygon", alpha=ALPHA) +
    # 6- proba >= 0.99
    stat_ellipse(Points_ellipses_over[Points_ellipses_over$Proba_overlap>=0.99 ,],
                 mapping = aes(d13C, d15N), fill = FILL_overlap,  color = COL_overlap, level=0.95, geom="polygon", alpha=ALPHA) +
    
    ### Observed individuals with pred values
    geom_errorbarh(data = PREY_iso_TEF, mapping=aes(y = d15N_mean, xmax = d13C_mean + d13C_SD, xmin = d13C_mean - d13C_SD, color = Prey_functionnal , height = 0.0), size=0.8) +
    geom_errorbar(data = PREY_iso_TEF, mapping=aes(x = d13C_mean, ymax = d15N_mean + d15N_SD, ymin = d15N_mean - d15N_SD, color = Prey_functionnal , width = 0.0), size=0.8) +
    geom_point(data = PREY_iso_TEF, mapping=aes(x = d13C_mean, y = d15N_mean, color = Prey_functionnal ), shape=15, size=7)+
    
    scale_color_manual(name = "Potential prey", values = c("goldenrod",  "darkorchid2", "green3", "royalblue4","deepskyblue"),
                       breaks = c( "Benthic flatfish",  "Benthic non flatfish", "Demersal fish", "Pelagic fish", "Pelagic squids"),
                       labels = c( "Benthic flatfish",  "Benthic non flatfish", "Demersal fish", "Pelagic fish", "Pelagic squids"), guide =F
    )+
    scale_fill_manual(name = "Potential prey", values = c( "goldenrod",  "darkorchid2", "green3", "royalblue4","deepskyblue"),
                      breaks = c( "Benthic flatfish",  "Benthic non flatfish", "Demersal fish", "Pelagic fish", "Pelagic squids"),
                      labels = c( "Benthic flatfish",  "Benthic non flatfish", "Demersal fish", "Pelagic fish", "Pelagic squids"), guide =F
    )+
    xlab(quote(delta^13*C)) + ylab(quote(delta^15*N))+
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(size = 14, color="black"),
          axis.title.x = element_text(size = 18)) + 
    scale_x_continuous(limits = c(-18.6, -11.9), breaks=Breaks_X) +
    scale_y_continuous(limits = c( Ranges_Y[1], 21), breaks=Breaks_Y) +
    geom_text(data = PREY_iso_TEF, mapping=aes(x = d13C_mean, y = d15N_mean, label= Code, fontface=2),
              size=3, color="white") +
    labs(caption = "Source of isotopic data for potential prey: Kopp et al. (2015)\nApplied TEF from prey to predator: d13C = +2.4 ± 1.3‰, d15N = +2.6 ± 1.4‰ (Lerner et al., 2018)") 
  
  Plot_05_potential_prey_overlap
  
  WID <- 7
  HEI <- 6
  ggsave(Plot_05_potential_prey_overlap,
         filename = paste(Direction, "Plot", "05_SI_overlap_potential_prey.png", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )
  
  ggsave(Plot_05_potential_prey_overlap,
         filename = paste(Direction, "Plot", "05_SI_overlap_potential_prey_PDF.pdf", sep = "/"), dpi = 300,
         width = WID, height = HEI
  )
  
  
save.image(paste(Direction, "Output", "hmodel_iso_Overlap1.RData", sep = "/"), safe = TRUE,
           compress = "xz")
#load(paste(Direction, "Output", "hmodel_iso_Overlap1.RData", sep = "/"))

###########################################################################################################  


###########################################################################################################
##################################### END #################################################################
###########################################################################################################