## ===========================
## QJRMS Paper - Bayesian verification of climate trends
## Bayesian trends + pooled posterior + optional verification
## ===========================

## ---- 0) Setup & utilities ----
options(stringsAsFactors = FALSE)

# Gaussian CRPS (Gneiting & Raftery 2007):
crps_norm <- function(y, mu, sig) {
  if (!is.finite(sig) || sig <= 0) return(NA_real_)
  z <- (y - mu)/sig
  # φ = dnorm, Φ = pnorm
  sig * (1/sqrt(pi) - 2*dnorm(z) - z*(2*pnorm(z) - 1))
}

# Conjugate Bayesian regression (Normal likelihood; Normal-Inverse-Gamma prior)
# y ~ N(X b, sigma^2 I)
# Prior: b|sigma^2 ~ N(b0, sigma^2 V0), sigma^2 ~ InvGam(a0, d0)
# Returns posterior draws of slope (b[2]) and predictive dist at new x (optional)
bayes_reg_conjugate <- function(x, y, x_new = NULL,
                                b0 = c(0,0), V0 = diag(2)*1e6,
                                a0 = 1e-2, d0 = 1e-2,
                                draws = 20000, seed = 42) {
  set.seed(seed)
  
  # keep finite
  idx <- is.finite(x) & is.finite(y)
  x <- x[idx]; y <- y[idx]
  n <- length(y)
  if (n < 24) return(NULL)
  
  # Posterior parameters
  X      <- cbind(1, x)
  V0i    <- solve(V0)
  XtX    <- crossprod(X)
  XtY    <- crossprod(X, y)
  Vn_inv <- V0i + XtX
  Vn     <- solve(Vn_inv)
  bn     <- Vn %*% (V0i %*% b0 + XtY)
  
  # sigma^2 ~ InvGamma(a_n, d_n)
  a_n <- a0 + n/2
  d_n <- d0 + 0.5*(sum(y^2) + t(b0)%*%V0i%*%b0 - t(bn)%*%Vn_inv%*%bn)
  
  # Draw sigma^2 and zero-mean Normal part ~ N(0, Vn)
  sigma2 <- 1 / rgamma(draws, shape = a_n, rate = d_n)
  
  Z <- matrix(rnorm(draws * 2), ncol = 2)   # draws × 2
  L <- chol(Vn)                              # 2 × 2, upper-triangular
  eps <- Z %*% L                             # draws × 2 ~ N(0, Vn)
  
  # Row-wise scale of eps by sqrt(sigma2_i), then add mean bn
  eps         <- sweep(eps, 1, sqrt(sigma2), `*`)           # (eps_i * sqrt(sigma2_i))
  beta_draws  <- sweep(eps, 2, as.numeric(bn), `+`)         # add bn to each column
  
  out <- list(beta_draws = beta_draws, sigma2 = sigma2, bn = bn, Vn = Vn)
  
  # Optional one-step-ahead predictive (plug-in using E[sigma^2])
  if (!is.null(x_new)) {
    Xnew   <- cbind(1, x_new)
    mu_new <- as.numeric(Xnew %*% bn)
    sig_new <- sqrt(mean(sigma2) * (1 + Xnew %*% Vn %*% t(Xnew)))
    out$pred_mu  <- mu_new
    out$pred_sig <- as.numeric(sig_new)
  }
  out
}

# OLS slope & p-value (two-sided)
ols_slope_p <- function(x, y) {
  idx <- is.finite(x) & is.finite(y)
  x <- x[idx]; y <- y[idx]
  if (length(y) < 24) return(c(NA_real_, NA_real_))
  fit <- lm(y ~ x)
  slope <- coef(summary(fit))[2,1]
  pval  <- coef(summary(fit))[2,4]
  c(slope, pval)
}

# Decimal years since a base date
year_frac <- function(d, base) as.numeric(d - base) / 365.25

############################################
## ---- 0) Stylised diagrams ----

# Figure 1 stylised overview of Bayesian approach 
library(ggplot2)
library(dplyr)
library(patchwork)   # for multi-panel layout

set.seed(123)

n <- 120  # e.g., 10 years of monthly data
time <- seq.Date(from = as.Date("2010-01-01"),
                 by   = "month",
                 length.out = n)

# Simulated slowly drifting mean + noise
true_trend <- 0.015
y <- 35 + true_trend * seq_len(n) + rnorm(n, sd = 0.6)

dat <- tibble(
  time = time,
  y    = y
)

# Panel (a): Data + Bayesian posterior summary (evolving predictive)
# Placeholder posterior summaries
dat <- dat %>%
  mutate(
    post_mean = 35 + 0.014 * seq_len(n),
    post_lwr  = post_mean - 1.1,
    post_upr  = post_mean + 1.1
  )

p_a <- ggplot(dat, aes(time, y)) +
  geom_point(size = 1.5, colour = "grey30") +
  geom_ribbon(aes(ymin = post_lwr, ymax = post_upr),
              fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = post_mean),
            colour = "steelblue", linewidth = 1) +
  labs(
    title = "(a) Observations and evolving predictive model",
    y = "Monthly temp. (°C)",
    x = NULL
  ) +
  theme_bw()

# Panel (b): Frequentist verification framing (fixed null)
ols <- lm(y ~ as.numeric(time), data = dat)

p_val <- summary(ols)$coefficients[2, 4]

dat <- dat %>%
  mutate(
    ols_fit = fitted(ols)
  )

p_b <- ggplot(dat, aes(time, y)) +
  geom_point(size = 1.5, colour = "grey30") +
  geom_line(aes(y = ols_fit),
            colour = "firebrick", linewidth = 1) +
  annotate(
    "text",
    x = dat$time[15],
    y = max(dat$y),
    label = paste0("OLS p-value = ", signif(p_val, 2)),
    hjust = 0,
    size = 3.5
  ) +
  labs(
    title = "(b) Frequentist inference against a fixed null",
    y = "Monthly temp. (°C)",
    x = NULL
  ) +
  theme_bw()

# Panel (c): Bayesian verification – evidence accumulation
# Simulated log-score differences (nonstationary vs baseline)
dat <- dat %>%
  mutate(
    log_score_diff = rnorm(n, mean = 0.03, sd = 0.08),
    cum_log_bf = cumsum(log_score_diff)
  )

p_c <- ggplot(dat, aes(time, cum_log_bf)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(colour = "darkblue", linewidth = 1) +
  labs(
    title = "(c) Posterior-predictive verification (evidence accumulation)",
    y = "Cum. log predictive score diff",
    x = "Time"
  ) +
  theme_bw()

final_figure <- p_a / p_b / p_c +
  plot_layout(heights = c(1, 1, 1))

final_figure # Figure 1

### Map of station locations
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

# Load Australia polygon
aus <- ne_countries(
  scale = "medium",
  country = "Australia",
  returnclass = "sf"
)

temp_stations <- tibble(
  station = c("Darwin", "Alice Springs", "Oodnadatta",
              "Sydney", "Melbourne", "Brisbane", "Perth"),
  lon = c(130.85, 133.88, 135.45, 151.21, 144.96, 153.03, 115.86),
  lat = c(-12.46, -23.70, -27.55, -33.87, -37.81, -27.47, -31.95),
  type = "Temperature"
)

sea_stations <- tibble(
  station = c("Sydney", "Darwin",
              "Newcastle", "Eden", "Broome",
              "Brisbane", "Port Pirie"),
  lon = c(151.23, 130.85, 151.78, 149.90, 122.24, 153.03, 137.94),
  lat = c(-33.85, -12.46, -32.93, -37.06, -17.96, -27.47, -33.18),
  type = "Sea level"
)

# Combine into a single data frame
stations <- bind_rows(temp_stations, sea_stations)

# Convert to sf object:
stations_sf <- st_as_sf(
  stations,
  coords = c("lon", "lat"),
  crs = 4326
)

# Plot map
p_map <- ggplot() +
  geom_sf(data = aus, fill = "grey95", colour = "grey40") +
  geom_sf(
    data = stations_sf,
    aes(colour = type, shape = type),
    size = 2.8
  ) +
  scale_colour_manual(
    values = c(
      "Temperature" = "firebrick",
      "Sea level"   = "steelblue"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Temperature" = 16,
      "Sea level"   = 17
    )
  ) +
  coord_sf(
    xlim = c(110, 155),
    ylim = c(-45, -10),
    expand = FALSE
  ) +
  labs(
    title = "Locations of temperature and sea-level stations",
    colour = "Station type",
    shape  = "Station type"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    legend.position = "inside",
    legend.position.inside = c(0.25,0.5),
    plot.title = element_text(size = 12, face = "bold")
  ) + xlab("Longitude") + ylab("Latitude")

p_map

# Figure 2
p_map +
  geom_sf_text(
    data = stations_sf %>%
      filter(station %in% c("Sydney", "Darwin",
                            "Newcastle", "Eden", "Broome",
                            "Brisbane", "Port Pirie", "Alice Springs", 
                            "Oodnadatta","Melbourne","Perth")),
    aes(label = station),
    size = 3,
    nudge_y = 0.8
  )

##########################################
## ---- 1) Load data ----
setwd("C:/Users/jwest.INTERNAL/OneDrive - Bureau of Meteorology/Documents/Data/Sea Levels")
f <- "MSL_Raw_Month.csv"
dat <- read.csv(f, check.names = FALSE)
names(dat)[1] <- "Date"
dat$Date <- as.Date(dat$Date, format = "%d/%m/%Y")  # day-first format in your file
# Coerce numeric
for (j in 2:ncol(dat)) dat[[j]] <- suppressWarnings(as.numeric(dat[[j]]))

stations <- setdiff(names(dat), "Date")

## ---- 2) Windows and station-wise trends (OLS + Bayesian) ----
analyze_window <- function(start_date, end_date, min_months = 60, seed = 123) {
  W <- subset(dat, Date >= as.Date(start_date) & Date <= as.Date(end_date))
  base <- as.Date(start_date)
  t_years <- year_frac(W$Date, base)
  res <- list(rows = list())
  for (st in stations) {
    y <- W[[st]]
    if (sum(is.finite(y)) < min_months) next
    ols <- ols_slope_p(t_years, y)
    fitB <- bayes_reg_conjugate(t_years, y, b0 = c(0,0), V0 = diag(2)*1e6,
                                a0 = 1e-2, d0 = 1e-2,
                                draws = 20000, seed = seed)
    if (is.null(fitB)) next
    slope_draws <- fitB$beta_draws[,2]
    mean_slope <- mean(slope_draws)
    ci <- quantile(slope_draws, c(.025,.975))
    ppos <- mean(slope_draws > 0)
    res$rows[[length(res$rows)+1]] <- data.frame(
      Station = st,
      N = sum(is.finite(y)),
      OLS_slope = ols[1],
      OLS_p = ols[2],
      Bayes_mean = mean_slope,
      Bayes_CI_low = ci[1],
      Bayes_CI_high = ci[2],
      Bayes_sd = sd(slope_draws),
      Bayes_Ppos = ppos
    )
  }
  do.call(rbind, res$rows)
}

# Run both windows (as in the manuscript)
res_2012 <- analyze_window("2012-01-01", "2022-12-31", min_months = 60, seed = 123)
res_2015 <- analyze_window("2015-01-01", "2024-12-31", min_months = 48, seed = 124)

write.csv(res_2012, "table_S1_sea_level_trends_2012_2022.csv", row.names = FALSE)
write.csv(res_2015, "table_S1b_sea_level_trends_2015_2024.csv", row.names = FALSE)

## ---- 3) Random-effects pooled posterior (meta-analytic) + shrinkage ----
pooled_random_effects <- function(df) {
  # df needs Bayes_mean, Bayes_sd
  d <- subset(df, is.finite(Bayes_mean) & is.finite(Bayes_sd) & Bayes_sd > 0)
  if (nrow(d) < 2) return(NULL)
  wi <- 1/(d$Bayes_sd^2)
  Yi <- d$Bayes_mean
  wbar <- sum(wi)
  mu_FE <- sum(wi*Yi) / wbar
  Q <- sum(wi*(Yi - mu_FE)^2)
  C <- wbar - sum(wi^2)/wbar
  tau2 <- max(0, (Q - (nrow(d)-1)) / C)
  wi_star <- 1/(d$Bayes_sd^2 + tau2)
  mu_RE <- sum(wi_star*Yi)/sum(wi_star)
  sd_RE <- sqrt(1/sum(wi_star))
  list(mu = mu_RE, sd = sd_RE, tau2 = tau2,
       shrunk = if (tau2 > 0) {
         wt <- 1/tau2
         ( (1/d$Bayes_sd^2)*d$Bayes_mean + wt*mu_RE ) / ( (1/d$Bayes_sd^2) + wt )
       } else rep(mu_RE, nrow(d)),
       df = d)
}

pool <- pooled_random_effects(res_2012)
if (!is.null(pool)) {
  cat(sprintf("Regional mean slope (2012–2022): %.6f, sd: %.6f, tau^2: %.6f, P>0 ≈ %.3f\n",
              pool$mu, pool$sd, pool$tau2,
              1 - pnorm(0, mean = pool$mu, sd = pool$sd)))
}

## ---- 4) Figures (base plotting; save to PNG) ----
#png("fig_S1_station_examples.png", width = 1000, height = 700)
par(mfrow = c(2,2), mar = c(4,4,3,1))
example_stations <- intersect(c("Sydney (Fort Denison)","Darwin","Newcastle","Eden"), stations)
W <- subset(dat, Date >= as.Date("2012-01-01") & Date <= as.Date("2022-12-31"))
t_years <- year_frac(W$Date, as.Date("2012-01-01"))
for (st in example_stations) {
  y <- W[[st]]
  plot(W$Date, y, type = "l", col = "#1f77b4", xlab = "Date", ylab = "Sea-level anomaly",
       main = st)
  ols <- ols_slope_p(t_years, y)
  if (is.finite(ols[1])) {
    yhat <- mean(y, na.rm=TRUE) + ols[1]*(t_years - mean(t_years))
    lines(W$Date, yhat, col = "#d62728", lwd = 2)
  }
  fitB <- bayes_reg_conjugate(t_years, y, draws = 30000, seed = 2024)
  if (!is.null(fitB)) {
    slope_draws <- fitB$beta_draws[,2]
    ppos <- mean(slope_draws>0)
    # inset density
    usr <- par("usr")
    xleft <- usr[1] + 0.67*(usr[2]-usr[1]); xright <- usr[2] - 0.02*(usr[2]-usr[1])
    ybottom <- usr[3] + 0.08*(usr[4]-usr[3]); ytop <- usr[3] + 0.45*(usr[4]-usr[3])
    par(xpd=NA)
    rect(xleft, ybottom, xright, ytop, col = "white", border = "gray70")
    par(xpd=NA)
    dens <- density(slope_draws)
    xs <- seq(xleft, xright, length.out = length(dens$x))
    ys <- ybottom + (dens$y/max(dens$y))*(ytop - ybottom)
    lines(xs, ys, col = "black")
    # zero line in inset (approximate mapping of x=0 to inset x-domain)
    zpos <- approx(dens$x, xs, xout = 0)$y
    segments(zpos, ybottom, zpos, ytop, col = "gray30", lty = 2)
    text(xleft + 0.03*(xright-xleft), ytop - 0.06*(ytop-ybottom),
         labels = sprintf("P(slope>0)=%.2f", ppos), adj = c(0,1), cex = 0.9)
  }
}
#dev.off()

# Fig S2: regional pooled posterior (normal approx)
if (!is.null(pool)) {
  xs <- seq(pool$mu - 5*pool$sd, pool$mu + 5*pool$sd, length.out = 400)
  ys <- dnorm(xs, mean = pool$mu, sd = pool$sd)
  #png("fig_S2_pooled_posterior.png", width = 800, height = 500)
  plot(xs, ys, type = "l", col = "purple", lwd = 2,
       xlab = "Regional slope (per year), 2012–2022",
       ylab = "Posterior density (approx)",
       main = "Regional pooled posterior of trend (random effects)")
  abline(v = pool$mu, col = "purple", lwd = 1.5)
  xx <- xs[xs >= pool$mu - 1.96*pool$sd & xs <= pool$mu + 1.96*pool$sd]
  yy <- ys[xs >= pool$mu - 1.96*pool$sd & xs <= pool$mu + 1.96*pool$sd]
  polygon(c(xx, rev(xx)), c(rep(0,length(xx)), rev(yy)), border = NA, col = rgb(0.5,0,0.5,0.2))
  abline(v = 0, col = "gray40", lty = 2)
  txt <- sprintf("mean=%.4f, 95%% CI=[%.4f, %.4f]\nP(slope>0)≈%.3f",
                 pool$mu, pool$mu - 1.96*pool$sd, pool$mu + 1.96*pool$sd,
                 1 - pnorm(0, mean = pool$mu, sd = pool$sd))
  mtext(txt, side = 3, adj = 0.02, cex = 0.9)
  #dev.off()
}

#Improve Figure S2:
# Fig S2: Regional pooled posterior (Normal approximation)
  if (!is.null(pool)) {
    
    xs <- seq(pool$mu - 5 * pool$sd,
              pool$mu + 5 * pool$sd,
              length.out = 500)
    
    ys <- dnorm(xs, mean = pool$mu, sd = pool$sd)
    
    plot(
      xs, ys,
      type = "l",
      lwd = 2.5,
      col = "#6A3D9A",
      xlab = "Regional sea-level trend (per year), 2012–2022",
      ylab = "Posterior density",
      main = "Regionally pooled posterior distribution of sea-level trend"
    )
    
    # 95% credible interval
    ci_lo <- pool$mu - 1.96 * pool$sd
    ci_hi <- pool$mu + 1.96 * pool$sd
    
    xx <- xs[xs >= ci_lo & xs <= ci_hi]
    yy <- ys[xs >= ci_lo & xs <= ci_hi]
    
    polygon(
      c(xx, rev(xx)),
      c(rep(0, length(xx)), rev(yy)),
      col = adjustcolor("#6A3D9A", alpha.f = 0.25),
      border = NA
    )
    
    # Zero-slope reference
    abline(v = 0, col = "grey40", lty = 2, lwd = 1.5)
    
    # Posterior mean
    abline(v = pool$mu, col = "#6A3D9A", lwd = 3)
    
    # Annotation (placed in margin, not data region)
    p_pos <- 1 - pnorm(0, mean = pool$mu, sd = pool$sd)
    
    legend(
      "topright",
      inset = 0.02,
      bty = "n",
      legend = c(
        sprintf("Mean = %.4f", pool$mu),
        sprintf("95%% CI = [%.4f, %.4f]", ci_lo, ci_hi),
        sprintf("P(slope > 0) = %.3f", p_pos)
      ),
      text.col = "#333333"
    )
  }

# Fig S3: window sensitivity
brks_2012 <- seq(
  floor(min(res_2012$Bayes_Ppos, na.rm = TRUE) / 0.02) * 0.02,
  ceiling(max(res_2012$Bayes_Ppos, na.rm = TRUE) / 0.02) * 0.02,
  by = 0.02
)

brks_2015 <- seq(
  floor(min(res_2015$Bayes_Ppos, na.rm = TRUE) / 0.02) * 0.02,
  ceiling(max(res_2015$Bayes_Ppos, na.rm = TRUE) / 0.02) * 0.02,
  by = 0.02
)

par(mfrow = c(1,2), mar = c(4,4,3,1))

hist(res_2012$Bayes_Ppos,
     breaks = brks_2012,
     col = "steelblue",
     border = "white",
     main = "P(slope > 0) across stations\n2012–2022",
     xlab = "Posterior probability")

hist(res_2015$Bayes_Ppos,
     breaks = brks_2015,
     col = "seagreen",
     border = "white",
     main = "P(slope > 0) across stations\n2015–2024",
     xlab = "Posterior probability")

par(mfrow = c(1,2), mar = c(4,4,3,1))

hist(res_2012$Bayes_Ppos,
     breaks = brks_2012,
     col = "steelblue",
     border = "white",
     xlim = c(0.4, 1.0),
     main = "Posterior evidence across stations\n2012–2022",
     xlab = "P(slope > 0)")

abline(v = 0.5, lty = 2, lwd = 2, col = "grey30")
abline(v = 0.9, lty = 3, col = "grey40")

hist(res_2015$Bayes_Ppos,
     breaks = brks_2015,
     col = "seagreen",
     border = "white",
     xlim = c(0.4, 1.0),
     main = "Posterior evidence across stations\n2015–2024",
     xlab = "P(slope > 0)")

abline(v = 0.5, lty = 2, lwd = 2, col = "grey30")
abline(v = 0.9, lty = 3, col = "grey40")


# #png("fig_S3_window_sensitivity.png", width = 1000, height = 450)
# par(mfrow = c(1,2), mar = c(4,4,3,1))
# hist(res_2012$Bayes_Ppos, breaks = seq(0.5,1.0,by=0.02), col = "steelblue",
#      main = "P(slope>0) across stations\n2012–2022", xlab = "Posterior probability")
# hist(res_2015$Bayes_Ppos, breaks = seq(0.5,1.0,by=0.02), col = "seagreen",
#      main = "P(slope>0) across stations\n2015–2024", xlab = "Posterior probability")
# #dev.off()

# Fig S4: hierarchical shrinkage scatter
if (!is.null(pool)) {
  d <- pool$df; d$beta_shrunk <- pool$shrunk
  rng <- range(c(d$Bayes_mean, d$beta_shrunk), na.rm = TRUE)
#  png("fig_S4_hierarchical_shrinkage.png", width = 800, height = 600)
  plot(d$Bayes_mean, d$beta_shrunk, pch = 21, bg = "orange", col = "black",
       xlab = "Station-wise posterior mean slope (per yr)",
       ylab = "Hierarchical shrunken slope (per yr)",
       main = "Hierarchical robustness: shrinkage of station trends",
       xlim = rng, ylim = rng)
  abline(0,1,lty=2,col="gray40")
  # annotate a few
  tag <- intersect(c("Sydney (Fort Denison)","Darwin","Newcastle","Eden","Brisbane","Broome"),
                   d$Station)
  for (st in tag) {
    i <- which(d$Station == st)[1]
    text(d$Bayes_mean[i], d$beta_shrunk[i], labels = sub("\\(.*","",st), pos = 4, cex=0.8)
  }
#  dev.off()
}

# Re-run the analysis to produce a better figure
d <- pool$df
d$beta_shrunk <- pool$shrunk
d$shrinkage <- d$beta_shrunk - d$Bayes_mean

rng <- range(c(d$Bayes_mean, d$beta_shrunk), na.rm = TRUE)

plot(
  d$Bayes_mean,
  d$beta_shrunk,
  pch = 21,
  bg  = "steelblue",
  col = "grey20",
  cex = 1.2,
  xlim = rng,
  ylim = rng,
  xlab = "Station-wise posterior mean slope (per year)",
  ylab = "Hierarchically pooled slope (per year)",
  main = "Effect of hierarchical pooling on station-level sea-level trends"
)

# 1:1 reference
abline(0, 1, lwd = 2, lty = 2, col = "grey50")

# arrows showing shrinkage
segments(
  x0 = d$Bayes_mean,
  y0 = d$Bayes_mean,
  x1 = d$Bayes_mean,
  y1 = d$beta_shrunk,
  col = adjustcolor("darkorange", alpha.f = 0.6),
  lwd = 1
)

# highlight key stations only
key <- intersect(
  c("Sydney (Fort Denison)", "Darwin", "Newcastle", "Eden"),
  d$Station
)

for (st in key) {
  i <- which(d$Station == st)[1]
  points(d$Bayes_mean[i], d$beta_shrunk[i],
         pch = 21, bg = "orange", col = "black", cex = 1.5)
  # text(d$Bayes_mean[i], d$beta_shrunk[i],
  #      labels = sub("\\(.*", "", st),
  #      pos = 4, cex = 0.8, srt = 30)
}

legend(
  "topleft",
  bty = "n",
  legend = c(
    "Elastic shrinkage toward regional mean",
    "1:1 line (no shrinkage)"
  ),
  lty = c(1, 2),
  col = c("darkorange", "grey50")
)


## ---- 5) OPTIONAL: Posterior-predictive verification with proper scores ----
## Rolling-origin one-step-ahead for a single station (extend to all stations):
verify_station <- function(st, start_date="2012-01-01", end_date="2022-12-31",
                           warm_up_months = 24) {
  W <- subset(dat, Date >= as.Date(start_date) & Date <= as.Date(end_date))
  base <- as.Date(start_date)
  t_years <- year_frac(W$Date, base)
  y <- W[[st]]
  n <- length(y)
  # storage
  out <- data.frame(Date=W$Date, y= y,
                    mu = NA_real_, sig = NA_real_,
                    logS = NA_real_, CRPS = NA_real_,
                    mu_base = NA_real_, sig_base = NA_real_,
                    logS_base = NA_real_, CRPS_base = NA_real_,
                    PIT = NA_real_)
  for (k in seq_len(n)) {
    if (k <= warm_up_months || sum(is.finite(y[1:(k-1)])) < 24) next
    x_tr <- t_years[1:(k-1)]; y_tr <- y[1:(k-1)]
    # Model A: trend (x)
    fitA <- bayes_reg_conjugate(x_tr, y_tr, x_new = t_years[k], draws = 20000, seed = 99+k)
    # Model B (baseline): Bayesian climatology (intercept only)
    fitB <- bayes_reg_conjugate(rep(0, length(x_tr)), y_tr, x_new = 0, draws = 20000, seed = 77+k)
    if (!is.null(fitA) && !is.null(fitB) && is.finite(y[k])) {
      muA <- fitA$pred_mu; sigA <- fitA$pred_sig
      muB <- fitB$pred_mu; sigB <- fitB$pred_sig
      out$mu[k] <- muA; out$sig[k] <- sigA
      out$mu_base[k] <- muB; out$sig_base[k] <- sigB
      out$logS[k] <- dnorm(y[k], mean = muA, sd = sigA, log = TRUE)         # Log Score
      out$CRPS[k] <- crps_norm(y[k], muA, sigA)                              # CRPS
      out$logS_base[k] <- dnorm(y[k], mean = muB, sd = sigB, log = TRUE)
      out$CRPS_base[k] <- crps_norm(y[k], muB, sigB)
      out$PIT[k] <- pnorm((y[k]-muA)/sigA)                                  # PIT for Normal
    }
  }
  out
}

# Example verification for Sydney; loop for others if desired
ver_syd <- verify_station("Sydney (Fort Denison)")
# Aggregate verification metrics
with(subset(ver_syd, is.finite(CRPS) & is.finite(CRPS_base)), {
  cat(sprintf("Sydney CRPS (model): %.4f | baseline: %.4f | skill = %.2f%%\n",
              mean(CRPS), mean(CRPS_base),
              100*(1 - mean(CRPS)/mean(CRPS_base))))
  cat(sprintf("Sydney mean logS gain (model - baseline): %.4f\n",
              mean(logS - logS_base, na.rm=TRUE)))
})
# Can also inspect histogram(ver_syd$PIT, breaks=20) for calibration.

#####################################################################
#### TEMPERATURE CASE STUDY
#####################################################################

# --- CONFIG -------------------------------------------------------------------
year_min <- 1990
year_max <- 2025

core_stations <- c(
  # Capitals
  "SYDNEY AIRPORT", "MELBOURNE AIRPORT", "BRISBANE", "ADELAIDE AIRPORT",
  "PERTH AIRPORT", "HOBART (ELLERSLIE RD)", "DARWIN AIRPORT", "CANBERRA AIRPORT",
  # Climate-sensitive exemplars (adjust to exact names in tempdat.csv)
  "MARBLE BAR", "PORT HEDLAND", "OODNADATTA AIRPORT", "ALICE SPRINGS AIRPORT",
  "CAIRNS", "BROOME", "MILDURA AIRPORT", "BIRDSVILLE AIRPORT"
)

# --- LOAD DAILY TEMPS ---------------------------------------------------------
library(data.table); library(dplyr); library(lubridate); library(readr); library(janitor)
library(rstanarm); library(tidyverse); library(scoringRules); library(posterior); library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("C:/Users/jwest.INTERNAL/OneDrive - Bureau of Meteorology/Documents/Data/Climate/Attribution/Bayesian Verification")

dt <- fread("tempdat.csv", select = c("Date","Stname","TMax","Tmin"))    # daily format & columns
# (file metadata and structure as per user-provided tempdat.csv)              # [1](https://bom365-my.sharepoint.com/personal/jason_west_bom_gov_au/_layouts/15/Doc.aspx?sourcedoc=%7BAAC08CC3-790A-4F0E-9A8A-662B514B8C37%7D&file=tempdat.csv&action=default&mobileredirect=true)
dt[, Date := as.Date(Date)]
dt[, Year := year(Date)]
dt[, Month := month(Date)]
dt <- dt[Year >= year_min & Year <= year_max]
dt <- dt[Stname %in% core_stations]

# --- MONTHLY BLOCK MAXIMA (Txx), require >=25 days ---------------------------
tx_monthly <- dt[!is.na(TMax),
                 .(n_days = sum(!is.na(TMax)),
                   Txx = ifelse(sum(!is.na(TMax)) >= 25, max(TMax, na.rm = TRUE), NA_real_)),
                 by = .(station_id = Stname, Year, Month)] %>%
  as_tibble() %>% filter(!is.na(Txx)) %>%
  mutate(date = ymd(paste(Year, Month, "15")))

# --- LOAD AND PREPARE INDICES (incl. SAM wide->long) -------------------------
enso <- read_csv("enso_anom_mthly.csv", show_col_types = FALSE) %>%
  transmute(Year, Month, ENSO = stANOM)                                   # [1](https://bom365-my.sharepoint.com/personal/jason_west_bom_gov_au/_layouts/15/Doc.aspx?sourcedoc=%7BAAC08CC3-790A-4F0E-9A8A-662B514B8C37%7D&file=tempdat.csv&action=default&mobileredirect=true)
iod  <- read_csv("iod.csv", show_col_types = FALSE) %>%
  transmute(Year, Month, IOD = iod)                                       # [2](https://bom365-my.sharepoint.com/personal/jason_west_bom_gov_au/_layouts/15/Doc.aspx?sourcedoc=%7BC9D8713A-436F-4E14-B666-F82F8D7EDE46%7D&file=Bayesian%20Verification%20Changing%20Climate%20v2.docx&action=default&mobileredirect=true)
sst  <- read_csv("nAusSST.csv", show_col_types = FALSE) %>%
  transmute(Year, Month, nAusSST)                                         # [3](https://bom365-my.sharepoint.com/personal/jason_west_bom_gov_au/_layouts/15/Doc.aspx?sourcedoc=%7B3CD9E7E6-04AA-45E4-8038-1E5A1FA13B3C%7D&file=SAM_wd.csv&action=default&mobileredirect=true)

sam <- read_csv("SAM.csv",show_col_types = FALSE) %>%
  pivot_longer(
    cols = -Year,
    names_to = "Month",
    values_to = "Value"
  ) %>%
  mutate(
    MonthNum = match(
      Month,
      c("JAN","FEB","MAR","APR","MAY","JUN",
        "JUL","AUG","SEP","OCT","NOV","DEC")
    )
  ) %>%
  arrange(Year, MonthNum) %>% 
  select(Year, MonthNum, Value)
colnames(sam)<-c("Year", "Month", "SAM")

indices <- enso %>% full_join(iod, by=c("Year","Month")) %>%
  full_join(sam, by=c("Year","Month")) %>%
  full_join(sst, by=c("Year","Month")) %>%
  mutate(date = ymd(paste(Year, Month, "15"))) %>%
  arrange(date) %>%
  mutate(ENSO_z = as.numeric(scale(ENSO)),
         IOD_z  = as.numeric(scale(IOD)),
         SAM_z  = as.numeric(scale(SAM)),
         SST_z  = as.numeric(scale(nAusSST)))

dat <- tx_monthly %>%
  left_join(indices, by = c("Year","Month","date")) %>%
  mutate(m = month(date),
         s1 = sin(2*pi*m/12), c1 = cos(2*pi*m/12),
         s2 = sin(4*pi*m/12), c2 = cos(4*pi*m/12),
         t  = as.numeric((date - min(date)) / lubridate::ddays(30))) %>%
  arrange(station_id, date)

dat <- dat %>%
  mutate(
    station_id = factor(station_id),
    date = as.Date(date)
  )

# Next steps: fit baseline & regime-aware BAYES models and run rolling-origin
# scoring (CRPS/logS) exactly as in Section 3 workflow 

# write.csv(dat,"C:/Users/jwest.INTERNAL/OneDrive - Bureau of Meteorology/Documents/Data/Climate/Attribution/Bayesian Verification/dat.csv")
# Factorise

# --- 2.1 Baseline (trend + seasonality, no regimes) -------------------------------------
# Bayesian climatology‑plus‑trend reference

form_base <- Txx ~
  t + s1 + c1 + s2 + c2 + 
  (1 | station_id) 

# --- 2.2 Regime‑aware model (hierarchical teleconnections) ------------------------------
# Physically motivated alternative.

form_regime <- Txx ~
  t + s1 + c1 + s2 + c2 +
  ENSO_z + IOD_z + SAM_z + SST_z +
  (1 + ENSO_z + IOD_z + SAM_z + SST_z | station_id)

# --- 3. Prequential (rolling‑origin) verification -------------------------------------
# Non‑negotiable for this paper
issue_dates <- sort(unique(dat$date))
issue_dates <- issue_dates[seq(1, length(issue_dates), by = 2)] # sample every second month
min_train   <- 120   # e.g. 10 years

score_store <- list()

for (i in seq(min_train + 1, length(issue_dates))) {
  
  train_date  <- issue_dates[i - 1]
  verify_date <- issue_dates[i]
  
  train_dat <- dat %>% filter(date <= train_date)
  test_dat  <- dat %>% filter(date == verify_date)
  
  fit_base <- stan_glmer(
    form_base,
    data   = train_dat,
    family = gaussian(),  # ✅ FIXED
    prior_intercept = normal(0, 5),
    prior = normal(0, 1),
    prior_aux = student_t(3, 0, 5),  # ✅ robust variance
    chains = 2, iter = 1000, refresh = 0
  )
  
  fit_regime <- stan_glmer(
    form_regime,
    data   = train_dat,
    family = gaussian(),  # ✅ FIXED
    prior_intercept = normal(0, 5),
    prior = normal(0, 1),
    prior_aux = student_t(3, 0, 5),
    chains = 2, iter = 1000, refresh = 0
  )
  
  pp_base   <- posterior_predict(fit_base, newdata = test_dat)
  pp_regime <- posterior_predict(fit_regime, newdata = test_dat)
  
  for (j in seq_len(nrow(test_dat))) {
    
    y <- test_dat$Txx[j]
    
    score_store[[length(score_store) + 1]] <- tibble(
      date = verify_date,
      station_id = test_dat$station_id[j],
      model = "Baseline",
      crps = crps_sample(y, pp_base[, j]),
      logS = logs_sample(y, pp_base[, j])
    )
    
    score_store[[length(score_store) + 1]] <- tibble(
      date = verify_date,
      station_id = test_dat$station_id[j],
      model = "Regime",
      crps = crps_sample(y, pp_regime[, j]),
      logS = logs_sample(y, pp_regime[, j])
    )
  }
}

# posterior predictive distributions
# proper scores
# station‑wise and time‑wise results
scores <- bind_rows(score_store)
# write.csv(scores,"C:/Users/jwest.INTERNAL/OneDrive - Bureau of Meteorology/Documents/Research/Papers/Bayesian Beliefs/Climate Attribution QJMS/score.csv")

# --- 4. Bayes‑factor evidence on the predictive scale ----------------------------------
bf_time <- scores %>%
  group_by(date, model) %>%
  summarise(mean_logS = mean(logS), .groups = "drop") %>%
  pivot_wider(names_from = model, values_from = mean_logS) %>%
  mutate(
    dlogS = Baseline - Regime,
    cum_dlogS = cumsum(dlogS),
    BayesFactor = exp(cum_dlogS)
  )


# --- 5. Calibration (PIT) from rstanarm ----------------------------------
# Example using final regime model
pp <- posterior_predict(fit_regime)
y  <- train_dat$Txx

pit <- apply(pp, 2, function(x) mean(x <= y))

pit_df <- tibble(pit = pit)

################################################
## Figures
library(ggplot2)
library(dplyr)

# --- Figure 1 — Station‑wise predictive skill (CRPS skill)
# To show that the Bayesian regime‑aware model improves probabilistic skill 
# across stations — not just in aggregate.
crps_skill <- scores %>%
  group_by(station_id, model) %>%
  summarise(mean_crps = mean(crps), .groups = "drop") %>%
  pivot_wider(names_from = model, values_from = mean_crps) %>%
  mutate(CRPS_skill = 1 - Regime / Baseline)

fig1 <- ggplot(crps_skill,
               aes(x = reorder(station_id, CRPS_skill),
                   y = CRPS_skill)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2) +
  coord_flip() +
  labs(
    x = NULL,
    y = "CRPS skill (Regime vs Baseline)",
    title = "Improvement in probabilistic skill using regime-aware Bayesian model"
  ) +
  theme_bw()

fig1

# --- Figure 2 — Distribution of skill across stations
# To show that improvements are systematic, not driven by a few outliers.
fig2 <- ggplot(crps_skill, aes(x = CRPS_skill)) +
  geom_histogram(bins = 20, fill = "grey70", colour = "white") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "CRPS skill (Regime vs Baseline)",
    y = "Number of stations",
    title = "Distribution of probabilistic skill improvements"
  ) +
  theme_bw()

fig2

# --- Figure 3 — Predictive evidence through time (Bayes factors)
# To demonstrate that evidence accumulates gradually and asymmetrically, rather 
# than flipping at an arbitrary p‑value threshold.
fig3 <- ggplot(bf_time, aes(x = date, y = cum_dlogS)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(size = 0.8) +
  labs(
    x = NULL,
    y = expression(sum(Delta*log*S)),
    title = "Cumulative predictive evidence favouring the regime-aware model",
    subtitle = "Positive values indicate increasing Bayes-factor support"
  ) +
  theme_bw()

fig3

# And add a secondary axis for the Bayes factor itself:
fig3 + scale_y_continuous(
  sec.axis = sec_axis(~ exp(.), name = "Bayes factor")
)

# --- Figure 4 — Calibration of posterior predictives (PIT)
# To show that increased resolution does not come at the cost of miscalibration.
pit_df <- tibble(pit = pit)

fig4 <- ggplot(pit_df, aes(x = pit)) +
  geom_histogram(binwidth = 0.05, boundary = 0,
                 fill = "grey70", colour = "white") +
  labs(
    x = "Probability Integral Transform",
    y = "Count",
    title = "Calibration of posterior predictive distributions"
  ) +
  theme_bw()

fig4

# PIT by region, a faceted version is better
ggplot(pit_df, aes(x = pit)) +
  geom_histogram(binwidth = 0.05, fill = "grey70") +
  facet_wrap(~ station_group)


# --- Figure 5 — Information–Noise (optional but powerful)

# To construct an Info-Noise graph there are several steps needed:
library(rstanarm)
library(tibble)

# Gives a Normal approximation to the full posterior predictive
# Suitable for CRPS, PIT, and variance recalibration
# Fast and stable
extract_mu_sd <- function(fit, newdata) {
  
  # Posterior samples of the linear predictor (on data scale)
  epred <- posterior_epred(
    fit,
    newdata = newdata
  )  # matrix: iterations × observations
  
  # Posterior samples of residual SD (sigma)
  sigma_draws <- as.matrix(fit, pars = "sigma")
  
  mu <- colMeans(epred)
  
  sd <- sqrt(
    apply(epred, 2, var) +
      mean(sigma_draws^2)
  )
  
  tibble(
    mu = mu,
    sd = sd
  )
}

# After fitting the model
pred_stats <- extract_mu_sd(fit_regime, newdata = test_dat)

for (j in seq_len(nrow(test_dat))) {
  
  score_store[[length(score_store) + 1]] <- tibble(
    date       = verify_date,
    station_id = test_dat$station_id[j],
    model      = "Regime",
    
    y  = test_dat$Txx[j],
    mu = pred_stats$mu[j],
    sd = pred_stats$sd[j]
  )
}

library(scoringRules)
library(tibble)

scores_norm <- tibble(
  date       = as.Date(character()),
  station_id = character(),
  model      = character(),
  y          = numeric(),
  mu         = numeric(),
  sd         = numeric()
)

scores_norm <- scores_norm %>%
  mutate(
    crps = crps_norm(y, mean = mu, sd = sd),
    pit  = pnorm(y, mean = mu, sd = sd)
  )

str(scores_norm)

# refit one representative model
fit_diag <- stan_glmer(
  form_regime,
  data   = dat %>% filter(date <= as.Date("2015-12-31")),
  family = gaussian(),
  prior_intercept = normal(0, 5),
  prior = normal(0, 1),
  prior_aux = student_t(3, 0, 5),
  chains = 2,
  iter = 1000,
  refresh = 0
)

pred_diag <- extract_mu_sd(
  fit_diag,
  newdata = dat %>% filter(date > as.Date("2015-12-31"))
)

diag_data <- dat %>% 
  filter(date > as.Date("2015-12-31"))

pred_diag <- pred_diag %>%
  mutate(
    y          = diag_data$Txx,
    station_id = diag_data$station_id,
    date       = diag_data$date
  )

str(pred_diag) # Valid forecast–verification table

# Compute CRPS for the regime model (Normal approximation)
library(scoringRules)

# raw predictive CRPS, before recalibration.
pred_diag <- pred_diag %>%
  mutate(
    crps_raw = crps_norm(y, mean = mu, sd = sd)
  )

# Estimate and apply variance recalibration
# Step 1: estimate optimal scale factor
estimate_alpha <- function(y, mu, sd) {
  objective <- function(alpha) {
    mean(crps_norm(y, mean = mu, sd = alpha * sd))
  }
  optimize(objective, interval = c(0.5, 2))$minimum
}

alpha_hat <- estimate_alpha(
  y  = pred_diag$y,
  mu = pred_diag$mu,
  sd = pred_diag$sd
)

alpha_hat

# Step 2: compute recalibrated CRPS
pred_diag <- pred_diag %>%
  mutate(
    crps_recal = crps_norm(y, mean = mu, sd = alpha_hat * sd)
  )

# Bring in the baseline CRPS
baseline_diag <- scores %>%
  filter(model == "Baseline") %>%
  select(date, station_id, crps_baseline = crps)

# Join to pred_diag:
pred_diag <- pred_diag %>%
  left_join(baseline_diag, by = c("date", "station_id"))

pred_diag<-na.omit(pred_diag)

# Compute Information and Noise
pred_diag <- pred_diag %>%
  mutate(
    Noise       = crps_raw   - crps_recal,
    Information = crps_recal - crps_baseline
  )
# Negative CRPS differences = improvement
# Therefore "more negative" = better

# Aggregate for plotting - stationwise
info_noise_station <- pred_diag %>%
  group_by(station_id) %>%
  summarise(
    Information = mean(Information),
    Noise       = mean(Noise),
    .groups = "drop"
  )

# Information-Noise Plot
library(ggplot2)

fig5 <- ggplot(
  info_noise_station,
  aes(x = reorder(station_id, Information), y = Information)
) +
  geom_col(fill = "steelblue") +
  geom_errorbar(
    aes(
      ymin = Information - Noise,
      ymax = Information
    ),
    width = 0,
    colour = "orange",
    size = 1
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "CRPS difference",
    title = "Predictive information dominates removable calibration noise",
    subtitle = "Blue: Information (after recalibration); Orange: Noise (removed by variance scaling)"
  ) +
  theme_bw()

fig5

plot_df <- info_noise_station %>%
  mutate(
    Total = Information + Noise
  ) %>%
  select(station_id, Information, Noise) %>%
  pivot_longer(
    cols = c(Information, Noise),
    names_to = "Component",
    values_to = "Value"
  )

fig5 <- ggplot(
  plot_df,
  aes(
    x = reorder(station_id, Value),
    y = Value,
    fill = Component
  )
) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      Information = "steelblue",
      Noise       = "darkorange"
    ),
    labels = c(
      Information = "Information (retained after recalibration)",
      Noise       = "Noise (removed by recalibration)"
    )
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "CRPS difference (negative = improvement)",
    title = "Predictive skill decomposed into information and noise",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

fig5

# Why this method is correct
# posterior_epred() captures parameter uncertainty
# sigma captures irreducible observation noise
# Summing variances is correct for Normal mixtures
# CRPS and PIT have closed‑form expressions
# No look‑ahead or leakage is introduced
