# 02_baseline_acf.R -- baseline trend fit + residual ACF for HAC lag selection
# Inputs : data/clean/daily_panel.csv
# Outputs: output/figs/fig03_acf.pdf
#          output/figs/fig04_outcome_dist.pdf
#          output/tables/tab_pre_baseline.md
#          data/clean/panel_transformed.csv  (level + log-transformed taxa)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(splines)
})

ROOT <- "/Users/haofeixu/Desktop/Causal Inference/final_project"
CLN  <- file.path(ROOT, "data/clean")
FIG  <- file.path(ROOT, "output/figs")
TAB  <- file.path(ROOT, "output/tables")
dir.create(TAB, showWarnings=FALSE, recursive=TRUE)

UPW <- 215:221; HUR <- 247:249
panel <- read.csv(file.path(CLN, "daily_panel.csv"))

# ---- 1. log-transform taxa ----
TAXA <- c("vibrio", "flavobact", "cyano", "sar11", "rhodo")
for (tx in TAXA) {
  x <- panel[[tx]]
  c_eps <- min(x[x > 0], na.rm = TRUE) / 2   # half the min positive value
  panel[[paste0("log_", tx)]] <- log(x + c_eps)
  cat(sprintf("  %-10s: min=%.2e  pseudocount=%.2e  log_range=[%.2f, %.2f]\n",
              tx, min(x, na.rm=TRUE), c_eps,
              min(log(x + c_eps), na.rm=TRUE), max(log(x + c_eps), na.rm=TRUE)))
}

write.csv(panel, file.path(CLN, "panel_transformed.csv"), row.names = FALSE)

# ---- 2. outcomes to diagnose ----
OUT <- c(Shannon="Shannon", Richness="Richness", Simpson="Simpson",
         log_vibrio="log_vibrio", log_flavobact="log_flavobact",
         log_cyano="log_cyano",  log_sar11="log_sar11", log_rhodo="log_rhodo")

# ---- 3. fit trend-only on FULL series, exclude shock days for trend identification ----
# Use natural spline with df=6 (middle of spec curve)
panel$not_shock <- as.integer(!(panel$day %in% c(UPW, HUR)))
fit_trend_resid <- function(y, day, mask) {
  keep <- mask == 1 & !is.na(y)
  m <- lm(y[keep] ~ ns(day[keep], df = 6))
  resid_full <- rep(NA_real_, length(y))
  resid_full[keep] <- residuals(m)
  list(m = m, resid = resid_full, r2 = summary(m)$r.squared)
}

baseline_tab <- data.frame(outcome=character(), r2=numeric(), n_pre=integer(),
                           sigma_resid=numeric(), ljungbox_p=numeric(),
                           ar1_acf=numeric(), suggested_lag=integer(),
                           stringsAsFactors = FALSE)
resid_list <- list()

for (nm in names(OUT)) {
  y <- panel[[OUT[nm]]]
  fit <- fit_trend_resid(y, panel$day, panel$not_shock)
  resid_list[[nm]] <- fit$resid
  r <- fit$resid[!is.na(fit$resid)]
  lb <- Box.test(r, lag = 10, type = "Ljung-Box")$p.value
  a1 <- as.numeric(acf(r, plot=FALSE, lag.max=1)$acf[2])
  # Andrews/NW rule of thumb for HAC lag: floor(4*(T/100)^(2/9)) ~= 3 at T=93;
  # bump to 5 if AR(1) is strong.
  sug <- if (abs(a1) > 0.3) 5L else 3L
  baseline_tab <- rbind(baseline_tab, data.frame(
    outcome=nm, r2=round(fit$r2,3), n_pre=length(r),
    sigma_resid=round(sd(r),4), ljungbox_p=signif(lb,3),
    ar1_acf=round(a1,3), suggested_lag=sug
  ))
}

# ---- 4. ACF panel ----
acf_df <- do.call(rbind, lapply(names(resid_list), function(nm){
  r <- resid_list[[nm]]; r <- r[!is.na(r)]
  a <- acf(r, plot=FALSE, lag.max=20)
  data.frame(outcome=nm, lag=as.integer(a$lag), acf=as.numeric(a$acf))
}))
ci95 <- 1.96 / sqrt(sum(panel$not_shock == 1))
p_acf <- ggplot(acf_df, aes(lag, acf)) +
  geom_hline(yintercept = c(-ci95, ci95), lty = 2, color="grey50") +
  geom_hline(yintercept = 0, color="black") +
  geom_segment(aes(xend = lag, yend = 0)) +
  facet_wrap(~ outcome, ncol = 2) +
  labs(title = "Residual ACF after df=6 natural-spline trend (shock days excluded)",
       x = "Lag (days)", y = "ACF") +
  theme_bw(9)
ggsave(file.path(FIG, "fig03_acf.pdf"), p_acf, width = 9, height = 9)

# ---- 5. outcome distributions (level vs log) ----
dist_long <- panel |>
  select(day, all_of(TAXA), paste0("log_", TAXA)) |>
  pivot_longer(-day) |>
  mutate(
    taxon  = sub("^log_", "", name),
    scale  = ifelse(grepl("^log_", name), "log", "level")
  )
p_dist <- ggplot(dist_long, aes(value)) +
  geom_histogram(bins = 25) +
  facet_grid(scale ~ taxon, scales = "free") +
  labs(title = "Taxa outcome distributions: raw relative abundance vs log",
       x = NULL, y = "Count") +
  theme_bw(9)
ggsave(file.path(FIG, "fig04_outcome_dist.pdf"), p_dist, width = 10, height = 5)

# ---- 6. write baseline table ----
md_lines <- c(
  "# Pre-period baseline (df=6 natural spline, shock days excluded)",
  "",
  paste(capture.output(print(baseline_tab, row.names = FALSE)), collapse = "\n")
)
writeLines(md_lines, file.path(TAB, "tab_pre_baseline.md"))
print(baseline_tab, row.names = FALSE)
cat("\nOutputs:\n",
    "  ", file.path(FIG, "fig03_acf.pdf"), "\n",
    "  ", file.path(FIG, "fig04_outcome_dist.pdf"), "\n",
    "  ", file.path(TAB, "tab_pre_baseline.md"), "\n",
    "  ", file.path(CLN, "panel_transformed.csv"), "\n")
