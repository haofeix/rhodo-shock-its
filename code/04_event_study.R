# 04_event_study.R -- dynamic event studies for upwelling and hurricane
# Model (per shock S):
#   Y_t = alpha + ns(day, 6) + sum_{k != -1} beta_k * 1{day - T_S = k} + (other shock dummy) + eps_t
# k ranges over [-7, +14]; k = -1 is reference (beta_{-1} := 0).
# Pre-trend test: joint Wald H0: beta_k = 0 for k in -7..-2 using HAC vcov.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(ggtext)
  library(splines); library(sandwich); library(lmtest)
})

ROOT <- "/Users/haofeixu/Desktop/Causal Inference/final_project"
CLN  <- file.path(ROOT, "data/clean")
FIG  <- file.path(ROOT, "output/figs")
TAB  <- file.path(ROOT, "output/tables")
panel <- read.csv(file.path(CLN, "panel_transformed.csv"))

UPW_T0 <- 215L; HUR_T0 <- 247L
UPW_SPAN <- c(0L, 6L)   # treatment window in event time (7 days)
HUR_SPAN <- c(0L, 2L)   # treatment window in event time (3 days)
KLO <- -7L; KHI <- 14L
ks <- KLO:KHI
ks_use <- setdiff(ks, -1L)

out_spec <- tribble(
  ~var,            ~display,                ~type,        ~lag,
  "Shannon",       "Shannon",               "primary",    3L,
  "Richness",      "Richness",              "primary",    3L,
  "Simpson",       "Simpson",               "primary",    3L,
  "log_vibrio",    "log Vibrionaceae",      "secondary",  7L,
  "log_flavobact", "log Flavobacteriaceae", "secondary",  3L,
  "log_cyano",     "log Cyanobacteria",     "secondary",  7L,
  "log_sar11",     "log SAR11",             "secondary",  3L,
  "log_rhodo",     "log Rhodobacteraceae",  "secondary", 14L
)

# Event-time dummies. Names use m/p prefix (minus/plus) to avoid formula parsing of "-"
k_to_name <- function(k) ifelse(k < 0, paste0("m", abs(k)),
                                paste0("p", sprintf("%02d", k)))
name_to_k <- function(s) {
  kk <- as.integer(sub("^[mp]", "", s))
  ifelse(startsWith(s, "m"), -kk, kk)
}
make_ev_dummies <- function(day_vec, T0, ks_use, prefix) {
  out <- lapply(ks_use, function(k) as.integer((day_vec - T0) == k))
  names(out) <- paste0(prefix, "_", k_to_name(ks_use))
  as.data.frame(out)
}
ev_U <- make_ev_dummies(panel$day, UPW_T0, ks_use, "upw")
ev_H <- make_ev_dummies(panel$day, HUR_T0, ks_use, "hur")
pdat <- cbind(panel, ev_U, ev_H)

# Fit one shock's event study
fit_event <- function(var, lag, which_shock) {
  ev_cols <- if (which_shock == "U") names(ev_U) else names(ev_H)
  nuis    <- if (which_shock == "U") "hurricane" else "upwelling"
  rhs <- c("ns(day, 6)", ev_cols, nuis)
  f <- as.formula(paste0(var, " ~ ", paste(rhs, collapse = " + ")))
  m <- lm(f, data = pdat)
  V <- NeweyWest(m, lag = lag, prewhite = FALSE, adjust = TRUE)
  ct <- coeftest(m, vcov. = V)
  prefix <- if (which_shock == "U") "upw_" else "hur_"
  rows <- grep(paste0("^", prefix, "[mp]"), rownames(ct), value = TRUE)
  ks_fitted <- name_to_k(sub(prefix, "", rows))
  est_df <- data.frame(
    var = var, shock = which_shock, term = rows, k = ks_fitted,
    est = ct[rows, "Estimate"], se = ct[rows, "Std. Error"],
    p   = ct[rows, "Pr(>|t|)"], row.names = NULL
  )
  est_df$lo <- est_df$est - 1.96 * est_df$se
  est_df$hi <- est_df$est + 1.96 * est_df$se
  ref <- data.frame(var = var, shock = which_shock,
                    term = paste0(prefix, "m1"), k = -1L,
                    est = 0, se = NA, p = NA, lo = NA, hi = NA)
  est_df <- rbind(est_df, ref) |> arrange(k)

  # Pre-trend summary: count of |t| > 1.96 among k in -7..-2 and max |t|.
  # (HAC vcov on one-day indicator dummies is structurally rank-deficient,
  # so the standard joint F/Wald test is not well-defined here.)
  pre_rows <- rows[ks_fitted %in% -7:-2]
  pre_t <- if (length(pre_rows) > 0) ct[pre_rows, "Estimate"] / ct[pre_rows, "Std. Error"] else numeric(0)
  list(est = est_df,
       pre_n_sig = sum(abs(pre_t) > 1.96),
       pre_max_abs_t = if (length(pre_t)) max(abs(pre_t)) else NA_real_)
}

# Run all 16 event studies
all_est <- list(); pt <- data.frame()
for (i in seq_len(nrow(out_spec))) {
  for (s in c("U","H")) {
    res <- fit_event(out_spec$var[i], out_spec$lag[i], s)
    all_est[[paste(out_spec$var[i], s, sep="_")]] <- res$est
    pt <- rbind(pt, data.frame(var = out_spec$var[i], shock = s,
                               pre_n_sig = res$pre_n_sig,
                               pre_max_abs_t = res$pre_max_abs_t))
  }
}
all_est <- bind_rows(all_est) |>
  left_join(out_spec |> select(var, display, type), by = "var") |>
  mutate(display = factor(display, levels = out_spec$display))

# Standardize by each outcome's SD
all_est <- all_est |>
  group_by(var) |>
  mutate(sd_y  = sd(pdat[[var[1]]], na.rm = TRUE),
         est_s = est / sd_y,
         lo_s  = lo  / sd_y,
         hi_s  = hi  / sd_y) |>
  ungroup()

# Plotter
plot_shock <- function(which_shock, shock_title, fname, span) {
  d <- all_est |> filter(shock == which_shock)
  pt_sub <- pt |> filter(shock == which_shock) |>
    left_join(out_spec |> select(var, display), by = "var") |>
    mutate(display = factor(display, levels = out_spec$display),
           annot = sprintf("pre-period |t|>1.96: %d / 6  (max |t|=%.2f)",
                           pre_n_sig, pre_max_abs_t))

  p <- ggplot(d, aes(k, est_s)) +
    annotate("rect", xmin = span[1]-0.5, xmax = span[2]+0.5,
             ymin = -Inf, ymax = Inf, alpha = 0.18, fill = "#d7301f") +
    geom_hline(yintercept = 0, color = "grey40") +
    geom_vline(xintercept = -1, lty = 3, color = "grey60") +
    geom_errorbar(aes(ymin = lo_s, ymax = hi_s), width = 0.35,
                  color = "grey40", na.rm = TRUE) +
    geom_line(color = "grey30", linewidth = 0.35, na.rm = TRUE) +
    geom_point(aes(color = (k >= span[1] & k <= span[2])), size = 1.6) +
    scale_color_manual(values = c(`TRUE` = "#a50f15", `FALSE` = "black"),
                       guide = "none") +
    geom_text(data = pt_sub,
              aes(x = KLO, y = Inf, label = annot),
              hjust = 0, vjust = 1.3, size = 2.8, inherit.aes = FALSE,
              color = "grey30") +
    facet_wrap(~ display, ncol = 4, scales = "free_y") +
    labs(x = "Event time k (days relative to shock onset)",
         y = "Standardized effect (beta_k / SD of outcome)",
         title = shock_title,
         subtitle = "Red band = treatment window; dotted line at k = -1 (reference)") +
    theme_bw(9) +
    theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(FIG, fname), p, width = 11.5, height = 5.8)
}

plot_shock("U", "Event study: Upwelling (T0 = day 215)",
           "fig06_event_upwelling.pdf", UPW_SPAN)
plot_shock("H", "Event study: Hurricane Earl (T0 = day 247)",
           "fig07_event_hurricane.pdf", HUR_SPAN)

# ---- Pre-trend summary table ----
pt_wide <- pt |>
  mutate(shock = recode(shock, U = "Upwelling", H = "Hurricane")) |>
  left_join(out_spec |> select(var, display), by = "var") |>
  mutate(display = factor(display, levels = out_spec$display)) |>
  arrange(display, shock)

pt_md <- c(
  "# Pre-trend diagnostic (event time k in -7..-2)",
  "",
  "Note: HAC vcov on one-day event-time dummies is rank-deficient in this design,",
  "so a joint F/Wald test is not well-defined. We report:",
  "  - `pre_n_sig`: count of pre-period coefficients with individual HAC |t| > 1.96 (out of 6)",
  "  - `pre_max_abs_t`: max |t| over the six pre-period coefficients.",
  "Under the null of no pre-trend, we expect pre_n_sig <= 1 by chance.",
  "",
  "| Outcome | Shock | pre_n_sig | max |t| |",
  "|---|---|---|---|"
)
for (i in seq_len(nrow(pt_wide))) {
  r <- pt_wide[i, ]
  pt_md <- c(pt_md, sprintf("| %s | %s | %d | %.2f |",
                            r$display, r$shock, r$pre_n_sig, r$pre_max_abs_t))
}
writeLines(pt_md, file.path(TAB, "tab_pretrend.md"))

cat("\n=== Pre-trend diagnostic (|t|>1.96 count out of 6, max |t|) ===\n")
print(pt_wide |> select(display, shock, pre_n_sig, pre_max_abs_t), row.names = FALSE)

# Highlight which post-shock event-time coefficients are individually significant
cat("\n=== Individually significant post-shock coefficients (k >= 0, p < 0.05) ===\n")
sig_post <- all_est |>
  filter(k >= 0, !is.na(p), p < 0.05) |>
  mutate(shock = recode(shock, U = "Upwelling", H = "Hurricane")) |>
  select(shock, display, k, est = est_s, p) |>
  arrange(shock, display, k)
print(sig_post, row.names = FALSE)

cat("\nOutputs:\n  ", file.path(FIG, "fig06_event_upwelling.pdf"),
    "\n  ", file.path(FIG, "fig07_event_hurricane.pdf"),
    "\n  ", file.path(TAB, "tab_pretrend.md"), "\n")
