# 03_its_main.R -- main ITS regression + HAC SE + primary/secondary tables
# Model: Y_t = alpha + ns(day, 6) + beta_U*upwelling + beta_H*hurricane + eps_t
# No environmental controls (avoid bad-control bias from post-treatment mediators).

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(ggtext)
  library(splines); library(sandwich); library(lmtest)
})

ROOT <- "/Users/haofeixu/Desktop/Causal Inference/final_project"
CLN  <- file.path(ROOT, "data/clean")
FIG  <- file.path(ROOT, "output/figs")
TAB  <- file.path(ROOT, "output/tables")
panel <- read.csv(file.path(CLN, "panel_transformed.csv"))

# Outcome specification
out_spec <- tribble(
  ~var,            ~display,                ~type,        ~lag, ~pred_U, ~pred_H,
  "Shannon",       "Shannon",               "primary",    3L,   "-",     "-",
  "Richness",      "Richness",              "primary",    3L,   "?",     "?",
  "Simpson",       "Simpson",               "primary",    3L,   "-",     "-",
  "log_vibrio",    "log Vibrionaceae",      "secondary",  7L,   "?",     "-",
  "log_flavobact", "log Flavobacteriaceae", "secondary",  3L,   "?",     "+",
  "log_cyano",     "log Cyanobacteria",     "secondary",  7L,   "-",     "?",
  "log_sar11",     "log SAR11",             "secondary",  3L,   "-",     "?",
  "log_rhodo",     "log Rhodobacteraceae",  "secondary", 14L,   "+",     "?"
)

# Fit + HAC SE for one outcome
fit_one <- function(var, lag) {
  f <- as.formula(paste0(var, " ~ ns(day, 6) + upwelling + hurricane"))
  m <- lm(f, data = panel)
  V <- NeweyWest(m, lag = lag, prewhite = FALSE, adjust = TRUE)
  ct <- coeftest(m, vcov. = V)
  sigma_y <- sd(panel[[var]], na.rm = TRUE)
  get_row <- function(name) {
    b <- ct[name, ]
    list(est=b["Estimate"], se=b["Std. Error"], p=b["Pr(>|t|)"],
         lo=b["Estimate"]-1.96*b["Std. Error"],
         hi=b["Estimate"]+1.96*b["Std. Error"])
  }
  U <- get_row("upwelling"); H <- get_row("hurricane")
  data.frame(
    var=var, lag=lag, sigma_y=sigma_y,
    est_U=U$est, se_U=U$se, p_U=U$p, lo_U=U$lo, hi_U=U$hi,
    est_H=H$est, se_H=H$se, p_H=H$p, lo_H=H$lo, hi_H=H$hi,
    row.names = NULL
  )
}

results <- do.call(rbind, Map(fit_one, out_spec$var, out_spec$lag)) |>
  left_join(out_spec |> select(-lag), by = "var")

# BH-FDR across 5 secondary taxa, per shock
sec_idx <- results$type == "secondary"
results$qU <- NA_real_; results$qH <- NA_real_
results$qU[sec_idx] <- p.adjust(results$p_U[sec_idx], method = "BH")
results$qH[sec_idx] <- p.adjust(results$p_H[sec_idx], method = "BH")

# Pretty-print helper
fmt_pval <- function(p, q=NULL){
  if (is.na(p)) return("NA")
  s <- if (p < 0.001) "<0.001" else sprintf("%.3f", p)
  if (!is.null(q) && !is.na(q)) s <- paste0(s, " (q=", sprintf("%.2f", q), ")")
  s
}
fmt_est <- function(e, lo, hi) sprintf("%.3f [%.3f, %.3f]", e, lo, hi)
dir_match <- function(est, pred){
  if (pred == "?") return("")
  s <- if (est > 0) "+" else "-"
  if (s == pred) "✓" else "✗"
}

# ---- Primary table (Shannon, Richness, Simpson) ----
prim <- results |> filter(type == "primary")
prim_md <- c(
  "# Primary outcomes: ITS main spec",
  "",
  "Model: `Y_t = alpha + ns(day, 6) + beta_U*upwelling + beta_H*hurricane + eps_t`",
  "No environmental controls. HAC Newey-West SE (lag per outcome). Pre-registered predictions in parentheses.",
  "",
  "| Outcome | HAC lag | beta_U (pred) | 95% CI | p | beta_H (pred) | 95% CI | p |",
  "|---|---|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(prim))) {
  r <- prim[i, ]
  prim_md <- c(prim_md, sprintf(
    "| %s | %d | %.3f (%s%s) | [%.3f, %.3f] | %s | %.3f (%s%s) | [%.3f, %.3f] | %s |",
    r$display, r$lag,
    r$est_U, r$pred_U, dir_match(r$est_U, r$pred_U),
    r$lo_U, r$hi_U, fmt_pval(r$p_U),
    r$est_H, r$pred_H, dir_match(r$est_H, r$pred_H),
    r$lo_H, r$hi_H, fmt_pval(r$p_H)
  ))
}
writeLines(prim_md, file.path(TAB, "tab_main_primary.md"))

# ---- Secondary table (5 taxa with BH-FDR) ----
sec <- results |> filter(type == "secondary")
sec_md <- c(
  "# Secondary outcomes (log taxa): ITS main spec",
  "",
  "BH-FDR applied across 5 taxa separately within each shock.",
  "Log-fold-change interpretation: beta = ln(post/pre).",
  "",
  "| Taxon | HAC lag | beta_U (pred) | 95% CI | p (q) | beta_H (pred) | 95% CI | p (q) |",
  "|---|---|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(sec))) {
  r <- sec[i, ]
  sec_md <- c(sec_md, sprintf(
    "| %s | %d | %.3f (%s%s) | [%.3f, %.3f] | %s | %.3f (%s%s) | [%.3f, %.3f] | %s |",
    r$display, r$lag,
    r$est_U, r$pred_U, dir_match(r$est_U, r$pred_U),
    r$lo_U, r$hi_U, fmt_pval(r$p_U, r$qU),
    r$est_H, r$pred_H, dir_match(r$est_H, r$pred_H),
    r$lo_H, r$hi_H, fmt_pval(r$p_H, r$qH)
  ))
}
writeLines(sec_md, file.path(TAB, "tab_main_secondary.md"))

# ---- Forest plot (standardized effect sizes) ----
long_U <- results |> mutate(
  shock = "Upwelling (days 215-221)", pred = pred_U,
  est_s = est_U / sigma_y, lo_s = lo_U / sigma_y, hi_s = hi_U / sigma_y,
  p = p_U, q = qU
) |> select(display, type, shock, pred, est_s, lo_s, hi_s, p, q)

long_H <- results |> mutate(
  shock = "Hurricane Earl (days 247-249)", pred = pred_H,
  est_s = est_H / sigma_y, lo_s = lo_H / sigma_y, hi_s = hi_H / sigma_y,
  p = p_H, q = qH
) |> select(display, type, shock, pred, est_s, lo_s, hi_s, p, q)

# Build y-axis labels: bold primary via markdown
out_levels_rev <- rev(out_spec$display)
out_labels <- setNames(
  ifelse(out_spec$type == "primary",
         paste0("**", out_spec$display, "**"),
         out_spec$display),
  out_spec$display
)

long <- bind_rows(long_U, long_H) |>
  mutate(
    display = factor(display, levels = out_levels_rev),
    shock   = factor(shock, levels = c("Upwelling (days 215-221)",
                                       "Hurricane Earl (days 247-249)")),
    direction_ok = dplyr::case_when(
      pred == "?"                                              ~ "no prediction",
      (pred == "+" & est_s > 0) | (pred == "-" & est_s < 0)    ~ "matches prediction",
      TRUE                                                     ~ "against prediction"
    ),
    is_sig = !is.na(p) & p < 0.05,
    # combined status: hue = prediction verdict; saturation = significance
    status = dplyr::case_when(
      is_sig  & direction_ok == "matches prediction" ~ "Sig. - matches",
      is_sig  & direction_ok == "against prediction" ~ "Sig. - against",
      is_sig  & direction_ok == "no prediction"      ~ "Sig. - no prior",
      !is_sig & direction_ok == "matches prediction" ~ "n.s. - matches",
      !is_sig & direction_ok == "against prediction" ~ "n.s. - against",
      !is_sig & direction_ok == "no prediction"      ~ "n.s. - no prior"
    ),
    status = factor(status, levels = c(
      "Sig. - matches","Sig. - against","Sig. - no prior",
      "n.s. - matches","n.s. - against","n.s. - no prior"
    )),
    point_size = ifelse(type == "primary", 4.5, 3.2)
  )

status_cols <- c(
  "Sig. - matches"  = "#08519c",   # dark blue
  "Sig. - against"  = "#a50f15",   # dark red
  "Sig. - no prior" = "#252525",   # near black
  "n.s. - matches"  = "#9ecae1",   # pale blue
  "n.s. - against"  = "#fcae91",   # pale red
  "n.s. - no prior" = "#cccccc"    # pale grey
)

# Shared-x range: cover all CI ends with small padding
xlims <- range(c(long$lo_s, long$hi_s), na.rm = TRUE)
xlims <- xlims + c(-0.1, 0.1) * diff(xlims)

# Horizontal separator line between primary (top 3) and secondary (bottom 5) in the rev ordering
n_sec <- sum(out_spec$type == "secondary")
sep_y <- n_sec + 0.5  # between secondary and primary blocks on the rev-ordered y

p_fore <- ggplot(long, aes(est_s, display)) +
  geom_hline(yintercept = sep_y, lty = 3, color = "grey60") +
  geom_vline(xintercept = 0, lty = 2, color = "grey50") +
  geom_errorbar(aes(xmin = lo_s, xmax = hi_s), width = 0.22,
                color = "grey35", linewidth = 0.45, orientation = "y") +
  geom_point(aes(shape = type, fill = status, size = point_size),
             color = "black", stroke = 0.5) +
  scale_shape_manual(values = c(primary = 21, secondary = 22),
                     guide = guide_legend(override.aes = list(size = 4, fill = "grey60"))) +
  scale_fill_manual(values = status_cols, drop = FALSE,
                    guide = guide_legend(override.aes = list(shape = 21, size = 4.5))) +
  scale_size_identity() +
  scale_y_discrete(labels = out_labels) +
  coord_cartesian(xlim = xlims) +
  facet_wrap(~ shock, nrow = 1) +
  labs(x = "Standardized effect (beta / SD of outcome)", y = NULL,
       shape = "Outcome tier", fill = "Estimate status",
       title = "ITS main specification: standardized shock effects, 95% HAC CI",
       subtitle = "Primary outcomes in bold (dotted separator); saturation = significance, hue = direction vs. prior") +
  theme_bw(10) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        axis.text.y = element_markdown(size = 10),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG, "fig05_coef_forest.pdf"), p_fore, width = 10, height = 6.5)

# ---- Console summary ----
cat("\n=== PRIMARY ===\n"); print(prim[, c("display","lag","est_U","p_U","est_H","p_H")], row.names=FALSE)
cat("\n=== SECONDARY (with BH-FDR) ===\n")
print(sec[, c("display","lag","est_U","p_U","qU","est_H","p_H","qH")], row.names=FALSE)
cat("\nTables:\n  ", file.path(TAB, "tab_main_primary.md"),
    "\n  ", file.path(TAB, "tab_main_secondary.md"),
    "\nFigure:\n  ", file.path(FIG, "fig05_coef_forest.pdf"), "\n")
