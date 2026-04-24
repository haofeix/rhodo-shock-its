# 05_headline_fig.R -- polished poster-ready figures focused on the core story
# Headline:  Rhodobacteraceae shows opposing responses to the two shocks
# Supporting: Vibrio x Hurricane (bloom termination); Cyano x Upwelling (regime-shift trigger)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(ggtext)
  library(splines); library(sandwich); library(lmtest); library(patchwork)
})

ROOT <- "/Users/haofeixu/Desktop/Causal Inference/final_project"
CLN  <- file.path(ROOT, "data/clean")
FIG  <- file.path(ROOT, "output/figs")
panel <- read.csv(file.path(CLN, "panel_transformed.csv"))

UPW_T0 <- 215L; HUR_T0 <- 247L
UPW_SPAN <- c(0L, 6L); HUR_SPAN <- c(0L, 2L)
KLO <- -7L; KHI <- 14L
ks_use <- setdiff(KLO:KHI, -1L)

k_to_name <- function(k) ifelse(k < 0, paste0("m", abs(k)), paste0("p", sprintf("%02d", k)))
name_to_k <- function(s) {
  kk <- as.integer(sub("^[mp]", "", s))
  ifelse(startsWith(s, "m"), -kk, kk)
}
make_ev_dummies <- function(day_vec, T0, prefix) {
  out <- lapply(ks_use, function(k) as.integer((day_vec - T0) == k))
  names(out) <- paste0(prefix, "_", k_to_name(ks_use))
  as.data.frame(out)
}
ev_U <- make_ev_dummies(panel$day, UPW_T0, "upw")
ev_H <- make_ev_dummies(panel$day, HUR_T0, "hur")
pdat <- cbind(panel, ev_U, ev_H)

# Fit event study for one outcome x shock, return tidy coef table
get_es <- function(var, lag, shock) {
  ev_cols <- if (shock == "U") names(ev_U) else names(ev_H)
  nuis    <- if (shock == "U") "hurricane" else "upwelling"
  prefix  <- if (shock == "U") "upw_" else "hur_"
  f <- as.formula(paste0(var, " ~ ns(day, 6) + ",
                         paste(ev_cols, collapse=" + "), " + ", nuis))
  m <- lm(f, data = pdat)
  V <- NeweyWest(m, lag = lag, prewhite = FALSE, adjust = TRUE)
  ct <- coeftest(m, vcov. = V)
  rows <- grep(paste0("^", prefix, "[mp]"), rownames(ct), value = TRUE)
  ks_fitted <- name_to_k(sub(prefix, "", rows))
  sd_y <- sd(pdat[[var]], na.rm = TRUE)
  df <- data.frame(k = ks_fitted,
                   est = ct[rows, "Estimate"],
                   se  = ct[rows, "Std. Error"],
                   p   = ct[rows, "Pr(>|t|)"])
  df$lo <- df$est - 1.96 * df$se
  df$hi <- df$est + 1.96 * df$se
  df$est_s <- df$est / sd_y
  df$lo_s  <- df$lo  / sd_y
  df$hi_s  <- df$hi  / sd_y
  ref <- data.frame(k = -1L, est = 0, se = NA, p = NA,
                    lo = NA, hi = NA, est_s = 0, lo_s = NA, hi_s = NA)
  df <- rbind(df, ref) |> arrange(k)
  # Pre-trend summary
  pre <- df |> filter(k %in% -7:-2) |> mutate(tt = est/se)
  df_attrs <- list(
    pre_n_sig    = sum(abs(pre$tt) > 1.96, na.rm = TRUE),
    pre_max_abs_t = max(abs(pre$tt), na.rm = TRUE),
    # post-shock mean effect (inside treatment window)
    post_mean_s = mean(df$est_s[df$k >= (if (shock=="U") UPW_SPAN else HUR_SPAN)[1] &
                                df$k <= (if (shock=="U") UPW_SPAN else HUR_SPAN)[2]], na.rm = TRUE),
    post_peak_s = {
      vals <- df$est_s[df$k >= 0 & df$k <= 14]
      vals[which.max(abs(vals))]
    }
  )
  attr(df, "summary") <- df_attrs
  df
}

# Polished event-study panel plotter
plot_es_panel <- function(es, span, title, subtitle_text,
                          shade_color = "#d7301f", panel_letter = NULL,
                          ylims = NULL) {
  s <- attr(es, "summary")
  annot <- sprintf("pre-trend |t|>1.96: %d / 6   (max |t| = %.2f)",
                   s$pre_n_sig, s$pre_max_abs_t)

  if (is.null(ylims)) {
    ylims <- range(c(es$lo_s, es$hi_s), na.rm = TRUE)
    ylims <- ylims + c(-0.15, 0.25) * diff(ylims)
  }

  p <- ggplot(es, aes(k, est_s)) +
    annotate("rect", xmin = span[1]-0.5, xmax = span[2]+0.5,
             ymin = -Inf, ymax = Inf, fill = shade_color, alpha = 0.14) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = -1, lty = 3, color = "grey55") +
    geom_errorbar(aes(ymin = lo_s, ymax = hi_s), width = 0.35,
                  color = "grey40", linewidth = 0.5, na.rm = TRUE) +
    geom_line(color = "grey25", linewidth = 0.5) +
    geom_point(aes(fill = (k >= span[1] & k <= span[2])),
               shape = 21, size = 2.4, stroke = 0.45, color = "black") +
    scale_fill_manual(values = c(`TRUE` = shade_color, `FALSE` = "white"),
                      guide = "none") +
    annotate("text", x = KLO, y = ylims[2], label = annot,
             hjust = 0, vjust = 1.2, size = 3.1, color = "grey25") +
    coord_cartesian(ylim = ylims) +
    labs(title = title, subtitle = subtitle_text,
         x = "Event time k (days)",
         y = "Standardized effect (beta_k / SD)") +
    theme_bw(11) +
    theme(
      plot.title = element_text(face = "bold", size = 12.5),
      plot.subtitle = element_text(color = "grey30", size = 10),
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 4, 2, 2, "pt")
    )
  if (!is.null(panel_letter)) {
    p <- p + annotate("text", x = KLO - 1.5, y = ylims[2], label = panel_letter,
                       fontface = "bold", size = 5, hjust = 0, vjust = 1.1)
  }
  p
}

# ===== Headline figure: Rhodobacteraceae ×2 =====
es_rhodo_U <- get_es("log_rhodo", 14, "U")
es_rhodo_H <- get_es("log_rhodo", 14, "H")

# Shared y-axis, symmetric around zero so both panels are visually balanced
max_abs <- max(abs(c(es_rhodo_U$lo_s, es_rhodo_U$hi_s,
                     es_rhodo_H$lo_s, es_rhodo_H$hi_s)), na.rm = TRUE)
shared_ylim <- c(-1.1, 1.15) * max_abs

p_r_U <- plot_es_panel(
  es_rhodo_U, UPW_SPAN,
  "Upwelling (days 215-221)",
  "Cold nutrient-rich water replaces warm-water community",
  shade_color = "#215732", panel_letter = "A", ylims = shared_ylim
)
p_r_H <- plot_es_panel(
  es_rhodo_H, HUR_SPAN,
  "Hurricane Earl (days 247-249)",
  "Mixing & sediment resuspension release organic matter",
  shade_color = "#BA0C2F", panel_letter = "B", ylims = shared_ylim
)

headline <- (p_r_U | p_r_H) +
  plot_annotation(
    title    = "Rhodobacteraceae responds in opposite directions to two exogenous shocks",
    subtitle = "Event-study coefficients relative to day before shock (k = −1); 95% HAC CI; pre-trend clean (0/6 violations on both sides)",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14.5),
      plot.subtitle = element_text(color = "grey30", size = 10.5),
      plot.margin   = margin(2, 4, 2, 2, "pt")
    )
  ) &
  theme(plot.margin = margin(2, 4, 2, 2, "pt"))
ggsave(file.path(FIG, "fig08_rhodo_headline.pdf"), headline, width = 12, height = 4.8)

# ===== Console summary =====
cat("\n=== Core result: Rhodobacteraceae event studies ===\n\n")
summarize_one <- function(es, name) {
  s <- attr(es, "summary")
  cat(sprintf("%-30s pre-trend: %d/6 (max|t|=%.2f)   in-window mean: %+.2f SD   peak: %+.2f SD\n",
              name, s$pre_n_sig, s$pre_max_abs_t, s$post_mean_s, s$post_peak_s))
}
summarize_one(es_rhodo_U, "log_rhodo x Upwelling")
summarize_one(es_rhodo_H, "log_rhodo x Hurricane")

cat("\nFigure: ", file.path(FIG, "fig08_rhodo_headline.pdf"), "\n")
