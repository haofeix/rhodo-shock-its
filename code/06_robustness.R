# 06_robustness.R -- simple specification curve for the main finding
# Only Rhodo x Upwelling (the headline ITS result).
# Variations: trend df (4/6/8/10), HAC lag (7/10/14/21), shock window (+/-1 day).

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
  library(splines); library(sandwich); library(lmtest)
})

ROOT <- "/Users/haofeixu/Desktop/Causal Inference/final_project"
CLN  <- file.path(ROOT, "data/clean")
FIG  <- file.path(ROOT, "output/figs")
panel <- read.csv(file.path(CLN, "panel_transformed.csv"))

HUR_DEFAULT <- 247:249
VAR <- "log_rhodo"

fit_its <- function(df, lag, upw_win) {
  d <- panel
  d$upw_var <- as.integer(d$day %in% upw_win)
  d$hur_var <- as.integer(d$day %in% HUR_DEFAULT)
  rhs <- c(paste0("ns(day, ", df, ")"), "upw_var", "hur_var")
  f <- as.formula(paste0(VAR, " ~ ", paste(rhs, collapse = " + ")))
  m <- lm(f, data = d)
  V <- NeweyWest(m, lag = lag, prewhite = FALSE, adjust = TRUE)
  ct <- coeftest(m, vcov. = V)
  sd_y <- sd(d[[VAR]], na.rm = TRUE)
  e <- ct["upw_var", "Estimate"]; se <- ct["upw_var", "Std. Error"]
  p <- ct["upw_var", "Pr(>|t|)"]
  data.frame(est = e, se = se, p = p,
             est_s = e/sd_y,
             lo_s  = (e - 1.96*se)/sd_y,
             hi_s  = (e + 1.96*se)/sd_y)
}

specs <- list(
  "Main (df=6, lag=14, 215-221)" = list(df=6,  lag=14, win=215:221),
  "df = 4"       = list(df=4,  lag=14, win=215:221),
  "df = 8"       = list(df=8,  lag=14, win=215:221),
  "df = 10"      = list(df=10, lag=14, win=215:221),
  "HAC lag = 7"  = list(df=6,  lag=7,  win=215:221),
  "HAC lag = 10" = list(df=6,  lag=10, win=215:221),
  "HAC lag = 21" = list(df=6,  lag=21, win=215:221),
  "Window 214-222" = list(df=6, lag=14, win=214:222),
  "Window 216-220" = list(df=6, lag=14, win=216:220)
)

res <- bind_rows(lapply(names(specs), function(nm) {
  s <- specs[[nm]]
  r <- fit_its(s$df, s$lag, s$win)
  r$spec <- nm
  r$category <- dplyr::case_when(
    nm == names(specs)[1]    ~ "Main",
    grepl("^df",      nm)    ~ "Trend df",
    grepl("^HAC lag", nm)    ~ "HAC lag",
    grepl("^Window",  nm)    ~ "Shock window"
  )
  r
}))

res <- res |>
  mutate(spec = factor(spec, levels = names(specs)),
         category = factor(category,
                           levels = c("Main", "Trend df", "HAC lag", "Shock window")))

cat_cols <- c("Main" = "#BA0C2F", "Trend df" = "#007D8A",
              "HAC lag" = "#215732", "Shock window" = "#F1B434")

p <- ggplot(res, aes(est_s, spec)) +
  geom_vline(xintercept = 0, lty = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = lo_s, xmax = hi_s, color = category),
                 height = 0.25, linewidth = 0.55) +
  geom_point(aes(color = category, size = category == "Main")) +
  scale_color_manual(values = cat_cols, name = "Variation") +
  scale_size_manual(values = c(`TRUE` = 4.5, `FALSE` = 2.8), guide = "none") +
  scale_y_discrete(limits = rev(levels(res$spec))) +
  labs(x = "Standardized effect (beta / SD)", y = NULL,
       title = "Specification curve: Rhodobacteraceae x Upwelling",
       subtitle = "9 specs across trend flexibility, HAC lag, and shock window; 95% HAC CI") +
  theme_bw(11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(color = "grey30", size = 10),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG, "fig11_spec_curve.pdf"), p, width = 8, height = 5)

# Console summary
cat("\n=== Rhodo x Upwelling: specification curve (9 specs) ===\n")
print(res |> mutate(across(where(is.numeric), \(x) round(x, 3))) |>
        select(spec, category, est_s, p), row.names = FALSE)

main_e <- res$est_s[res$category == "Main"]
cat(sprintf("\n  Main estimate        : %+.2f SD  (p = %.3f)\n",
            main_e, res$p[res$category == "Main"]))
cat(sprintf("  Estimate range       : [%+.2f, %+.2f] SD across 9 specs\n",
            min(res$est_s), max(res$est_s)))
cat(sprintf("  All same sign        : %s\n",
            all(sign(res$est_s) == sign(main_e))))
cat(sprintf("  Significant at p<0.05: %d / 9\n", sum(res$p < 0.05)))

cat("\nOutput: ", file.path(FIG, "fig11_spec_curve.pdf"), "\n")
