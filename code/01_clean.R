# 01_clean.R -- merge raw data into daily_panel.csv
# Outputs:
#   data/clean/daily_panel.csv
#   output/logs/data_audit.md
#   output/figs/fig01_env_timeseries.pdf
#   output/figs/fig02_outcomes_timeseries.pdf

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

ROOT <- "/Users/haofeixu/Desktop/Causal Inference/final_project"
RAW  <- file.path(ROOT, "data/raw")
CLN  <- file.path(ROOT, "data/clean")
FIG  <- file.path(ROOT, "output/figs")
LOG  <- file.path(ROOT, "output/logs")
dir.create(CLN, showWarnings=FALSE, recursive=TRUE)
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
dir.create(LOG, showWarnings=FALSE, recursive=TRUE)

UPW <- 215:221    # upwelling window
HUR <- 247:249    # Hurricane Earl window

# ---- 1. alpha diversity ----
alpha <- read_excel(file.path(RAW, "daily_alpha_diversity-1.0.xlsx")) |>
  mutate(day = as.integer(sub("Day_", "", Day))) |>
  select(day, Richness, Shannon, Simpson)

# ---- 2. environment ----
env_raw <- read_excel(file.path(RAW, "环境数据.xlsx"))
# "nd" -> NA, coerce to numeric
env <- env_raw |>
  rename(day = OrdinalDay) |>
  mutate(across(-day, ~ suppressWarnings(as.numeric(ifelse(.x == "nd", NA, .x)))))
# rename to short names
names(env) <- c("day", "air_temp", "ammonium", "atm_pressure", "chl_a",
                "dom_wave_pd", "macroalgae", "nitrite_nitrate", "phosphate",
                "salinity", "silicate", "tidal_dir", "water_level",
                "water_temp", "wave_height", "wind_speed")

# ---- 3. relative abundance -> indicator taxa ----
ab <- read_excel(file.path(RAW, "daily_mean_relative_abundance.xlsx"))
tax <- ab$taxonomy
day_cols <- grep("^Day_", names(ab), value = TRUE)

# taxonomic pattern matches (Greengenes format)
is_vibrio  <- grepl("f__Vibrionaceae",      tax, fixed = TRUE)
is_flavo   <- grepl("f__Flavobacteriaceae", tax, fixed = TRUE)
is_cyano   <- grepl("p__Cyanobacteria",     tax, fixed = TRUE)
is_sar11   <- grepl("f__SAR11",              tax, fixed = TRUE)  # Greengenes labels SAR11 / Pelagibacter here
is_rhodo   <- grepl("f__Rhodobacteraceae",   tax, fixed = TRUE)

cat("Indicator taxa OTU counts:\n")
cat("  Vibrionaceae:       ", sum(is_vibrio), "\n")
cat("  Flavobacteriaceae:  ", sum(is_flavo), "\n")
cat("  Cyanobacteria (p):  ", sum(is_cyano), "\n")
cat("  Pelagibacteraceae:  ", sum(is_sar11), "\n")
cat("  Rhodobacteraceae:   ", sum(is_rhodo), "\n")

sum_group <- function(mask) {
  if (sum(mask) == 0) return(rep(NA_real_, length(day_cols)))
  colSums(as.matrix(ab[mask, day_cols, drop = FALSE]), na.rm = TRUE)
}
taxa_df <- tibble(
  day            = as.integer(sub("Day_", "", day_cols)),
  vibrio         = sum_group(is_vibrio),
  flavobact      = sum_group(is_flavo),
  cyano          = sum_group(is_cyano),
  sar11          = sum_group(is_sar11),
  rhodo          = sum_group(is_rhodo)
)

# ---- 4. merge ----
panel <- alpha |>
  full_join(taxa_df, by = "day") |>
  full_join(env,     by = "day") |>
  arrange(day) |>
  mutate(
    upwelling = as.integer(day %in% UPW),
    hurricane = as.integer(day %in% HUR)
  )

write.csv(panel, file.path(CLN, "daily_panel.csv"), row.names = FALSE)
cat("\nWrote", file.path(CLN, "daily_panel.csv"), "rows =", nrow(panel), "\n")

# ---- 5. data audit ----
miss <- sapply(panel, function(x) sum(is.na(x)))
audit_lines <- c(
  "# Data audit",
  paste0("- Rows: ", nrow(panel), " (days ", min(panel$day), "-", max(panel$day), ")"),
  paste0("- Upwelling window: ", min(UPW), "-", max(UPW)),
  paste0("- Hurricane window: ", min(HUR), "-", max(HUR)),
  "",
  "## Missing counts by column",
  paste0("- `", names(miss), "`: ", miss)
)
writeLines(audit_lines, file.path(LOG, "data_audit.md"))

# ---- 6. diagnostic figures ----
shade <- list(
  annotate("rect", xmin = min(UPW)-0.5, xmax = max(UPW)+0.5, ymin = -Inf, ymax = Inf,
           alpha = 0.18, fill = "#2c7fb8"),
  annotate("rect", xmin = min(HUR)-0.5, xmax = max(HUR)+0.5, ymin = -Inf, ymax = Inf,
           alpha = 0.22, fill = "#d7301f")
)

env_long <- panel |>
  select(day, water_temp, salinity, atm_pressure, wind_speed, wave_height, nitrite_nitrate) |>
  pivot_longer(-day)
p_env <- ggplot(env_long, aes(day, value)) + shade +
  geom_line(na.rm = TRUE) + geom_point(size = 0.6, na.rm = TRUE) +
  facet_wrap(~ name, scales = "free_y", ncol = 2) +
  theme_bw(9) + labs(title = "Environmental variables; blue = upwelling, red = Hurricane Earl",
                     x = "Ordinal day", y = NULL)
ggsave(file.path(FIG, "fig01_env_timeseries.pdf"), p_env, width = 9, height = 6)

out_long <- panel |>
  select(day, Shannon, Richness, Simpson, vibrio, flavobact, cyano, sar11, rhodo) |>
  pivot_longer(-day)
p_out <- ggplot(out_long, aes(day, value)) + shade +
  geom_line(na.rm = TRUE) + geom_point(size = 0.6, na.rm = TRUE) +
  facet_wrap(~ name, scales = "free_y", ncol = 2) +
  theme_bw(9) + labs(title = "Outcomes (alpha diversity + indicator taxa)",
                     x = "Ordinal day", y = NULL)
ggsave(file.path(FIG, "fig02_outcomes_timeseries.pdf"), p_out, width = 9, height = 8)

cat("Figures written to", FIG, "\n")
