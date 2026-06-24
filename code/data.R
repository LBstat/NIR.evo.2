
build_data_main <- function() {

  # ── 1. Import data ────────────────────────────────────────────────────────
  data <- read_excel("data/raw/Dati filetti.xlsx")
  data_lightness <- read_excel("data/raw/Dati luminosità.xlsx")
  spectral_measures <- read_excel("data/raw/Spettri.xlsx")

  setDT(data)
  setDT(data_lightness)
  setDT(spectral_measures)

  setnames(data, c("ipoxantina", "Inosina", "K value", "mg di N/100g"), c("Hypoxanthine", "Inosine", "Kvalue", "mg N/100g"))

  # ── 2. Aggregation ────────────────────────────────────────────────────────
  data[, colnames(data) := lapply(.SD, function(x) {
    x[x < 0] <- 0
    x
  }), .SDcols = colnames(data)][, day := ((seq_len(.N) - 1) %/% 5) + 1]

  setnames(data_lightness, "campione", "sample")
  data_lightness <- data_lightness[, .(
    `L*` = median(`L*`, na.rm = TRUE),
    `a*` = median(`a*`, na.rm = TRUE),
    `b*` = median(`b*`, na.rm = TRUE),
    `C*` = median(`C*`, na.rm = TRUE),
    h = median(h, na.rm = TRUE)
  ), by = sample][, sample := as.numeric(gsub("campione_", "", sample))]

  data_analysis <- cbindlist(list(data, data_lightness))
  setcolorder(data_analysis, "sample")

  nir_vars <- grep("^[[:digit:]]*$", colnames(spectral_measures), value = TRUE)

  spectral_measures <- spectral_measures[spettro %between% c(1, 5)][, ..nir_vars]
  spectral_measures <- spectral_measures[
    , sample := ((seq_len(.N) - 1) %/% 5) + 1
  ][, lapply(.SD, function(x) median(x, na.rm = TRUE)), by = sample]

  data_analysis <- merge(data_analysis, spectral_measures)

  saveRDS(data_analysis, file = "data/intermediate/data analysis.rds")
}
