
source("settings.R")

lapply(list.files("code",
  pattern = "\\.R$",
  ignore.case = TRUE,
  full.names = TRUE,
  recursive = TRUE),
  source
)
