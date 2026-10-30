
preprocess_nir <- function(data, window = 11, poly = 2, m = 0) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertDataTable(data, types = "numeric")
  assertNumber(window)
  assertNumber(poly, lower = 1)
  assertNumber(m, lower = 0)

  if (window %% 2 != 1) stop("Window must be odd")
  if (poly >= window) stop("Polynom must be smaller than the window")

  data <- as.matrix(data)

  # ── 2. Smoothing ───────────────────────────────────────────────────────────
  smoothed <- t(apply(data, 1, function(x) {
    sgolayfilt(x, p = poly, n = window, m = m)
    }
  ))

  # ── 3. SNV normalization ───────────────────────────────────────────────────
  normalized <- t(apply(smoothed, 1, function(x) {
      s <- sd(x)
      if (s == 0) return(rep(0, length(x)))
      (x - mean(x)) / s
    }
  ))

  normalized_dt <- as.data.table(normalized)
  colnames(normalized_dt) <- colnames(data)

  normalized_dt
}

assign_nir_region <- function(wl_nm) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertNumeric(wl_nm)

  vapply(wl_nm, function(w) {
    hit <- which(w >= nir_regions$lo & w <= nir_regions$hi)
    if (length(hit)) nir_regions$region[hit[1]] else "other"
    }, character(1))
}

regional_features <- function(Xmat, regions) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertMatrix(Xmat, mode = "numeric")
  assertCharacter(regions)

  wl <- as.numeric(sub("^W", "", colnames(Xmat)))
  rof <- assign_nir_region(wl)

  feats <- lapply(regions, function(r) {
    cols <- colnames(Xmat)[rof == r]
    sub <- Xmat[, cols, drop = FALSE]

    w <- wl[rof == r]
    slope <- apply(sub, 1, function(y) if (length(w) < 2) 0 else coef(lm(y ~ w))[2])

    data.table(
      mean = rowMeans(sub),
      slope = slope,
      sd = apply(sub, 1, sd)
    )
  })

  out <- cbindlist(feats)
  colnames(out) <- paste(rep(regions, each = 3), c("mean","slope","sd"), sep = "_")

  out
}

build_band_region_map <- function(band_names) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertCharacter(band_names)

  wl <- as.numeric(sub("^W", "", band_names))
  stopifnot(!any(is.na(wl)))
  dt <- data.table(
    band = band_names,
    wavelength = wl,
    region = assign_nir_region(wl),
    stringsAsFactors = FALSE
  )

  dt[order(wavelength)]
}

make_learner <- function(ln) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertCharacter(ln)

  switch(ln,
    NB = lrn("classif.naive_bayes", predict_type = "prob"),
    LDA = lrn("classif.lda", predict_type = "prob"),
    QDA = lrn("classif.qda", predict_type = "prob"),
    KNN = lrn("classif.kknn", predict_type = "prob"),
    SVM = lrn("classif.svm", predict_type = "prob", kernel = "radial", type = "C-classification"),
    RF = lrn("classif.ranger", predict_type = "prob", num.trees = 500),
    XGB = lrn("classif.xgboost", predict_type = "prob", nrounds = 40, eval_metric = "logloss"))
}

# Tuned search spaces (mtry.ratio, node.size, ...)
search_space <- function(ln) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertCharacter(ln)

  switch(ln,
    RF = ps(
      classif.ranger.mtry.ratio = p_dbl(0.2, 0.9),
      classif.ranger.min.node.size = p_int(1, 5)
    ),
    XGB = ps(
      classif.xgboost.eta = p_dbl(1e-2, 0.3, logscale = TRUE),
      classif.xgboost.max_depth = p_int(1, 3),
      classif.xgboost.nrounds   = p_int(20, 80)
    ),
    SVM = ps(
      classif.svm.cost  = p_dbl(1e-2, 1e2, logscale = TRUE),
      classif.svm.gamma = p_dbl(1e-3, 1e1, logscale = TRUE)
    ),
    KNN = ps(
      classif.kknn.k = p_int(3, 11)
    ),NULL)
}

# GraphLearner: oversample (minority, train only) -> learner. No in-fold
# feature selection (avoids ANOVA/contrasts failure under single-class days
# and the leakage of selecting features from labels).
build_learner <- function(ln, resample = "cv", folds = 3, n_evals = 25) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertCharacter(ln)

  gl <- GraphLearner$new(
    po("classbalancing", id = "oversample", adjust = "minor", reference = "major", shuffle = FALSE) %>>%
      po("learner", make_learner(ln)), id = ln
  )

  ss <- search_space(ln)
  if (is.null(ss)) return(gl)

  AutoTuner$new(
    learner = gl,
    resampling = rsmp(resample, folds = folds),
    measure = msr("classif.bacc"),
    search_space = ss,
    terminator = trm("evals", n_evals = n_evals),
    tuner = tnr("random_search")
  )
}

# One task from a feature matrix; day = group, never a feature
make_lodo_task <- function(Xmat, cls, day) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertDataTable(Xmat)
  assertFactor(cls, levels = CLASS_LEVELS)
  assertNumeric(day)

  dt <- as.data.table(Xmat)
  setnames(dt, make.names(colnames(dt)))

  dt[, c("class", "day") := .(cls, day)]

  task <- TaskClassif$new("nir", backend = dt, target = "class", positive = POS_LABEL)
  task$set_col_roles("day", roles = "group")
  task
}

# Evaluate all learners under LODO with pooled metrics
eval_lodo <- function(Xmat, cls, day, learners, tuned = TRUE) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertDataTable(Xmat)
  assertFactor(cls, levels = CLASS_LEVELS)
  assertNumeric(day)
  assertCharacter(learners)
  assertFlag(tuned)

  task  <- make_lodo_task(Xmat, cls, day)
  outer <- rsmp("loo")

  rbindlist(lapply(learners, function(ln) {
    lrnr <- if (tuned) build_learner(ln) else
      GraphLearner$new(
        po("classbalancing", id="oversample", adjust="minor", reference="major", shuffle=FALSE) %>>% 
          po("learner", make_learner(ln)), id = ln
      )

    rr <- tryCatch(resample(task, lrnr, outer, store_models = FALSE), error = function(e) NULL)
    if (is.null(rr)) return(data.table(
      learner = ln,
      auc = NA_real_,
      bacc = NA_real_,
      sens = NA_real_,
      spec = NA_real_)
    )

    pred <- rr$prediction()

    g <- function(m) tryCatch(as.numeric(pred$score(msr(m))), error=function(e) NA_real_)

    data.table(
      learner = ln,
      auc  = round(g("classif.auc"), 3),
      bacc = round(g("classif.bacc"), 3),
      sens = round(g("classif.sensitivity"), 3),
      spec = round(g("classif.specificity"), 3)
    )
  }))
}

seed_intervals <- function(dt, metric = "auc", conf = 0.95) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertDataTable(dt)
  assertString(metric)
  assertNumber(conf)

  m <- dt[!is.na(get(metric)), .(
    n = .N,
    mean = mean(get(metric)),
    sd = sd(get(metric)),
    lo_emp = quantile(get(metric), (1 - conf) / 2),
    hi_emp = quantile(get(metric), 1 - (1 - conf) / 2)
  ), by = learner]

  tcrit <- qt(1 - (1 - conf) / 2, df = pmax(m$n - 1, 1))

  m[, `:=`(
    se = sd / sqrt(n),
    # CI of the mean (narrow): expected value of the metric
    ci_mean_lo = mean - tcrit * sd / sqrt(n),
    ci_mean_hi = mean + tcrit * sd / sqrt(n),
    # spread of a single run (wide, honest): mean ± t*sd
    run_lo = mean - tcrit * sd,
    run_hi = mean + tcrit * sd
  )]

  m[, metric := metric][order(-mean)]
}

# Default GraphLearner (oversample %>>% learner) — identical to eval_lodo's
# tuned = FALSE branch.
default_glrn <- function(ln) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertCharacter(ln)

  GraphLearner$new(
    po("classbalancing", id = "oversample", adjust = "minor",
       reference = "major", shuffle = FALSE) %>>%
      po("learner", make_learner(ln)), id = ln
  )
}

# Returns truth + P(POS) by running LODO with the default GraphLearner, exactly
# as eval_lodo does, but keeping the prediction object instead of scoring.
pooled_pred <- function(ln, cols, days) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertCharacter(ln)
  assertNames(cols)
  assertNumeric(days)

  if (any(duplicated(days))) days <- unique(days)

  truth <- character(0)
  prob <- numeric(0)

  for (d in days) {
    te <- which(day_real == d)
    tr <- which(day_real %in% setdiff(days, d))

    if (length(unique(cls_real[tr])) < 2 || length(te) == 0) next

    task_tr <- make_lodo_task(Xr_real[tr, ..cols], cls_real[tr], day_real[tr])
    lrnr <- default_glrn(ln)

    ok <- tryCatch({ suppressWarnings(lrnr$train(task_tr)); TRUE }, error = function(e) FALSE)

    if (!ok) next

    dd <- as.data.table(Xr_real[te, ..cols])
    setnames(dd, make.names(colnames(dd)))

    p <- tryCatch(suppressWarnings(lrnr$predict_newdata(dd)$prob[, POS_LABEL]), error = function(e) rep(NA_real_, nrow(dd)))

    if (length(p) == length(te)) {
      prob  <- c(prob, as.numeric(p))
      truth <- c(truth, as.character(cls_real[te]))
    }
  }

  list(truth = factor(truth, levels = c(NEG_LABEL, POS_LABEL)), prob = prob)
}

auc_of <- function(pr) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertNumeric(pr$prob)

  ok <- !is.na(pr$prob)
  if (length(unique(pr$truth[ok])) < 2) return(NA_real_)

  as.numeric(pROC::auc(pROC::roc(pr$truth[ok], pr$prob[ok], levels = c(NEG_LABEL, POS_LABEL), direction = "<", quiet = TRUE)))
}

delong_ci <- function(pr) {

  # ── 1. Assertions ──────────────────────────────────────────────────────────
  assertNumeric(pr$prob)

  r <- pROC::roc(pr$truth, pr$prob, levels = c(NEG_LABEL, POS_LABEL), direction = "<", quiet = TRUE)
  as.numeric(pROC::ci.auc(r, method = "delong"))
}

save_tab <- function(df, name) {
  if (save_tables) write.csv2(df, file.path("data/results", paste0(name, ".csv")), row.names = FALSE)
}

save_fig <- function(plot, name, w, h) {
  if (save_figures) ggsave(file.path("plot", name), plot, width = w, height = h, dpi = 300, bg = "white")
}
