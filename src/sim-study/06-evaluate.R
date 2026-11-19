# Evaluate extracted fits against the structural predictor, independently of DGP.
# Source 01-dgp.R and 03-extract.R first (or use the already loaded worker).
# All slopes use original covariate units. MAD means mean absolute error.

# Truth uses the ORIGINAL training sample for centering and projection, even
# on a grid. Spatial scalar truth gives each region one value, not n repeats.
evaluation_truth <- function(scenario, grid = seq(-.5, .5, length.out = 50L)) {
  dat <- scenario$dat
  parts <- lapply(1:2, function(j) {
    effect <- function(x) f_true(scenario$predictor, x, x)[[paste0("f", j)]]
    split_true_effect(effect, dat[[paste0("x", j)]], new_x = grid)
  })
  names(parts) <- c("f1", "f2")
  offsets <- truth_centering(scenario)
  train <- training_truth(scenario, offsets)
  gamma <- if (scenario$predictor == "X3") scenario$gamma[levels(dat$region)] else NULL
  on_grid <- grid_truth(scenario$predictor, grid, gamma, offsets)
  train$fcont_full <- train$int + train$fcont
  on_grid$fcont_full <- on_grid$int + on_grid$fcont
  for (j in 1:2) {
    # Center the affine projection by its training mean. Its remainder is
    # then orthogonal to both the intercept and linear term, as in DR.
    key <- paste0("f", j)
    train[[paste0(key, "L")]] <- parts[[j]]$slope * (dat[[paste0("x", j)]] - parts[[j]]$mean_x)
    on_grid[[paste0(key, "L")]] <- parts[[j]]$slope * (grid - parts[[j]]$mean_x)
    train[[paste0(key, "NL")]] <- train[[key]] - train[[paste0(key, "L")]]
    on_grid[[paste0(key, "NL")]] <- on_grid[[key]] - on_grid[[paste0(key, "L")]]
  }
  parameters <- list(int = train$int, slope1 = parts$f1$slope, slope2 = parts$f2$slope,
                     spat = if (!is.null(gamma)) gamma - mean(gamma) else NULL)
  list(parameters = parameters, grid = on_grid, training = train, decomposition = parts)
}

# Preserve vector names and grid dimensions when computing scalar/location metrics.
evaluation_shape <- function(values, template) {
  attributes(values) <- attributes(template)
  values
}

evaluation_error <- function(estimate, truth) {
  if (is.list(truth)) return(setNames(lapply(seq_along(truth), function(i) {
    region <- if (is.list(estimate)) estimate[[if (is.null(names(estimate))) i else names(truth)[i]]] else estimate
    evaluation_error(region, truth[[i]])
  }), names(truth)))
  stopifnot(length(estimate) == length(truth), identical(dim(estimate), dim(truth)))
  evaluation_shape(as.numeric(estimate) - as.numeric(truth), truth)
}

evaluation_interval <- function(bounds, truth) {
  if (is.list(truth)) return(setNames(lapply(seq_along(truth), function(i) {
    region <- if (is.list(bounds)) bounds[[if (is.null(names(bounds))) i else names(truth)[i]]] else bounds
    evaluation_interval(region, truth[[i]])
  }), names(truth)))
  # The last interval dimension is lower/upper, for scalars, curves and surfaces.
  stopifnot(length(bounds) == 2L * length(truth))
  bounds <- matrix(bounds, ncol = 2L)
  list(covered = evaluation_shape(bounds[, 1L] <= as.numeric(truth) &
                                   as.numeric(truth) <= bounds[, 2L], truth),
       width = evaluation_shape(bounds[, 2L] - bounds[, 1L], truth))
}

# One cell supplies location-wise errors/coverage/width and per-replicate
# training MAD/RMSE. Pass fitted to reuse extraction; otherwise ... is passed
# to results_on_train(), e.g. B_M3=100. Its grid must match this fixed grid.
evaluate_replicate <- function(s, fitted = NULL, grid_N = 50L, ...) {
  grid <- seq(-.5, .5, length.out = grid_N)
  if (is.null(fitted)) fitted <- results_on_train(s, xval_N = grid_N, ...)
  truth <- evaluation_truth(s$scenario, grid)
  effects <- c("f1", "f2", "fcont", "fcont_full", "f", "f_full")
  target <- c(truth$parameters, truth$grid[effects])
  models <- setNames(vector("list", length(fitted)), names(fitted))
  interval_probs <- NULL
  for (model in names(fitted)) {
    fit <- fitted[[model]]
    if (is.null(fit)) next
    pe <- fit$point_est; pg <- fit$point_est_grid
    if (is.null(pg$fcont_full)) stop("Re-extract fits to obtain fcont_full and scalar slopes.")
    stopifnot(identical(fit$grid, grid))
    if (!is.null(interval_probs)) stopifnot(identical(fit$interval_probs, interval_probs))
    interval_probs <- fit$interval_probs
    estimate <- c(list(int = pe$int, slope1 = pe$slope1, slope2 = pe$slope2,
                       spat = pg$spat), pg[effects])
    estimate <- estimate[!vapply(estimate, is.null, logical(1))]
    errors <- setNames(lapply(names(estimate), function(eff)
      evaluation_error(estimate[[eff]], target[[eff]])), names(estimate))
    interval_types <- intersect(c("CI_quant", "CI_hpd"), names(fit))
    intervals <- setNames(lapply(interval_types, function(type)
      setNames(lapply(names(estimate), function(eff)
        evaluation_interval(fit[[type]][[eff]], target[[eff]])), names(estimate))), interval_types)

    accuracy <- function(error) c(MAD = mean(abs(error)), RMSE = sqrt(mean(error^2)))
    training <- setNames(lapply(effects, function(eff) {
      error <- evaluation_error(pe[[eff]], truth$training[[eff]])
      if (s$scenario$predictor == "X3" && eff %in% c("f", "f_full")) {
        # Region r uses ONLY its own observed covariate pairs, not all n pairs.
        return(lapply(split(error, s$scenario$dat$region), accuracy))
      }
      accuracy(error)
    }), effects)
    models[[model]] <- list(errors = errors, intervals = intervals, training = training)
  }
  list(scenario = s$scenario[c("dgp", "predictor", "n", "replicate", "seed")],
       grid = grid, interval_probs = interval_probs, settings = s$settings,
       region_counts = table(s$scenario$dat$region), models = models)
}

# Aggregate ONE (DGP, predictor, n) setting across any number of replicates.
# Failed fits retain NA training rows and are counted. Missing intervals (e.g.
# disabled M3 bootstrap) have zero valid counts, rather than false noncoverage.
aggregate_evaluations <- function(evaluations) {
  first <- evaluations[[1L]]
  ids <- vapply(evaluations, function(z) as.numeric(z$scenario$replicate), numeric(1))
  setting <- first$scenario[c("dgp", "predictor", "n")]
  known_probs <- Filter(Negate(is.null), lapply(evaluations, `[[`, "interval_probs"))
  if (length(known_probs)) stopifnot(all(vapply(known_probs, identical, logical(1), known_probs[[1L]])))
  stopifnot(!anyDuplicated(ids), all(vapply(evaluations, function(z)
    identical(unlist(z$scenario[c("dgp", "predictor", "n")]), unlist(setting)) &&
      identical(z$grid, first$grid) && identical(names(z$models), names(first$models)), logical(1))))
  reduce_tree <- function(nodes, summarize, leaf = function(x) !is.list(x)) {
    if (leaf(nodes[[1L]])) return(summarize(nodes))
    keys <- names(nodes[[1L]])
    setNames(lapply(keys, function(key)
      reduce_tree(lapply(nodes, `[[`, key), summarize, leaf)), keys)
  }
  # Columns are replicates. Each row is one scalar, grid location or region.
  collect <- function(nodes) matrix(vapply(nodes, as.numeric, numeric(length(nodes[[1L]]))),
                                     nrow = length(nodes[[1L]]))
  divide <- function(total, n) ifelse(n > 0, total / n, NA_real_)
  point_summary <- function(nodes) {
    values <- collect(nodes)
    valid <- is.finite(values); n <- rowSums(valid)
    values[!valid] <- 0
    shape <- function(x) evaluation_shape(x, nodes[[1L]])
    list(bias = shape(divide(rowSums(values), n)), MAD = shape(divide(rowSums(abs(values)), n)),
         RMSE = shape(sqrt(divide(rowSums(values^2), n))), n = shape(n))
  }
  interval_summary <- function(nodes) {
    covered <- collect(lapply(nodes, `[[`, "covered"))
    width <- collect(lapply(nodes, `[[`, "width"))
    valid <- is.finite(covered) & is.finite(width)
    n <- rowSums(valid); covered[!valid] <- width[!valid] <- 0
    coverage <- divide(rowSums(covered), n)
    shape <- function(x) evaluation_shape(x, nodes[[1L]]$width)
    list(coverage = shape(coverage), width = shape(divide(rowSums(width), n)),
         mcse = shape(sqrt(coverage * (1 - coverage) / n)), n = shape(n))
  }
  models <- setNames(lapply(names(first$models), function(model) {
    nodes <- lapply(evaluations, function(z) z$models[[model]])
    available <- which(!vapply(nodes, is.null, logical(1)))
    if (!length(available)) return(list(n_fits = 0L, n_failed = length(ids)))
    nodes <- nodes[available]
    training_summary <- function(records) {
      out <- data.frame(replicate = ids, MAD = NA_real_, RMSE = NA_real_)
      out[available, c("MAD", "RMSE")] <- t(vapply(records, as.numeric, numeric(2)))
      out
    }
    list(n_fits = length(available), n_failed = length(ids) - length(available),
         point = reduce_tree(lapply(nodes, `[[`, "errors"), point_summary),
         uncertainty = reduce_tree(lapply(nodes, `[[`, "intervals"), interval_summary,
           leaf = function(x) is.list(x) && identical(names(x), c("covered", "width"))),
         training = reduce_tree(lapply(nodes, `[[`, "training"), training_summary))
  }), names(first$models))
  list(scenario = setting, grid = first$grid, replicates = ids,
       interval_probs = if (length(known_probs)) known_probs[[1L]] else NULL, models = models)
}
