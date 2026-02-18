.sam_validate_p <- function(p) {
  if (length(p) != 1L || !is.numeric(p) || is.na(p) || p < 1 || p %% 1 != 0) {
    stop("`p` must be a positive integer.", call. = FALSE)
  }
  as.integer(p)
}

.sam_validate_nlambda <- function(nlambda, default_nlambda) {
  if (is.null(nlambda)) {
    return(as.integer(default_nlambda))
  }

  if (length(nlambda) != 1L || !is.numeric(nlambda) || is.na(nlambda) ||
      nlambda < 1 || nlambda %% 1 != 0) {
    stop("`nlambda` must be a positive integer.", call. = FALSE)
  }

  as.integer(nlambda)
}

.sam_validate_lambda <- function(lambda, nlambda, lambda.min.ratio, default_nlambda) {
  if (is.null(lambda)) {
    if (length(lambda.min.ratio) != 1L || !is.numeric(lambda.min.ratio) ||
        is.na(lambda.min.ratio) || lambda.min.ratio <= 0 || lambda.min.ratio > 1) {
      stop("`lambda.min.ratio` must be in (0, 1].", call. = FALSE)
    }

    return(list(
      lambda = NULL,
      nlambda = .sam_validate_nlambda(nlambda, default_nlambda),
      supplied = FALSE
    ))
  }

  lambda <- as.numeric(lambda)
  if (length(lambda) < 1L || any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("`lambda` must contain finite positive values.", call. = FALSE)
  }
  if (is.unsorted(-lambda)) {
    warning("`lambda` should be a decreasing sequence for efficient warm starts.", call. = FALSE)
  }

  list(lambda = lambda, nlambda = length(lambda), supplied = TRUE)
}

.sam_validate_xy <- function(X, y, fn_name) {
  X <- as.matrix(X)
  y <- as.vector(y)

  if (!is.numeric(X) || !is.numeric(y)) {
    stop(sprintf("`%s` requires numeric `X` and `y`.", fn_name), call. = FALSE)
  }
  if (nrow(X) < 1L || ncol(X) < 1L) {
    stop("`X` must contain at least one row and one column.", call. = FALSE)
  }
  if (length(y) != nrow(X)) {
    stop("`y` length must match `nrow(X)`.", call. = FALSE)
  }
  if (any(!is.finite(X)) || any(!is.finite(y))) {
    stop("`X` and `y` must not contain NA/NaN/Inf values.", call. = FALSE)
  }

  list(X = X, y = y, n = nrow(X), d = ncol(X))
}

.sam_validate_weights <- function(w, n) {
  if (is.null(w)) {
    return(rep(1, n))
  }

  w <- as.numeric(w)
  if (length(w) != n || any(!is.finite(w)) || any(w <= 0)) {
    stop("`w` must be a positive numeric vector with length `nrow(X)`.", call. = FALSE)
  }
  w
}

.sam_scale_training <- function(X) {
  X.min <- apply(X, 2, min)
  X.max <- apply(X, 2, max)
  X.ran <- X.max - X.min
  zero_range <- (X.ran == 0)
  if (any(zero_range)) {
    X.ran[zero_range] <- 1
  }

  X.scaled <- sweep(X, 2, X.min, "-")
  X.scaled <- sweep(X.scaled, 2, X.ran, "/")

  list(X = X.scaled, X.min = X.min, X.ran = X.ran)
}

.sam_scale_newdata <- function(object, newdata) {
  newdata <- as.matrix(newdata)
  if (!is.numeric(newdata)) {
    stop("`newdata` must be numeric.", call. = FALSE)
  }
  if (ncol(newdata) != length(object$X.min)) {
    stop("`newdata` column count must match training data.", call. = FALSE)
  }
  if (any(!is.finite(newdata))) {
    stop("`newdata` must not contain NA/NaN/Inf values.", call. = FALSE)
  }

  X.ran <- object$X.ran
  zero_range <- (X.ran == 0)
  if (any(zero_range)) {
    X.ran[zero_range] <- 1
  }

  scaled <- sweep(newdata, 2, object$X.min, "-")
  scaled <- sweep(scaled, 2, X.ran, "/")
  scaled <- pmax(scaled, 0)
  pmin(scaled, 1)
}

.sam_build_basis <- function(X, p) {
  n <- nrow(X)
  d <- ncol(X)
  m <- d * p

  Z <- matrix(0, n, m)
  knots <- if (p > 1L) matrix(0, p - 1L, d) else matrix(numeric(0), 0, d)
  boundary.knots <- matrix(0, 2, d)

  for (j in seq_len(d)) {
    idx <- (j - 1L) * p + seq_len(p)
    basis_j <- splines::ns(X[, j], df = p)
    Z[, idx] <- basis_j
    if (p > 1L) {
      knot_j <- attr(basis_j, "knots")
      if (!is.null(knot_j) && length(knot_j) == (p - 1L)) {
        knots[, j] <- knot_j
      }
    }
    boundary.knots[, j] <- attr(basis_j, "Boundary.knots")
  }

  list(Z = Z, knots = knots, Boundary.knots = boundary.knots)
}

.sam_resolve_knots <- function(object) {
  if (!is.null(object$knots)) {
    return(object$knots)
  }
  if (!is.null(object$nkots)) {
    return(object$nkots)
  }

  stop("Model object does not contain spline knots.", call. = FALSE)
}

.sam_build_basis_newdata <- function(object, newdata) {
  p <- object$p
  d <- ncol(newdata)
  m <- d * p
  nt <- nrow(newdata)

  knots <- .sam_resolve_knots(object)
  boundary.knots <- object$Boundary.knots

  if (is.null(boundary.knots) || ncol(boundary.knots) != d) {
    stop("Model object has incompatible boundary knots.", call. = FALSE)
  }
  if (p > 1L && (is.null(knots) || ncol(knots) != d)) {
    stop("Model object has incompatible spline knots.", call. = FALSE)
  }

  Zt <- matrix(0, nt, m)
  for (j in seq_len(d)) {
    idx <- (j - 1L) * p + seq_len(p)
    if (p > 1L) {
      Zt[, idx] <- splines::ns(
        newdata[, j],
        df = p,
        knots = knots[, j],
        Boundary.knots = boundary.knots[, j]
      )
    } else {
      Zt[, idx] <- splines::ns(
        newdata[, j],
        df = p,
        Boundary.knots = boundary.knots[, j]
      )
    }
  }

  Zt
}
