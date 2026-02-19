#-----------------------------------------------------------------------#
# Package: SAM                                                          #
# Method: Sparse Additive Modelling using Hinge Loss                    #
#-----------------------------------------------------------------------#

#' Training function of Sparse Additive Hinge-Loss Classifier
#'
#' Fit a sparse additive classifier with hinge loss.
#'
#' The solver combines block coordinate descent, fast iterative soft-thresholding,
#' and Newton updates. Computation is accelerated by warm starts and active-set
#' screening.
#'
#' @param X Numeric training matrix with \code{n} rows (samples) and \code{d}
#'   columns (features).
#' @param y Training labels of length \code{n}. Labels must be coded as
#'   \code{-1} and \code{1}.
#' @param p The number of basis spline functions. The default value is 3.
#' @param lambda Optional user-supplied regularization sequence. If provided,
#'   use a decreasing sequence; warm starts are used along the path and are
#'   usually much faster than fitting a single value.
#' @param nlambda The number of lambda values. The default value is 20.
#' @param lambda.min.ratio Smallest lambda as a fraction of \code{lambda.max}
#'   (the smallest value that keeps all component functions at zero). The
#'   default is \code{0.4}.
#' @param thol Stopping tolerance. The default value is \code{1e-5}.
#' @param mu Smoothing parameter used to approximate hinge loss. The default
#'   value is \code{0.05}.
#' @param max.ite Maximum number of iterations. The default value is \code{1e5}.
#' @param w Optional positive observation weights of length \code{n}. The
#'   default is \code{1} for all observations.
#' @return
#' \item{p}{
#'   The number of basis spline functions used in training.
#' }
#' \item{X.min}{
#'   Per-feature minimums from training data (used to rescale test data).
#' }
#' \item{X.ran}{
#'   Per-feature ranges from training data (used to rescale test data).
#' }
#' \item{lambda}{
#'   Sequence of regularization parameters used in training.
#' }
#' \item{w}{
#'   Solution path matrix with size \code{d*p+1} by \code{length(lambda)}; each column corresponds to one regularization parameter.
#' }
#' \item{df}{
#'   Degrees of freedom along the solution path (number of non-zero component
#'   functions).
#' }
#' \item{knots}{
#'   The \code{p-1} by \code{d} matrix. Each column contains the knots applied to the corresponding variable.
#' }
#' \item{Boundary.knots}{
#'   The \code{2} by \code{d} matrix. Each column contains the boundary points applied to the corresponding variable.
#' }
#' \item{func_norm}{
#'   Functional norm matrix (\code{d} by \code{length(lambda)}); each column
#'   corresponds to one regularization parameter.
#' }
#' @seealso \code{\link{SAM}},\code{\link{plot.samHL},\link{print.samHL},\link{predict.samHL}}
#' @examples
#'
#' ## generating training data
#' n = 200
#' d = 100
#' X = 0.5*matrix(runif(n*d),n,d) + matrix(rep(0.5*runif(n),d),n,d)
#' y = sign(((X[,1]-0.5)^2 + (X[,2]-0.5)^2)-0.06)
#'
#' ## flipping about 5 percent of y
#' y = y*sign(runif(n)-0.05)
#'
#' ## Training
#' out.trn = samHL(X,y)
#' out.trn
#'
#' ## plotting solution path
#' plot(out.trn)
#'
#' ## generating testing data
#' nt = 1000
#' Xt = 0.5*matrix(runif(nt*d),nt,d) + matrix(rep(0.5*runif(nt),d),nt,d)
#'
#' yt = sign(((Xt[,1]-0.5)^2 + (Xt[,2]-0.5)^2)-0.06)
#'
#' ## flipping about 5 percent of y
#' yt = yt*sign(runif(nt)-0.05)
#'
#' ## predicting response
#' out.tst = predict(out.trn,Xt)
#' @export
samHL = function(X, y, p=3, lambda = NULL, nlambda = NULL, lambda.min.ratio = 0.4, thol=1e-5, mu = 5e-2, max.ite = 1e5, w = NULL){

  p = .sam_validate_p(p)
  checked = .sam_validate_xy(X, y, "samHL")
  X = checked$X
  y = checked$y
  n = checked$n
  d = checked$d
  m = d * p

  w = .sam_validate_weights(w, n)

  np = sum(y==1)
  nn = sum(y==-1)

  if (np == 0 || nn == 0) {
    stop("Please provide both 1 and -1 labels.", call. = FALSE)
  }

  if((np+nn)!=n){
    stop("Please check the labels. (Must be coded in 1 and -1)", call. = FALSE)
  }

  lambda.info = .sam_validate_lambda(lambda, nlambda, lambda.min.ratio, default_nlambda = 20)
  if (lambda.info$supplied) {
    lambda = lambda.info$lambda
    nlambda = lambda.info$nlambda
  } else {
    nlambda = lambda.info$nlambda
  }

  if(np>nn) a0 = 1-nn/np*mu else a0 = np/nn*mu - 1

  scaled = .sam_scale_training(X)
  X = scaled$X
  fit = list(p = p, X.min = scaled$X.min, X.ran = scaled$X.ran)

  basis = .sam_build_basis(X, p)
  Z = basis$Z
  fit$knots = basis$knots
  fit$nkots = basis$knots
  fit$Boundary.knots = basis$Boundary.knots

  Z = cbind(matrix(rep(y,m),n,m)*Z,y)

  if(!lambda.info$supplied){
    u = cbind((rep(1,n) - a0*y)/mu,rep(0,n),rep(1,n))
    u = apply(u,1,median)

    lambda_max = max(sqrt(colSums(matrix(t(Z[,1:(p*d)])%*%u,p,d)^2)))
    lambda = exp(seq(log(1),log(lambda.min.ratio),length=nlambda))*lambda_max
  }

  L0 = norm(Z,'f')^2/mu

  out = .C("grpSVM", Z = as.double(Z), w = as.double(w), lambda = as.double(lambda), nnlambda = as.integer(nlambda), LL0 = as.double(L0), nn = as.integer(n), dd = as.integer(d), pp = as.integer(p),aa0 = as.double(a0), xx = as.double(matrix(0,m+1,nlambda)), mmu = as.double(mu), mmax_ite = as.integer(max.ite), tthol = as.double(thol),aalpha = as.double(0.5),df=as.double(rep(0,nlambda)),func_norm=as.double(matrix(0,d,nlambda)),PACKAGE="SAM")

  fit$lambda = out$lambda
  fit$w = matrix(out$xx,ncol=nlambda)
  fit$df = out$df
  fit$func_norm = matrix(out$func_norm,ncol=nlambda)

  class(fit) = "samHL"
  return(fit)
}

#' Printing function for S3 class \code{"samHL"}
#'
#' Print a summary of an object of class \code{"samHL"}.
#'
#' The output includes the regularization path length and its degrees of
#' freedom.
#'
#' @param x An object with S3 class \code{"samHL"}
#' @param \dots Additional arguments passed to methods; currently unused.
#' @seealso \code{\link{samHL}}
#' @export
print.samHL = function(x,...){
  cat("Path length:",length(x$df),"\n")
  cat("d.f.:",x$df[1],"--->",x$df[length(x$df)],"\n")
}

#' Plot function for S3 class \code{"samHL"}
#'
#' Plot the regularization path (regularization parameter versus functional
#' norm).
#'
#' The x-axis shows regularization parameters on a log scale. The y-axis shows
#' the functional norm of each component function.
#'
#' @param x An object with S3 class \code{"samHL"}
#' @param \dots Additional arguments passed to methods; currently unused.
#' @seealso \code{\link{samHL}}
#' @export
plot.samHL = function(x,...){
  par = par(omi = c(0.0, 0.0, 0, 0), mai = c(1, 1, 0.1, 0.1))
  matplot(x$lambda,t(x$func_norm),type="l",xlab="Regularization Parameters",ylab = "Functional Norms",cex.lab=2,log="x",lwd=2)
}

#' Prediction function for S3 class \code{"samHL"}
#'
#' Predict decision values and class labels for test data.
#'
#' The test matrix is rescaled using the training \code{X.min}/\code{X.ran},
#' truncated to \code{[0, 1]}, and expanded with the same spline basis used
#' during training.
#'
#' @param object An object with S3 class \code{"samHL"}.
#' @param newdata Numeric test matrix with \code{n} rows and \code{d} columns.
#' @param thol Decision-value threshold used to convert scores to labels. The
#'   default value is \code{0}.
#' @param \dots Additional arguments passed to methods; currently unused.
#' @return
#' \item{values}{
#'   Predicted decision values as an \code{n} by \code{length(lambda)} matrix.
#' }
#' \item{labels}{
#'   Predicted class labels (\code{-1}/\code{1}) as an \code{n} by
#'   \code{length(lambda)} matrix.
#' }
#' @seealso \code{\link{samHL}}
#' @export
predict.samHL = function(object, newdata, thol = 0, ...){
  newdata = .sam_scale_newdata(object, newdata)
  nt = nrow(newdata)
  Zt = .sam_build_basis_newdata(object, newdata)

  out = list()
  out$values = cbind(Zt, rep(1, nt)) %*% object$w
  out$labels = (out$values > thol) * 2 - 1

  return(out)
}
