#-----------------------------------------------------------------------#
# Package: SAM                                                          #
# Method: Sparse Additive Modelling using Quadratic Loss                #
#-----------------------------------------------------------------------#

#' Training function of Sparse Additive Regression with Quadratic Loss
#'
#' Fit a sparse additive regression model with quadratic loss.
#'
#' The solver combines block coordinate descent, fast iterative
#' soft-thresholding, and Newton updates. Computation is accelerated by warm
#' starts and active-set screening.
#'
#' @param X Numeric training matrix with \code{n} rows (samples) and
#'   \code{d} columns (features).
#' @param y Numeric response vector of length \code{n}.
#' @param p The number of basis spline functions. The default value is 3.
#' @param lambda Optional user-supplied regularization sequence. If provided,
#'   use a decreasing sequence; warm starts are used along the path and are
#'   usually much faster than fitting a single value.
#' @param nlambda The number of lambda values. The default value is 30.
#' @param lambda.min.ratio Smallest lambda as a fraction of \code{lambda.max}
#'   (the smallest value that keeps all component functions at zero). The
#'   default is \code{5e-3}.
#' @param thol Stopping tolerance. The default value is \code{1e-5}.
#' @param max.ite Maximum number of iterations. The default value is \code{1e5}.
#' @param regfunc A string indicating the regularizer. The default value is "L1". You can also assign "MCP" or "SCAD" to it.
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
#'   Solution path matrix with size \code{d*p} by \code{length(lambda)}; each
#'   column corresponds to one regularization parameter.
#' }
#' \item{intercept}{
#'   The solution path of the intercept.
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
#' \item{sse}{
#'   Sums of square errors of the solution path.
#' }
#' @seealso \code{\link{SAM}},\code{\link{plot.samQL},\link{print.samQL},\link{predict.samQL}}
#' @examples
#'
#' ## generating training data
#' n = 100
#' d = 500
#' X = 0.5*matrix(runif(n*d),n,d) + matrix(rep(0.5*runif(n),d),n,d)
#'
#' ## generating response
#' y = -2*sin(X[,1]) + X[,2]^2-1/3 + X[,3]-1/2 + exp(-X[,4])+exp(-1)-1
#'
#' ## Training
#' out.trn = samQL(X,y)
#' out.trn
#'
#' ## plotting solution path
#' plot(out.trn)
#'
#' ## generating testing data
#' nt = 1000
#' Xt = 0.5*matrix(runif(nt*d),nt,d) + matrix(rep(0.5*runif(nt),d),nt,d)
#'
#' yt = -2*sin(Xt[,1]) + Xt[,2]^2-1/3 + Xt[,3]-1/2 + exp(-Xt[,4])+exp(-1)-1
#'
#' ## predicting response
#' out.tst = predict(out.trn,Xt)
#' @export
samQL = function(X, y, p=3, lambda = NULL, nlambda = NULL, lambda.min.ratio = 5e-3, thol=1e-5, max.ite = 1e5, regfunc="L1"){

  p = .sam_validate_p(p)
  checked = .sam_validate_xy(X, y, "samQL")
  X = checked$X
  y = checked$y
  n = checked$n
  d = checked$d
  m = d * p

  lambda.info = .sam_validate_lambda(lambda, nlambda, lambda.min.ratio, default_nlambda = 30)
  if (lambda.info$supplied) {
    lambda = lambda.info$lambda
    nlambda = lambda.info$nlambda
  } else {
    nlambda = lambda.info$nlambda
  }

  scaled = .sam_scale_training(X)
  X = scaled$X

  fit = list(p = p, X.min = scaled$X.min, X.ran = scaled$X.ran)

  basis = .sam_build_basis(X, p)
  Z = basis$Z
  fit$knots = basis$knots
  fit$nkots = basis$knots
  fit$Boundary.knots = basis$Boundary.knots

  Z.mean = apply(Z,2,mean)
  Z = sweep(Z, 2, Z.mean, "-")
  y.mean = mean(y)
  y = y - y.mean

  lambda_input = 1
  if(!lambda.info$supplied){
    lambda_input = 0
    lambda = exp(seq(log(1),log(lambda.min.ratio),length = nlambda))
  }


  out = .C("grplasso",y = as.double(y), X = as.double(Z), lambda = as.double(lambda), nnlambda = as.integer(nlambda), nn = as.integer(n), dd = as.integer(d), pp = as.integer(p), ww = as.double(matrix(0,m,nlambda)), mmax_ite = as.integer(max.ite), tthol = as.double(thol), regfunc = as.character(regfunc), iinput = as.integer(lambda_input), df=as.integer(rep(0,nlambda)), sse=as.double(rep(0,nlambda)), func_norm = as.double(matrix(0,d,nlambda)), PACKAGE="SAM")

  fit$lambda = out$lambda
  fit$w = matrix(out$w,ncol=nlambda)
  fit$df = out$df
  fit$sse = out$sse
  fit$func_norm = matrix(out$func_norm,ncol=nlambda)
  fit$intercept = rep(y.mean,nlambda) - t(Z.mean)%*%fit$w
  fit$XX = out$X

  class(fit) = "samQL"
  return(fit)
}

#' Printing function for S3 class \code{"samQL"}
#'
#' Print a summary of an object of class \code{"samQL"}.
#'
#' The output includes the regularization path length and its degrees of
#' freedom.
#'
#' @param x An object with S3 class \code{"samQL"}
#' @param \dots Additional arguments passed to methods; currently unused.
#' @seealso \code{\link{samQL}}
#' @export
print.samQL = function(x,...){
  cat("Path length:",length(x$df),"\n")
  cat("d.f.:",x$df[1],"--->",x$df[length(x$df)],"\n")
}

#' Plot function for S3 class \code{"samQL"}
#'
#' Plot the regularization path (regularization parameter versus functional
#' norm).
#'
#' The x-axis shows regularization parameters on a log scale. The y-axis shows
#' the functional norm of each component function.
#'
#' @param x An object with S3 class \code{"samQL"}
#' @param \dots Additional arguments passed to methods; currently unused.
#' @seealso \code{\link{samQL}}
#' @export
plot.samQL = function(x,...){
  par = par(omi = c(0.0, 0.0, 0, 0), mai = c(1, 1, 0.1, 0.1))
  matplot(x$lambda,t(x$func_norm),type="l",xlab="Regularization Parameters",ylab = "Functional Norms",cex.lab=2,log="x",lwd=2)
}

#' Prediction function for S3 class \code{"samQL"}
#'
#' Predict responses for test data.
#'
#' The test matrix is rescaled using the training \code{X.min}/\code{X.ran},
#' truncated to \code{[0, 1]}, and expanded with the same spline basis used
#' during training.
#'
#' @param object An object with S3 class \code{"samQL"}.
#' @param newdata Numeric test matrix with \code{n} rows and \code{d} columns.
#' @param \dots Additional arguments passed to methods; currently unused.
#' @return
#' \item{values}{
#'   Predicted responses as an \code{n} by \code{length(lambda)} matrix.
#' }
#' @seealso \code{\link{samQL}}
#' @export
predict.samQL = function(object, newdata,...){
  newdata = .sam_scale_newdata(object, newdata)
  nt = nrow(newdata)
  Zt = .sam_build_basis_newdata(object, newdata)

  out = list()
  out$values = cbind(Zt, rep(1, nt)) %*% rbind(object$w, object$intercept)

  return(out)
}
