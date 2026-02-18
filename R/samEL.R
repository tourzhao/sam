#-----------------------------------------------------------------------#
# Package: SAM                                                          #
# Method: Sparse Additive Modelling using Exponential Loss              #
#-----------------------------------------------------------------------#

#' Training function of Sparse Additive Possion Regression
#'
#' The log-linear model is learned using training data.
#'
#' We adopt various computational algorithms including the block coordinate descent, fast iterative soft-thresholding algorithm, and newton method. The computation is further accelerated by "warm-start" and "active-set" tricks.
#'
#' @param X The \code{n} by \code{d} design matrix of the training set, where \code{n} is sample size and \code{d} is dimension.
#' @param y The \code{n}-dimensional response vector of the training set, where \code{n} is sample size. Responses must be non-negative integers.
#' @param p The number of basis spline functions. The default value is 3.
#' @param lambda A user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on nlambda and lambda.min.ratio. Supplying a value of lambda overrides this. WARNING: use with care. Do not supply a single value for lambda. Supply instead a decreasing sequence of lambda values. samEL relies on its warms starts for speed, and its often faster to fit a whole path than compute a single fit.
#' @param nlambda The number of lambda values. The default value is 20.
#' @param lambda.min.ratio Smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero). The default is 0.1.
#' @param thol Stopping precision. The default value is 1e-5.
#' @param max.ite The number of maximum iterations. The default value is 1e5.
#' @param regfunc A string indicating the regularizer. The default value is "L1". You can also assign "MCP" or "SCAD" to it.
#' @return
#' \item{p}{
#'   The number of basis spline functions used in training.
#' }
#' \item{X.min}{
#'   A vector with each entry corresponding to the minimum of each input variable. (Used for rescaling in testing)
#' }
#' \item{X.ran}{
#'   A vector with each entry corresponding to the range of each input variable. (Used for rescaling in testing)
#' }
#' \item{lambda}{
#'   A sequence of regularization parameter used in training.
#' }
#' \item{w}{
#'   The solution path matrix (\code{d*p+1} by length of \code{lambda}) with each column corresponding to a regularization parameter. Since we use the basis expansion with the intercept, the length of each column is \code{d*p+1}.
#' }
#' \item{df}{
#'   The degree of freedom of the solution path (The number of non-zero component function)
#' }
#' \item{knots}{
#'   The \code{p-1} by \code{d} matrix. Each column contains the knots applied to the corresponding variable.
#' }
#' \item{Boundary.knots}{
#'   The \code{2} by \code{d} matrix. Each column contains the boundary points applied to the corresponding variable.
#' }
#' \item{func_norm}{
#'   The functional norm matrix (\code{d} by length of \code{lambda}) with each column corresponds to a regularization parameter. Since we have \code{d} input variables, the length of each column is \code{d}.
#' }
#' @seealso \code{\link{SAM}},\code{\link{plot.samEL},\link{print.samEL},\link{predict.samEL}}
#' @examples
#'
#' ## generating training data
#' n = 200
#' d = 100
#' X = 0.5*matrix(runif(n*d),n,d) + matrix(rep(0.5*runif(n),d),n,d)
#' u = exp(-2*sin(X[,1]) + X[,2]^2-1/3 + X[,3]-1/2 + exp(-X[,4])+exp(-1)-1+1)
#' y = rep(0,n)
#' for(i in 1:n) y[i] = rpois(1,u[i])
#'
#' ## Training
#' out.trn = samEL(X,y)
#' out.trn
#'
#' ## plotting solution path
#' plot(out.trn)
#'
#' ## generating testing data
#' nt = 1000
#' Xt = 0.5*matrix(runif(nt*d),nt,d) + matrix(rep(0.5*runif(nt),d),nt,d)
#' ut = exp(-2*sin(Xt[,1]) + Xt[,2]^2-1/3 + Xt[,3]-1/2 + exp(-Xt[,4])+exp(-1)-1+1)
#' yt = rep(0,nt)
#' for(i in 1:nt) yt[i] = rpois(1,ut[i])
#'
#' ## predicting response
#' out.tst = predict(out.trn,Xt)
#' @export
samEL = function(X, y, p=3, lambda = NULL, nlambda = NULL, lambda.min.ratio = 0.25, thol=1e-5, max.ite = 1e5, regfunc="L1"){

  p = .sam_validate_p(p)
  checked = .sam_validate_xy(X, y, "samEL")
  X = checked$X
  y = checked$y
  n = checked$n
  d = checked$d
  m = d * p

  if(sum(y<0)>0){
    stop("Please check the responses. (Must be non-negative)", call. = FALSE)
  }
  if (mean(y) <= 0) {
    stop("At least one positive response is required.", call. = FALSE)
  }

  lambda.info = .sam_validate_lambda(lambda, nlambda, lambda.min.ratio, default_nlambda = 20)
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

  L0 = norm(Z,"f")^2

  #if(n>d) L0 = max(eigen(Z%*%t(Z), symmetric=T, only.values = T)$values) else L0 = max(eigen(t(Z)%*%Z, symmetric=T, only.values = T)$values)


  Z = cbind(Z,rep(1,n))

  a0 = log(mean(y))

  z = colSums(matrix(rep(y,m+1),n,m+1)*Z)

  if(!lambda.info$supplied){
    g = -z + colSums(matrix(rep(exp(a0),m+1),n,m+1)*Z)

    lambda_max=max(sqrt(colSums(matrix(g[1:(p*d)],p,d)^2)))
    lambda = exp(seq(log(1),log(lambda.min.ratio),length=nlambda))*lambda_max
  }

  out = .C("grpPR", A = as.double(Z), y = as.double(y), lambda = as.double(lambda), nlambda = as.integer(nlambda),
           LL0 = as.double(L0), nn = as.integer(n), dd = as.integer(d), pp = as.integer(p),
           xx = as.double(matrix(0,m+1,nlambda)), aa0 = as.double(a0), mmax_ite = as.integer(max.ite),
           tthol = as.double(thol), regfunc = as.character(regfunc), aalpha = as.double(0.5),
           z = as.double(z),df = as.integer(rep(0,nlambda)),func_norm = as.double(matrix(0,d,nlambda)), PACKAGE="SAM")

  fit$lambda = out$lambda
  fit$w = matrix(out$xx,ncol=nlambda)
  fit$df = out$df
  fit$func_norm = matrix(out$func_norm,ncol=nlambda)


  class(fit) = "samEL"
  return(fit)
}

#' Printing function for S3 class \code{"samEL"}
#'
#' Summarize the information of the object with S3 class \code{samEL}.
#'
#' The output includes length and d.f. of the regularization path.
#'
#' @param x An object with S3 class \code{"samEL"}
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{samEL}}
#' @export
print.samEL = function(x,...){
  cat("Path length:",length(x$df),"\n")
  cat("d.f.:",x$df[1],"--->",x$df[length(x$df)],"\n")
}

#' Plot function for S3 class \code{"samEL"}
#'
#' This function plots the regularization path (regularization parameter versus functional norm)
#'
#' The horizontal axis is for the regularization parameters in log scale. The vertical axis is for the functional norm of each component.
#'
#' @param x An object with S3 class \code{"samEL"}
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{samEL}}
#' @export
plot.samEL = function(x,...){
  par = par(omi = c(0.0, 0.0, 0, 0), mai = c(1, 1, 0.1, 0.1))
  matplot(x$lambda,t(x$func_norm),type="l",xlab="Regularization Parameters",ylab = "Functional Norms",cex.lab=2,log="x",lwd=2)
}

#' Prediction function for S3 class \code{"samEL"}
#'
#' Predict the labels for testing data.
#'
#' The testing dataset is rescale to the samELe range, and expanded by the samELe spline basis functions as the training data.
#'
#' @param object An object with S3 class \code{"samEL"}.
#' @param newdata The testing dataset represented in a \code{n} by \code{d} matrix, where \code{n} is testing sample size and \code{d} is dimension.
#' @param \dots System reserved (No specific usage)
#' @return
#'   \item{expectations}{
#'     Estimated expected counts also represented in a \code{n} by the length of \code{lambda} matrix, where \code{n} is testing sample size.
#'   }
#' @seealso \code{\link{samEL}}
#' @export
predict.samEL = function(object, newdata,...){
  newdata = .sam_scale_newdata(object, newdata)
  nt = nrow(newdata)
  Zt = .sam_build_basis_newdata(object, newdata)

  out = list()
  out$expectation = exp(cbind(Zt, rep(1, nt)) %*% object$w)
  out$expectations = out$expectation

  return(out)
}
