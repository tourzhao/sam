#' @useDynLib SAM, .registration = TRUE
#' @importFrom splines ns
#' @importFrom stats median
#' @importFrom graphics matplot
NULL


#' Sparse Additive Modelling
#'
#' SAM provides sparse additive models for high-dimensional prediction tasks
#' (regression and classification). It uses spline basis expansion and efficient
#' optimization routines to compute full regularization paths.
#'
#' The package exposes four model families:
#' \itemize{
#' \item \code{\link{samQL}}: quadratic-loss sparse additive regression.
#' \item \code{\link{samLL}}: logistic-loss sparse additive classification.
#' \item \code{\link{samHL}}: hinge-loss sparse additive classification.
#' \item \code{\link{samEL}}: Poisson-loss sparse additive regression.
#' }
#' All models share a common spline representation and return regularization
#' paths, allowing model selection after one fit.
#' @docType package
#' @author Tuo Zhao, Xingguo Li, Haoming Jiang, Han Liu, and Kathryn Roeder\cr
#' Maintainer: Tuo Zhao <tourzhao@gatech.edu>
#' @references 
#' P. Ravikumar, J. Lafferty, H.Liu and L. Wasserman. "Sparse Additive Models", \emph{Journal of Royal Statistical Society: Series B}, 2009.\cr
#' T. Zhao and H.Liu. "Sparse Additive Machine", \emph{International Conference on Artificial Intelligence and Statistics}, 2012.\cr
#' @seealso \code{\link{samQL}},\code{\link{samHL}},\code{\link{samLL}},\code{\link{samEL}}
"_PACKAGE"
#> [1] "_PACKAGE"
