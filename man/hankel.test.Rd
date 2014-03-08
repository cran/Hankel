\name{hankel.test}
\alias{hankel.test}
\alias{print.hankeltest}
\title{Hankel-Test for the univariate two-sample problem.}
\description{

Computes the Hankel test statistic and its empirically standardized version in the case that \eqn{\nu} equals the exponential distribution (see details below). Optionally, simple computations of the p-value and the critical value for a selectable significance level and number of replicates are provided. 

}
\usage{
hankel.test(x,y,standardized=FALSE,replicates=500,
calcCritVal=FALSE,calcPVal=TRUE,sigLevel=0.95,
probMeasure="exp",params=list(lambda=1.0))
}
\arguments{
\item{x}{
  First set of observations as vector.
}

\item{y}{
  Second set of observations as vector.
}

\item{standardized}{
  Boolean variable which indicates whether the empirically standardized version of the test statistic is supposed to be applied. 
  
  The default is \code{standardized=FALSE}.
}

\item{probMeasure}{
  The family of probability measures which is supposed to be applied. However, only the exponential distribution has been implemented so far. 
}

\item{params}{
  List of parameter which specify the probability measure \eqn{\nu} which is applied (e.g. the rate parameter \eqn{\lambda} of the exponential distribution). The default is \code{list(lambda=1.0)}.
}

\item{sigLevel}{
  Significance level of the test. The default is \code{sigLevel=0.95}.
}

\item{calcCritVal}{
  Boolean variable which indicates whether the critical value of the test statistic for the given significance level is supposed to be computed. 
The default is \code{calcCritVal=FALSE}.
}


\item{calcPVal}{
  Boolean variable which indicates whether the p-value of the statistic corresponding to the given observations is supposed to be computed. The default is \code{calcPVal=TRUE}.
}


\item{replicates}{
  Number of bootstrap replicates for the computation of the critical value and the p-value. The default is \code{replicates=500}.
}



}
\details{
The Hankel test statistic is of the Cramer- von Mises type and was proposed by L. Baringhaus, University of Hanover. It is given by
\deqn{
T_{m,n}=\frac{mn}{m+n} \int_{R_{\ge 0}}(\frac{1}{m} \sum_{i=1}^m J_0(2\sqrt{tX_i})-\frac{1}{n} \sum_{i=1}^n J_0(2\sqrt{tY_i}))^2 \,d\nu(t),
}
where \eqn{X_1,\dots,X_m,Y_1,\dots,Y_n} are real nonnegative random variables, \eqn{ J_0} the Bessel function of the first kind of order zero and \eqn{\nu} some suitable probability measure on the nonnegative half-line. 

In case that \eqn{\nu} equals the exponential distribution with parameter \eqn{\lambda >0}, an alternative expression for \eqn{T_{m,n}}, which is obtained by applying formula 6.615 
of Gradstein and Ryshik (1981), is used:
\deqn{
T_{m,n}(\lambda)=\frac{mn}{m+n} \int_{0}^{\infty} \left (\frac{1}{m} \sum_{i=1}^m J_0(2\sqrt{tX_i})-\frac{1}{n} \sum_{j=1}^n J_0(2\sqrt{tY_j})\right )^2 \lambda \exp(-\lambda t)\,dt
}
\deqn{
= \frac{mn}{m+n} \bigg [\frac{1}{m^2} \sum_{i=1}^m \sum_{j=1}^m I_0(2\sqrt{X_iX_j}/\lambda)\hspace{0.5mm} \exp(-(X_i+X_j)/\lambda) 
}
\deqn{
+  \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n I_0(2\sqrt{Y_iY_j}/\lambda)\hspace{0.5mm} \exp(-(Y_i+Y_j)/\lambda) 
}
\deqn{
- \frac{2}{mn} \sum_{i=1}^m \sum_{j=1}^n I_0(2\sqrt{X_iY_j}/\lambda)\hspace{0.5mm} \exp(-(X_i+Y_j)/\lambda)\bigg ],
}
where \eqn{I_0} denotes the Bessel function of the third kind of order 0. This representation is used for the computation of the test statistic.

The empirically standardized version of the test is obtained by replacing the  \eqn{X_1,\dots,X_m} and \eqn{Y_1,\dots,Y_m} by the empirically standardized variables \eqn{U_1=X_1/\eta_{m,n},\ldots,U_m=X_m/\eta_{m,n}} and \eqn{V_1=Y_1/\eta_{m,n},\ldots,V_n=Y_n/\eta_{m,n}}, where \eqn{\eta_{m,n}=\frac{1}{m+n}\left ( \sum_{i=1}^{m}X_i+ \sum_{j=1}^{n}Y_j\right)}, in the representation of \eqn{T_{m,n}(\lambda)}. For further details see Baringhaus and Kolbe (2014).
 
}
\value{
An object of class \code{"hankeltest"} is returned. The following components are retrievable:

\item{samplesizes}{
A list containing the samplesizes.
}
\item{statistic}{
Value of the Hankel test statistic for the given observations.
}
\item{standardized}{
Boolean value which indicates whether the standardized version of the test was applied.
}
\item{probMeasure}{
String value which indicates which family of probability measures was applied.
}
\item{params}{
List of parameter which specify the applied probability measure \eqn{\nu}. 
}
\item{sigLevel}{
Significance level of the test. 
}
\item{pValue}{
Boostrap estimation of the p-value of the test (if computed).
}
\item{critValue}{
Boostrap estimation of the critical value (if computed).  
}
\item{replicates}{
Number of bootstrap replicates taken for the computation of the critical value and the p-value.  
}
}
\references{

Baringhaus, L. and Kolbe, D. (2014). \emph{Two-sample tests based on empirical Hankel transforms.} Statistical Papers, 10.1007/s00362-014-0599-1.

Baringhaus, L. and Taherizadeh, F. (2010).\emph{ Empirical Hankel transforms and its applications to
goodness-of-fit tests.} J. Multivariate Anal., 101:1445-14457. 

Gradstein, I. and Ryshik, I. (1981). \emph{Tables.} Harri Deutsch, Frankfurt.
}
\author{Daniel Kolbe} 


\examples{
# comparison of an uniform distribution with an weibull distribution 
x<-runif(30)
y<-rweibull(50,5,0.75)
hankel.test(x,y)

# comparison of an uniform distribution with an weibull distribution 
# in standardized case with parameter lambda=0.1 and replicates=1000

x<-runif(30)
y<-rweibull(50,5,0.75)
hankel.test(x,y,standardized=TRUE,params=list(lambda=0.1),replicates=1000)
}

\keyword{htest}
\keyword{nonparametric}
