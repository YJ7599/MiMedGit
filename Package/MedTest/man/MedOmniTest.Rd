\name{MedOmniTest}
\alias{MedOmniTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Distance-based Omnibus Test for Mediation Effects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A distance-based omnibus test for mediation effects of a high-dimensional structured mediator.
}
\usage{
MedOmniTest(x, y, m.list, z = NULL, nperm = 9999)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
 independent variable (univariable)
}
  \item{y}{
%%     ~~Describe \code{y} here~~
 response variable (continuous univariable)
}
  \item{m.list}{
%%     ~~Describe \code{m.list} here~~
 a list of distance matrices that summarize the high-dimensional structured mediator in pair-wise distances.
}
  \item{z}{
%%     ~~Describe \code{z} here~~
 variables to be adjusted (e.g. potential confounders). Can be a vector (univariable) or a matrix  (multivariable) of numerical values. \code{model.matrix} should be used if it is a data frame. Intercept (first column) should be dropped. 
}
  \item{nperm}{
%%     ~~Describe \code{nperm} here~~
 number of permutations.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{margPs }{a list of single distance-based mediation test p-values.}
  \item{permP }{the omnibus test p-value calculated via permutation tests.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Jie Zhang, Zhi Wei and Jun Chen. A distance-based approach for testing the mediation effects of the human microbiome. 2018 Jun 1;34(11):1875-1883
}
\author{
%%  ~~who you are~~
Jun Chen, Jie Zhang and Zhi Wei
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(matrixStats)
library(vegan)
library(GUniFrac)
library(MedTest)

data(aofib)
data(bmi)
data(tree.rooted)
data(otu.tab)

unifracs <- GUniFrac(otu.tab, tree.rooted)$unifracs
m.list <- list(BC=vegdist(otu.tab, method="bray"), 
		JAC=as.matrix(vegdist(otu.tab, 'jaccard', binary=TRUE)),
		UniFrac=unifracs[, , c('d_UW')],
		GUniFrac=unifracs[, , c('d_0.5')], 
		WUniFrac=unifracs[, , c('d_1')])

set.seed(12345)
MedOmniTest(aofib, bmi, m.list)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Metagenomics }
\keyword{ Mediation }% __ONLY ONE__ keyword per line
