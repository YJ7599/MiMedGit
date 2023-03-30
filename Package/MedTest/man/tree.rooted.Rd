\name{tree.rooted}
\alias{tree.rooted}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Phylogenetic Tree
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("tree.rooted")}
\format{
  The format is:
List of 5
 $ edge       : num [1:13307, 1:2] 12694 12694 12693 12693 12692 ...
 $ Nnode      : num 6634
 $ tip.label  : chr [1:6674] "16687" "14531" "5260" "12208" ...
 $ edge.length: num [1:13307] 0.0146 0.0637 0.02627 0.0179 0.00014 ...
 $ node.label : chr [1:6633] "" "0.026" "0.007" "0.000" ...
 - attr(*, "class")= chr "phylo"
 - attr(*, "order")= chr "postorder"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(tree.rooted)
## maybe str(tree.rooted) ; plot(tree.rooted) ...
}
\keyword{datasets}
