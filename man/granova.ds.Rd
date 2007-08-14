\name{granova.ds}
\alias{granova.ds}
\title{Granova for display of Dependent Samples}
\description{Plots dependent samples as a scatterplot, together with various summary information: effect size, distribution of X-Y differences,
means for X and Y, a 95\% confidence interval for effect size.  Standard t-statistic and other numerical results also returned.
}
\usage{
granova.ds(xdata, revc = FALSE, sw = 0.4, ne = 0.5, labcex = 1, ident = FALSE, pt.lab = NULL,
xlab = NULL, ylab = NULL, main = NULL, sub = NULL, par.orig = TRUE)
}
\arguments{
  \item{xdata}{is an n X 2 dataframe. First column defines X (intially for horzontal axis), the second defines Y. }
  \item{revc}{reverses X,Y specifications.}
  \item{sw}{extends axes toward lower left, effectively moving data points to the southwest.}
  \item{ne}{extends axes toward upper right, effectively moving data points to northeast. Making both sw and ne smaller moves points farther
    apart, while making both larger moves data points closer together.}
  \item{labcex}{controls size of axes labels.}
  \item{ident}{logical, default FALSE. Allows user to identify individual points.}
  \item{pt.lab}{optional character vector defining labels for points.  Only used if ident is TRUE.
                If NULL, rownames(xdata) are used if available; if not 1:n is used.}
  \item{xlab}{optional label (as character) for horizontal axis. If not defined, axis labels are taken from colnames of xdata.}
  \item{ylab}{optional label (as character) for vertical axis.}
  \item{main}{optional main title (as character); if not supplied by user generic title is provided.}
  \item{sub}{optional subtile (as character).}
  \item{par.orig}{returns par to original settings; if multipanel plots it is advisable to specify FALSE.}
}
\details{Paired X & Y values are plotted as scatterplot. The identity reference line (for Y=X) is drawn. Since the better data view often entails
 having X's > Y's the revc argument facilitates reversal of the X, Y specifications. Because all data points are plotted relative to the identity
 line, and summary results are shown graphically, clusters, data trends, outliers, and possible uses of transformations are readily seen, possibly
 to be accommodated. Parallel projections of data points to (a lower-left) line segment show how each point relates to its X-Y = D difference; blue
 'crosses' are used to display the distribution of difference scores and the mean difference is displayed as a heavy dashed (red) line, parallel to
 the identity reference line. Means for X and Y are also plotted (as thin dashed vertical and horizontal lines), and rug plots are shown for the
 distributions of X (at the top of graphic) and Y (on the right side). Several summary statistics are plotted as well, to facilitate both description
 and inference; see below. The 95\% confidence interval for the population mean difference is also shown graphically.
  }
\value{A list is returned with the following components:
  \item{mean(X)}{Mean of X values}
   \item{mean(Y)}{Mean of Y values}
    \item{mean(D=X-Y)}{Mean of differences D = X - Y}
     \item{SD(D)}{Standard deviation of differences D}
      \item{ES(D)}{Effect Size for differences D: mean(D)/SD(D)}
       \item{r(X,Y)}{Correlation based on X,Y pairs}
        \item{r(x+y,D)}{Correlation based on X+Y,D pairs}
         \item{LL 95\%CI}{Lower bound for 95\% confidence interval for population mean(D)}
          \item{UL 95\%CI}{Upper bound for 95\% confidence interval for population mean(D)}
           \item{t(D-bar)}{t-statistic associated w/ test of hypothesis that population mean(D) = 0.0}
            \item{df.t}{Degrees of freedom for the t-statistic}
             \item{pval.t}{P-value for two sided t-test of null hypothesis that population mean(D) does not equal zero.}
}
\references{Exploratory Plots for Paired Data, Rosenbaum P., The American Statistician, May 1989, vol. 43, no. 2, pp. 108-9.}
\author{Robert M. Pruzek \email{RMPruzek@yahoo.com}

James E. Helmreich \email{ James.Helmreich@Marist.edu}
 }

\examples{
library(DAAG) #Data analysis and graphics, by John McDonald
data(pair65) #length of stretched rubber bands; scores show length after being heated (X), as well as at ambient temperature (Y).
granova.ds(pair65, main = "Dependent sample assessment plot for pair65 data, n = 9")

}
\keyword{hplot}
\keyword{htest}