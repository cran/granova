granova.2w <- function(data.A.B, fit = "linear", ident = FALSE){

# Function depicts two-way ANOVA using data-based contrasts.
# Argument 'data.A.B' should be an n X 3 dataframe.  If the rows are named uniquely, then points
# will be able to be identified with those labels, otherwise the row number is used.
# The first column should be the response values, the second and third columns should
# be factors defining the levels.  The factors are ordered
# by the data (row/col means); they are are not taken to be ordered at the outset.
# Function 'scatter3d' be must be available in R, which is easiest to do by loading its source, Rcmdr (thanks,John Fox).
# The 'fit' is defaulted to 'linear' where interaction is not fit (i.e. flat
# surface based on predicting cell means using row and column marginal effects).
# Replace 'linear' with, say, quadratic to produce a curved surface.
# This version uses a product contrast predictor in scatter3d, in effect depicting graphically Tukey's one
# degree of freedom for non-additivity; the last line of the model summary for the 2nd ANOVA shows the Tukey
# 1 df.non-additivity.effect
# Note: right click on the scatterplot to terminate 'identify' and return the output from the function.

mtx <- is.data.frame(data.A.B)
if(!mtx)data.A.B <- data.frame(data.A.B)

vA <- length(unique(data.A.B[,2]))
vB <- length(unique(data.A.B[,3]))
N  <- dim(data.A.B)[1]
rnd2 <- function(x)round(x,2)

A <-  data.A.B[,2]
B <-  data.A.B[,3]
yy <- data.A.B[,1]

#Means in each cell
mns <- tapply(yy, list(A,B), mean)
mns.vec <- as.vector(mns)
mns.matrx<-matrix(mns.vec,ncol=vB)

#Next two lines are of class 'array'; must be changed to class numeric
cell.mnsA <- tapply(yy, list(A), mean)
cell.mnsB <- tapply(yy, list(B), mean)

#For scat3d; points to group means as separate group and thus colors them differently.
group.factor <- factor(c(rep(0,N), rep(1, vA*vB)))

#mnsA,B are lengths of A,B vectors with appropriate means
mnsA <- as.numeric(cell.mnsA)[A]
mnsB <- as.numeric(cell.mnsB)[B]
ordA<-order(cell.mnsA)
ordB<-order(cell.mnsB)

mns.matrx<-mns.matrx[ordA,ordB]
cell.cnts<-table(A,B)[ordA,ordB]

#grandmean
grndmean <- mean(yy)

#these two vectors are effects cum data-based contrasts for A & B factors
facA.mn.cntrst <- mnsA - grndmean
facB.mn.cntrst <- mnsB - grndmean

aov.yy <- aov(yy ~ factor(A)*factor(B))

#Trying to put means in: something of a hack.  Adding points at the means for each cell,
#then using the group feature of scatter3d to give them a different color.   Scatter3d
#calculates the separate regression plane for the means, but it is the same as the plane
#calculated for the data.  

mnsAA <- rep(as.numeric(cell.mnsA), length(unique(B))) - grndmean
mnsBB <- rep(as.numeric(cell.mnsB), ea = length(unique(A))) - grndmean
facA.mn.cntrst <- c(facA.mn.cntrst, mnsAA)
facB.mn.cntrst <- c(facB.mn.cntrst, mnsBB)

yy <- c(yy,mns.vec)
aaa<- paste(1:length(unique(A)), rep(paste(1:length(unique(B)), "mean", sep=""), ea = length(unique(A))), sep = "")

out <- list('A.effects' = signif(sort(cell.mnsA-grndmean),3), 'B.effects' = signif(sort(cell.mnsB-grndmean),3),
             CellCounts.Reordered = signif(cell.cnts,3), CellMeans.Reordered = signif(mns.matrx,3), aov.summary = summary(aov.yy))

if(is.null(rownames(data.A.B))){rownames(data.A.B) <- 1:N}

scatter3d(facA.mn.cntrst, yy, facB.mn.cntrst, xlab = colnames(data.A.B)[2], ylab = 'Response', 
    zlab = colnames(data.A.B)[3], group = group.factor, fogtype='exp2',fov=55, surface = TRUE, fit = fit, surface.col = c(4,8))
if(ident){identify3d(facA.mn.cntrst, yy, facB.mn.cntrst, labels = c(rownames(data.A.B), aaa))}

return(out)
}
