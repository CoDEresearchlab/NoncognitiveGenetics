#### Multivariate twin model-fitting (NGS 151115)
####
#### Model: multivariate saturated and Cholesky ACE decomposition
#### Script based on BivACE.R, 1Saturated.R and 2ACE_Chol.R from 2015 SGDP summer school.


### ---------------------------------------------------------------------------
### Setup / retrieve data

## clean workspace
rm(list=ls())

## import dependencies
Sys.setenv(OMP_NUM_THREADS=parallel::detectCores()) 
library(OpenMx)
library(psych)
library(foreign)

mxOption(NULL,"Default optimizer","SLSQP")
mxOption(key='Number of Threads', value=parallel::detectCores())
## set working directory
setwd('/Users/margheritamalanchini/Dropbox/2019 Queen Mary/MM papers 2019/MM_AGA_Genomic SEM Achievement/Twin modelling OpenMx')

## input/output configuration

inFile <- './data/130417tt&al.dat'				# prepared data to analyse
SPSSFile <- './data/WIP_100417tt&al.sav'
outFileStem <- './results/Test-teach-Alevel-July2018/mult_'

## open data
rawDat <- read.spss(SPSSFile, use.value.labels=F, to.data.frame=T)
labels <- attr(rawDat, 'variable.labels')

## open cleaned data
dat <- read.table(
	inFile,
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)



## select variables

 # define variables. For Cholesky analysis, decomposition is sequential, predicting
 # the last variable specified.
VarsList <- list(

	#3 variables combinations 
	c( 'EngKS1testarsr','NPDTeachEngKS1arsr','rcqalmrars'),
	c( 'EngKS2testarsr','NPDTeachEngKS2arsr','rcqalmrars'),
	c( 'EngKS3testarsr','NPDTeachEngKS3arsr','rcqalmrars'),
	c( 'MatKS1testarsr','NPDTeachMatKS1arsr','rcqalmrars'),
	c( 'MatKS2testarsr','NPDTeachMatKS2arsr','rcqalmrars'),
	c( 'MatKS3testarsr','NPDTeachMatKS3arsr','rcqalmrars'),
	c( 'SciKS2testarsr','NPDTeachSciKS2arsr','rcqalmrars'),
	c( 'SciKS3testarsr','NPDTeachSciKS3arsr','rcqalmrars')
)



### ---------------------------------------------------------------------------
### Start loop (repeat for each vector of variables in VarsList)

for (i in 1:length(VarsList)) {


## select variables (all variables for twin 1, then all for twin 2)
Vars <- VarsList[[i]]
selVars <- c(paste(Vars, "1", sep=""), paste(Vars, "2", sep=""))

## number of variables
nv		<- length(Vars)		# number of variables for a twin
ntv		<- 2*nv				# number of variables for a pair
ncor	<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

## report progress
cat("\n", paste("--- Starting run for", Vars[1], ", ", Vars[2], "... - group", i, "of", length(VarsList)), "---\n")

## select data
mzData <- subset(dat, zygos==1, selVars)
dzData <- subset(dat, zygos==2, selVars)


## -----------------------------------------------------------------------
## Print descriptive statistics by zygosity group
#
#describe(mzData)
#describe(dzData)
#colMeans(mzData,na.rm=TRUE)
#cov(mzData,use="complete")
#cor(mzData,use="complete")
#colMeans(dzData,na.rm=TRUE)
#cov(dzData,use="complete")
#cor(dzData,use="complete")


## -----------------------------------------------------------------------
## Create starting values
##
## TODO: is this satisfactory?  I don't like randomness for this, as it leads
## to results even the original researcher can't reproduce precisely. Currently
## adding a fairly large fixed increment here, as smaller ones (< c. 0.5) are
## producing status RED (code 6) warnings, probably indicating relatively low
## confidence that the optimiser hasn't got stuck at a local maximum/minimum.
##
## (a, c, and e path coefficients in Cholesky ACE decomposition are currently
## using arbitrary default values).

#jiggle		<-rnorm(nlower, mean = 0, sd = .1)
StMZmean	<-colMeans(mzData,na.rm=TRUE) + 1  # (i.e., add fixed increment instead of jiggle)
StDZmean	<-colMeans(dzData,na.rm=TRUE) + 1
StMZcov		<-vech(t(chol(cov(mzData,use="complete")))) + 1
StDZcov		<-vech(t(chol(cov(dzData,use="complete")))) + 1

Stmean <-colMeans(mzData[,1:nv],na.rm=TRUE)  # for the ACE model


## -----------------------------------------------------------------------
## 1. Specify and run a fully saturated model (Cholesky decomposition)

## matrices for the MZ and DZ data
meansMZ		<-mxMatrix("Full", 1, ntv, free=TRUE, values=StMZmean, name="expMeanMZ") 
pathsMZ		<-mxMatrix("Lower", ntv, ntv, free = TRUE, values = StMZcov, name="CholMZ")
covMZ		<-mxAlgebra(CholMZ %*% t(CholMZ), name="expCovMZ")

meansDZ		<-mxMatrix("Full", 1, ntv, free=TRUE, values=StDZmean, name="expMeanDZ") 
pathsDZ		<-mxMatrix("Lower", ntv, ntv, free=TRUE, values=StDZcov, name="CholDZ")
covDZ		<-mxAlgebra(CholDZ %*% t(CholDZ), name="expCovDZ")

## data objects for multiple groups
dataMZ		<- mxData( observed=mzData, type="raw" )
dataDZ		<- mxData( observed=dzData, type="raw" )

## objective objects for multiple groups
objMZ		<- mxExpectationNormal( covariance="expCovMZ", means="expMeanMZ", dimnames=selVars )
objDZ		<- mxExpectationNormal( covariance="expCovDZ", means="expMeanDZ", dimnames=selVars )

fitFunction	<- mxFitFunctionML()

## combine groups
modelMZ		<- mxModel( meansMZ,pathsMZ, covMZ, dataMZ, objMZ, fitFunction, name="MZ" )
modelDZ		<- mxModel( meansDZ,pathsDZ, covDZ, dataDZ, objDZ, fitFunction, name="DZ" )
minus2ll	<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj			<- mxFitFunctionAlgebra( "m2LL" )
SatModel	<- mxModel( "Sat", modelMZ, modelDZ, minus2ll, obj )

## run saturated model
SatFit		<- mxRun(SatModel)
(SatSum		<- summary(SatFit))

## generate saturated output
#round(SatFit@output$estimate,4)
#
#mxEval(MZ.expMeanMZ, SatFit)
#mxEval(MZ.expCovMZ, SatFit)
#mxEval(DZ.expMeanDZ, SatFit)
#mxEval(DZ.expCovDZ, SatFit)


## -----------------------------------------------------------------------
## 2. Specify and run Cholesky ACE decomposition 
## (Assumes equal means/variances across birth order and zygosity)

## create labels for lower triangular matrices
aLabs		<- paste("a", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep="")
cLabs		<- paste("c", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep="")
eLabs		<- paste("e", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep="")

## create labels for column and diagonal matrices
mLabs		<- paste("m",1:nv,sep="")

## matrices declared to store a, c, and e path coefficients
pathA		<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=.4, labels=aLabs, name="a" )
pathC		<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=.1, labels=cLabs, name="c" )
pathE		<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=.2, labels=eLabs, name="e" )
	
## matrices generated to hold A, C, and E computed variance components
covA		<- mxAlgebra( expression=a %*% t(a), name="A" )
covC		<- mxAlgebra( expression=c %*% t(c), name="C" ) 
covE		<- mxAlgebra( expression=e %*% t(e), name="E" )
covP		<- mxAlgebra( expression=A+C+E, name="V" )
StA			<- mxAlgebra( expression=A/V, name="h2" )
StC			<- mxAlgebra( expression=C/V, name="c2" )
StE			<- mxAlgebra( expression=E/V, name="e2" )

## algebra to compute correlations
matI		<- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
Rph			<- mxAlgebra( expression=solve(sqrt(I*V)) %&% V , name="Phcor")
Rg			<- mxAlgebra( expression=solve(sqrt(I*A)) %&% A , name="Acor")
Rc			<- mxAlgebra( expression=solve(sqrt(I*C)) %&% C , name="Ccor")
Re			<- mxAlgebra( expression=solve(sqrt(I*E)) %&% E , name="Ecor")

# algebra to compute standardised paths
invSD		<- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")
sta			<- mxAlgebra( expression=iSD %&% a, name="sta")
stc			<- mxAlgebra( expression=iSD %&% c, name='stc')
ste			<- mxAlgebra( expression=iSD %&% e, name='ste')

# algebra to compute squared standardised paths
sta2		<- mxAlgebra( expression=sta * sta, name='sta2')
stc2		<- mxAlgebra( expression=stc * stc, name='stc2')
ste2		<- mxAlgebra( expression=ste * ste, name='ste2')

## algebra for expected means and variance/covariance matrices in MZ & DZ twins
Mean		<- mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=Stmean, labels=c(mLabs,mLabs), name="ExpMean" )

covMZ		<- mxAlgebra( expression= rbind( cbind(A+C+E , A+C),
                                           cbind(A+C , A+C+E)),		name="ExpCovMZ" )
covDZ		<- mxAlgebra( expression= rbind( cbind(A+C+E       , 0.5%x%A+C),
                                           cbind(0.5%x%A+C , A+C+E)),	name="ExpCovDZ" )

## data objects for multiple groups
dataMZ		<- mxData( observed=mzData, type="raw" )
dataDZ		<- mxData( observed=dzData, type="raw" )

## objective objects for multiple groups
objMZ		<- mxExpectationNormal( covariance="ExpCovMZ", means="ExpMean", dimnames=selVars )
objDZ		<- mxExpectationNormal( covariance="ExpCovDZ", means="ExpMean", dimnames=selVars )

fitFunction	<- mxFitFunctionML()

## combine groups
params		<- list( pathA, pathC, pathE, covA, covC, covE, covP, StA, StC, StE, matI, Rph, Rg, Rc, Re, invSD, sta, stc, ste, sta2, stc2, ste2) 
modelMZ		<- mxModel( params, covMZ, Mean, dataMZ, objMZ, fitFunction, name="MZ" )
modelDZ		<- mxModel( params, covDZ, Mean, dataDZ, objDZ, fitFunction, name="DZ" )
minus2ll	<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj			<- mxFitFunctionAlgebra( "m2LL" )

## define required confidence intervals

 # for 4 variables:
#conf1		<- mxCI (c ('MZ.h2[1,1]', 'MZ.h2[2,2]', 'MZ.h2[3,3]', 'MZ.h2[4,4]') )  # h2
#conf2		<- mxCI (c ('MZ.c2[1,1]', 'MZ.c2[2,2]', 'MZ.c2[3,3]', 'MZ.c2[4,4]') )  # c2
#conf3		<- mxCI (c ('MZ.e2[1,1]', 'MZ.e2[2,2]', 'MZ.e2[3,3]', 'MZ.e2[4,4]') )  # e2
#conf4		<- mxCI (c ('MZ.Acor[2,1]','MZ.Acor[3,1]','MZ.Acor[4,1]','MZ.Acor[3,2]','MZ.Acor[4,2]','MZ.Acor[4,3]'))  # Rg
#conf5		<- mxCI (c ('MZ.Ccor[2,1]','MZ.Ccor[3,1]','MZ.Ccor[4,1]','MZ.Ccor[3,2]','MZ.Ccor[4,2]','MZ.Ccor[4,3]'))  # Rc
#conf6		<- mxCI (c ('MZ.Ecor[2,1]','MZ.Ecor[3,1]','MZ.Ecor[4,1]','MZ.Ecor[3,2]','MZ.Ecor[4,2]','MZ.Ecor[4,3]'))  # Re
#conf7		<- mxCI( c(
#				'sta2[1,1]', 'sta2[2,1]', 'sta2[3,1]', 'sta2[4,1]',
#				'sta2[2,2]', 'sta2[3,2]', 'sta2[4,2]',
#				'sta2[3,3]', 'sta2[4,3]',
#				'sta2[4,4]'
#				) )  # standardised, squared genetic paths

 # for 3 variables:
conf1		<- mxCI (c ('MZ.h2[1,1]', 'MZ.h2[2,2]', 'MZ.h2[3,3]') )  # h2
conf2		<- mxCI (c ('MZ.c2[1,1]', 'MZ.c2[2,2]', 'MZ.c2[3,3]') )  # c2
conf3		<- mxCI (c ('MZ.e2[1,1]', 'MZ.e2[2,2]', 'MZ.e2[3,3]') )  # e2
conf4		<- mxCI (c ('MZ.Acor[2,1]', 'MZ.Acor[3,1]', 'MZ.Acor[3,2]') )  # Rg
conf5		<- mxCI (c ('MZ.Ccor[2,1]', 'MZ.Ccor[3,1]', 'MZ.Ccor[3,2]') )  # Rc
conf6		<- mxCI (c ('MZ.Ecor[2,1]', 'MZ.Ecor[3,1]', 'MZ.Ecor[3,2]') )  # Re
conf7		<- mxCI( c(
			'sta2[1,1]', 'sta2[2,1]', 'sta2[3,1]',
			'sta2[2,2]', 'sta2[3,2]',
			'sta2[3,3]'
			) )  # standardised, squared genetic paths

## define model
CholAceModel<- mxModel( "CholACE", params, modelMZ, modelDZ, minus2ll, obj, conf1, conf2, conf3, conf4, conf5, conf6, conf7)

## run model
CholAceFit	<- mxRun(CholAceModel, intervals=T)
(CholAceSum	<- summary(CholAceFit))

## generate ACE output
#round(CholAceFit@output$estimate,4)
#
#mxEval(CholACE.h2, CholAceFit)
#mxEval(CholACE.c2, CholAceFit)
#mxEval(CholACE.e2, CholAceFit)
#
#mxEval(CholACE.Acor, CholAceFit)
#mxEval(CholACE.Ccor, CholAceFit)
#mxEval(CholACE.Ecor, CholAceFit)
#mxEval(CholACE.Phcor, CholAceFit)


## -----------------------------------------------------------------------
## Compare models: fit statistics
#
#mxCompare(SatFit, CholAceFit)


## -----------------------------------------------------------------------
## Output results to file

## open file
fileName <- paste(outFileStem, paste(Vars, collapse="-"), sep="")
sink( fileName, append=T)

## title (variables)
varStr <- ""
for (i in 1:(length(selVars)/2)) {
	varStr = paste(varStr, "  ", selVars[i], ": ", labels[selVars[i]], "\n", sep="")
}

cat( paste("Multivariate modelling for:\n", varStr, "\n", sep="") )

## test model assumptions. TODO: currently unused
#cat("\n--Fit statistics: compare saturated Gaussian decomposition to constrained submodel--\n")
#cat("(Equal means and variances between twin1/2 and zygosity group; equal within-individual\ncross-trait correlations; symmetric cross-twin cross-trait correlations in MZ and DZ groups.)\n")
#print( mxCompare(SatGFit, Sub1Fit) )

## test ACE model fit
cat("\n--Fit statistics: compare saturated Cholesky decomposition to ACE model--\n")
print( mxCompare(SatFit, CholAceFit) )

## path estimates
cat("\n--Path estimates (unstandardised)--\n")
print(round(CholAceFit@output$estimate,4))

## standardised paths
cat("\n--Path estimates (standardised, not squared)--\n")
stPaths <- cbind(CholAceFit$CholACE.sta@result, CholAceFit$CholACE.stc@result, CholAceFit$CholACE.ste@result)
colnames(stPaths) <- paste(rep( c("a", "c", "e"), each=length(selVars)/2), 1:nv, sep="")
rownames(stPaths) <- paste("v", 1:nrow(stPaths), sep="")
print(round(stPaths,4))

## standardised, squared paths
cat("\n--Path estimates (standardised, squared)--\n")
stSqPaths <- cbind(CholAceFit$CholACE.sta2@result, CholAceFit$CholACE.stc2@result, CholAceFit$CholACE.ste2@result)
colnames(stSqPaths) <- paste(rep( c("a", "c", "e"), each=length(selVars)/2), 1:nv, sep="")
rownames(stSqPaths) <- paste("v", 1:nrow(stSqPaths), sep="")
print(round(stSqPaths,4))

## components of rPh
#cat("\n--Component contributions to rPh--\n")
#rACE <- CholAceFit$CholACE.RphACE@result
#colnames(rACE) <- c("A", "C", "E")
#rownames(rACE) <- "contrib"
#print(round(rACE,4))

## multivariate ACE estimates
cat("\n--Multivariate ACE estimates--\n")
cat("(Diagonals: univariate results. Off-diagonals: proportions of rPh.)\n")
aceEsts <- cbind(CholAceFit$CholACE.h2@result, CholAceFit$CholACE.c2@result, CholAceFit$CholACE.e2@result)
colnames(aceEsts) <- rep( c("h2", "c2", "e2"), each=length(selVars)/2 )
rownames(aceEsts) <- paste("v", 1:nrow(aceEsts), sep="")
print(round(aceEsts,4))

## correlations
cat("\n--Correlations--\n")
corrs <- cbind(CholAceFit$CholACE.Phcor@result[2,1], CholAceFit$CholACE.Acor@result[2,1], CholAceFit$CholACE.Ccor@result[2,1], CholAceFit$CholACE.Ecor@result[2,1])
colnames(corrs) <- c("rPh", "rA", "rC", "rE")
print(round(corrs,4))

## full ACE model summary (including CIs)
cat("\n--ACE model summary--\n")
print(CholAceSum)

sink()


### ---------------------------------------------------------------------------
### End loop
}
