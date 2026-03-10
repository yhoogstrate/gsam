

if(!require(minfi)){
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("minfi")
}

if(!require(devtools)) install.packages("devtools")   
library(devtools)
# prevent warnings from beeing converted to errors when calling install_github
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
install_github('mwsill/RFpurify')


if(!require(minfiData)){
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("minfiData")
}



library(RFpurify)
library(minfiData)

# load example data
data(MsetEx)
MsetEx



# predict purity
absolute <- predict_purity(MsetEx,method="ABSOLUTE")
estimate <- predict_purity(MsetEx,method="ESTIMATE")

# plot predicted ABSOLUTE purity against ESTIMATE purity
plot(x=estimate,y=absolute,pch=19,ylim=c(0,1),xlim=c(0,1))
abline(a=0,b=1,col="red")


