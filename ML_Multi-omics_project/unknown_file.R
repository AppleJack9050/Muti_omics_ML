library(iClusterPlus) 
library(GenomicRanges) 
library(gplots) 
library(lattice)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%

setwd("/home/taserface/Documents/Precision Medicine ML project")
data <- read.csv(file = 'brca_data_w_subtypes.csv')

GE = as.matrix(data[, 0:604])
CNA = as.matrix(data[, 605:1464])
SM = as.matrix(data[, 1465:1713])
PE = as.matrix(data[, 1714:1936])
# 
# GE = GE[, 0:10]
# CNA = CNA[, 0:10]
# SM = SM[, 0:10]
# PE = PE[, 0:10]


#SM = SM[,-(which(colSums(SM)==0))]
#SM %>% as.data.frame %>% select(where(~ sum(.) != 0)) 

# head(CNA)
# CNA[CNA == -1] <- 2
# CNA[CNA == -2] <- 3

# for(k in 1:5){ 
#   cv.fit = tune.iClusterPlus(cpus=12, dt1=GE, dt2=PE,
#                              type=c("gaussian", "gaussian"), K=k,
#                              scale.lambda=c(1,1), maxiter=20) 
#   save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep="")) 
# }

# data(variation.hg18.v10.nov.2010) 
# gbm.cn=CNregions(seg=gbm.seg,epsilon=0,adaptive=FALSE,rmCNV=TRUE, 
#                  cnv=variation.hg18.v10.nov.2010[,3:5], 
#                  frac.overlap=0.5, rmSmallseg=TRUE,nProbes=5)
# 
# gbm.exp = gbm.exp[, 0 : 5]
# gbm.cn = gbm.cn[, 0 : 5]
# 

# fit.single = iClusterPlus(dt1=GE, dt2=CNA, dt3=PE, dt4=SM,
#                            type=c("gaussian","multinomial","gaussian", 'binomial'),
#                            lambda=c(0.04, 0.61, 1, 1),K=2,maxiter=10)

# fit.single = iClusterPlus(dt1=GE, dt2=CNA, dt3=PE,
#                           type=c("gaussian","multinomial","gaussian"),
#                           lambda=c(0.04, 0.61, 1),K=2,maxiter=10)

for(k in 1:5){ 
  cv2.fit = tune.iClusterPlus(cpus=8,dt1=GE, dt2=CNA, dt3=SM,dt4=PE, 
                              type=c("gaussian","multinomial","binomial","gaussian"), K=k,  
                              scale.lambda=c(1,1), maxiter=20) 
  save(cv2.fit, file=paste("cv2.fit.k",k,".Rdata",sep=""))
}

output2=alist() 
files=grep("cv2.fit",dir()) 
for(i in 1:length(files)){ 
  load(dir()[files[i]])
  output2[[i]]=cv2.fit 
} 
nLambda = nrow(output2[[1]]$lambda) 
nK = length(output2) 
BIC = getBIC(output2) 
devR = getDevR(output2)

minBICid = apply(BIC,2,which.min) 
devRatMinBIC = rep(NA,nK) 
for(i in 1:nK){ 
  devRatMinBIC[i] = devR[minBICid[i],i] 
}

plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)", 
     ylab="%Explained Variation")
