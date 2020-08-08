# Author: Julie BOGOIN
# Modified by: Jinu Han

library(cn.mops)

all <- list.files(path=".", pattern=".analysisready.bam$")


#######################################################################################
print('Working on all...')

setwd(dir=".")


bamDataRanges <- getReadCountsFromBAM(all, WL = 5000, parallel=20)

resCNMOPS <- cn.mops(bamDataRanges)

#plot(resCNMOPS, which=5)

setwd(dir="./cn.mops_output")

segm <- as.data.frame(segmentation(resCNMOPS))
write.csv(segm, file="segmentation_all.csv")

CNVs <- as.data.frame(cnvs(resCNMOPS))
write.csv(CNVs, file="cnvs_all.csv")

CNVRegions <- as.data.frame(cnvr(resCNMOPS))
write.csv(CNVRegions, file="cnvr_all.csv")
