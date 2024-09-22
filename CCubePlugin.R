rm(list = ls())
library(dplyr)
library(ccube)
library(doParallel)

registerDoParallel(cores=3)

set.seed(1234)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")


input <- function(inputfile) {
        pfix = prefix()
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
}

run <- function() {}

output <- function(outputfile) {
numSv <- as.integer(parameters["numSv", 2])#100
maxiter <- as.integer(parameters["maxiter", 2])#200
init =as.integer(parameters["init", 2])#5  # number of clusters
numOfClusterPool = as.integer(parameters["clusterA", 2]):as.integer(parameters["clusterB", 2])#1:6
numOfRepeat = as.integer(parameters["numrepeats", 2])#5
baseDepth = as.integer(parameters["baseDepth", 2])#40
ccfSet <- c(as.numeric(parameters["ccfX", 2]), as.numeric(parameters["ccfY", 2]), as.numeric(parameters["ccfZ", 2]))#c(1, 0.3, 0.7) # true ccf pool
ccfTrue <- sample(ccfSet, numSv, c(0.5,0.3,0.2), replace = T)

purity <- 0.8


cnPoolMaj <- c(1,2,3,4)
cnPoolMin <- c(0,1,2)


# 1st break point
cnProfile1 <- sample(cnPoolMaj, numSv, c(0.25, 0.25, 0.25,0.25), replace =T)
cnProfile2 <- sample(cnPoolMin, numSv, c(1/3, 1/3, 1/3), replace =T)

cnProfileTot <- cnProfile1 + cnProfile2

cnProfile <- cbind(cnProfile1, cnProfile2, cnProfileTot )
cnProfile <- t( apply(cnProfile, 1, sort) )


mydata <- data.frame(mutation_id = paste0("ss", seq_len(numSv)) ,
                     ccf_true = ccfTrue,
                     minor_cn1 = cnProfile[,1],
                     major_cn1 = cnProfile[,2],
                     total_cn1 = cnProfile[,3],
                     stringsAsFactors = F)


mydata$purity <- purity
mydata$normal_cn <- 2

mydata <- mutate(rowwise(mydata),
                 mult1 = sample(c(1,if (major_cn1 ==1) { 1 } else {major_cn1}), 1),
                 true_vaf1 = cp2ap(ccf_true, purity, normal_cn, total_cn1, total_cn1, mult1),
                 total_counts1 = rpois(1, total_cn1/2 * baseDepth),
                 var_counts1 = rbinom(1, total_counts1, true_vaf1),
                 ref_counts1 = total_counts1 - var_counts1)

# 2nd break point

cnProfile1 <- sample(cnPoolMaj, numSv, c(0.5, 0.2, 0.2,0.1), replace =T)
cnProfile2 <- sample(cnPoolMin, numSv, c(0.3, 0.5, 0.2), replace =T)

cnProfileTot <- cnProfile1 + cnProfile2

cnProfile <- cbind(cnProfile1, cnProfile2, cnProfileTot )
cnProfile <- t( apply(cnProfile, 1, sort) )

mydata$minor_cn2 = cnProfile[,1]
mydata$major_cn2 = cnProfile[,2]
mydata$total_cn2 = cnProfile[,3]

mydata <- mutate(rowwise(mydata),
                 mult2 = sample(c(1,if (major_cn2 ==1) { 1 } else {major_cn2}), 1),
                 true_vaf2 = cp2ap(ccf_true, purity, normal_cn, total_cn2, total_cn2, mult2),
                 total_counts2 = rpois(1, total_cn1/2 * baseDepth ),
                 var_counts2 = rbinom(1, total_counts2, true_vaf2),
                 ref_counts2 = total_counts2 - var_counts2)




# SV prototype
doubleBreakPtsRes <- RunCcubePipeline(ssm = mydata, modelSV = T,
                                      numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
                                      runAnalysis = T, runQC = T, , multiCore = T)

fn1 = paste(outputfile, "double_break_points_results.pdf", sep=".")
MakeCcubeStdPlot_sv(res = doubleBreakPtsRes$res, ssm = doubleBreakPtsRes$ssm, printPlot = T, fn = fn1)

# Run ccube seperately

mydata1 = mydata[, c("mutation_id", "var_counts1", "ref_counts1", "major_cn1", "minor_cn1", "purity", "normal_cn")]
mydata1 = dplyr::rename(mydata1, var_counts = var_counts1, ref_counts = ref_counts1,
                        major_cn = major_cn1, minor_cn = minor_cn1)

breakPt1Res <- RunCcubePipeline(ssm = mydata1,
                                numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
                                runAnalysis = T, runQC = T, multiCore = T)


fn2 = paste(outputfile, "ccube_1st_breakpoint.pdf", sep=".")
MakeCcubeStdPlot(ssm = breakPt1Res$ssm, res = breakPt1Res$res, printPlot = T, fn = fn2)

mydata2= mydata[, c("mutation_id", "var_counts2", "ref_counts2", "major_cn2", "minor_cn2", "purity", "normal_cn")]
mydata2 = dplyr::rename(mydata2, var_counts = var_counts2, ref_counts = ref_counts2, major_cn = major_cn2, minor_cn = minor_cn2)

breakPt2Res <- RunCcubePipeline(ssm = mydata2,
                                numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
                                runAnalysis = T, runQC = T, multiCore = T)

fn3 = paste(outputfile, "ccube_2nd_breakpoint.pdf", sep=".")
MakeCcubeStdPlot(ssm = breakPt2Res$ssm, res = breakPt2Res$res, printPlot = T, fn = fn3)

# compare event ccf

label1 = doubleBreakPtsRes$res$label
label2 = breakPt1Res$res$label
label3 = breakPt2Res$res$label

mydata$ccube_double_mult1 = doubleBreakPtsRes$res$full.model$bv1
mydata$ccube_double_mult2 = doubleBreakPtsRes$res$full.model$bv2

mydata$ccube_single_mult1 = breakPt1Res$res$full.model$bv
mydata$ccube_single_mult2 = breakPt2Res$res$full.model$bv

mydata <- mutate(rowwise(mydata),
                 vaf1 = var_counts1/(var_counts1+ref_counts1),
                 ccube_double_ccf1 = MapVaf2CcfPyClone(vaf1,
                                                purity,
                                                normal_cn,
                                                major_cn1 + minor_cn1,
                                                major_cn1 + minor_cn1,
                                                ccube_double_mult1,
                                                constraint=F),
                 ccube_single_ccf1 = MapVaf2CcfPyClone(vaf1,
                                                       purity,
                                                       normal_cn,
                                                       major_cn1 + minor_cn1,
                                                       major_cn1 + minor_cn1,
                                                       ccube_single_mult1,
                                                       constraint=F),
                 true_obs_ccf1 = MapVaf2CcfPyClone(vaf1,
                                                   purity,
                                                   normal_cn,
                                                   major_cn1 + minor_cn1,
                                                   major_cn1 + minor_cn1,
                                                   mult1,
                                                   constraint=F),
                 vaf2 = var_counts2/(var_counts2+ref_counts2),
                 ccube_double_ccf2 = MapVaf2CcfPyClone(vaf2,
                                                purity,
                                                normal_cn,
                                                major_cn2 + minor_cn2,
                                                major_cn2 + minor_cn2,
                                                ccube_double_mult2,
                                                constraint=F),
                 ccube_single_ccf2 = MapVaf2CcfPyClone(vaf2,
                                                       purity,
                                                       normal_cn,
                                                       major_cn2 + minor_cn2,
                                                       major_cn2 + minor_cn2,
                                                       ccube_single_mult2,
                                                       constraint=F),
                 true_obs_ccf2 = MapVaf2CcfPyClone(vaf2,
                                                   purity,
                                                   normal_cn,
                                                   major_cn2 + minor_cn2,
                                                   major_cn2 + minor_cn2,
                                                   mult2,
                                                   constraint=F)
)





myColors=gg_color_hue(10)

fn = paste(outputfile, "event_ccf_comparsions.pdf", sep=".")
pdf(fn, width=8, height=8)

par(mfrow=c(2,2))
plot(mydata$true_obs_ccf1, mydata$ccube_double_ccf1, col = myColors[label1],
     xlim = c(0, max( c(mydata$true_obs_ccf1, mydata$ccube_double_ccf1) ) ),
     ylim = c(0, max( c(mydata$true_obs_ccf1, mydata$ccube_double_ccf1) ) ),
     xlab = "true ccf", ylab = "estimated ccf", main = "double model: 1st break point")
points( seq(0, max( c(mydata$true_obs_ccf1, mydata$ccube_double_ccf1) ), length.out = 100 ),
        seq(0, max( c(mydata$true_obs_ccf1, mydata$ccube_double_ccf1) ), length.out = 100 ),
        type = "l" )

plot(mydata$true_obs_ccf2, mydata$ccube_double_ccf2, col = myColors[label1],
     xlim = c(0, max( c(mydata$true_obs_ccf2, mydata$ccube_double_ccf2) ) ),
     ylim = c(0, max( c(mydata$true_obs_ccf2, mydata$ccube_double_ccf2) ) ),
     xlab = "true ccf", ylab = "estimated ccf", main = "double model: 2nd break point"
)

points( seq(0, max( c(mydata$true_obs_ccf2, mydata$ccube_double_ccf2) ), length.out = 100 ),
        seq(0, max( c(mydata$true_obs_ccf2, mydata$ccube_double_ccf2) ), length.out = 100 ),
        type = "l" )

plot(mydata$true_obs_ccf1, mydata$ccube_single_ccf1, col = myColors[label2],
     xlim = c(0, max( c(mydata$true_obs_ccf1, mydata$ccube_single_ccf1) ) ),
     ylim = c(0, max( c(mydata$true_obs_ccf1, mydata$ccube_single_ccf1) ) ),
     xlab = "true ccf", ylab = "estimated ccf", main = "single model: 1st break point")
points( seq(0, max( c(mydata$true_obs_ccf1, mydata$ccube_single_ccf1) ), length.out = 100 ),
        seq(0, max( c(mydata$true_obs_ccf1, mydata$ccube_single_ccf1) ), length.out = 100 ),
        type = "l" )

plot(mydata$true_obs_ccf2, mydata$ccube_single_ccf2, col = myColors[label3],
     xlim = c(0, max( c(mydata$true_obs_ccf2, mydata$ccube_single_ccf2) ) ),
     ylim = c(0, max( c(mydata$true_obs_ccf2, mydata$ccube_single_ccf2) ) ),
     xlab = "true ccf", ylab = "estimated ccf", main = "single model: 2nd break point"
)

points( seq(0, max( c(mydata$true_obs_ccf2, mydata$ccube_single_ccf2) ), length.out = 100 ),
        seq(0, max( c(mydata$true_obs_ccf2, mydata$ccube_single_ccf2) ), length.out = 100 ),
        type = "l" )

dev.off()
}





