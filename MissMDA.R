library(rmatio)
library(missMDA)
#library(mice)

cstr <- read.mat('datasets/cstr.mat'); cstr.data <- cstr$fea
webACE <- read.mat('datasets/WebACE.mat'); webACE.data <- webACE$fea
classic3 <- read.mat('datasets/classic3.mat'); classic3.data <- classic3$A
#amputeWebACE <- ampute(webACE$fea, 0.3)

## Function to generate missing Data
genMissingData <- function(datas, tauxMissing = 0.3) {
  dimx = dim(datas)[1]
  dimy = dim(datas)[2]
  samples <- as.integer(length(datas) * tauxMissing)
  vals <- sample.int(dimx*dimy, samples, replace = F)
  
  #dim(unique(x))
  datas[vals] = NA
  print(sum(is.na(datas))==samples)
  return(datas)
}

# -------- Generate missing in range (5->35, by 5)
## Cstr
cstr5  <- genMissingData(cstr.data, 0.05); NA_cstr5  <- which(is.na(cstr5), arr.ind=TRUE); cstr5_I <- imputeCA(cstr5)
cstr10 <- genMissingData(cstr.data, 0.1);  NA_cstr10 <- which(is.na(cstr10), arr.ind=TRUE); cstr10_I <- imputeCA(cstr10)
cstr15 <- genMissingData(cstr.data, 0.15); NA_cstr15 <- which(is.na(cstr15), arr.ind=TRUE); cstr15_I <- imputeCA(cstr15)
cstr20 <- genMissingData(cstr.data, 0.2);  NA_cstr20 <- which(is.na(cstr20), arr.ind=TRUE); cstr20_I <- imputeCA(cstr20)
cstr25 <- genMissingData(cstr.data, 0.25); NA_cstr25 <- which(is.na(cstr25), arr.ind=TRUE); cstr25_I <- imputeCA(cstr25)
cstr30 <- genMissingData(cstr.data, 0.3);  NA_cstr30 <- which(is.na(cstr30), arr.ind=TRUE); cstr30_I <- imputeCA(cstr30)
cstr35 <- genMissingData(cstr.data, 0.35); NA_cstr35 <- which(is.na(cstr35), arr.ind=TRUE); cstr35_I <- imputeCA(cstr35)

# WebACE
webACE5  <- genMissingData(webACE.data, 0.05); NA_webACE5  <- which(is.na(webACE5), arr.ind=TRUE);  webACE5_I <- imputeCA(webACE5)
webACE10 <- genMissingData(webACE.data, 0.1);  NA_webACE10 <- which(is.na(webACE10), arr.ind=TRUE); webACE10_I <- imputeCA(webACE10)
webACE15 <- genMissingData(webACE.data, 0.15); NA_webACE15 <- which(is.na(webACE15), arr.ind=TRUE); webACE15_I <- imputeCA(webACE15)
webACE20 <- genMissingData(webACE.data, 0.2);  NA_webACE20 <- which(is.na(webACE20), arr.ind=TRUE); webACE20_I <- imputeCA(webACE20)
webACE25 <- genMissingData(webACE.data, 0.25); NA_webACE25 <- which(is.na(webACE25), arr.ind=TRUE); webACE25_I <- imputeCA(webACE25)
webACE30 <- genMissingData(webACE.data, 0.3);  NA_webACE30 <- which(is.na(webACE30), arr.ind=TRUE); webACE30_I <- imputeCA(webACE30)
webACE35 <- genMissingData(webACE.data, 0.35); NA_webACE35 <- which(is.na(webACE35), arr.ind=TRUE); webACE35_I <- imputeCA(webACE35)

# Classic3
classic5  <- genMissingData(classic3.data, 0.05); NA_classic5  <- which(is.na(as.matrix(classic5)), arr.ind=TRUE);  classic5_I  <- imputeCA(classic5)
#classic10 <- genMissingData(classic3.data, 0.1);  NA_classic10 <- which(is.na(as.matrix(classic10)), arr.ind=TRUE); classic10_I <- imputeCA(classic10)
classic15 <- genMissingData(classic3.data, 0.15); NA_classic15 <- which(is.na(as.matrix(classic15)), arr.ind=TRUE); classic15_I <- imputeCA(classic15)
#classic20 <- genMissingData(classic3.data, 0.2);  NA_classic20 <- which(is.na(as.matrix(classic20)), arr.ind=TRUE); classic20_I <- imputeCA(classic20)
classic25 <- genMissingData(classic3.data, 0.25); NA_classic25 <- which(is.na(as.matrix(classic25)), arr.ind=TRUE); classic25_I <- imputeCA(classic25)
#classic30 <- genMissingData(classic3.data, 0.3);  NA_classic30 <- which(is.na(as.matrix(classic30)), arr.ind=TRUE); classic30_I <- imputeCA(classic30)
#classic35 <- genMissingData(classic3.data, 0.35); NA_classic35 <- which(is.na(as.matrix(classic35)), arr.ind=TRUE); classic35_I <- imputeCA(classic35)

## ·· Write Files ··
# Cstr
write.mat(as.data.frame(NA_cstr5),'datasets/cstr/NA_cstr5.mat')
write.csv(as.data.frame(cstr5_I),'datasets/cstr/cstr5_I.csv', row.names = F)

write.mat(as.data.frame(NA_cstr10),'datasets/cstr/NA_cstr10.mat')
write.csv(as.data.frame(cstr10_I),'datasets/cstr/cstr10_I.csv', row.names = F)

write.mat(as.data.frame(NA_cstr15),'datasets/cstr/NA_cstr15.mat')
write.csv(as.data.frame(cstr15_I),'datasets/cstr/cstr15_I.csv', row.names = F)

write.mat(as.data.frame(NA_cstr20),'datasets/cstr/NA_cstr20.mat')
write.csv(as.data.frame(cstr20_I),'datasets/cstr/cstr20_I.csv', row.names = F)

write.mat(as.data.frame(NA_cstr25),'datasets/cstr/NA_cstr25.mat')
write.csv(as.data.frame(cstr25_I),'datasets/cstr/cstr25_I.csv', row.names = F)

write.mat(as.data.frame(NA_cstr30),'datasets/cstr/NA_cstr30.mat')
write.csv(as.data.frame(cstr30_I),'datasets/cstr/cstr30_I.csv', row.names = F)

write.mat(as.data.frame(NA_cstr35),'datasets/cstr/NA_cstr35.mat')
write.csv(as.data.frame(cstr35_I),'datasets/cstr/cstr35_I.csv', row.names = F)

# WebACE
write.mat(as.data.frame(NA_webACE5),'datasets/webACE/NA_webACE5.mat')
write.csv(as.data.frame(webACE5_I),'datasets/webACE/webACE5_I.csv', row.names = F)

write.mat(as.data.frame(NA_webACE10),'datasets/webACE/NA_webACE10.mat')
write.csv(as.data.frame(webACE10_I),'datasets/webACE/webACE10_I.csv', row.names = F)

write.mat(as.data.frame(NA_webACE15),'datasets/webACE/NA_webACE15.mat')
write.csv(as.data.frame(webACE15_I),'datasets/webACE/webACE15_I.csv', row.names = F)

write.mat(as.data.frame(NA_webACE20),'datasets/webACE/NA_webACE20.mat')
write.csv(as.data.frame(webACE20_I),'datasets/webACE/webACE20_I.csv', row.names = F)

write.mat(as.data.frame(NA_webACE25),'datasets/webACE/NA_webACE25.mat')
write.csv(as.data.frame(webACE25_I),'datasets/webACE/webACE25_I.csv', row.names = F)

write.mat(as.data.frame(NA_webACE30),'datasets/webACE/NA_webACE30.mat')
write.csv(as.data.frame(webACE30_I),'datasets/webACE/webACE30_I.csv', row.names = F)

write.mat(as.data.frame(NA_webACE35),'datasets/webACE/NA_webACE35.mat')
write.csv(as.data.frame(webACE35_I),'datasets/webACE/webACE35_I.csv', row.names = F)

# Classic3
write.mat(as.data.frame(NA_classic5),'datasets/classic3/NA_classic5.mat')
write.csv(as.data.frame(classic5_I),'datasets/classic3/classic5_I.csv', row.names = F)

write.mat(as.data.frame(NA_classic15),'datasets/classic3/NA_classic15.mat')
write.csv(as.data.frame(classic15_I),'datasets/classic3/classic15_I.csv', row.names = F)

write.mat(as.data.frame(NA_classic25),'datasets/classic3/NA_classic25.mat')
write.csv(as.data.frame(classic25_I),'datasets/classic3/classic25_I.csv', row.names = F)



