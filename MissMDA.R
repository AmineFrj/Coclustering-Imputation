library(rmatio)
library(missMDA)
#library(mice)

webACE <- read.mat('datasets/WebACE.mat')
#webACE$fea
#amputeWebACE <- ampute(webACE$fea, 0.3)

# Data : Données
# tauxMissing : % de données missing (valeurs entre 0 et 1 avec 0.5 = 50%)
genMissingData <- function(datas, tauxMissing = 0.3) {
  print("just check if 2 next numbers are equals...")
  dimx = dim(datas)[1]
  dimy = dim(datas)[2]
  samples <- as.integer(length(datas) * tauxMissing)
  print(samples)
  vals <- sample.int(dimx*dimy, samples, replace = F)
  
  x = cbind(vals %% dimx + 1, vals %% dimy + 1)
  #dim(unique(x))
  datas[x] = NA
  print(sum(is.na(datas)))
  return(datas)
}


missWebACE <- genMissingData(webACE$fea, 0.2)
sum(is.na(missWebACE))
NAs <- which(is.na(missWebACE), arr.ind=TRUE)

res <- imputeCA(missWebACE)
which(is.na(amputeWebACE$amp), arr.ind=TRUE)

write.mat(as.data.frame(NAs),'NaNs.mat')
write.csv(as.data.frame(res),'MissWebACE.csv', row.names = F)

















