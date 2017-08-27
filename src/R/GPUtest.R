library(gputools)

#_______________________________________________________________________________
magnitude <- 10
dimA <- 2*magnitude;dimB <- 3*magnitude;dimC <- 4*magnitude
matA <- matrix(runif(dimA*dimB), dimA, dimB)
matB <- matrix(runif(dimB*dimC), dimB, dimC)

system.time(matA%*%matB);
system.time(gpuMatMult(matA, matB))

#_______________________________________________________________________________
magnitude <- 1000
dimA <- 2*magnitude;dimB <- 3*magnitude;dimC <- 4*magnitude
matA <- matrix(runif(dimA*dimB), dimA, dimB)
matB <- matrix(runif(dimB*dimC), dimB, dimC)

system.time(matA%*%matB);
system.time(gpuMatMult(matA, matB))
