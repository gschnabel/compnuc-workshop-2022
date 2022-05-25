###################################################################
#
# Apply the Generalized Least Squares method to reconstruct
# faces from pictures where some pixels intensities are missing.
#
###################################################################

library(pixmap)
library(Matrix)
library(mvtnorm)

picWidth <- 112
picHeight <- 92
numPts <- picWidth*picHeight

plotFace <- function(faceVec) {
  dim(faceVec) <- c(picWidth,picHeight)
  facePic <- pixmapGrey(faceVec)
  plot(facePic)
}

# read all pictures and convert it to vector
basePath <- "../data/att_faces"
picFiles <- list.files(basePath,pattern="\\.pgm",recursive=TRUE,full.names=TRUE)
picList <- lapply(picFiles,function(x) read.pnm(x))
# build a covariance matrix
picVecList <- t(sapply(picList,function(x) {res <- getChannels(x); dim(res) <- NULL; res}))
covStruc <- cov.wt(picVecList)

# show faces
#png("face_database.png",width=800,height=500)
faceSel <- 1:36
par(mfrow=c(6,6),mar=c(0,0,0,0),oma=c(0,0,0,0))
for (nr in faceSel) {
  #baseScale <- range(covStruc$center) - mean(covStruc$center)
  plot(picList[[(nr-1)*10+1]])
}
#dev.off()

#show one face
#png("face_sample03.png",width=112*3,height=92*3)
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(picList[[(1-1)*10+1]])
#dev.off()

#png("face_mean.png",width=112*3,height=92*3)
# show the mean face
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
workVec <- covStruc$center
plotFace(workVec)
#dev.off()

# create mapping from complete face to partial observation
createS <- function(idx) {
  sparseMatrix(i=seq_along(idx),j=idx,x=1,dims=c(length(idx),numPts))
}
obsIdx <- sort(sample(numPts,5e3,replace=FALSE))
S <- createS(obsIdx)
B <- diag(rep(0.01,length(obsIdx)))

# remove information from face

#png("reduced_face.png",width=112*5,height=92*5)
par(mfrow=c(1,2),mar=c(5,0,0,5),oma=c(0,0,0,1))
#faceNr <- 28*10-3
faceNr <- 382
plot(picList[[faceNr]])
workVec <- as.vector(getChannels(picList[[faceNr]]))
workVec[-obsIdx] <- 0
plotFace(workVec)
#dev.off()

# apply Bayesian statistics to complete the picture
A <- covStruc$cov
SAST <- as.matrix(S%*%A%*%t(S))
Xmat <- chol2inv(chol(SAST + B))
newMean <- covStruc$center + A%*%t(S)%*%Xmat%*%(workVec[obsIdx]-S%*%covStruc$center)
#newCov <- A - A%*%t(S)%*%Xmat%*%S%*%A

# show the predicted face
#png("myface_reconstruction_05.png",width=6*3,height=6*112/92,units="cm",res=150)
par(mfrow=c(1,3),mar=c(0,0,0,0),oma=c(0,0,0,0))
#plotFace(covStruc$center)
#plot(myFace)
plotFace(workVec)
plotFace(newMean)
plot(picList[[faceNr]])
#dev.off()

# sample from the prior distribution
# or sample from the posterior distribution 
#png("face_priorsamples.png",width=92*8,height=112*5)
# use this line to sample from prior
# the eigenvector decomposition may take a while
# maybe something in the ballpark of ten minutes
system.time(eigStruc <- eigen(A,symmetric = TRUE))
# or that line to sample from posterior
#system.time(eigStruc <- eigen(newCov,symmetric=TRUE))
numFaces <- 8*5
par(mfrow=c(5,8),mar=c(0,0,0,0),oma=c(0,0,0,0))
for (i in seq(numFaces)) {
  sel <- 1:400
  zsmpl <- rnorm(length(sel))*sqrt(eigStruc$values[1:400])
  priorFaceSmpl <- covStruc$center + eigStruc$vectors[,sel] %*% zsmpl
  plotFace(priorFaceSmpl)
}
#dev.off()


# show some eigenfaces
faceSel <- 1:36
par(mfrow=c(6,6),mar=c(0,0,0,0),oma=c(0,0,0,0))
for (nr in faceSel) {
  #baseScale <- range(covStruc$center) - mean(covStruc$center)
  baseScale <- range(newMean) - mean(newMean)
  eigFaceScale <- max(baseScale/range(eigStruc$vectors[,nr]))
  #eigFace <- covStruc$center + eigStruc$vectors[,nr]*eigFaceScale
  eigFace <- newMean + eigStruc$vectors[,nr]*eigFaceScale
  print(range(eigFace))
  plotFace(eigFace)
}
