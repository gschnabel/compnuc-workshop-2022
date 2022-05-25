# kalmann example

library(animation)
library(mixtools)
library(png)
library(mvtnorm)
library(Matrix)

plotCar <- function(xc,dx,img) {
  
  l <- 10; b <- l * 159/296
  dx <- dx / sqrt(sum(dx^2))
  rx <- c(dx[2],-dx[1])
  xf <- xc + l*dx*0.5
  xb <- xc - l*dx*0.5
  xfl <- xf + b*rx/2
  xfr <- xf - b*rx/2
  xbl <- xb - b*rx/2
  xbr <- xb + b*rx/2
  
  rasterImage(img,xbr[1],xbr[2],
              xbr[1]+l,xbr[2]+b,
              atan2(dx[2],dx[1])*180/pi)
  #lines(c(xfl[1],xfr[1],xbr[1],xbl[1],xfl[1]),
  #      c(xfl[2],xfr[2],xbr[2],xbl[2],xfl[2]),col="red")
}

followPath <- function(curstate,path) {
  maxphi <- 10
  tardist <- 7
  # current state
  #curstate: pos, dir, idx
  idx <- curstate$est$idx
  speed <- sqrt(sum(curstate$est$speed^2))
  dir <- curstate$est$speed / speed
  rx <- curstate$est$pos
  # checkpoint reached?
  curdist <- 0; idx <- idx-1
  while(curdist < tardist) {
    idx <- idx+1
    if (idx > length(path$x)) idx <- 1
    nx <- c(path$x[idx],path$y[idx])
    curdist <- sqrt(sum((nx-rx)^2)) 
  }
  # do path correction
  vecang <- nx - rx 
  vecang <- vecang/sqrt(sum(vecang^2))
  if (acos(min(sum(dir*vecang),1))*180/pi < maxphi) {
    newdir <- vecang
  } else {
    vdir <- c(dir[2],-dir[1])
    s <- sign(sum(vdir*vecang))
    newdir <- dir + s*vdir * tan(maxphi*pi/180)
    newdir <- newdir / sqrt(sum(newdir^2))
  }
  speed <- as.vector(newdir*speed)
  truespeed <- as.vector(speed + rmvnorm(1,c(0,0),curstate$est$uncspeed)) 
  curstate$est$pos <- rx+speed
  curstate$est$speed <- speed
  curstate$est$uncpos <- propagatePriorCov(curstate)
  curstate$est$idx <- idx
  curstate$true$pos <- curstate$true$pos + truespeed
  curstate$true$speed <- truespeed
  curstate
}

makeplot <- function() {
  
  curstate <- list(
    timestep = 0,
    est = list(pos=c(90,35),speed=c(0,2),idx=1, 
               uncpos=diag(c(50^2,50^2)),uncspeed=diag(c(0.25^2,0.25^2))),
    true = list(pos=c(90,35),speed=c(0,2)))
  smoothspeed <- matrix(curstate$est$speed,2,10)

  tseq <- seq(0,1,length=500)
  for (t in tseq) {
    curstate$timestep <- curstate$timestep + 1
    par(pty="s",oma=c(0,0,0,0),mai=c(0,0,0,0),
        xaxs="i",yaxs="i",bg='NA')
    plot(0,0,xlim=c(0,100),ylim=c(0,100),type="n",axes=FALSE,xlab="",ylab="")
    rasterImage(marsTexture,0,0,100,100,0)
    curstate <- followPath(curstate,myPath)
    curstate <- updateEst(curstate)
    #smoothspeed <- cbind(smoothspeed[,1:9],curstate$est$speed)
    # plot the situation
    with(myPath,lines(x,y))
    with(curstate$est,plotCar(pos,speed,ferrariPic2))
    with(curstate$true,plotCar(pos,speed,ferrariPic1))
    with(curstate$true,points(pos[1],pos[2],pch=16,col="white"))
    with(curstate$est,points(pos[1],pos[2],pch=16,col="white"))
    with(curstate$est,ellipse(pos,uncpos,col="white",lwd=3))
  }
}

# statistics part

gpsObs <- function(truestate) {
  B <- diag(rep(400,2))
  list(y = as.vector(truestate$pos + rmvnorm(1,c(0,0),B)),
       B = B)
}

laserObs <- function(truestate) {
  pos <- truestate$pos
  if ((pos[2] > 53 && pos[2] < 57) || 
      (pos[2] > 88 && pos[2] < 92)) {
    B <- diag(c(1e3^2,0.2^2))
    y <- as.vector(pos + rmvnorm(1,c(0,0),B))
    list(y=y,B=B)
  } else if ((pos[1] > 38 && pos[1] < 42) || 
             (pos[1] > 78 && pos[1] < 82)) {
    B <- diag(c(0.2^2,1e3^2))
    y <- as.vector(pos + rmvnorm(1,c(0,0),B))
    list(y=y,B=B)
  } else
    NULL
}

propagatePriorCov <- function(curstate) {
  T <- rbind(c(1,0,1,0),
             c(0,1,0,1))
  x0 <- with(curstate$est,T%*%c(pos,speed))
  A0 <- with(curstate$est,T%*%bdiag(uncpos,uncspeed)%*%t(T))
  print(diag(A0))
  A0
}

updateEst <- function(curstate) {

  useLaser <- TRUE
  useGPS <- TRUE
  
  x0 <- curstate$est$pos
  A0 <- curstate$est$uncpos
  obsGPS <- gpsObs(curstate$true)
  obsLaser <- laserObs(curstate$true)
  #with(obs,ellipse(y,B,col="violet",lwd=3))
  x1 <- x0
  A1 <- A0
  # gps update
  if (isTRUE(useGPS)) {
    #with(obsGPS,points(y[1],y[2],pch=16,col="violet",cex=3))
    X <- chol2inv(chol(A0 + obsGPS$B))
    x1 <- x0 + A0%*%X%*%(obsGPS$y - x0)
    A1 <- A0 - A0%*%X%*%A0
  }
  if (isTRUE(useLaser)) {
    abline(h=c(55,90),col="red",lwd=3)
    abline(v=c(40,80),col="red",lwd=3)  
    # laser update
    if (!is.null(obsLaser)) {
      
      X <- chol2inv(chol(A1 + obsLaser$B))
      x1 <- x1 + A1%*%X%*%(obsLaser$y - x1)
      A1 <- A1 - A1%*%X%*%A1
    }
  }
  # updating done
  curstate$est$pos <- as.vector(x1)
  curstate$est$uncpos <- A1
  curstate
}

# load the images

myPath <- readRDS("../data/example_path.RData")
ferrariPic1 <- readPNG("../images/ferrari_red.png")
ferrariPic2 <- readPNG("../images/ferrari_yellow.png")
marsTexture <- readPNG("../images/mars_texture.png")

# produce animations

set.seed(30)
system.time(saveVideo(makeplot(),video.name="test.avi",interval=1/30,
          other.opts="-vcodec ffvhuff -threads 4"))

png("mars_laser_beam.jpg",width=7,height=7,res=150, units="cm")

par(pty="s",oma=c(0,0,0,0),mai=c(0,0,0,0),
    xaxs="i",yaxs="i",bg='NA')
plot(0,0,xlim=c(0,100),ylim=c(0,100),type="n",axes=FALSE,xlab="",ylab="")
rasterImage(marsTexture,0,0,100,100,0)
# plot the situation
with(myPath,lines(x,y))
abline(h=c(55,90),col="red",lwd=3)
abline(v=c(40,80),col="red",lwd=3)  
dev.off()
