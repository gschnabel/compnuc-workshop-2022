library(ggplot2)
library(gganimate)
library(data.table)
library(mvtnorm)

# for reproducibility of the random sampling
set.seed(5)

##################################################
# examples of covariance functions
##################################################

covfun_sqr <- function(x, y, delta=10, lambda=5, nugget=1e-7) {
  delta^2 * exp(-outer(x,y,`-`)^2/(2*lambda^2)) + outer(x,y,`==`)*nugget^2
}

covfun_mat12 <- function(x, y, delta=10, rho=5, nugget=1e-7) {
  delta^2 * exp(-1/rho * abs(outer(x,y,`-`))) + outer(x,y,`==`)*nugget^2
}

covfun_mat32 <- function(x, y, delta=10, rho=5, nugget=1e-7) {
  d <- abs(outer(x,y,`-`))
  delta^2 * (1+sqrt(3)*d/rho) * exp(-sqrt(3)*d/rho) + outer(x,y,`==`)*nugget^2
}

covfun_lin <- function(x, y, dk=1, dd=10, nugget=1e-7) {
  outer(x,y) * dk^2 + dd^2 + nugget^2
}

##################################################
# convenience functions
##################################################

order_samples <- function(smpl, take_num=nrow(smpl)) {
  
  ord <- rep(0, nrow(smpl))
  sample_used <- logical(nrow(smpl))
  ord[1] <- 1
  sample_used[1] <- TRUE
  distmat <- as.matrix(dist(smpl, method="maximum"))
  num_used <- 1
  cur_idx <- 1
  while (num_used <= take_num) {
    ord[num_used] <- cur_idx
    sample_used[cur_idx] <- TRUE
    next_idx <- which(!sample_used)[which.min(distmat[!sample_used,cur_idx])]
    if (length(next_idx) == 0) break
    num_used <- num_used + 1
    cur_idx <- next_idx
  }
  smpl <- smpl[ord,]
}

##################################################
# implementation of GLS (Gaussian process style)
##################################################

gls <- function(xgrid, xexp, yexp, Cexp, covfun, selcomp=NULL) {
  if (!is.list(covfun)) { covfun <- list(covfun) }
  compound_covfun <- function(x, y) {
    res <- 0
    for (curcovfun in covfun) {
      res <- res + curcovfun(x, y)
    }
    res
  }
  selected_covfun <- if (is.null(selcomp)) compound_covfun else covfun[[selcomp]]
  # assume zero mean function for simplicity
  y1 <- rep(0, length(xgrid))
  y2 <- rep(0, length(xexp))
  A11 <- selected_covfun(xgrid, xgrid)
  A12 <- selected_covfun(xgrid, xexp) 
  A22 <- compound_covfun(xexp, xexp) + Cexp 
  A11post <- A11 - A12 %*% solve(A22) %*% t(A12)
  y1post <- y1 + A12 %*% solve(A22) %*% (yexp - y2) 
  # quick hack to restore symmetry of A11post due to small numerical errors
  A11post <- 0.5 * (A11post + t(A11post))
  list(x = xgrid, y = as.vector(y1post), A = A11post)
}

##################################################
# definition of true function
# and experimental data points
##################################################

truefun <- function(x) {
  exp(0.2*x) + exp(0.1*x) * sin(2*x)
}

create_expdata <- function(n=5, fun, equidist=FALSE) {
  Cexp <- diag(x = 1, nrow=n)
  x <- if (equidist) seq(0,10,length=n) else runif(n, 0, 10)
  y <- fun(x)
  yexp <- as.vector(rmvnorm(1, y, Cexp))
  list(x = x, y = yexp, Cexp = Cexp)
}

# prepare a data table with the true function
xgrid <- seq(0, 10, length=400)
truedf <- data.table(x = xgrid, y = truefun(xgrid)) 

##################################################
# comparison of Gaussian processes (animations)
##################################################

y1 <- rep(0, length(xgrid))
A11 <- covfun_mat12(xgrid, xgrid)
smpl <- order_samples(rmvnorm(400, y1, A11), take_num=400)
smpldt1 <- data.table(
  type = 'matern12',
  frame = rep(1:nrow(smpl), each=length(xgrid)),
  x=xgrid, y=as.vector(t(smpl)))

A11 <- covfun_mat32(xgrid, xgrid)
smpl <- order_samples(rmvnorm(400, y1, A11), take_num=400)
smpldt2 <- data.table(
  type = 'matern32',
  frame = rep(1:nrow(smpl), each=length(xgrid)),
  x=xgrid, y=as.vector(t(smpl)))

A11 <- covfun_sqr(xgrid, xgrid)
smpl <- order_samples(rmvnorm(400, y1, A11), take_num=400)
smpldt3 <- data.table(
  type = 'sqrexp',
  frame = rep(1:nrow(smpl), each=length(xgrid)),
  x=xgrid, y=as.vector(t(smpl)))

A11 <- covfun_lin(xgrid, xgrid)
smpl <- order_samples(rmvnorm(400, y1, A11), take_num=400)
smpldt4 <- data.table(
  type = 'linear',
  frame = rep(1:nrow(smpl), each=length(xgrid)),
  x=xgrid, y=as.vector(t(smpl)))

smpldt <- rbind(smpldt1, smpldt2, smpldt3, smpldt4)

ggp <- ggplot(data=smpldt, aes(x=x,y=y))
ggp <- ggp + theme(panel.grid.minor=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.background=element_rect(fill='white', colour='white'))
ggp <- ggp + geom_line(col="black")
ggp <- ggp + transition_time(frame)
ggp <- ggp + ease_aes('linear')
ggp <- ggp + facet_wrap(vars(type), ncol=2)

animate(ggp, detail=1, nframes=400, width=400, height=400, bg="white")
#anim_save('gp_prior_animation.gif')

##################################################
# comparison of combined Gaussian processes (animations)
##################################################

y1 <- rep(0, length(xgrid))
my_covfun_lin <- function(x, y) {
  covfun_lin(x, y, dk=2, dd=0)
}
A11 <- my_covfun_lin(xgrid, xgrid)
smpl <- order_samples(rmvnorm(50, y1, A11))
smpldt1 <- data.table(
  type = 'linear',
  frame = rep(seq(1,400,length=nrow(smpl)), each=length(xgrid)),
  x=xgrid, y=as.vector(t(smpl)))

my_covfun_sqrexp <- function(x, y) {
  covfun_sqr(x, y, delta=2, lambda=2)
}
A11 <- my_covfun_sqrexp(xgrid, xgrid)
smpl <- order_samples(rmvnorm(400, y1, A11), take_num=400)
smpldt2 <- data.table(
  type = 'sqrexp',
  frame = rep(1:nrow(smpl), each=length(xgrid)),
  x=xgrid, y=as.vector(t(smpl)))

covfun_lin_plus_sqrexp <- function(x, y) {
  covfun_lin(x, y, dk=3, dd=0) + covfun_sqr(x, y, delta=2, lambda=2)
}
A11 <- covfun_lin_plus_sqrexp(xgrid, xgrid)
smpl <- order_samples(rmvnorm(400, y1, A11), take_num=400)
smpldt3 <- data.table(
  type = 'linear+sqrexp',
  frame = rep(1:nrow(smpl), each=length(xgrid)),
  x=xgrid, y=as.vector(t(smpl)))

covfun_lin_times_sqrexp <- function(x, y) {
  covfun_lin(x, y, dk=3, dd=0) * covfun_sqr(x, y, delta=2, lambda=2)
}
A11 <- covfun_lin_times_sqrexp(xgrid, xgrid)
smpl <- order_samples(rmvnorm(400, y1, A11), take_num=400)
smpldt4 <- data.table(
  type = 'linear*sqrexp',
  frame = rep(1:nrow(smpl), each=length(xgrid)),
  x=xgrid, y=as.vector(t(smpl)))

smpldt <- rbind(smpldt1, smpldt2, smpldt3, smpldt4)

smpldt$type <- factor(smpldt$type, ordered=TRUE,
                      levels=c('linear','sqrexp','linear+sqrexp','linear*sqrexp'))

ggp <- ggplot(data=smpldt, aes(x=x,y=y))
ggp <- ggp + theme(panel.grid.minor=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.background=element_rect(fill='white', colour='white'))
ggp <- ggp + geom_line(col="black")
ggp <- ggp + transition_time(frame)
ggp <- ggp + ease_aes('linear')
ggp <- ggp + facet_wrap(vars(type), ncol=2, scales='free')

animate(ggp, detail=1, nframes=400, width=400, height=400, bg="white")
#anim_save('gp_prior_combination_animation.gif')

##################################################
# animations of Gaussian processes constrained with observations
##################################################

# prepare the data
explist <- create_expdata(10, fun=truefun)
expdf <- data.table(x = explist$x, y = explist$y, dy = sqrt(diag(explist$Cexp)))

# fit the Gaussian process to the measurements
my_covfun_sqr <- function(x, y) { covfun_sqr(x, y, delta=10, lambda=2) }
covfun_list <- list(covfun_lin, covfun_mat12, covfun_mat32, my_covfun_sqr)
covfun_names <- c('linear', 'matern12', 'matern32', 'sqrexp')

smpldt <- data.table()
resdt <- data.table()

for (idx in seq_along(covfun_list)) { 
  print(paste0('fitting data points using ', covfun_names[idx]))
  covfun <- covfun_list[[idx]]
  glsres <- gls(xgrid, explist$x, explist$y, explist$Cexp, covfun)
  curresdf <- data.table(type=covfun_names[[idx]],
                         x = glsres$x, 
                         y = glsres$y, 
                         dy = sqrt(diag(glsres$A)))
  smpl <- rmvnorm(400, glsres$y, glsres$A)
  smpl <- order_samples(smpl)
  cursmpldt <- data.table(type=covfun_names[idx],
                          frame=rep(1:nrow(smpl), each=length(xgrid)),
                          x=xgrid,
                          y=as.vector(t(smpl)))
  smpldt <- rbind(smpldt, cursmpldt)
  resdt <- rbind(resdt, curresdf)
}

ggp <- ggplot(data=smpldt, aes(x=x,y=y))
ggp <- ggp + theme(panel.grid.minor=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.background=element_rect(fill='white', colour='white'))
ggp <- ggp + geom_line(col="black")

ggp <- ggp + geom_errorbar(aes(x=x, ymin=y-dy, ymax=y+dy), data=expdf, col='red')
ggp <- ggp + geom_point(aes(x=x, y=y), data=expdf, col='red')

ggp <- ggp + transition_time(frame)
ggp <- ggp + ease_aes('linear')
ggp <- ggp + facet_wrap(vars(type), ncol=2, scales='free')

animate(ggp, detail=1, nframes=400, width=600, height=600, bg="white")
#anim_save('gp_posterior_animation.gif')

##################################################
# plot the posterior uncertainties
##################################################

ggp <- ggplot(data=resdt, aes(x=x,y=y))
ggp <- ggp + theme(panel.grid.minor=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.background=element_rect(fill='white', colour='white'))
ggp <- ggp + geom_line(col="black")
ggp <- ggp + geom_ribbon(aes(x=x, ymin=y-dy, ymax=y+dy), col='gray', alpha=0.2)
ggp <- ggp + geom_errorbar(aes(x=x, ymin=y-dy, ymax=y+dy), data=expdf, col='red')
ggp <- ggp + geom_point(aes(x=x, y=y), data=expdf, col='red')
ggp <- ggp + facet_wrap(vars(type), ncol=2, scales='free')
ggp
#ggsave('gp_posterior_static.png', units='cm', dpi=150, width=20, height=15)


##################################################
# final example: use a compound covariance function
# model + defect and disentangle the 
# posterior of the model contribution
##################################################

parabola_fun <- function(x) { (x-5)^2 }

explist <- create_expdata(50, parabola_fun, equidist=TRUE)
expdt <- data.table(x = explist$x, y = explist$y, dy = sqrt(diag(explist$Cexp)))

my_covfun_mat32 <- function(x, y) { covfun_mat32(x, y, delta=10, rho=3) }
my_covfun_sqr <- function(x, y) { covfun_sqr(x, y, delta=10, lambda=3) }
my_covfun_lin <- function(x, y) { covfun_lin(x, y, dk = 20, dd=50) }

glsres <- gls(xgrid, explist$x, explist$y, explist$Cexp, covfun = list(my_covfun_lin, my_covfun_mat32), selcomp = NULL)
resdf <- data.table(x = glsres$x, 
                    y = glsres$y, 
                    dy = sqrt(diag(glsres$A)))

glsres_mod <- gls(xgrid, explist$x, explist$y, explist$Cexp, covfun = list(my_covfun_lin, my_covfun_mat32), selcomp = 1)
resdf_mod <- data.table(x = glsres_mod$x, 
                        y = glsres_mod$y, 
                        dy = sqrt(diag(glsres_mod$A)))


ggp <- ggplot() + theme_bw() 
# plot experimental data
ggp <- ggp + geom_errorbar(aes(x=x, ymin=y-dy, ymax=y+dy), data=expdt)
ggp <- ggp + geom_point(aes(x=x, y=y), data=expdt)
# plot overall prediction
ggp <- ggp + geom_ribbon(aes(x=x, ymin=y-dy, ymax=y+dy), data=resdf, alpha=0.3, fill="blue")
ggp <- ggp + geom_line(aes(x=x, y=y), data=resdf, col="blue")
# plot model component prediction
ggp <- ggp + geom_ribbon(aes(x=x, ymin=y-dy, ymax=y+dy), data=resdf_mod, alpha=0.3, fill="green")
ggp <- ggp + geom_line(aes(x=x, y=y), data=resdf_mod, col="green")
ggp
#ggsave("disentangled_posterior_lin_sqrexp_50points.png", units="cm", width=12, height=10)
