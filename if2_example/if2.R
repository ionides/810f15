## ----setup,include=FALSE, results="hide"---------------------------------

# TEST = TRUE is quick version, to check the code runs
TEST <- TRUE

# knitr options
require("knitr")
opts_chunk$set(
               progress=T,prompt=F,tidy=F,highlight=T,
               warning=F,message=F,error=F,
               results='hide',echo=F,dev='png',
               size='scriptsize',
               fig.path='figure/',fig.lp="fig:",
               fig.align='left',
               fig.show='asis',
               fig.height=3.35,fig.width=3.35,
               out.width="3.35in",
               dpi=300,
               dev=c('png','tiff'),
               dev.args=list(
                 png=list(bg='transparent'),
                 tiff=list(compression='lzw')
                 )
               )

options(
        stringsAsFactors=FALSE,
        help_type="html",
        scipen=-1
        )

scinot <- function (x, digits = 2, type = c("expression","latex")) {
  type <- match.arg(type)
  x <- signif(x,digits=digits)
  ch <- floor(log10(abs(x)))
  mn <- x/10^ch
  switch(type,
         expression={
           bquote(.(mn)%*%10^.(ch))
         },
         latex={
           paste0("\\scinot{",mn,"}{",ch,"}")
         }
         )
}


## ----load-packages,include=FALSE-----------------------------------------

## the following installs 'pomp' if it isn't already installed.
if (!("pomp" %in% rownames(installed.packages()))) {
  install.packages("pomp")
}

require(pomp)

stopifnot(packageVersion("pomp")>="0.50-9")

require(xtable)

options(
        xtable.caption.placement="top",
        xtable.include.rownames=FALSE
        )


## ----cluster,include=FALSE, results="hide"-------------------------------

nodefile <- Sys.getenv("PBS_NODEFILE")
## are we on a cluster? 
CLUSTER <- nchar(nodefile)>0
## if CLUSTER=FALSE, assume we are on a multicore machine
  

## ----toy,results='hide'--------------------------------------------------

binary.file <- "toy.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {

CORES <- 30    ## number of parallel processes
JOBS <- 30     ## number of estimation runs
NMIF <- 100     ## number of IF iterations per estimation run
NMCMC <- 10000  ## number of filtering iterations per estimation run
NP <- 100     ## number of particles in pfilter operations

if(TEST){
  CORES <- 10    ## number of parallel processes
  JOBS <- 10     ## number of estimation runs
  NMIF <- 20     ## number of IF iterations per estimation run
  NMCMC <- 20  ## number of filtering iterations per estimation run
  NP <- 50     ## number of particles in pfilter operations
}


if (CLUSTER) {
  require(doSNOW)
  hostlist <- tail(scan(nodefile,what=character(0)),-1)
  cl <- makeCluster(hostlist,type='SOCK')
  registerDoSNOW(cl)
} else {
  require(doParallel)
  registerDoParallel(CORES)
}

require(pomp)

theta.true <- c(th1=1,th2=1,x1.0=0,x2.0=0)
sd.y.true <- c(10,1)

toy.proc.sim <- function (x, t, params, ...) {
 th1 <- params["th1"]
 th2 <- params["th2"]
 xx <- c(exp(th1),exp(th1)*th2)
 names(xx) <- c("x1","x2")
 return(xx)
}

toy.meas.sim <- function (x, t, params, sd.y, ...) {
 yy <- rnorm(n=2,mean=x,sd=sd.y)
 c(y1=yy[1],y2=yy[2])
}

toy.meas.dens <- function (y, x, t, params, sd.y, ..., log) {
f <- dnorm(x=y,mean=x,sd=sd.y,log=log)
if (log) sum(f) else prod(f)
}

Times <- 100
dat <- matrix(1,nrow=2,ncol=Times,dimnames=list(c("y1","y2"),NULL))

toy <- pomp(data=dat,times=1:Times, t0=0,
            rprocess=discrete.time.sim(step.fun=toy.proc.sim,delta.t=1),
            rmeasure=toy.meas.sim,dmeasure=toy.meas.dens,
            sd.y=sd.y.true)

loglik.exact <- function(po,params){
 th1 <- params["th1"]
 th2 <- params["th2"]
 xx <-  c(x1=exp(th1),x2=exp(th1)*th2)
 sum(apply(obs(po),2,toy.meas.dens,x=xx,t=0,params=params,sd.y=sd.y.true,log=T))
}

set.seed(555)
po <- simulate(toy,params=theta.true) 
mean.dat <- apply(obs(po),1, mean)
mle.exact <- c(th1=unname(log(mean.dat[1])),
               th2=unname(mean.dat[2]/mean.dat[1]))
max.exact <- loglik.exact(po,mle.exact)

## RUN MIF1 TO ESTIMATE PARAMETERS

rwsd <- c(th1=0.1,th2=0.1)

th1.range <- c(-2,2)
th2.range <- c(0,10)

tic <- Sys.time()
mpar1 <- foreach(i=1:JOBS,
                 .packages='pomp',
                 .inorder=FALSE) %dopar% 
{
  set.seed(123+i)
  th.draw <-c(th1=runif(1,min=th1.range[1],max=th1.range[2]),
              th2=runif(1,min=th2.range[1],max=th2.range[2]))
  m <- mif(
           po,
           Nmif=NMIF,
           start=c(th.draw,x1.0=0,x2.0=0),
           pars=c("th1","th2"),
           Np=NP,
           ic.lag=100,
           method="mif",
           rw.sd=rwsd,
           cooling.type="geometric",
           cooling.fraction=sqrt(0.1),
           var.factor=2
           )
  list(pomp=m,start=th.draw)
}
toc <- Sys.time()
etime1 <- toc-tic

m1.out <- rbind(
                sapply(mpar1,function(x)coef(x$pomp,c("th1","th2"))),
                lik.exact = sapply(mpar1,function(x) loglik.exact(x$pomp,coef(x$pomp)))
                )

## RUN MIF2 TO ESTIMATE PARAMETERS

tic <- Sys.time()
mpar2 <- foreach(i=1:JOBS,
                 .packages='pomp',
                 .inorder=FALSE) %dopar% 
{
  set.seed(123+i)
  th.draw <-c(th1=runif(1,min=th1.range[1],max=th1.range[2]),
              th2=runif(1,min=th2.range[1],max=th2.range[2]))
  m <- mif(
           po,
           Nmif=NMIF,
           start=c(th.draw,x1.0=0,x2.0=0),
           pars=c("th1","th2"),
           Np=NP,
           ic.lag=100,
           method="mif2",
           rw.sd=rwsd,
           cooling.type="geometric",
           cooling.fraction=sqrt(0.1),
           var.factor=2
           )
  list(pomp=m,start=th.draw)
}
toc <- Sys.time()
etime2 <- toc-tic

m2.out <- rbind(
                sapply(mpar2,function(x)coef(x$pomp,c("th1","th2"))),
                lik.exact = sapply(mpar2,function(x) loglik.exact(x$pomp,coef(x$pomp)))
                )

## RUN PMCMC TO ESTIMATE PARAMETERS

th.min <- c(th1=-2,th2=0)
th.max <- c(th1=2,th2=10)
toy.hyperparams <- list(min=th.min,max=th.max)
toy.dprior <- function(params, hyperparams, ..., log)
{
  f <- sum(dunif(params,
                 min=hyperparams$min,
                 max=hyperparams$max,
                 log=TRUE))
  if (log) f else exp(f)
}

toy.rprior <- function(params,hyperparams, ...)
{
  params[c("th1","th2")] <- runif(n=2,min=hyperparams$min,max=hyperparams$max)
  params
}

po <- pomp(po,
           rprior=toy.rprior,
           dprior=toy.dprior,
           hyperparams=toy.hyperparams)

tic <- Sys.time()
mpar3 <- foreach (i=1:30,
                  .packages='pomp',
                  .inorder=FALSE) %dopar% 
{
  set.seed(567+i)
  th.draw <- rprior(po,coef(po))[,1]
  m <- pmcmc(po,
             Nmcmc=NMCMC,
             start=th.draw,
             Np=NP,
             rw.sd=rwsd,
             max.fail=1000
             )
  list(pomp=m,start=th.draw)
}
toc <- Sys.time()
etime3 <- toc-tic

m3.out <- rbind(
                sapply(mpar3,function(x)coef(x$pomp,c("th1","th2"))),
                lik.exact = sapply(mpar3,function(x) loglik.exact(x$pomp,coef(x$pomp)))
                )

m3.trace <- lapply(mpar3,function(x)conv.rec(x$pomp,c("th1","th2")))

### COMPUTE EXACT LIKELIHOODS

L <- loglik.exact(po,coef(po))
L1 <- dnorm(obs(po)[1,],mean=exp(coef(po)["th1"]),sd=sd.y.true[1],log=T)
L2 <- dnorm(obs(po)[2,],mean=coef(po)["th2"]*exp(coef(po)["th1"]),sd=sd.y.true[2],log=T)

N1 <- 200
N2 <- 200
th1.vec <- seq(from=-2,to=2,length=N1)
th2.vec <- seq(from=0,to=10,length=N2)
lik.array <- matrix(NA,nr=N1,nc=N2)
for (i1 in 1:N1) {
 for (i2 in 1:N2) {
   th <- c(th1=th1.vec[i1],th2=th2.vec[i2])
   lik.array[i1,i2] <- loglik.exact(po,th)
 }
}

save(m1.out,etime1,
     m2.out,etime2,
     m3.out,m3.trace,etime3,
     mle.exact,lik.array,
     th1.vec,th2.vec,
     file="toy.rda",compress='xz')

if (CLUSTER) stopCluster(cl)

}


## ----fig-toy,cache=TRUE--------------------------------------------------

XLIM <- c(-2,2)
YLIM <- c(0,10)
LINE.YAXIS <- 2
LINE.XAXIS <- 2.5
LEVELS <- c(0,3,10,100,10000000)
X.LABEL <- 0.87
Y.LABEL <- 0.87	

op <- par(mfrow=c(2,2),mai=c(0,0.2,0.2,0),omi=c(0.7,0.4,0,0.1))
CEX.TRIANGLE <- CEX.SQUARE <- 1.5
CEX.POINTS <- 1.5
CEX.LAB <- 1.5
CEX.AXIS <- 1.2

CEX.TRIANGLE <- CEX.SQUARE <- 1
CEX.POINTS <- 1
CEX.LAB <- 1.2
CEX.AXIS <- 0.8
CEX.AXIS.NUMBERS <- 1

########### if1 ###############

plot(x=m1.out["th1",],y=m1.out["th2",],xlim=XLIM,ylim=YLIM,
     xlab='',ylab='',axes=F,type='n')
box()
axis(side=2,cex.axis=CEX.AXIS.NUMBERS)
axis(side=1,labels=F,cex.axis=CEX.AXIS.NUMBERS)
mtext(side=2,bquote(theta[2]),line=LINE.YAXIS,cex=CEX.AXIS)

.filled.contour(x=th1.vec,y=th2.vec,z=(max(lik.array)-lik.array),
                levels=LEVELS,col=c('white','red','orange','yellow'))

points(x=m1.out["th1",],y=m1.out["th2",],cex=CEX.POINTS)
points(x=mle.exact["th1"],y=mle.exact["th2"],pch=17,col="green",cex=CEX.TRIANGLE)
points(x=mle.exact["th1"],y=mle.exact["th2"],pch=2,cex=CEX.TRIANGLE)

plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"A",cex=CEX.LAB)

########### if2 ###############

plot(x=m2.out["th1",],y=m2.out["th2",],xlim=XLIM,ylim=YLIM,
     xlab='',ylab='',axes=F,type='n')
box()
axis(side=2,labels=F,cex.axis=CEX.AXIS.NUMBERS)
axis(side=1,labels=F,cex.axis=CEX.AXIS.NUMBERS)

.filled.contour(x=th1.vec,y=th2.vec,z=(max(lik.array)-lik.array),
                levels=LEVELS,col=c('white','red','orange','yellow'))

points(x=m2.out["th1",],y=m2.out["th2",],cex=CEX.POINTS)
#points(x=m2.out["th1.start",],y=m2.out["th2.start",],pch=3,cex=0.3)
points(x=mle.exact["th1"],y=mle.exact["th2"],pch=17,col="green",cex=CEX.TRIANGLE)
points(x=mle.exact["th1"],y=mle.exact["th2"],pch=2,cex=CEX.TRIANGLE)

plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"B",cex=CEX.LAB)

########### pmcmc final points ######

plot(x=m3.out["th1",],y=m3.out["th2",],xlim=XLIM,ylim=YLIM,xlab='',ylab='',axes=F,type='n')
box()
axis(side=2,cex.axis=CEX.AXIS.NUMBERS)
axis(side=1,cex.axis=CEX.AXIS.NUMBERS)
mtext(side=2,bquote(theta[2]),line=LINE.YAXIS,cex=CEX.AXIS)
mtext(side=1,bquote(theta[1]),line=LINE.XAXIS,cex=CEX.AXIS)

.filled.contour(x=th1.vec,y=th2.vec,z=(max(lik.array)-lik.array),
                levels=LEVELS,col=c('white','red','orange','yellow'))

points(x=m3.out["th1",],y=m3.out["th2",],cex=CEX.POINTS)
#points(x=m.out["th1.start",],y=m.out["th2.start",],pch=3,cex=0.3)
points(x=mle.exact["th1"],y=mle.exact["th2"],pch=17,col="green",cex=CEX.TRIANGLE)
points(x=mle.exact["th1"],y=mle.exact["th2"],pch=2,cex=CEX.TRIANGLE)
#points(x=post.mean["th1"],y=post.mean["th2"],pch=15,col="blue",cex=CEX.SQUARE)
#points(x=post.mean["th1"],y=post.mean["th2"],pch=0,cex=CEX.SQUARE)

plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"C",cex=CEX.LAB)

######### pmcmc marginals #########

m3.trace <- lapply(m3.trace,tail,0.9*nrow(m3.trace[[1]]))

lik.dev <- lik.array - max(lik.array)
post.array <- exp(lik.dev)/sum(exp(lik.dev))
post.th1 <- apply(post.array,1,sum)
post.th2 <- apply(post.array,2,sum)
post.mean <- c(th1=sum(post.th1*th1.vec),th2=sum(post.th2*th2.vec))

bw <- 0.1
k1 <- lapply(m3.trace,
             function(x){
               density(x[,1],bw=bw)
             })

k1x <- sapply(k1,getElement,"x")
k1y <- sapply(k1,getElement,"y")
matplot(k1x[,1:8],k1y[,1:8],lty=1,type='l',main="",axes=F,xlab='',ylab='',xlim=XLIM)

box()
axis(side=1,cex=CEX.AXIS.NUMBERS)
mtext(side=1,bquote(theta[1]),line=LINE.XAXIS,cex=CEX.AXIS)

th1.post <- ksmooth(y=post.th1/(th1.vec[2]-th1.vec[1]),x=th1.vec,bandwidth=0.2)
lines(th1.post,lty="dotted",lwd=2)
abline(h=0)

plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"D",cex=CEX.LAB)

par(op)


## ----parallel-mif,eval=F,echo=F------------------------------------------
## 
## mpar <- foreach(
##                 i=1:JOBS,
##                 .packages='pomp',
##                 .inorder=FALSE
##                 ) %dopar% {
## 
##  set.seed(7777+i)
##  th.draw <- rprior(dacca,coef(dacca))[,1]
##  m <- try(
##           mif(
##               dacca,
##               Nmif=NMIF,
##               Np=NP,
##               start=th.draw,
##               pars=dacca.pars,
##               ivps=dacca.ivps,
##               ic.lag=IC.LAG,
##               method=METHOD,
##               rw.sd=dacca.rw.sd,
##               cooling.type="geometric",
##               cooling.fraction=sqrt(0.1),
##               var.factor=2,
##               transform=TRUE
##               )
##           )
##  ll <- replicate(n=NLIK,logLik(pfilter(m,Np=NPLIK)))
##  list(pomp=m,start=th.draw,ll=ll)
## }
## 

## ----cholera-mif1-mif2,eval=T,echo=F,results='hide'----------------------

binary.file <- "cholera-mif1-mif2.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {


CORES <- 100    ## number of parallel processes
JOBS <- 100     ## number of estimation runs
NMIF <- 100     ## number of IF iterations per estimation run
NP <- 10000     ## number of particles in pfilter operations
NLIK <- 10      ## number of likelihood evaluations
NPLIK <- 20000    ## number of particles in likelihood evaluation
## for mif1 computations:
METHOD1 <- "mif" ## use IF2
IC.LAG1 <- 60 ## fixed smoothing lag for initial conditions
## for mif2 computations:
METHOD2 <- "mif2" ## use IF2
IC.LAG2 <- 1000 ## fixed smoothing lag for initial conditions

if (TEST) {
  CORES <- 10
  JOBS <- 10
  NMIF <- 2
  NP <- 100
  NLIK <- 3
  NPLIK <- 1000
}

if (CLUSTER) {
  require(doSNOW)
  hostlist <- tail(scan(nodefile,what=character(0)),-1)
  cl <- makeCluster(hostlist,type='SOCK')
  registerDoSNOW(cl)
} else {
  require(doParallel)
  registerDoParallel(CORES)
}

pompExample(dacca)

param.tab <- as.matrix(read.table(row.names=1,header=TRUE,text="
                  mle1 box_min box_max
gamma      20.800000000   10.00   40.00
eps        19.100000000    0.20   30.00
rho         0.000000000    0.00    0.00
delta       0.020000000    0.02    0.02
deltaI      0.060000000    0.03    0.60
clin        1.000000000    1.00    1.00
alpha       1.000000000    1.00    1.00
beta.trend -0.004980000   -0.01    0.00
log.beta1   0.747000000   -4.00    4.00
log.beta2   6.380000000    0.00    8.00
log.beta3  -3.440000000   -4.00    4.00
log.beta4   4.230000000    0.00    8.00
log.beta5   3.330000000    0.00    8.00
log.beta6   4.550000000    0.00    8.00
log.omega1 -1.692819521  -10.00    0.00
log.omega2 -2.543383579  -10.00    0.00
log.omega3 -2.840439389  -10.00    0.00
log.omega4 -4.691817993  -10.00    0.00
log.omega5 -8.477972478  -10.00    0.00
log.omega6 -4.390058806  -10.00    0.00
sd.beta     3.130000000    1.00    5.00
tau         0.230000000    0.10    0.50
S.0         0.621000000    0.00    1.00
I.0         0.378000000    0.00    1.00
Rs.0        0.000000000    0.00    0.00
R1.0        0.000843000    0.00    1.00
R2.0        0.000972000    0.00    1.00
R3.0        0.000000116    0.00    1.00
nbasis      6.000000000    6.00    6.00
nrstage     3.000000000    3.00    3.00
"))

dacca.pars <- c("gamma","eps","deltaI","beta.trend","log.beta1","log.beta2", 
                "log.beta3","log.beta4", "log.beta5", "log.beta6", "log.omega1",
                "log.omega2","log.omega3","log.omega4","log.omega5","log.omega6",
                "sd.beta",   "tau")
dacca.ivps <- c("S.0","I.0","R1.0","R2.0","R3.0")
dacca.rw.sd <- c(rep(0.1,length(dacca.pars)),rep(0.2,length(dacca.ivps)))
names(dacca.rw.sd) <- c(dacca.pars,dacca.ivps)

dacca.hyperparams <- list(min=param.tab[,"box_min"],
                          max=param.tab[,"box_max"])

dacca.rprior <- function(params, hyperparams, ...)
{
  r <- runif(
             n=length(hyperparams$min),
             min=hyperparams$min,
             max=hyperparams$max
             )
  names(r) <- names(hyperparams$min)
  return(r)
}

dacca.dprior <- function(params, hyperparams, ..., log = FALSE)
{
  d <- sum(
           dunif(
                 x=params,
                 min=hyperparams$min,
                 max=hyperparams$max,
                 log=TRUE
                 )
           )
  if (log) d else exp(d)
}

dacca <- pomp(dacca,
              rprior=dacca.rprior,
              dprior=dacca.dprior,
              hyperparams=dacca.hyperparams)

METHOD <- METHOD1
IC.LAG <- IC.LAG1
tic <- Sys.time()

mpar <- foreach(
                i=1:JOBS,
                .packages='pomp',
                .inorder=FALSE
                ) %dopar% {
                  
 set.seed(7777+i)
 th.draw <- rprior(dacca,coef(dacca))[,1]
 m <- try(
          mif(
              dacca,
              Nmif=NMIF,
              Np=NP,
              start=th.draw,
              pars=dacca.pars,
              ivps=dacca.ivps,
              ic.lag=IC.LAG,
              method=METHOD,
              rw.sd=dacca.rw.sd,
              cooling.type="geometric",
              cooling.fraction=sqrt(0.1),
              var.factor=2,
              transform=TRUE
              )
          )
 ll <- replicate(n=NLIK,logLik(pfilter(m,Np=NPLIK)))
 list(pomp=m,start=th.draw,ll=ll)
}

toc <- Sys.time()
etime1 <- toc-tic

m1.in <- rbind(sapply(mpar,function(x)x$start))
m1.out <- rbind(sapply(mpar,function(x)coef(x$pomp)))
m1.lik <- rbind(sapply(mpar,function(x)x$ll))

METHOD <- METHOD2
IC.LAG <- IC.LAG2
tic <- Sys.time()

mpar <- foreach(
                i=1:JOBS,
                .packages='pomp',
                .inorder=FALSE
                ) %dopar% {
                  
 set.seed(7777+i)
 th.draw <- rprior(dacca,coef(dacca))[,1]
 m <- try(
          mif(
              dacca,
              Nmif=NMIF,
              Np=NP,
              start=th.draw,
              pars=dacca.pars,
              ivps=dacca.ivps,
              ic.lag=IC.LAG,
              method=METHOD,
              rw.sd=dacca.rw.sd,
              cooling.type="geometric",
              cooling.fraction=sqrt(0.1),
              var.factor=2,
              transform=TRUE
              )
          )
 ll <- replicate(n=NLIK,logLik(pfilter(m,Np=NPLIK)))
 list(pomp=m,start=th.draw,ll=ll)
}

toc <- Sys.time()
etime2 <- toc-tic

m2.in <- rbind(sapply(mpar,function(x)x$start))
m2.out <- rbind(sapply(mpar,function(x)coef(x$pomp)))
m2.lik <- rbind(sapply(mpar,function(x)x$ll))

save(m1.in,m1.out,m1.lik,etime1,
     m2.in,m2.out,m2.lik,etime2,
     file=binary.file,compress='xz')

if (CLUSTER) stopCluster(cl)

}


## ----cholera-compare-----------------------------------------------------
# COMPARES IF1 AND IF2 FOR DACCA
load("cholera-mif1-mif2.rda")

if2.sd <- apply(m2.lik,1,sd)
if2.median <- apply(m2.lik,2,median)

if1.sd <- apply(m1.lik,1,sd)
if1.median <- apply(m1.lik,2,median)

## op <- par(mai=c(0.75,0.75,0.1,0.5))
op <- par(mar=c(3.5,3.5,0.5,0.5))
y <- if2.median
x <- if1.median
#lo=min(c(x,y))
lo <- -6000
hi <- max(c(x,y,-3760))
plot(y=y,x=x,ylab="",xlab="",xlim=c(lo,hi),ylim=c(lo,hi))
CEX.AXIS.LABELS <- 1
mtext("log likelihood for IF2",2,line=2.5,cex=CEX.AXIS.LABELS)
mtext("log likelihood for IF1",1,line=2.5,cex=CEX.AXIS.LABELS)
abline(a=0,b=1)
abline(h=-3763.8,lty="dotted",lwd=2)
abline(v=-3763.8,lty="dotted",lwd=2)
par(op)



