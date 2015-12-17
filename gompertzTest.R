rm(list=ls())
require(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
trials <- 100
ptime <- system.time(
 {
  r <- foreach(icount(trials)) %dopar% {
   require(pomp)
   data(gompertz)
   sim <- simulate(gompertz)
  }
 }
)

save(r,file="sims.Rda")






