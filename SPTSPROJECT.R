rm(list=ls()) #clearing the worksheet
detach(data)
data <- read.csv("C:/Users/andre/Downloads/ratWeightUpdated.csv")
n_params <- 158 * 3 + 3 * 2 + 1
n_rats <- 158

attach(data)

Y<-matrix(weight,nrow=n_rats,byrow=TRUE)
lY<-log(Y)
t<-matrix(week,nrow=n_rats,byrow=TRUE)
rats<-matrix(id,nrow=n_rats,byrow=TRUE)

post.samp.size<-1000
thin<-10

par.num<- n_params

sd.proposal<-rep(1,par.num)
sd.proposal[1:158]<-0.065

sd.proposal[159:316]<-0.19

sd.proposal[168]<-0.2
sd.proposal[180]<-0.202
sd.proposal[195:196]<-0.205
sd.proposal[201]<-0.2
sd.proposal[210]<-0.2
sd.proposal[238]<-0.265
sd.proposal[239:240]<-0.215
sd.proposal[241]<-0.275
sd.proposal[244]<-0.25
sd.proposal[246:250]<-0.27
sd.proposal[249]<-0.29
sd.proposal[252]<-0.275
sd.proposal[253:254]<-0.25
sd.proposal[255:256]<-0.205
sd.proposal[258:262]<-0.25
sd.proposal[264]<-0.22
sd.proposal[265:268]<-0.25
sd.proposal[269]<-0.21
sd.proposal[270]<-0.23
sd.proposal[271]<-0.21
sd.proposal[272]<-0.225
sd.proposal[273:276]<-0.215
sd.proposal[277]<-0.265
sd.proposal[278]<-0.21
sd.proposal[279:284]<-0.25
sd.proposal[281]<-0.21
sd.proposal[285]<-0.21
sd.proposal[286]<-0.26
sd.proposal[287]<-0.262
sd.proposal[289:290]<-0.23
sd.proposal[292:296]<-0.24
sd.proposal[298]<-0.25
sd.proposal[299]<-0.2
sd.proposal[300]<-0.25
sd.proposal[303]<-0.25
sd.proposal[305]<-0.22
sd.proposal[306:307]<-0.24
sd.proposal[309]<-0.225
sd.proposal[310]<-0.205
sd.proposal[312:315]<-0.25

sd.proposal[317:474]<-15
sd.proposal[319]<-12.5
sd.proposal[320]<-4.85
sd.proposal[323]<-5.1
sd.proposal[325]<-6.45
sd.proposal[337]<-2
sd.proposal[340]<-0.68
sd.proposal[342]<-9.65
sd.proposal[344]<-7.25
sd.proposal[352]<-7.3
sd.proposal[360]<-3.75
sd.proposal[364]<-1.32
sd.proposal[365]<-1.18
sd.proposal[372]<-0.9
sd.proposal[374]<-4.9
sd.proposal[380]<-1.25
sd.proposal[475]<-0.1
sd.proposal[476]<-0.12
sd.proposal[477]<-1.4
sd.proposal[478:480]<-0.5
sd.proposal[n_params]<-0.12

chain<-matrix(rep(0,par.num*post.samp.size), byrow= TRUE, nrow=post.samp.size)

lw <- function(t, la, lb, lc) la - exp(-exp(lb)*(t-exp(lc))) # logarithm of weight



likelihood <- function(lY, t, lAi, lbi, lci, ltau) sum(dnorm(lY, mean=lw(t, lAi, lbi, lci),sd= sqrt(1/exp(ltau)),log=TRUE))
condPriorAi <- function(lAi,la,logalpha) sum(dnorm(lAi,mean=la,sd=sqrt(1/exp(logalpha)),log=TRUE))
condPriorbi <- function(lbi,lb,logbeta) sum(dnorm(lbi,mean=lb,sd=sqrt(1/exp(logbeta)),log=TRUE))
condPriorci <- function(lci,lc,loggamma) sum(dnorm(lci,mean=lc,sd=sqrt(1/exp(loggamma)),log=TRUE))
priorla<-function(la) dnorm(la,mean=0,sd=sqrt(100),log=TRUE)
priorlb<-function(lb) dnorm(lb,mean=0,sd=sqrt(100),log=TRUE)
priorlc<-function(lc) dnorm(lc,mean=0,sd=sqrt(100),log=TRUE)
logpriorloggamma<-function(precision) precision*0.1-0.1*exp(precision) #this is the log of the prior of the log precisions, lalpha, lbeta, lgamma,ltau

posterior<-function(parameters){
  lAi<-matrix(parameters[1:158], byrow=FALSE, ncol=14, nrow=n_rats)
  lbi<-matrix(parameters[159:316], byrow=FALSE, ncol=14, nrow=n_rats)
  lci<-matrix(parameters[317:474], byrow=FALSE, ncol=14, nrow=n_rats)
  lk<-likelihood(lY, t, lAi,lbi,lci, parameters[n_params])
  cAi<-condPriorAi(parameters[1:158], parameters[158 * 3 + 1], parameters[158 * 3 + 4])
  cbi <-condPriorbi(parameters[159:316], parameters[158 * 3 + 2], parameters[158 * 3 + 5])
  cci<-condPriorci(parameters[317:474], parameters[158 * 3 + 3], parameters[158 * 3 + 6])
  return(lk+cAi+cbi+cci+priorla(parameters[158 * 3 + 1])+priorlb(parameters[158 * 3 + 2])+priorlc(parameters[158 * 3 + 3])+logpriorloggamma(parameters[158 * 3 + 4])+logpriorloggamma(parameters[158 * 3 + 5])+logpriorloggamma(parameters[158 * 3 + 6])+logpriorloggamma(parameters[n_params]))
}

######################################

print("beginning main loop")

acceptances<-rep(0, n_params)


old<-new<-chain[1,]
for (iter in 2:post.samp.size){
  print(iter)
  for (inner in 1:thin){
    old<-new	
    for (i in 1: n_params){
      updated.pars<-new
      candidate<-old[i]+rnorm(1, mean=0, sd=sd.proposal[i])
      updated.pars[i]<-candidate 
      log.acc.ratio<- posterior(updated.pars)-posterior(new)
      if (log(runif(1))<log.acc.ratio){new[i]<-candidate
      acceptances[i]<-acceptances[i]+1
      } else new[i]<-old[i]
    }
    chain[iter,]<-new
  }    }   
final_chain<-chain
acceptances/(post.samp.size*thin)
post.means<-apply(chain,2,mean)


##############  PLOT THE CHAINS

dev.new()
par(mfrow=c(4, 3))
for (i in 1:12){
  plot(chain[,i],type="l",ylab="Individual Ai")
}
dev.new()
par(mfrow=c(4, 3))
for (i in 159:(159+11)){
  plot(chain[,i],type="l",ylab="Individual bi")
}
dev.new()
par(mfrow=c(4,3))
for (i in 345:(345+11)){
  plot(chain[,i],type="l",ylab="Individual ci")
}
dev.new()
par(mfrow=c(3,3))
for (i in (158 * 3 + 1): n_params){
  plot(chain[,i],type="l",ylab="Population parameters")
}

sequence<-201:post.samp.size
llss<-length(sequence)

trajectories<-array(0,dim=c(n_rats,14,llss))
for (iter in 1:llss) {
  lA<-matrix(chain[sequence[iter],1:158],byrow=FALSE,ncol=14,nrow=n_rats)
  lb<-matrix(chain[sequence[iter],159:316],byrow=FALSE,ncol=14,nrow=n_rats)
  lc<-matrix(chain[sequence[iter],317:474],byrow=FALSE,ncol=14,nrow=n_rats)
  trajectories[,,iter]<-exp(lw(t, lA, lb, lc))
}
medians<-apply(trajectories,c(1,2),median)
par(mfrow=c(4,3))
for (rat in 1:12){
  plot(x=t[rat,], y=Y[rat,])
  lines(x=t[rat,], y=medians[rat,])
}

lAi_bar<-matrix(apply(chain[201:1000, 1:158], 2, mean), byrow=FALSE, ncol=14, nrow=n_rats)
lbi_bar<-matrix(apply(chain[201:1000, 159:316], 2, mean), byrow=FALSE, ncol=14, nrow=n_rats)
lci_bar<-matrix(apply(chain[201:1000, 317:474], 2, mean), byrow=FALSE, ncol=14, nrow=n_rats)
ltau_bar<-mean(chain[201:1000, n_params])
lk_bar<-likelihood(lY, t, lAi_bar, lbi_bar, lci_bar, ltau_bar)
avg_lk<-0
for (iter in 1:llss){
  lA_i<-matrix(chain[sequence[iter],1:158],byrow=FALSE,ncol=14,nrow=n_rats)
  lb_i<-matrix(chain[sequence[iter],159:316],byrow=FALSE,ncol=14,nrow=n_rats)
  lc_i<-matrix(chain[sequence[iter],317:474],byrow=FALSE,ncol=14,nrow=n_rats)
  ltau<-chain[sequence[iter], n_params]
  avg_lk<-avg_lk + likelihood(lY, t, lA_i, lb_i, lc_i, ltau)
}
avg_lk<-avg_lk/llss
pd<- -2*(avg_lk)+2*lk_bar
DIC<- pd-2*(avg_lk)




