# set working directory to 'wd', using 'setwd'
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(ggplot2)
library(dplyr)
library(gridExtra)
theme_set(theme_minimal())

### Read the iterates for three randomly chosen psi and zeta coefficients 
# (all coefficients, including burnin)
psi <- read.csv('./out/psi.csv')
zeta <- read.csv('./out/zeta.csv')
IT <- nrow(psi) # nr of iterations

# ------------------ Trace, acf, hist plots for psi and zeta ----------------------------------------
# We take the first and third column from psi; and the first columns from zeta

# trace plots (all iterates)
psizeta <- data.frame(cbind(psi[-IT,1],psi[-IT,3],zeta[,1]))
names(psizeta) <- c(names(psi)[c(1,3)], names(zeta)[1])
psizeta_st <- stack(psizeta)
psizeta_st$iterate <- rep(1:(IT-1),3)
names(psizeta_st) <- c('value','coefficient','iterate')
p_psizeta_trace <- ggplot(data=psizeta_st, aes(x=iterate,y=value)) + geom_line() + facet_grid(coefficient~.,scales='free')+
  theme(text = element_text(size=15), axis.text.x = element_text(angle=90, hjust=1))+ylab("")

# histogram plots (first half left out, considered burnin)
indi <- round(IT/2):(IT-1)
psizeta_wb <- psizeta[indi,]  # wb = without burnin
psizeta_wb_st <- stack(psizeta_wb)
psizeta_wb_st$iterate <- rep(indi,3)
names(psizeta_wb_st) <- c('value','coefficient','iterate')
p_psizeta_hist <-   ggplot(data=psizeta_wb_st, aes(x=value))+ 
  facet_wrap(~ coefficient,scales='free') +geom_histogram(aes(y=..density..),fill='white',colour='Black')+
     geom_density(alpha=.2, fill="#E69F00")

# autocorrelation plots (first half left out, considered burnin)
LM <- 50 # value for maximal lag
aa1 <- acf(psizeta_wb[,1],plot=FALSE,lag.max=LM)
a1 <- aa1$acf
a2 <- acf(psizeta_wb[,2],plot=FALSE,lag.max=LM)$acf
a3 <- acf(psizeta_wb[,3],plot=FALSE,lag.max=LM)$acf
dAcf <- data.frame(acf=c(a1,a2,a3), lag=rep(aa1$lag,3), 
                   coordinate=rep(c(names(psi)[c(1,3)],names(zeta)[1]),each=LM+1))
p_psizeta_acf <- ggplot(data=dAcf,aes(x=lag,y=acf))+
  geom_hline(aes(yintercept = 0)) +
  geom_line()+geom_area(fill="#E69F00",alpha=0.2)+
  #  geom_segment(mapping = aes(xend = lag, yend = 0))+
  facet_grid(coordinate~.)

### Read the iterates for alpha (all coefficients, including burnin)
alpha <- read.csv('./out/alpha.csv',header=FALSE)
dAlpha <- data.frame(alpha=alpha[,1],iterate=1:nrow(alpha))
pa1 <- ggplot(data=dAlpha, aes(x=iterate,y=alpha))+geom_line()
pa2 <- ggplot(data=dAlpha[indi,], aes(x=alpha))+geom_histogram(aes(y=..density..),fill='white',colour='Black')+ geom_density(alpha=.2, fill="lightsteelblue1") 

# Read breaks (these are in total N points, so there are N bins)
breaks <- as.numeric(read.csv('./out/breaks.csv',header=FALSE)[,1])
N <- length(breaks)-1

# Read observations (to be added as rug plot)
observations <- read.csv('./out/observations.csv',header=FALSE)
names(observations) <- 'x'

# Read parameters fo the posterior gamma distribution on each bin 
# (assuming independent Gamma distr prior)
postInd <- read.csv('./out/post_ind.csv',header=FALSE) # first column is the shape, second is the rate
# Construct a dataframe that contains the posterior average, 2.5% percentile and 97.5% percentile
# xmin = all break points excepts the last one
# xmax = all break points excepts the first one
dInd <- data.frame(xmin=breaks[-(N+1)],xmax=breaks[-1],average=rep(0,N), lower=rep(0,N), upper=rep(0,N))
for (k in 1:N)
{  dInd[k,3] <- postInd[k,1]/postInd[k,2]
   dInd[k,4:5] <- qgamma(c(0.025,0.975),shape=postInd[k,1],rate=postInd[k,2])
}

# construct dataframe with posterior summary measures with the Markov chain prior
dMarkov <- read.csv('./out/summarypsi.csv')
dMarkov$xmin <- breaks[-(N+1)] # all break points excepts the last one
dMarkov$xmax <- breaks[-1]     # all break points excepts the first one

# Read true intensity function on fine grid
dTrue <- read.csv("./out/true_intensity.csv",header=FALSE)
names(dTrue) <- c("x","intensity")

# to thet the step function right, we need to add the rightmost points, 
tMarkov <- dMarkov[,c(5,1)] # first pick xmin, then average
tMarkov[N+1,] <- c(breaks[N+1],tMarkov[N,2])
tInd <- dInd[,c(1,3)]  # first pick xmin, then average
tInd[N+1,] <- c(breaks[N+1],tInd[N,2])

# read marginal likelihood results
dMll <- read.csv("./out/marglikelihood_N.csv",header=FALSE)
names(dMll) <- c("N", "log_marginal_likelihood")

p_dMll <- ggplot(data=dMll, aes(x=N, y=log_marginal_likelihood))+ 
  geom_point(size=1.0)+geom_line()+ylab('log marginal likelihood(N)') +
   scale_x_continuous(breaks = 1:20,limits=c(1,20))
   

#############################################################################################################
#                      plot results and save to 'out' directory  
#############################################################################################################

# psizeta - trace with burnin, acf without burnin

pdf('./out/psi_zeta_tr-acf.pdf',width=8,height=5)
grid.arrange(p_psizeta_trace, p_psizeta_acf, ncol=2)
dev.off()

# psizeta - histogram without burnin
pdf('./out/psi_zeta_hist.pdf',width=7,height=3)
show(p_psizeta_hist)
dev.off()

# alpha - trace with burnin, histogram without burnin
pdf('./out/alpha_tr_hist.pdf',width=6,height=5)
grid.arrange(pa1,pa2,ncol=1)
dev.off()

titel_en_rug <- 1
y_max <- 10.0 # maximal value on the y-axis 

if (titel_en_rug==1)
{
### Plot posterior mean, posterior quantiles + truth  
# Markov chain prior
pMarkov <- ggplot() + geom_rect(data=dMarkov, aes(xmin=xmin,xmax=xmax,ymin = lower, ymax = upper), fill = "lightsteelblue1")+
    geom_line(data=dTrue, aes(x=x,y=intensity),colour='red',linetype = "solid",size=1.5) +
  geom_step(data=tMarkov, aes(x=xmin,y=average),colour='black',size=1.3)+ 
   ggtitle("Markov chain prior")+  theme(plot.title = element_text(hjust = 0.5))  +ylim(0,y_max)+xlab("")+ylab("")+
  geom_rug(data=observations, mapping=aes(x=x), color="black",sides='t')

# Independent gamma prior
pInd <- ggplot() + geom_rect(data=dInd, aes(xmin=xmin,xmax=xmax,ymin = lower, ymax = upper), fill = "lightsteelblue1") +
  geom_line(data=dTrue, aes(x=x,y=intensity),colour='red',linetype = "solid",size=1.5) +
  geom_step(data=tInd, aes(x=xmin,y=average),colour='black',size=1.3)+ 
  ggtitle("Independent gamma prior")+  theme(plot.title = element_text(hjust = 0.5))+xlab("")+ylab("")+ylim(0,y_max)+
geom_rug(data=observations, mapping=aes(x=x), color="black",sides='t') 

} else {
  pMarkov <- ggplot() + geom_rect(data=dMarkov, aes(xmin=xmin,xmax=xmax,ymin = lower, ymax = upper), fill = "lightsteelblue1")+
    geom_line(data=dTrue, aes(x=x,y=intensity),colour='red',linetype = "solid",size=1.5) +
    geom_step(data=tMarkov, aes(x=xmin,y=average),colour='black',size=1.3)+ xlab("")+ylab("")+
      ylim(0,y_max)
  
  # Independent gamma prior
  pInd <- ggplot() + geom_rect(data=dInd, aes(xmin=xmin,xmax=xmax,ymin = lower, ymax = upper), fill = "lightsteelblue1") +
    geom_line(data=dTrue, aes(x=x,y=intensity),colour='red',linetype = "solid",size=1.5) +
    geom_step(data=tInd, aes(x=xmin,y=average),colour='black',size=1.3)+xlab("")+ylab("")+
   ylim(0,y_max)
}

pdf('./out/estimates.pdf',width=8,height=4)
grid.arrange(pInd,pMarkov,ncol=2)
dev.off()

# plot log of marginal likelihood
pdf('./out/marglikelihood.pdf',width=6,height=3)
show(p_dMll)
dev.off()

