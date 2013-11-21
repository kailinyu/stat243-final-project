fx <- function(x) dnorm(x,mean=0,sd=2)
# k=5; xLower=-2; xUpper=2; n=10; plotEnv <- FALSE
ARsamp(fx,plotEnv=TRUE)  ### set plotEnv=TRUE to see the envelope function at initilization and updates



###main function for adaptive-rejection sampling
ARsamp <- function(fx,n=100,k=5,xLower=-2,xUpper=2,plotEnv=FALSE){
### fx is the distribution we're sampling from
### n is the number of samples
## k is the number of abscissae to start with
## xLower and xUpper give bounds (define D in paper's notation)
samples <- array(NA,n)

## need to check of log-concave (probably make a function to do this) 
	### output will be true or false-and based on this we either continue with the 
	### function or we exit with an error message saying it needs to be log-concave

## also check that hp(xLower) > 0  &&  hp(xUpper) < 0


## initialize--> pick k values of x for abscissae (function to do this)
	## how many (k) abscissae to use for start? 
	
	out <- initialize(fx,k,xLower,xUpper,plotEnv=plotEnv)
	xa <- out$xa
	za <- out$za
	h <- out$h
	hp <- out$hp
 
## sample
	n.accept <- 0
	while(n.accept < n){
		w <- runif(1)
		xStar <- runif(1,xLower,xUpper) ## FIX THIS LINE! need to replace with call to function that can sample x* form sk(k)
		if(w <= exp(flk(xStar,xa,h)-fuk(xStar,xa,za,h,hp))){
			n.accept <- n.accept + 1
			samples[n.accept] <- xStar
		}else{
			out <- update(xStar,xa,za,h,hp,k,fx,plotEnv=plotEnv)
			xa <- out$xa
			za <- out$za
			h <- out$h
			hp <- out$hp
			k <- out$k
			if(w <= exp(h[out$kStar] - fuk(xStar,xa,za,h,hp)) ){
				n.accept <- n.accept + 1
				samples[n.accept] <- xStar
			}
		}
	}
	return(samples)
}


### auxilary functions

#check log-concavity
# LCcheck <- function(fx,xLower,xUpper) ##return either TRUE or FALSE 

update <- function(xStar,xa,za,h,hp,k,fx,plotEnv=FALSE){
	ins <- which(xStar<xa)[1] - 1 ## ins such that xa[ins] < xStar < xa[ins+1]
	xa <- c(xa[1:ins],xStar,xa[(ins+1):k])
	h <- c(h[1:ins],fh(xStar,fx),h[(ins+1):k])
	hp <- c(hp[1:ins],fhp(xStar,fx),h[(ins+1):k])
	k <- k+1
	za <- array(NA,k)  #### can make this faster by inserting new value
	za[k] <- xa[k] ##see note about this in intialize()
	for(i in 1:(k-1)){
		za[i] <- (h[i+1]-h[i]-xa[i+1]*hp[i+1]+xa[i]*hp[i])/(hp[i]-hp[i+1])
	}
	if(plotEnv){ plot.envelope(xa[1],xa[k],fx,xa,za,h,hp)}
	return(list(xa=xa,za=za,h=h,hp=hp,k=k,kStar=(ins+1)))
}


initialize <- function(fx, k, xLower, xUpper, plotEnv = FALSE){
	xa <- array(NA,k) ##x-coord of abscissae
	xa[1] <- xLower
	xa[k] <- xUpper
	za <- array(NA,k) ##x-coord of intersections of tangents LAST ELEMENT IS XUPPER, WHICH 
	####IS NOT A Z-VALUE. it is there to facilitate calculation of uk & sk when x>za[k-1]
	za[k] <- xUpper ##see above note
	for(i in 1:(k-2)){
		xa[i+1] <- xLower + i*(xUpper-xLower)/(k-1) #linearly interpolate intermediate points
	}
	h <- fh(xa,fx)
	hp <- fhp(xa,fx)
	for(i in 1:(k-1)){
		za[i] <- (h[i+1]-h[i]-xa[i+1]*hp[i+1]+xa[i]*hp[i])/(hp[i]-hp[i+1])
	}
	if(plotEnv){ plot.envelope(xLower,xUpper,fx,xa,za,h,hp)}
	return(list(xa=xa,za=za,h=h,hp=hp))
}

plot.envelope <- function(xLower,xUpper,fx,xa,za,h,hp){
	x.plot <- seq(xLower,xUpper,by=0.05)
	dev.new()
	plot(x.plot,log(fx(x.plot)),type='l',col='black')
	points(xa,log(fx(xa)))
	points(za,fuk(za,xa,za,h,hp))
	lines(x.plot,fuk(x.plot,xa,za,h,hp),col='blue')
	lines(x.plot,flk(x.plot,xa,h),col='red')
}

##function to calculate h = log(f)
fh <- function(x,fx){ log(fx(x)) }
##function to calculate dh/dx = d(log(fx(x)))/dx
fhp <- function(x,fx,eps=0.005){ (log(fx(x+eps))-log(fx(x)))/eps }

## function to calculate upper-hull as function of x (eqn 2 in paper)
fuk <- function(x,xa,za,h,hp){ 
	out <- array(NA,length(x))
	for(i in 1:length(x)){	
		j <- which(x[i]<za)[1] ## j such that za[j-1] < x < za[j]
		out[i] <- h[j] + (x[i]-xa[j])*hp[j]
	}
	return(out)
} 

### function (and helper) to calculate sk(x)  (eqn 3 in paper)
fsk <- function(x,xa,za,h,hp,fuk,s.aux) exp(fuk(x,xa,za,h,hp))/ integrate(s.aux,xLower,xUpper,xa,za,h,hp,fuk)[[1]]
s.aux <- function(x,xa,za,h,hp,fuk) exp(fuk(x,xa,za,h,hp))

## function to calculate lower-hull as function of x  (eqn 4 in paper)
flk <- function(x,xa,h){
	out <- array(NA,length(x))
	for(i in 1:length(x)){
		j <- which(x[i]<xa)[1] - 1 ## j such that xa[j] < x < xa[j+1]
		out[i] <- ((xa[j+1]-x[i])*h[j] + (x[i]-xa[j])*h[j+1]) / (xa[j+1] - xa[j])
	}
	return(out)
}

