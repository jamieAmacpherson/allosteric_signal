library(minpack.lm)

## Command line arguments: (1) data 1, (2) data 2, (3) data 3, (4) last time point (ns)
args=commandArgs(TRUE)
stopifnot(! is.na(args[1]))
print(args)

## User inputs
dat1=read.table(args[1])	#data 1
dat2=read.table(args[2])	#data 2
dat3=read.table(args[3])	#data 3
lastt = as.numeric(args[4])	#length of simulation (ns)
outprefix = args[5]		#ouput prefix
plateau = as.numeric(args[6])	#predicted plateau of exponential fit


## Calculate the mean and sd of the data sets
dat=as.data.frame(cbind(dat1,dat2,dat3))

datmeans=apply(dat, 1, mean)
datsds = apply(dat, 1, sd)


time = seq(from=0, to=lastt, length.out=nrow(dat1))

dat = as.data.frame(cbind(time, datmeans, datsds))
names(dat)=c("time", "datmeans", "datsds")


## Fit the data with a single-exponential regression
fit.nls = nlsLM(datmeans ~ a + (b-a) * (1-exp(-c*time)),
	data=dat,
	start=list(a=datmeans[1], b=plateau, c=0.02),
	lower = c(datmeans[1] - 0.01, plateau - (plateau*0.2), 0.001),
	upper = c(datmeans[1] + 0.01, plateau + (plateau * 0.2), 0.05),
	trace = T,
	algorithm = "port")

        fit.nls.summary <- summary(fit.nls)
        fit.nls.Kd      <- fit.nls.summary$param[1]
        fit.nls.predict <- predict(fit.nls)
        results <- as.data.frame(cbind(dat, fit.nls.predict))
	residuals <- residuals(fit.nls)

dat = cbind(dat, residuals)


# set the relative text size
cexval=2.2

## determine the yaxis label based on the output prefix defined by the user
if(grepl("entropy", outprefix) == TRUE) {
	yaxtit="Entropy (J/K mol)"
}	else {
		yaxtit="Covariance overlap"
}


## plot the curve with residuals as a subplot
pdf(file=paste(outprefix, "_fit.pdf", sep=""), width=7, height=7)

# arrange the plot and subplot
layout(matrix(1:2, ncol=1), widths=1, heights=c(2.5, 1.5), respect=FALSE)
par(mar = c(0, 5, 1, 1.2))

# plot the main curve
plot(dat$time, dat$datmeans, pch=16,
	cex=cexval-0.5, cex.lab=cexval, cex.axis=cexval,
	ylab=yaxtit, xaxt = 'n',
	ylim = c((min(datmeans) - (min(datmeans)*0.05)), max(datmeans) + (max(datmeans)*0.05)),
	panel.first=grid())

# error bars
arrows(dat$time, dat$datmeans - dat$datsds, dat$time, dat$datmeans + dat$datsds,
	length=0.025, angle=90, code=3)

# regression fit
lines(dat$time, results$fit.nls.predict , lty=2,col="red",lwd=7)

# subplot of the residuals
par(mar = c(5, 5, 0, 1.2))
plot(dat$time, residuals,
	pch=16, cex=cexval-0.5,
	cex.lab=cexval, cex.axis=cexval,
	xlab="Time (ns)", ylab="Residuals",
	ylim = c((min(residuals) + (min(residuals)*0.4)), max(residuals) + (max(residuals)*0.4)),
	panel.first=grid())

# line going through y=0 
abline(0,0, lty=2, lwd=1)


dev.off()

## Print out the fitting parameters
fit.nls.summary

print("t1= ")
print(log(2)/fit.nls.summary$coefficients[3,1])

print("k1=")
print(fit.nls.summary$coefficients[3,1])

