library(minpack.lm)

## Command line arguments: (1) data 1, (2) data 2, (3) data 3, (4) last time point (ns)
args=commandArgs(TRUE)
stopifnot(! is.na(args[1]))
print(args)

#
dat1=read.table(args[1])
dat2=read.table(args[2])
dat3=read.table(args[3])
lastt = as.numeric(args[4])
outprefix = args[5]

dat=as.data.frame(cbind(dat1,dat2,dat3))

datmeans=apply(dat, 1, mean)
datsds = apply(dat, 1, sd)


time = seq(from=0, to=lastt, length.out=nrow(dat1))

dat = as.data.frame(cbind(time, datmeans, datsds))
names(dat)=c("time", "datmeans", "datsds")

fit.nls = nlsLM(datmeans ~ a + (b-a) * (1-exp(-c*time)),
	data=dat,
	start=list(a=datmeans[1], b=1, c=0.02),
	lower = c(datmeans[1] - 0.01, 0.6, 0.001),
	upper = c(datmeans[1] + 0.01, 0.96, 0.05),
	trace = T,
	algorithm = "port")

        fit.nls.summary <- summary(fit.nls)
        fit.nls.Kd      <- fit.nls.summary$param[1]
        fit.nls.predict <- predict(fit.nls)
        results <- as.data.frame(cbind(dat, fit.nls.predict))
	residuals <- residuals(fit.nls)

dat = cbind(dat, residuals)

cexval=2.2

pdf(file=paste(outprefix, "_fit.pdf", sep=""), width=7, height=7)

layout(matrix(1:2, ncol=1), widths=1, heights=c(2.5, 1.5), respect=FALSE)
par(mar = c(0, 5, 1, 1.2))
plot(dat$time, dat$datmeans, pch=16,
	cex=cexval-0.5, cex.lab=cexval, cex.axis=cexval,
	ylab="Covariance overlap", xaxt = 'n',
	ylim=c(0.35, 1.05), panel.first=grid())

arrows(dat$time, dat$datmeans - dat$datsds, dat$time, dat$datmeans + dat$datsds,
	length=0.025, angle=90, code=3)

lines(dat$time, results$fit.nls.predict , lty=2,col="red",lwd=7)


par(mar = c(5, 5, 0, 1.2))
plot(dat$time, residuals,
	pch=16, cex=cexval-0.5,
	cex.lab=cexval, cex.axis=cexval,
	xlab="Time (ns)", ylab="Residuals",
	ylim = c((min(residuals) + (min(residuals)*0.4)), max(residuals) + (max(residuals)*0.4)),
	panel.first=grid())

abline(0,0, lty=2, lwd=1)


dev.off()

fit.nls.summary

print("t1= ")
print(log(2)/fit.nls.summary$coefficients[3,1])

print("k1=")
print(fit.nls.summary$coefficients[3,1])

