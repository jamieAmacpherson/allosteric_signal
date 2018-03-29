#! /usr/bin/R
library(foreach)
library(doMC)

## Register the number of CPUs
registerDoMC(4)


# log-axis minor ticks function
# (https://stackoverflow.com/questions/6955440/displaying-minor-logarithmic-ticks-in-x-axis-in-r)
minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx,...){

  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]

  major.ticks <- pretty(lims,n=5)
  if(missing(mn)) mn <- min(major.ticks)
  if(missing(mx)) mx <- max(major.ticks)

  major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]

  labels <- sapply(major.ticks,function(i)
            as.expression(bquote(10^ .(i)))
          )
  axis(ax,at=major.ticks,labels=labels, cex.axis=1.7,...)

  n <- n+2
  minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
  minors <- minors[-c(1,n)]

  minor.ticks = c(outer(minors,major.ticks,`+`))
  minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]


  axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}


# reconstruct a sequence given a sliding window average
reconstruct = function(sumseq, n.window, datmean, datsd, ...){

	# Initialise unknown sequence
	unknown.seq = c()

	# Hardcode sumseq for now
	known.sumseq = sumseq

	# Hardcode sum window length
	n = n.window

	# Assign positions n - 1
	assign.n = n - 1

	for (i in c(1:assign.n)){

		# initialisation value to start the sequence
		init.val = known.sumseq[1]/n

		# sample a variance from the initialised value
		init.val = init.val + rnorm(1, mean=0, sd=(init.val/5))

		unknown.seq = append(unknown.seq, init.val)
	}



	print(paste('First (n-1) random assignments:', unknown.seq, sep=' '))

	# solve subsequent positions

	for (positions in c(n:(length(known.sumseq) + assign.n) )){

		pos.val =  known.sumseq[positions - assign.n] -
			unknown.seq[positions - 1] -
			unknown.seq[positions - 2] - 
			unknown.seq[positions - 3]

		unknown.seq = append(unknown.seq, pos.val)

		print(pos.val)

	}

	return(unknown.seq)

}



mean.sliding.window = function(datain){
	# Function to compute the sliding-window mean of a vector of numbers
	# INPUT: vector of numbers
	# OUTPUT: sliding window average (n=4)

	# loop through the vector of numbers to calculate the running mean of 
	# the sliding window
	avout = c()

	for(i in c(1 : (length(datain) - 3))){

		# sum over length n=4
		tmp.mean = sum(c(datain[i], datain[i+1], datain[i+2], datain[i+3]))
		
		# append the mean value to `avout`
		avout = append(avout, tmp.mean)

	}

	return(avout)
}

#____________________________________________________________
# Compute the numerical convergence of the solution
#____________________________________________________________
convergence.plt = function(iteration.data){
	# compute convergence of the solution to the numerical problem
	# of assigning positional values, given sliding window sums
	# (see schematic for details of the numerical problem)
	# INPUT: dataframe output of the computed solution
	# OUTPUT: plot showing the convergence of the numerical solution

	n.iterations = ncol(iteration.data)

	# initialise a vector to catch the output of the convergence 
	# calculation
	convergence.out = c()

	# loop through the iterations of the numerical solution
	convergence.out = foreach(i=1:n.iterations, .combine=rbind) %dopar% {
		print(paste(i, n.iterations, sep='/'))
		# calculate the difference between the average solution and the next
		# iteration

		if(i == 1){
			prev.average = iteration.data[,1]
			curr.average = iteration.data[,1]
		}
		if(i == 2){
			prev.average = iteration.data[,1]
			curr.average = apply(iteration.data[,c(1:2)], 1, mean)
		}
		if(i > 2){
			prev.average = apply(iteration.data[,c(1:(i-1))], 1, mean)
			curr.average = apply(iteration.data[,c(1:(i))], 1, mean)
		}
		
		diff = (curr.average - prev.average)/prev.average

		# compute the mean difference
		average.diff = mean(diff)

		# if the difference is negative, convert to positive value
		if(average.diff < 0){
			average.diff = average.diff * -1
		}

		# append the result to the output vector
		#convergence.out = append(convergence.out, average.diff)

		

	}

	iterations = seq(from=1, to=length(convergence.out))
	par(mar=c(5,5,2,2))
	plot(log10(iterations),
		convergence.out,
		type = 'l',
		xlab = 'Iteration',
		ylab = 'Fractional difference',
		xaxt='n',
		cex.axis=1.7,
		cex.lab=1.7)

	minor.ticks.axis(1,9,mn=0, mx=round(log10(n.iterations)))
}

#____________________________________________________________
# Testing the algorithm
#____________________________________________________________

# Generate a random known sequence to benchmark the algorithm
benchmark.test = function(n.repeats, sequence.length, plt.switch, ...){
	## Benchmark test 

	known.sequence = rnorm(sequence.length, mean=0.001, sd=0.001)

	# Calculate the sliding sum of the known sequence
	known.sliding.sum = mean.sliding.window(known.sequence)



	reconstructed.sequence = reconstruct(known.sliding.sum,
		n.window = 4,
		datmean = 0.001,
		datsd = 0.001)

	if(plt.switch == 'YES'){
		par(mfrow = c(2,1),
			mar=c(5,5,2,2))
		plot(reconstructed.sequence,
			type='l',
			cex.axis = 1.7,
			cex.lab = 1.7,
			...)
	}
	

	# catch the results of each of the iterations
	iteration.results = c()

	# loop through a defined number of iterations, solving for the positional 
	# sums (see schematic for a description of the numerical problem)
	iteration.results = foreach(i=1:n.repeats, .combine=cbind) %dopar% {
		test.lines = reconstruct(known.sliding.sum, 4, datmean = 0.001, datsd = 0.001)
		
		# catch results
		iteration.results = cbind(iteration.results, test.lines)

		# plot the results of the iteration
		
	}

	if(plt.switch == 'YES'){

		apply(iteration.results, 2, lines)
	}

	# compute the mean reconstructed sequence
	average.result = apply(iteration.results, 1, mean)

	# compute the standard deviation of the result
	sd.result = apply(iteration.results, 1, sd)

	# plot a trace of the 'known sequence'
	if(plt.switch == 'YES'){
		lines(known.sequence, col='red')
	}

	# real vs. calculated dataframe
	#out.df = as.data.frame(cbind(known.sequence, average.result))
	lm.df = lm(average.result ~ known.sequence)

	if(plt.switch == 'YES'){
		# plot the correlation between known sequence and the average result
		plot(known.sequence, average.result,
			cex.axis = 1.7,
			cex.lab = 1.7)
	
		# plot the standard deviation of the computation as error bars
		arrows(known.sequence, average.result - sd.result,
			known.sequence, average.result + sd.result,
			length=0.05,
			angle=90,
			code=3)

		abline(lm.df)
	}	

	return(iteration.results)
}




pdf('benchmark_test_10.pdf')
test.out = benchmark.test(n.repeats = 10000,
	sequence.length = 10,
	plt.switch='YES',
	ylim=c(-0.005, 0.005))
#convergence.plt(test.out)
dev.off()

pdf('benchmark_test_100.pdf')
test.out = benchmark.test(n.repeats = 10000,
	sequence.length = 100,
	plt.switch='YES',
	ylim=c(-0.005, 0.005))
#convergence.plt(test.out)
dev.off()


pdf('benchmark_test_100_convergence.pdf')
test.out = benchmark.test(n.repeats = 50000, sequence.length = 100, plt.switch='NO')
convergence.plt(test.out)
dev.off()

pdf('benchmark_test_1000.pdf')
test.out = benchmark.test(n.repeats = 1000, sequence.length = 1000)
convergence.plt(test.out)
dev.off()

pdf('benchmark_test_2000.pdf')
test.out = benchmark.test(n.repeats = 1000, sequence.length = 2000)
convergence.plt(test.out)
dev.off()

