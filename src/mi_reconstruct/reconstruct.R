#! /usr/bin/R

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
# Testing the algorithm
#____________________________________________________________

# Generate a random known sequence to benchmark the algorithm
benchmark.test = function(n.repeats, sequence.length){

	known.sequence = rnorm(sequence.length, mean=0.001, sd=0.001)

	# Calculate the sliding sum of the known sequence
	known.sliding.sum = mean.sliding.window(known.sequence)



	reconstructed.sequence = reconstruct(known.sliding.sum,
		n.window = 4,
		datmean = 0.001,
		datsd = 0.001)


	par(mfrow = c(2,1))
	plot(reconstructed.sequence, type='l')


	for (i in c(1:n.repeats)){
		test.lines = reconstruct(known.sliding.sum, 4, datmean = 0.001, datsd = 0.001)

		lines(test.lines)
	}


	lines(known.sequence, col='red')

	## return (1) residuals and (2) dataframe containing real vs. calculated
	# residuals:
	#out.resids = known.sequence - reconstructed.sequence

	# real vs. calculated dataframe
	out.df = as.data.frame(cbind(known.sequence, reconstructed.sequence))
	lm.df = lm(out.df$reconstructed.sequence ~ out.df$known.sequence)

	plot(out.df)
	abline(lm.df)
}


pdf('benchmark_test_10.pdf')
test.out = benchmark.test(n.repeats = 100, sequence.length = 10)
dev.off()


pdf('benchmark_test_100.pdf')
test.out = benchmark.test(n.repeats = 100, sequence.length = 100)
dev.off()

pdf('benchmark_test_1000.pdf')
test.out = benchmark.test(n.repeats = 100, sequence.length = 1000)
dev.off()

pdf('benchmark_test_2000.pdf')
test.out = benchmark.test(n.repeats = 100, sequence.length = 2000)
dev.off()

