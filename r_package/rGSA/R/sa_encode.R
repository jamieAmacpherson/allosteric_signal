#' sa_encode.R
#'
#' Encode a protein structure with the M32K25 structural alphabet.
#' 
#' Given a PDB structure this function will fit the C-alpha PDB coordinates
#' with the M32K25 structural alphabet, selecting the optimal fragment for each
#' successive four C-alpha residue set. The function returns the string of fitted
#' structural alphabet fragments as a string
#'
#'
#' @param pdbfile
#'
#' @return structural alphabet string
#'
#' @export

sa_encode = function(traj.xyz, str.format = "pdb"){
	## fragment letters
	fragment_letters = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y')

	## fragment coordinates
	fragment_coordinates = list(
	matrix( c( c(  2.630,  11.087,-12.054), c(  2.357,  13.026,-15.290), c(  1.365,  16.691,-15.389), c(  0.262,  18.241,-18.694)), byrow = TRUE, nrow = 4, ncol=3), #A
        matrix( c( c(  9.284,  15.264, 44.980), c( 12.933,  14.193, 44.880), c( 14.898,  12.077, 47.307), c( 18.502,  10.955, 47.619)), byrow = TRUE, nrow = 4, ncol=3), #B
        matrix( c( c(-25.311,  23.402, 33.999), c(-23.168,  25.490, 36.333), c(-23.449,  24.762, 40.062), c( -23.266, 27.976, 42.095)), byrow = TRUE, nrow = 4, ncol=3), #C
        matrix( c( c( 23.078,   3.265, -5.609), c( 21.369,   6.342, -4.176), c( 20.292,   6.283, -0.487), c( 17.232,   7.962,  1.027)), byrow = TRUE, nrow = 4, ncol=3), #D
        matrix( c( c( 72.856,  22.785, 26.895), c( 70.161,  25.403, 27.115), c( 70.776,  28.306, 29.539), c( 69.276,  31.709, 30.364)), byrow = TRUE, nrow = 4, ncol=3), #E
        matrix( c( c( 41.080,  47.709, 33.614), c( 39.271,  44.390, 33.864), c( 36.049,  44.118, 31.865), c( 32.984,  43.527, 34.064)), byrow = TRUE, nrow = 4, ncol=3), #F
        matrix( c( c( 59.399,  59.100, 40.375), c( 57.041,  57.165, 38.105), c( 54.802,  54.093, 38.498), c( 54.237,  51.873, 35.502)), byrow = TRUE, nrow = 4, ncol=3), #G
        matrix( c( c( -1.297,  14.123,  7.733), c(  1.518,  14.786,  5.230), c(  1.301,  17.718,  2.871), c( -0.363,  16.930, -0.466)), byrow = TRUE, nrow = 4, ncol=3), #H
        matrix( c( c( 40.106,  24.098, 63.681), c( 40.195,  25.872, 60.382), c( 37.528,  27.160, 58.053), c( 37.489,  25.753, 54.503)), byrow = TRUE, nrow = 4, ncol=3), #I
        matrix( c( c( 25.589,   1.334, 11.216), c( 27.604,   1.905, 14.443), c( 30.853,  -0.042, 14.738), c( 30.051,  -1.020, 18.330)), byrow = TRUE, nrow = 4, ncol=3), #J
        matrix( c( c( 17.239,  71.346, 65.430), c( 16.722,  74.180, 67.850), c( 18.184,  77.576, 67.092), c( 20.897,  77.030, 69.754)), byrow = TRUE, nrow = 4, ncol=3), #K
        matrix( c( c( 82.032,  25.615,  4.316), c( 81.133,  23.686,  7.493), c( 83.903,  21.200,  8.341), c( 81.485,  19.142, 10.443)), byrow = TRUE, nrow = 4, ncol=3), #L
        matrix( c( c( 28.972,  -1.893, -7.013), c( 28.574,  -5.153, -5.103), c( 30.790,  -7.852, -6.647), c( 30.144, -10.746, -4.275)), byrow = TRUE, nrow = 4, ncol=3), #M
        matrix( c( c( -4.676,  72.183, 52.250), c( -2.345,  71.237, 55.105), c(  0.626,  71.396, 52.744), c(  1.491,  72.929, 49.374)), byrow = TRUE, nrow = 4, ncol=3), #N
        matrix( c( c(  0.593,  -3.290,  6.669), c(  2.032,  -2.882,  3.163), c(  4.148,  -6.042,  3.493), c(  7.276,  -4.148,  2.496)), byrow = TRUE, nrow = 4, ncol=3), #O
        matrix( c( c( 29.683,  47.318, 25.490), c( 26.781,  47.533, 27.949), c( 26.068,  51.138, 26.975), c( 27.539,  52.739, 30.088)), byrow = TRUE, nrow = 4, ncol=3), #P
        matrix( c( c( 34.652,  36.550, 18.964), c( 33.617,  37.112, 15.311), c( 32.821,  40.823, 15.695), c( 34.062,  43.193, 12.979)), byrow = TRUE, nrow = 4, ncol=3), #Q
        matrix( c( c(  8.082,  44.667, 15.947), c(  5.161,  46.576, 17.520), c(  5.855,  49.813, 15.603), c(  3.022,  50.724, 13.161)), byrow = TRUE, nrow = 4, ncol=3), #R
        matrix( c( c( 64.114,  65.465, 28.862), c( 63.773,  68.407, 26.422), c( 67.481,  69.227, 26.232), c( 67.851,  68.149, 22.610)), byrow = TRUE, nrow = 4, ncol=3), #S
        matrix( c( c(-18.708,-123.580,-46.136), c(-18.724,-126.113,-48.977), c(-18.606,-123.406,-51.661), c(-14.829,-123.053,-51.400)), byrow = TRUE, nrow = 4, ncol=3), #T
        matrix( c( c( 61.732,  49.657, 35.675), c( 62.601,  46.569, 33.613), c( 65.943,  46.199, 35.408), c( 64.205,  46.488, 38.806)), byrow = TRUE, nrow = 4, ncol=3), #U
        matrix( c( c( 88.350,  40.204, 52.963), c( 86.971,  39.540, 49.439), c( 85.732,  36.159, 50.328), c( 83.085,  37.796, 52.614)), byrow = TRUE, nrow = 4, ncol=3), #V
        matrix( c( c( 23.791,  23.069,  3.102), c( 26.051,  22.698,  6.166), c( 23.278,  21.203,  8.349), c( 21.071,  19.248,  5.952)), byrow = TRUE, nrow = 4, ncol=3), #W
        matrix( c( c(  1.199,   3.328, 36.383), c(  1.782,   3.032, 32.641), c(  1.158,   6.286, 30.903), c(  1.656,   8.424, 34.067)), byrow = TRUE, nrow = 4, ncol=3), #X
        matrix( c( c( 33.001,  12.054,  8.400), c( 35.837,  11.159, 10.749), c( 38.009,  10.428,  7.736), c( 35.586,   7.969,  6.163)), byrow = TRUE, nrow = 4, ncol=3)) #Y

	names(fragment_coordinates) = fragment_letters;


        if(str.format == 'pdb'){
	       ## parse pdb file coordinates into matrices of four sequential Calpha coordinates
	       ## index vector of Calpha atom
               ca.ix = bio3d::atom.select(pdbfile, "calpha");
	       ## matrix of x,y,z coordinates of all Calpha atoms
	       ca.xyz = pdbfile$atom[ca.ix$atom, c("x","y","z")];
               } else if(str.format == 'dcd'){
                ca.xyz = traj.xyz
               } else
        stop(paste("structure extension", str.format, "not supported"));

	## create index matrix to fragment-subset coordinate matrix
	## each matrix column corresponds to a 4-Calpha fragment of the input structure
        nFs = dim(ca.xyz)[1] - 3;
	fs.m = matrix(nrow = 4, ncol = nFs);
	## index list from 1 to (length-3)
	fs.ix = seq(1:(dim(ca.xyz)[1] - 3));
	## assign index list and the following incremental lists
	fs.m[1, ] = fs.ix;
	fs.m[2, ] = fs.ix + 1;
	fs.m[3, ] = fs.ix + 2;
	fs.m[4, ] = fs.ix + 3;

	## fit all input structure fragments to all alphabet fragments
	rmsd.m = apply(fs.m, 2, function(x) {
		sapply(1:length(fragment_coordinates), function(y) {
			kabsch(ca.xyz[x, ], fragment_coordinates[[y]]);
		});
	});
	## vector of minimal rmsd values
	rmsd_min.v = apply(rmsd.m, 2, min);
	## vector of row indices of minimal rmsd values
	rmsd_min.ix = sapply(1:length(rmsd_min.v), function(z) {
		which(rmsd.m[ , z] %in% rmsd_min.v[z]);
	})

	## fragment string
	fragstring = fragment_letters[rmsd_min.ix];
	# return string of SA fragments
	return(fragstring);
}


## encode a dcd trajectory file with the M32K25 structural alphabet
encode_dcd_trajectory = function(traj, num.atoms, parallel.calc = 'TRUE'){

        ## length of the trajectory
        nframes = length(traj[,1]);

        ## append xyz coordinates of each frame to an element in a list
        trajectory = list();

        for (i in seq(from = 1, to = nframes, by = 1)) {

                prog = (i/nframes)*100;
                if(prog %% 10 == 0){
                        print(paste(prog, ' %', sep=''))
                }

                ## assign the xyz coordinates to a data frame
                coords = as.data.frame(matrix(traj[i,], nrow = num.atoms, byrow=T))
                names(coords) = c('x', 'y', 'z')

                ## append the coordinates to a list
                trajectory[[i]] = coords
        }

        ## switch to execute the fragment encoding in parallel
        if(parallel.calc == 'TRUE'){
                # determine the number of cores on the machine
                n_cores = parallel::detectCores() - 1;

                ## initiate a cluster to encode in parallel
                cluster = parallel::makeCluster(n_cores)

                ## export sa_encode to cluster
                parallel::clusterEvalQ(cluster, devtools::load_all('/home/macphej/jm.software/development/allosteric_signal/r_package/rGSA/'))

                ## encode the trajectory with the M32K25 structural alphabet
                ## in parallel
                sa.trajectory = pbapply::pblapply(trajectory, cl = cluster, sa_encode, 'dcd')

        } else if(parallel.calc == 'FALSE'){
                
                ## encode the trajectory with the M32K25 structural alphabet
                sa.trajectory = pbapply::pblapply(trajectory, sa_encode, 'dcd');

        }


        ## number of fragments in the alignment
        num.frags = num.atoms - 3;

        ## format the structural alphabet strings into a matrix
        sa.trajectory.mat = matrix(unlist(sa.trajectory), ncol = num.frags, byrow = FALSE);

        ## write the structural alphabet-encoded trajectory to the disk
        saveRDS(sa.trajectory.mat, paste0(format(Sys.time(), "%Y%m%d_%H%M%S_"), "_sa_trajectory_matrix.rds", sep=""))

        ## return the structural alphabet-encoded trajectory
        return(sa.trajectory.mat);
}



## visualise the structural alphabet alignment as a sequence logo
vis_seq_logo = function(sa_traj){

        ## initialise a vector to contain the structural alphabet fragments
        frags = c()

        ## 
        for(i in seq(from = 1, to = ncol(sa_traj))) {
                frags[i] = paste(sa_traj[i,], collapse = '')
        }
        
        plt = ggseqlogo::ggseqlogo(frags, method = 'prob' )

        print(plt)

}


## split the structural alphabet alignment into regular blocks
split_sa_align = function(sa_traj, nblocks){

        ## if the number of rows in the alignment is a prime number, remove the first row
        rows = nrow(sa_traj);

        if(primes::is_prime(rows) == 'TRUE'){
                sa_traj = sa_traj[-1,]

                print('WARNING: number if rows in the alignment is a prime number. The first encoded snapshot in the trajectory was removed')
        };

        ## again determine the number of rows
        rows = nrow(sa_traj);

        ## determine whether the alignment can be equally spit into the defined
        ## number of blocks
        if(rows %% nblocks != 0) {

                for(i in seq(from=1, to=9, by=1)) {

                        ## search for a lower number of blocks that evenly divides the trajectory
                        new.nblocks = nblocks - i

                        ## break the block size search if the trajectory is evenly divisible by the 
                        ## new block size
                        if(rows %% new.nblocks == 0){
                                break
                        } 
                }

                print(paste('WARNING: fragment-encoded trajectory is not evenly divisible by ', nblocks, '. The trajectory will instead be divided into ', new.nblocks, ' even blocks.'))

        } else if(rows %% nblocks == 0){
                new.nblocks = nblocks
        }

        ## determine the number of rows in each block
        block.rows = rows/new.nblocks

        ## initialise list to append into
        sa.list = list()

        ## generate vector dereferencing the start row for each block in the alignment
        ref.seq = seq(from = 1, to = rows, by = block.rows)

        ## loop over the structural alphabet alignment 
        for(i in c(1: length(ref.seq))) {

                rowstart = ref.seq[i];

                rowend = ref.seq[i] + (block.rows-1);
                #print(paste(rowstart, rowend))
                ## assign block to temporary matrix 
                mat = sa_traj[c(rowstart : rowend), ]

                ## append the block to sa.list
                sa.list[[i]] = mat
        }

        return(sa.list);
}

























