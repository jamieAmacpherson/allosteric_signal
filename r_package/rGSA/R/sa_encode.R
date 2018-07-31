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


 sa_encode = function(pdbfile){

        # fragment letters
        fragment_letters = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y')

 	# fragment coordinates
	fragment_coordinates = list(
 	        matrix(  c( c(2.630,  11.087, -12.054),  c(2.357,  13.026, -15.290),  c(1.365,  16.691, -15.389),  c(0.262,  18.241, -18.694)), byrow = TRUE, nrow = 4, ncol=3), #/*A*/
                matrix(  c( c(9.284,  15.264,  44.980),  c(12.933,  14.193,  44.880),  c(14.898,  12.077,  47.307),  c(18.502,  10.955,  47.619)), byrow = TRUE, nrow = 4, ncol=3), #/*B*/
                matrix(  c( c(-25.311,  23.402,  33.999), c(-23.168,  25.490,  36.333), c( -23.449,  24.762,  40.062), c( -23.266,  27.976,  42.095)), byrow = TRUE, nrow = 4, ncol=3), #/*C*/
                matrix(  c( c(23.078,   3.265,  -5.609), c(21.369,   6.342,  -4.176), c(20.292,   6.283,  -0.487), c(17.232,   7.962,   1.027)), byrow = TRUE, nrow = 4, ncol=3), #/*D*/
                matrix(  c( c(72.856,  22.785,  26.895), c(70.161,  25.403,  27.115), c(70.776,  28.306,  29.539), c(69.276,  31.709,  30.364)), byrow = TRUE, nrow = 4, ncol=3), #/*E*/
                matrix(  c( c(41.080,  47.709,  33.614), c(39.271,  44.390,  33.864), c(36.049,  44.118,  31.865), c(32.984,  43.527,  34.064)), byrow = TRUE, nrow = 4, ncol=3), #/*F*/
                matrix(  c( c(59.399,  59.100,  40.375), c(57.041,  57.165,  38.105), c(54.802,  54.093,  38.498), c(54.237,  51.873,  35.502)), byrow = TRUE, nrow = 4, ncol=3), #/*G*/
                matrix(  c( c(-1.297,  14.123,   7.733), c(1.518,  14.786,   5.230), c(1.301,  17.718,   2.871), c(-0.363,  16.930,  -0.466)), byrow = TRUE, nrow = 4, ncol=3), #/*H*/
                matrix(  c( c(40.106,  24.098,  63.681), c(40.195,  25.872,  60.382), c(37.528,  27.160,  58.053), c(37.489,  25.753,  54.503)), byrow = TRUE, nrow = 4, ncol=3), #/*I*/
                matrix(  c( c(25.589,   1.334,  11.216), c(27.604,   1.905,  14.443), c(30.853,  -0.042,  14.738), c(30.051,  -1.020,  18.330)), byrow = TRUE, nrow = 4, ncol=3), #/*J*/
                matrix(  c( c(17.239,  71.346,  65.430), c(16.722,  74.180,  67.850), c(18.184,  77.576,  67.092), c(20.897,  77.030,  69.754)), byrow = TRUE, nrow = 4, ncol=3), #/*K*/
                matrix(  c( c(82.032,  25.615,   4.316), c(81.133,  23.686,   7.493), c(83.903,  21.200,   8.341), c(81.485,  19.142,  10.443)), byrow = TRUE, nrow = 4, ncol=3), #/*L*/
                matrix(  c( c(28.972,  -1.893,  -7.013), c(28.574,  -5.153,  -5.103), c(30.790,  -7.852,  -6.647), c(30.144, -10.746,  -4.275)), byrow = TRUE, nrow = 4, ncol=3), #/*M*/
                matrix(  c( c(-4.676,  72.183,  52.250), c(-2.345,  71.237,  55.105), c(0.626,  71.396,  52.744), c(1.491,  72.929,  49.374)), byrow = TRUE, nrow = 4, ncol=3), #/*N*/
                matrix(  c( c(0.593,  -3.290,   6.669), c(2.032,  -2.882,   3.163), c(4.148,  -6.042,   3.493), c(7.276,  -4.148,   2.496)), byrow = TRUE, nrow = 4, ncol=3), #/*O*/
                matrix(  c( c(29.683,  47.318,  25.490), c(26.781,  47.533,  27.949), c(26.068,  51.138,  26.975), c(27.539,  52.739,  30.088)), byrow = TRUE, nrow = 4, ncol=3), #/*P*/
                matrix(  c( c(34.652,  36.550,  18.964), c(33.617,  37.112,  15.311), c(32.821,  40.823,  15.695), c(34.062,  43.193,  12.979)), byrow = TRUE, nrow = 4, ncol=3), #/*Q*/
                matrix(  c( c(8.082,  44.667,  15.947), c(5.161,  46.576,  17.520), c(5.855,  49.813,  15.603), c(3.022,  50.724,  13.161)), byrow = TRUE, nrow = 4, ncol=3), #/*R*/
                matrix(  c( c(64.114,  65.465,  28.862), c(63.773,  68.407,  26.422), c(67.481,  69.227,  26.232), c(67.851,  68.149,  22.610)), byrow = TRUE, nrow = 4, ncol=3), #/*S*/
                matrix(  c( c(-18.708,-123.580, -46.136), c(-18.724,-126.113, -48.977), c(-18.606,-123.406, -51.661), c(-14.829,-123.053, -51.400)), byrow = TRUE, nrow = 4, ncol=3), #/*T*/
                matrix(  c( c(61.732,  49.657,  35.675), c(62.601,  46.569,  33.613), c(65.943,  46.199,  35.408), c(64.205,  46.488,  38.806)), byrow = TRUE, nrow = 4, ncol=3), #/*U*/
                matrix(  c( c(88.350,  40.204,  52.963), c(86.971,  39.540,  49.439), c(85.732,  36.159,  50.328), c(83.085,  37.796,  52.614)), byrow = TRUE, nrow = 4, ncol=3), #/*V*/
                matrix(  c( c(23.791,  23.069,   3.102), c(26.051,  22.698,   6.166), c(23.278,  21.203,   8.349), c(21.071,  19.248,   5.952)), byrow = TRUE, nrow = 4, ncol=3), #/*W*/
                matrix(  c( c(1.199,   3.328,  36.383), c(1.782,   3.032,  32.641), c(1.158,   6.286,  30.903), c(1.656,   8.424,  34.067)), byrow = TRUE, nrow = 4, ncol=3), #/*X*/
                matrix(  c( c(33.001,  12.054,   8.400), c(35.837,  11.159,  10.749), c(38.009,  10.428,   7.736), c(35.586,   7.969,   6.163)), byrow = TRUE, nrow = 4, ncol=3)  #/*Y*/
                )

        # parse pdb file coordinates into matrices of four sequential C-alpha coordinates

        # select the calpha atoms
        ca.inds = atom.select(pdbfile, "calpha")

        # split the pdb file into xyz coordinates
        #x_coords = pdb$atom[ca.inds$atom,9 ]
        #y_coords = pdb$atom[ca.inds$atom,10 ]
        #z_coords = pdb$atom[ca.inds$atom,11 ]

        x_coords = pdb$atom[,9]
        y_coords = pdb$atom[,10]
        z_coords = pdb$atom[,11]

        # convert pdb coordinates into a list of matrices
        pdb_coordinates = list()

        natoms = length(x_coords)
        nfrags = natoms-3

        for(i in seq(1, nfrags)) {

                pdb_coordinates[[i]] = matrix( c(
                        x_coords[i:(i+3)],
                        y_coords[i:(i+3)],
                        z_coords[i:(i+3)]),
                        nrow = 4, ncol=3
                )
        }


        # perform the fitting procedure 
        fragstring = c()

        for(i in seq(1, nfrags)) {

                # initialise a vector to contain all of the
                # rmsd values
                rmsdvec = c()

                # loop over all possible vectors
                for(w in seq(1, length(fragment_coordinates))) {

                        # determine the fitting for the fragment upto the pdb coordinates
                        krmsd = kabsch(pdb_coordinates[[i]], fragment_coordinates[[w]])
                        #print(krmsd)

                        # append the fitted coordinates into rmsdvec
                        rmsdvec = append(rmsdvec, krmsd)

                }

                # plot the rmsd for each fragment fit
                pdf(paste('resid_', i, '_fragment_fit.pdf', sep=''))
                plot(rmsdvec, type='o')
                dev.off()

                # find the best SA fragment and append it

                minstr = fragment_letters[which.min(rmsdvec)]

                print(paste(min(rmsdvec), minstr))

                fragstring = append(fragstring, minstr)

        }

        # return string of SA fragments
        return(fragstring)

 }





#######
# debugging
#######

debugf = function(){

        fragstring = c()

        for(i in seq(1, length(fragment_coordinates))) {

                fitted_fragments = list() 

                rmsdlist = list()

                rmsdvec = c()

                # loop over all possible vectors
                for(w in seq(1, length(fragment_coordinates))) {

                        # determine the fitting for the fragment upto the pdb coordinates
                        ktmp = kabsch(fragment_coordinates[[i]], fragment_coordinates[[w]])



                        # append the fitted coordinates into rmsdvec
                        rmsdvec = append(rmsdvec, ktmp)


                }

                minstr = fragment_letters[which.min(rmsdvec)]
                print(paste(min(rmsdvec), minstr))

                fragstring = append(fragstring, minstr)
        }

}















