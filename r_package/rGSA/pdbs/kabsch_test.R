# test

# initialise test matrices
a = matrix(c(c(1,2,3),c(4,5,6),c(7,8,9)),3,3)

b = matrix(c(c(4,6,9),c(5,6,9),c(7,8,9)),3,3)

kabsch_test = function(a, b){
  # compute the rmsd before fit

  av = as.vector(a)
  bv = as.vector(b)

  print('rmsd before kabsch: ')
  tmp1 = rmsd(av, bv)
  print(tmp1)


  print('rmsd after kabsch: ')
  k = kabsch(a,b)

  kv = as.vector(k)

  tmp2 = rmsd(av,kv)
  print(tmp2)
}




kabsch_test(a,b)