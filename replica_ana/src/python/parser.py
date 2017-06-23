import argparse

parser = argparse.ArgumentParser(description='Calculate the covariance overlap and cosine content of Mutual information matrices.')

parser.add_argument('nmodes', type=int, nargs=1, 
                     help='the number of eigenmodes used for the calculation')

args = parser.parse_args()

print args.nmodes[0]


