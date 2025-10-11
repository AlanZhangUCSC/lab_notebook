import random
import argparse

parser = argparse.ArgumentParser(description='Shuffle input values with a given seed')
parser.add_argument('--seed', type=int, help='random seed for shuffling')
parser.add_argument('input', nargs='+', help='values to shuffle')

args = parser.parse_args()

if args.seed:
  random.seed(args.seed)
inseq = args.input
random.shuffle(inseq)
print(" ".join(inseq))