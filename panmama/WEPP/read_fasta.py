import gzip
import random
import math
import sys


def read_record(filename):
  label = None
  seq = []
  fp = None
  should_close = True
  
  if hasattr(filename, 'read'):
    fp = filename
    should_close = False
  elif filename == '-':
    fp = sys.stdin
    should_close = False
  elif filename.endswith('.gz'):
    fp = gzip.open(filename, 'rt')
  else:
    fp = open(filename)
  
  while True:
    line = fp.readline()
    if line == '': break
    line = line.rstrip()
    if line.startswith('>'):
      if len(seq) > 0:
        seq = ''.join(seq)
        yield(label, seq)
        label = line[1:]
        seq = []
      else:
        label = line[1:]
    else:
      seq.append(line)
  
  if label is not None:
    yield(label, ''.join(seq))
  
  if should_close:
    fp.close()