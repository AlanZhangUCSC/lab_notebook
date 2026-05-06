import sys

in_hmm = sys.argv[1]
out_prefix = sys.argv[2]

passed_info = f'{out_prefix}.passed.info.tsv'
failed_info = f'{out_prefix}.failed.info.tsv'
hmm_max_tblout = f'{out_prefix}.hmm_max.tblout'

passed_info_fh = open(passed_info, "w")
failed_info_fh = open(failed_info, "w")
hmm_max_tblout_fh = open(hmm_max_tblout, "w")

curcopyid = None
curannfam = None
curmaxfam = set()
curmaxbit = -1
curannbit = -1
curmax_lines= []
with open(in_hmm, "r") as f:
  for line in f:
    fields = line.strip().split("\t")
    copyid, family, bitscore = fields[2], fields[15], float(fields[13])
    if copyid != curcopyid:
      if curcopyid is not None:
        if curannfam in curmaxfam:
          print(f'{curcopyid}\t{curannfam}\t{curmaxbit}', file=passed_info_fh)
        else:
          print(f'{curcopyid}\t{curannfam}\t{curannbit}\t{','.join(curmaxfam)}\t{curmaxbit}', file=failed_info_fh)
        for curmax_line in curmax_lines:
          print(curmax_line, file=hmm_max_tblout_fh)
      curcopyid = copyid
      curmaxfam = set([family])
      curmaxbit = bitscore
      curannfam = copyid.split('|')[0]
      curannbit = -1
      curmax_lines = [line.strip()]
      if family == curannfam: curannbit = bitscore
    else:
      if bitscore > curmaxbit:
        curmaxfam = set([family])
        curmax_lines = [line.strip()]
        curmaxbit = bitscore
      elif bitscore == curmaxbit:
        curmaxfam.add(family)
        curmax_lines.append(line.strip())
      if family == curannfam: curannbit = bitscore

if curcopyid is not None:
  if curannfam in curmaxfam:
    print(f'{curcopyid}\t{curannfam}\t{curmaxbit}', file=passed_info_fh)
  else:
    print(f'{curcopyid}\t{curannfam}\t{curmaxbit}\t{','.join(curmaxfam)}', file=failed_info_fh)
  for curmax_line in curmax_lines:
    print(curmax_line, file=hmm_max_tblout_fh)
    
passed_info_fh.close()
failed_info_fh.close()
hmm_max_tblout_fh.close()