import sys

with open(sys.argv[1], 'r') as f:
  for line in f:
    fields = line.strip().split('\t')
    maxhit_family_name, maxhit_family_acc, copy_name = fields[0], fields[1], fields[2]
    maxhit_family_name = maxhit_family_name.split('#')[0]
    maxhit_family_acc = maxhit_family_acc.split('.')[0]

    if copy_name.startswith(maxhit_family_name):
      print(line.strip())
    else:
      coord_info = copy_name.split('::')[1]
      new_copy_name = f'{maxhit_family_name}|{maxhit_family_acc}::{coord_info}'
      print('\t'.join(fields[:2] + [new_copy_name] + fields[3:]))



    