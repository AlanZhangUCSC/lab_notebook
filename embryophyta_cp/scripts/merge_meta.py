import sys


class Name:
  def __init__(self, name):
    self.name = name
    self.genus = None
    self.species = None
    self.uniden_id = None
    self.subspecies = None
    self.varietas = None

    self.parse_name()
  
  def parse_name(self):
    fields = self.name.split(' ')
    if len(fields) == 1: self.raise_error_and_exit()

    if ' x ' in self.name:
      self.raise_error_and_exit()
    else:
      subsp = ' subsp. ' in self.name
      var = ' var. ' in self.name
      uniden = ' sp. ' in self.name
      aff = ' aff. ' in self.name
      cf = ' cf. ' in self.name
      if len(fields) == 2:
        self.genus = fields[0].strip('[]')
        self.species = fields[1]
        self.check_genus_species_name(self.genus, self.species)
      else:
        if sum([subsp, var, uniden, aff, cf]) == 0:
          self.genus = fields[0].strip('[]')
          self.species = fields[1]
          self.check_genus_species_name(self.genus, self.species)
          print(f'Warning: {self.name} resolved to {self.genus} {self.species}', file=sys.stderr)
        elif sum([subsp, var, uniden, aff, cf]) == 1:
          if   subsp:  self.genus, self.species, self.subspecies = self.process_subsp(self.name)
          elif var:    self.genus, self.species, self.varietas = self.process_var(self.name)
          elif uniden: self.genus, self.species, self.uniden_id = self.process_uniden(self.name)
          elif aff:    self.genus, self.species = self.process_aff(self.name)
          elif cf:     self.genus, self.species = self.process_cf(self.name)
        elif sum([subsp, var, uniden, aff, cf]) > 1:
          self.raise_error_and_exit()
 

  def process_subsp(self, _name):
    fields = _name.split(' ')
    if fields[2] != 'subsp.': self.raise_error_and_exit()
    genus, species, _, subspecies = fields[0:4]
    genus = genus.strip('[]')
    self.check_genus_species_name(genus, species)
    return genus, species, subspecies
  
  def process_var(self, _name):
    fields = _name.split(' ')
    if fields[2] != 'var.': self.raise_error_and_exit()
    genus, species, _, varietas = fields[0:4]
    genus = genus.strip('[]')
    self.check_genus_species_name(genus, species)
    return genus, species, varietas
  
  def process_uniden(self, _name):
    fields = _name.split(' ')
    if fields[1] != 'sp.': self.raise_error_and_exit()
    genus = fields[0]
    genus = genus.strip('[]')
    uniden_id = ' '.join(fields[2:])
    self.check_genus_species_name(genus, 'sp.')
    return genus, 'sp.', uniden_id
  
  def process_aff(self, _name):
    fields = _name.split(' ')
    if fields[1] != 'aff.': self.raise_error_and_exit()
    genus = fields[0]
    genus = genus.strip('[]')
    self.check_genus_species_name(genus, 'sp.')
    return genus, 'sp.'
  
  def process_cf(self, _name):
    fields = _name.split(' ')
    if fields[1] != 'cf.': self.raise_error_and_exit()
    genus = fields[0]
    genus = genus.strip('[]')
    self.check_genus_species_name(genus, 'sp.')
    return genus, 'sp.'

  def check_genus_species_name(self, genus, species):
    if not (
      len(genus) >= 2 and
      len(species) >= 2 and
      genus[0].isupper() and
      species[0].islower()
    ): self.raise_error_and_exit()
  
  def raise_error_and_exit(self):
    print(f'Warning: {self.name} unresolved... Please handle at Name level', file=sys.stderr)
    exit(1)

class FullName:
  def __init__(self, name):
    self.name = name
    self.name1 = None
    self.name2 = None
    self.hybrid_type = None
    self.parse_name()

  def parse_name(self):
    hybrid_fields = self.name.split(' x ')
    if len(hybrid_fields) == 1:
      if 'hybrid cultivar' in hybrid_fields[0]:
        self.hybrid_type = 'hybrid_cultivar'
        genus = hybrid_fields[0].split(' ')[0]
        self.name1 = Name(f'{genus} sp.')
      else:
        self.name1 = Name(hybrid_fields[0])
    elif len(hybrid_fields) == 2:
      par1, par2 = hybrid_fields
      par1_fields = par1.split(' ')
      par2_fields = par2.split(' ')
      if len(par1_fields) == 1 and len(par2_fields) == 1:
        self.name1 = Name(f'{par1_fields[0]} {par2_fields[0]}')
        self.hybrid_type = 'nothospecies'
      else:
        self.name1 = Name(par1)
        self.name2 = Name(par2)
        if self.name1.genus == self.name2.genus:
          self.hybrid_type = 'intragenic'
        else:
          self.hybrid_type = 'intergenic'
      
    elif len(hybrid_fields) == 3:
      par1, par2, par3 = hybrid_fields
      par1_fields = par1.split(' ')
      par2_fields = par2.split(' ')
      par3_fields = par3.split(' ')
      if len(par1_fields) == 2 and len(par2_fields) == 1 and len(par3_fields) == 1:
        self.name1 = Name(par1)
        self.name2 = Name(f'{par2_fields[0]} {par3_fields[0]}')
        if self.name1.genus == self.name2.genus:
          self.hybrid_type = 'intragenic_nothospecies'
        else:
          self.hybrid_type = 'intergenic_nothospecies'
      elif len(par1_fields) == 1 and len(par2_fields) == 1 and len(par3_fields) == 2:
        self.name1 = Name(f'{par1_fields[0]} {par2_fields[0]}')
        self.name2 = Name(par3)
        if self.name1.genus == self.name2.genus:
          self.hybrid_type = 'intragenic_nothospecies'
        else:
          self.hybrid_type = 'intergenic_nothospecies'
      else:
        self.raise_error_and_exit()
    else:
      self.raise_error_and_exit()



def main():
  metadata_file = sys.argv[1]
  taxonomy_file = sys.argv[2]

  taxid_to_taxonomy = {}
  with open(taxonomy_file, "r") as f:
    for line in f:
      taxid, scientific_name, king, phyl, clss, ordr, fmly, gnus = line.strip().split("\t")
      name = FullName(scientific_name)

      taxid_to_taxonomy[taxid] = {
        "name": name,
        "king": king,
        "phyl": phyl,
        "clss": clss,
        "ordr": ordr,
        "fmly": fmly,
        "gnus": gnus,
      }

  with open(metadata_file, "r") as f:
    for line in f:
      accession, taxid, organism, slen, create_date, update_date, title = line.strip().split("\t")
      if taxid not in taxid_to_taxonomy:
        print(f'Warning: {taxid} not found in taxonomy', file=sys.stderr)
        continue
      name = taxid_to_taxonomy[taxid]["name"]
      king = taxid_to_taxonomy[taxid]["king"]
      phyl = taxid_to_taxonomy[taxid]["phyl"]
      clss = taxid_to_taxonomy[taxid]["clss"]
      ordr = taxid_to_taxonomy[taxid]["ordr"]
      fmly = taxid_to_taxonomy[taxid]["fmly"]
      gnus = taxid_to_taxonomy[taxid]["gnus"]
      assert(organism == name.name)
      print(f'{accession}\t{taxid}\t{slen}\t{king}\t{phyl}\t{clss}\t{ordr}\t{fmly}\t{gnus}\t{organism}\t', end='')
      print(f'{name.hybrid_type}\t{name.name1.genus}\t{name.name1.species}\t{name.name1.subspecies}\t{name.name1.varietas}\t{name.name1.uniden_id}\t', end='')
      if name.name2:
        print(f'{name.name2.genus}\t{name.name2.species}\t{name.name2.subspecies}\t{name.name2.varietas}\t{name.name2.uniden_id}\t{title}')
      else:
        print(f'None\tNone\tNone\tNone\tNone\t{title}')

if __name__ == "__main__":
  main()