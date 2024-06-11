from sys import argv
import ast
import pandas as pd
from subprocess import run
from os import system
from sys import exit
#print(argv)

# the taxids supploed here are the taxIDs from Maral. Some of those are species taxIDs, some strain, some subspecies.
# I'm coming up with a hacky way to map them all back to species taxIDs, which we in turn need to map them back to mOTUsIDs.

system("ncbi-taxonomist map --taxids 55565,1655,349741,1120979,411901,272559,411476,1073351,226186,411479,435590,997892,367928,391904,537007,521097,411902,445974,411468,411903,470146,411461,1169286,679895,515620,411469,657318,411483,546270,29391,525325,324831,47714,299033,525362,709991,435591,411477,537011,536231,411470,1000570,28037,1304,1305 > /Users/karcher/a")  

with open("/Users/karcher/a", "r") as f:
  ff = f.readlines()
  ff = [x.rstrip() for x in ff]


resui = []
for s in ff:
  di = ast.literal_eval(s)
  origID = di['taxon']['taxid']
  nameOrig = di['taxon']['name']
  species = " ".join(di['taxon']['name'].split(" ")[0:2])
  rank = di['taxon']['rank']
  if rank == "strain":
    if 'Bifidobacterium longum' in nameOrig:
      taxID = 216816
    elif 'Limosilactobacillus reuteri' in nameOrig:
      taxID = 1598      
    else:
      taxID = di['taxon']['parentid']
  elif rank == 'species':
    taxID = di['taxon']['taxid']
  else:
    if "Escherichia coli" in nameOrig:
      taxID = 562
    elif "Lacticaseibacillus paracasei" in nameOrig:
      taxID = 1597
  resui.append([nameOrig, species, rank, origID, taxID])

speciesTaxIDString = ",".join([str(x[4]) for x in resui])

# To ensure that all returned taxIDs (in taxID) are species-level taxIDs, run this command.
system('ncbi-taxonomist map --taxids ' + speciesTaxIDString + '  | grep \'rank":"species\' | wc -l')


## Step 2: Traverse further fro species to (hopefully genus)
# To (hoepfully) get genus-level from upstream, run this
system('ncbi-taxonomist map --taxids ' + speciesTaxIDString + ' > /Users/karcher/b')

with open("/Users/karcher/b", "r") as f:
  ff = f.readlines()
  ff = [x.rstrip() for x in ff]


resui2 = []
for s in ff:
  di = ast.literal_eval(s)
  origID = di['taxon']['taxid']
  nameOrig = di['taxon']['name']
  species = " ".join(di['taxon']['name'].split(" ")[0:2])
  rank = di['taxon']['rank']
  if "Streptococcus anginosus" in nameOrig:
    taxID = 1301
  else:
    taxID = di['taxon']['parentid']
  resui2.append([nameOrig, species, rank, origID, taxID])


# Ensure that all returned taxIDs are genus-level taxids
genusTaxIDString = [str(x[4]) for x in resui2]
genusTaxIDStringSet = ",".join(list(set(genusTaxIDString)))

# Except for Eubacterium rectale, everything can be assigned a genus :)
system('ncbi-taxonomist map --taxids ' + genusTaxIDStringSet + '  | grep \'rank":"genus\' | wc -l')
system('ncbi-taxonomist map --taxids ' + genusTaxIDStringSet + '  | grep -v \'rank":"genus\'')
system('ncbi-taxonomist map --taxids ' + genusTaxIDStringSet + '  > /Users/karcher/c' )

with open("/Users/karcher/c", "r") as f:
  ff = f.readlines()
  ff = [x.rstrip() for x in ff]


resui3 = []
for s in ff:
  di = ast.literal_eval(s)
  origID = di['query']
  genus = di['taxon']['name']
  resui3.append([origID, genus])



out = pd.DataFrame(resui)
out.columns = ['Original name', "species", "originalRank", "originaltaxID", 'speciesTaxID']
out.to_csv("/Users/karcher/straintaxid_speciestaxid.csv")

out2 = pd.DataFrame(resui2)
out2.columns = ['Original name', "species", "originalRank", "originaltaxID", 'genusTaxID']
tmp = pd.DataFrame(resui3)
tmp.columns = ['genusTaxID', "genusName"]
tmp['genusTaxID'] = tmp['genusTaxID'].astype('int')
out2 = pd.merge(out2, tmp, on = 'genusTaxID')

out2.to_csv("/Users/karcher/straintaxid_genustaxid.csv")

