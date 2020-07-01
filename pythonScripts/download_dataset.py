import urllib.request

pdbid = "3cx5"
path = "/home/katebrich/Documents/diplomka/PDBe_files/"

url=f'https://www.ebi.ac.uk/pdbe/entry-files/pdb{pdbid}.ent'

req = urllib.request.Request(url)
with urllib.request.urlopen(req) as f:
    responseBody = f.read()

with open(f'{path}/{pdbid}.pdb', 'wb') as file:
    file.write(responseBody)
