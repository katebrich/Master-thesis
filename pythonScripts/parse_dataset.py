

def parse_dataset(filepath):
    list = []
    with open(filepath) as f:
        for line in f:
            pdb_id, chain_ids = line.split()
            chain_ids = chain_ids.upper()
            list.append((pdb_id, chain_ids))
    return list

def parse_prank_dataset(filepath):
    pdb_files = []
    with open(filepath) as f:
        for line in f:
            #todo jaka je presna struktura .ds souboru v pranku? nechybi mi neco?
            if (line[0] != '#'): #skip comments
                words = line.split()
                if (len(words) > 0 and words[0] != 'HEADER:'): # skip empty lines and header
                    pdb_files.append(words[0])

    return pdb_files