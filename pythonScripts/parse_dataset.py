

def parse_dataset(filepath):
    list = []
    with open(filepath) as f:
        for line in f:
            pdb_id, chain_ids = line.split()
            list.append((pdb_id, chain_ids))
    return list