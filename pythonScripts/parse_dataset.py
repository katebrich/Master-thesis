

def parse_dataset(filepath):
    list = []
    with open(filepath) as f:
        for line in f:
            (key, value) = line.split()
            list.append((key, value))
    return list