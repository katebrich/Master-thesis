
hydropathy_kyte_doolitle = {
    'A': 1.8,   #alanine
    'C': 2.5,   #cystein
    'D': -3.5,  #aspartic acid
    'E': -3.5,  #glumatic acid
    'F': 2.8,   #phenylalanine
    'G': -0.4,  #glycine
    'H': -3.2,  #histidine
    'I': 4.5,   #isoleucine
    'K': -3.9,  #lysine
    'L': 3.8,   #leucine
    'M': 1.9,   #methionine
    'N': -3.5,  #asparagine
    'P': -1.6,  #proline
    'Q': -3.5,  #glutamine
    'R': -4.5,  #arginine
    'S': -0.8,  #serine
    'T': -0.7,  #threonine
    'V': 4.2,   #valine
    'W': -0.9,  #tryptophan
    'Y': -1.3,   #tyrosin
    'X': None #unknown AA
}

molecular_weight = {
    'A': 89.1,   #alanine
    'C': 121.2,   #cystein
    'D': 133.1,  #aspartic acid
    'E': 147.1,  #glumatic acid
    'F': 165.2,   #phenylalanine
    'G': 75.1,  #glycine
    'H': 155.2,  #histidine
    'I': 131.2,   #isoleucine
    'K': 146.2,  #lysine
    'L': 131.2,   #leucine
    'M': 149.2,   #methionine
    'N': 132.1,  #asparagine
    'P': 115.1,  #proline
    'Q': 146.2,  #glutamine
    'R': 174.2,  #arginine
    'S': 105.1,  #serine
    'T': 119.1,  #threonine
    'V': 117.1,   #valine
    'W': 204.2,  #tryptophan
    'Y': 181.2,   #tyrosin
    'X': None #unknown AA
}

random_prop = {
    'A': 54,   #alanine
    'C': 73,   #cystein
    'D': 84,  #aspartic acid
    'E': 72,  #glumatic acid
    'F': 71,   #phenylalanine
    'G': 75,  #glycine
    'H': 41,  #histidine
    'I': 56,   #isoleucine
    'K': 46,  #lysine
    'L': 67,   #leucine
    'M': 76,   #methionine
    'N': 92,  #asparagine
    'P': 86,  #proline
    'Q': 63,  #glutamine
    'R': 21,  #arginine
    'S': 19,  #serine
    'T': 49,  #threonine
    'V': 53,   #valine
    'W': 22,  #tryptophan
    'Y': 59,   #tyrosin
    'X': None #unknown AA
}

def get_AA_scores(scores_dict, fasta):
    result = []
    for i in range(1, len(fasta) + 1):
        AA_score = scores_dict[fasta[i-1]]
        if (AA_score is None):
            continue
        result.append((i, AA_score))
    return result
