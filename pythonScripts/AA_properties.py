
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
    'Y': -1.3   #tyrosin
}

def get_AA_scores(scores_dict, fasta):
    result = []
    for i in range(1, len(fasta) + 1):
        AA_score = scores_dict[fasta[i-1]]
        result.append((i, AA_score))
    return result
