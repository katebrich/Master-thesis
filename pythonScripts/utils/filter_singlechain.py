
dataset_name="mix"
new_name="mix_noDupl"
datasets_path="/home/katebrich/Documents/diplomka/GitHub/datasets"
dataset_file=f"{datasets_path}/{dataset_name}.txt"
output_file= f"{datasets_path}/{new_name}.txt"

structures = {}
with open(dataset_file, 'r') as input:
    with open(output_file, 'w') as out:
        for line in input:
            s = line.split()
            if (s[0] in structures):
                if (structures[s[0]] == s[1]):
                    continue
            else:
                structures[s[0]] = s[1]

with open(output_file, 'w') as out:
    for key, value in structures.items():
        out.write(f"{key}\t{value}\n")