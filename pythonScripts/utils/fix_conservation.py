import os

datasets_path="/home/katebrich/Documents/diplomka/P2Rank/datasets/"
dataset_name="mix"
conservation_dir_old=f"{datasets_path}{dataset_name}/features/pdbekb_conservation_old"
conservation_dir_new=f"{datasets_path}{dataset_name}/features/pdbekb_conservation"

os.rename(conservation_dir_new, conservation_dir_old)
os.mkdir(conservation_dir_new)

for f in os.listdir(conservation_dir_old):
    old_path = os.path.join(conservation_dir_old, f)
    new_path = os.path.join(conservation_dir_new, f)
    with open(old_path, 'r') as old:
        with open(new_path, 'w') as new:
            for line in old:
                columns = line.split()
                if columns[1] == "4+":
                    columns[1] = "4"
                new.write(f"{columns[0]} {columns[1]}\n")
