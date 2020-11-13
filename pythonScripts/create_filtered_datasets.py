from DatasetLigandsFilter import DatasetLigandsFilter

filter="p2rank"
dataset_name="coach420"
label="11_09_test"
datasets_path="/home/katebrich/Documents/diplomka/GitHub/datasets"
dataset_file=f"{datasets_path}/{dataset_name}.txt"
output_file= f"{datasets_path}/{dataset_name}_filter_{filter}.txt"
pdb_dir=f"/home/katebrich/Documents/diplomka/P2Rank/datasets/{label}/{dataset_name}/PDB"
threads=4 #todo


df = DatasetLigandsFilter(filter)
df.run(dataset_file, output_file, pdb_dir, threads, remove_empty_lines=True)