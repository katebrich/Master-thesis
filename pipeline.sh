dataset_file="/home/katebrich/Documents/diplomka/datasets/chen11.txt"
output_dir="/home/katebrich/Documents/diplomka/TEST_pipeline"
filter_ligands=none #todo none, water, small molecules, given IDs, MOAD,...
#force=false

bash download_dataset.sh -d $dataset_file -o $output_dir -l $filter_ligands
