dataset_file="/home/katebrich/Documents/diplomka/datasets/chen11.txt"
data_dir="/home/katebrich/Documents/diplomka/datasets/TEST_pipeline"
filter_ligands=none #todo none, water, small molecules, given IDs, MOAD,...
features=hydropathy,PTM
#force=false

#bash download_dataset.sh -d $dataset_file -o $data_dir -l $filter_ligands
#todo check jestli to pred tim probehlo bez chyby
bash compute_binding_sites.sh -d $dataset_file -i ${data_dir}/PDB -o ${data_dir}/lbs

#bash compute_features.sh -f $features -d $dataset_file -i $data_dir -o $data_dir/features
