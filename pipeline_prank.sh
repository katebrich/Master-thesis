set -e
python_scripts_path=./pythonScripts/

p2rank_location="/home/katebrich/Documents/diplomka/P2Rank_with_csv_feature"

dataset_train="/home/katebrich/Documents/diplomka/datasets/chen11.txt"
dataset_eval="/home/katebrich/Documents/diplomka/datasets/joined(mlig).txt"
data_dir_train="/home/katebrich/Documents/diplomka/datasets/TEST_1/chen11"
data_dir_eval="/home/katebrich/Documents/diplomka/datasets/TEST_1/joined(mlig)"
features_list=unp_PTM
filter_ligands=none #todo
threads=4
label=TEST_unp_PTM
#strict=true #todo

#todo default vals

filename=$(basename -- "$dataset_train")
filename="${filename%.*}"
dataset_prank_train="${data_dir_train}/${filename##*/}.ds"

filename=$(basename -- "$dataset_eval")
filename="${filename%.*}"
dataset_prank_eval="${data_dir_eval}/${filename##*/}.ds"

#bash download_dataset.sh -d $dataset_train -o $data_dir_train -l $filter_ligands -t $threads -s
#bash download_dataset.sh -d $dataset_eval -o $data_dir_eval -l $filter_ligands -t $threads -s

#python3 ${python_scripts_path}create_mappings_cache.py -d $dataset_train -o $data_dir_train/mappings -t $threads
#python3 ${python_scripts_path}create_mappings_cache.py -d $dataset_eval -o $data_dir_eval/mappings -t $threads

#bash compute_features.sh -f $features_list -d $dataset_train -i $data_dir_train -o $data_dir_train/features -t $threads -s
bash compute_features.sh -f $features_list -d $dataset_eval -i $data_dir_eval -o $data_dir_eval/features -t $threads -s

bash create_p2rank_custom_feature_files.sh -f $features_list -d $dataset_train -i $data_dir_train -o $data_dir_train/p2rank/$label -t $threads -s
bash create_p2rank_custom_feature_files.sh -f $features_list -d $dataset_eval -i $data_dir_eval -o $data_dir_eval/p2rank/$label -t $threads -s

cp -a $data_dir_train/p2rank/$label/. ${p2rank_location}/custom_feature/${label}
cp -a $data_dir_eval/p2rank/$label/. ${p2rank_location}/custom_feature/${label}

python3 ${python_scripts_path}create_prank_ds.py -d $data_dir_train/PDB -o $dataset_prank_train
python3 ${python_scripts_path}create_prank_ds.py -d $data_dir_eval/PDB -o ${dataset_prank_eval}

#todo vytvorit config, slozku custom_feature atd. pokud chybi
bash ${p2rank_location}/p2rank/prank traineval -t $dataset_prank_train -e $dataset_prank_eval \
    -csv_file_feature_directories ",${p2rank_location}/custom_feature/${label}," -label $label -c ${p2rank_location}/p2rank/config/custom_feature \
    -threads $threads -rf_trees 128 -rf_depth 6 -delete_models 0 -loop 1 -seed 42
