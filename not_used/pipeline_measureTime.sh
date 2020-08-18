set -e
python_scripts_path=./pythonScripts/

dataset_file="/home/katebrich/Documents/diplomka/datasets/chen11.txt"
data_dir="/home/katebrich/Documents/diplomka/datasets/TEST_pipeline"
filter_ligands=none #todo none, water, small molecules, given IDs, MOAD,...
features_list=pKa_COOH,pKa_NH3,molecular_weight,unp_variants,hydropathy
threads=4
#strict=true #todo

START_TIME=$SECONDS
#bash download_dataset.sh -d $dataset_file -o $data_dir -l $filter_ligands -t $threads
ELAPSED_TIME=$((SECONDS - START_TIME))
echo "download_dataset.sh - elapsed time:" $ELAPSED_TIME

START_TIME=$SECONDS
#bash compute_binding_sites.sh -d $dataset_file -i ${data_dir}/PDB -o ${data_dir}/lbs -t $threads
ELAPSED_TIME=$((SECONDS - START_TIME))
echo "compute_binding_sites.sh - elapsed time:" $ELAPSED_TIME

START_TIME=$SECONDS
#python3 ${python_scripts_path}create_mappings_cache.py -d $dataset_file -o $data_dir/mappings -t $threads
ELAPSED_TIME=$((SECONDS - START_TIME))
echo "create_mappings_cache.py - elapsed time:" $ELAPSED_TIME

START_TIME=$SECONDS
bash compute_features.sh -f $features_list -d $dataset_file -i $data_dir -o $data_dir/features -t $threads -s
ELAPSED_TIME=$((SECONDS - START_TIME))
echo "compute_features.sh - elapsed time:" $ELAPSED_TIME

mkdir -p $data_dir/analysis

START_TIME=$SECONDS
oldIFS=$IFS
IFS=','
for feature in $features_list; do
    bash run_analysis.sh -d $dataset_file -l $data_dir/lbs -v $data_dir/features/$feature -o $data_dir/analysis/$feature -f $feature -t $threads -s
done
ELAPSED_TIME=$((SECONDS - START_TIME))
echo "run_analysis.sh - elapsed time:"$ELAPSED_TIME
