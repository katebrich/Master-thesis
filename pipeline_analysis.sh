#todo log experiment

ds=joined

python_scripts_path=./pythonScripts/
P2Rank_path="/home/katebrich/Documents/diplomka/P2Rank"
main_folder="/home/katebrich/Documents/diplomka/GitHub"

experiment_name="final"
threads=4

tasks="F,A"
#features_list="unp_disulfid,unp_PTM,unp_glycosylation,unp_lipidation,unp_mod_res,unp_variation,unp_topology,unp_sec_str,unp_non_standard,unp_natural_variant,unp_compbias,aa,aa_pairs,hydropathy,polarity,polarity_binary,charged,aromaticity,mol_weight,H_bond_atoms,HSE_up,HSE_down,exposureCN,bfactor,bfactor_CA,pdbekb_sec_str,pdbekb_conservation,dynamine,efoldmine,depth,mobiDB,phi_angle,psi_angle,cis_peptide,lbs"
features_list="aa_ratio,unp_disulfid"

# 1. dataset without filter

dataset_train_name=${ds}
dataset_train_file="${main_folder}/datasets/${dataset_train_name}.txt"
data_train_dir="${P2Rank_path}/datasets/${experiment_name}/${dataset_train_name}"

echo "Running experiment $experiment_name..." #todo more details

##copy conservation data
#mkdir -p $data_train_dir/conservation
#cp -r ${main_folder}/conservation_results/$dataset_train_name/*.json $data_train_dir/conservation
#echo "Conservation data copied for $dataset_train_name"

#download, compute binding sites, features and statistical analysis
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list
echo "Data prepared successfully for $dataset_train_file"

exit
# 2. p2rank filter

dataset_train_name="${ds}_filter_p2rank"
dataset_train_file="${main_folder}/datasets/${dataset_train_name}.txt"
data_train_dir="${P2Rank_path}/datasets/${experiment_name}/${dataset_train_name}"

##copy conservation data
#mkdir -p $data_train_dir/conservation
#cp -r ${main_folder}/conservation_results/$dataset_train_name/*.json $data_train_dir/conservation
#echo "Conservation data copied for $dataset_train_name"

#download, compute binding sites, features and statistical analysis
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list
echo "Data prepared successfully for $dataset_train_file"

# 3. moad filter

dataset_train_name="${ds}_filter_MOAD"
dataset_train_file="${main_folder}/datasets/${dataset_train_name}.txt"
data_train_dir="${P2Rank_path}/datasets/${experiment_name}/${dataset_train_name}"

##copy conservation data
#mkdir -p $data_train_dir/conservation
#cp -r ${main_folder}/conservation_results/$dataset_train_name/*.json $data_train_dir/conservation
#echo "Conservation data copied for $dataset_train_name"

#download, compute binding sites, features and statistical analysis
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list
echo "Data prepared successfully for $dataset_train_file"
