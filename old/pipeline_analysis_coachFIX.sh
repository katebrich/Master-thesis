#todo log experiment

python_scripts_path=./pythonScripts/
P2Rank_path="/home/katebrich/Documents/diplomka/P2Rank"
main_folder="/home/katebrich/Documents/diplomka/GitHub"

threads=4

#features_list="unp_disulfid,aa_ratio,unp_PTM,unp_glycosylation,unp_lipidation,unp_mod_res,unp_variation,unp_topology,unp_sec_str,unp_non_standard,unp_natural_variant,unp_compbias,aa,aa_pairs,hydropathy,polarity,charged,aromaticity,mol_weight,H_bond_atoms,HSE_up,HSE_down,exposureCN,bfactor,pdbekb_sec_str,pdbekb_conservation,dynamine,efoldmine,depth,mobiDB,phi_angle,psi_angle,cis_peptide,lbs"
features_list="x"
tasks="A"

# 1. dataset without filter

dataset_name_short=joined_filter_p2rank

dataset_train_name="${dataset_name_short}"
experiment_name="${dataset_train_name}_analysis_FIX"
dataset_train_file="${main_folder}/datasets/${dataset_train_name}.txt"
data_train_dir="${P2Rank_path}/datasets/${dataset_train_name}"

echo "Running experiment $experiment_name..." #todo more details

#download, compute binding sites, features and statistical analysis
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t F -m $threads -f "dynamine,efoldmine,pdbekb_conservation"
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 0 -i 1 -b 0
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 500 -i 1000 -b 1
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 100 -i 1000 -b 1
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 1000 -i 1000 -b 1
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 200 -i 1000 -b 1
echo "Data prepared successfully for $dataset_train_file"
exit
# 2. p2rank filter

dataset_train_name="${dataset_name_short}_filter_p2rank"
dataset_train_file="${main_folder}/datasets/${dataset_train_name}.txt"
data_train_dir="${P2Rank_path}/datasets/${dataset_train_name}"

#download, compute binding sites, features and statistical analysis
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t F -m $threads -f pdbekb_conservation
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 0 -i 1 -b 0
tasks="A"
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 500 -i 1000 -b 1
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 100 -i 1000 -b 1
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 1000 -i 1000 -b 1
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 200 -i 1000 -b 1
echo "Data prepared successfully for $dataset_train_file"

# 3. moad filter

dataset_train_name="${dataset_name_short}_filter_MOAD"
dataset_train_file="${main_folder}/datasets/${dataset_train_name}.txt"
data_train_dir="${P2Rank_path}/datasets/${dataset_train_name}"

#download, compute binding sites, features and statistical analysis
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t F -m $threads -f pdbekb_conservation
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 0 -i 1 -b 0
tasks="A"
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 500 -i 1000 -b 1
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 100 -i 1000 -b 1
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 1000 -i 1000 -b 1
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list -s 200 -i 1000 -b 1
echo "Data prepared successfully for $dataset_train_file"
