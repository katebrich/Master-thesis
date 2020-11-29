set -e

#todo log experiment

python_scripts_path=./pythonScripts/
P2Rank_path="/home/katebrich/Documents/diplomka/P2Rank"
main_folder="/home/katebrich/Documents/diplomka/GitHub"

experiment_name="lbsOnly_avg"
threads=4
loops=1

filter=""

dataset_train_name="chen11"
dataset_train_name_filter=${dataset_train_name}${filter}
dataset_train_file="${main_folder}/datasets/${dataset_train_name_filter}.txt"
data_train_dir="${P2Rank_path}/datasets/${dataset_train_name_filter}"

dataset_test_name="joined"
dataset_test_name_filter=${dataset_test_name}${filter}
dataset_test_file="${main_folder}/datasets/${dataset_test_name_filter}.txt"
data_test_dir="${P2Rank_path}/datasets/${dataset_test_name_filter}"

dataset_eval_name="coach420"
dataset_eval_name_filter=${dataset_eval_name}${filter}
dataset_eval_file="${main_folder}/datasets/${dataset_eval_name_filter}.txt"
data_eval_dir="${P2Rank_path}/datasets/${dataset_eval_name_filter}"

#TODO smazat
dataset_eval2_name="holo4k"
dataset_eval2_name_filter=${dataset_eval2_name}${filter}
dataset_eval2_file="${main_folder}/datasets/${dataset_eval2_name_filter}.txt"
data_eval2_dir="${P2Rank_path}/datasets/${dataset_eval2_name_filter}"

##########################################
#            PREPARE DATA                #
##########################################

#features_list="helix,strand,turn,aa_ALA,aa_CYS,aa_ASP,aa_GLU,aa_PHE,aa_GLY,aa_HIS,aa_ILE,aa_LYS,aa_LEU,aa_MET,aa_ASN,aa_PRO,aa_GLN,aa_ARG,aa_SER,aa_THR,aa_VAL,aa_TRP,aa_TYR,hydropathy,charged,aromaticity,mol_weight,H_bond_atoms,HSE_up,HSE_down,exposure_CN,bfactor,pdbekb_conservation,dynamine,efoldmine,depth,mobiDB,phi_angle,psi_angle"
features_list="lbs"

echo "Running experiment $experiment_name..." #todo log

l=""
if [ ! -z $filter ]; then #filter is not empty string
	l="-l"
fi

#create dataset files for P2Rank
python3 ${python_scripts_path}create_p2rank_ds.py -p ${data_train_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_train_name_filter}.ds $l -d $dataset_train_file
echo "P2Rank dataset file created for $dataset_train_name_filter"
#python3 ${python_scripts_path}create_p2rank_ds.py -p ${data_test_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_test_name_filter}.ds $l -d $dataset_test_file
#echo "P2Rank dataset file created for $dataset_test_name_filter"
python3 ${python_scripts_path}create_p2rank_ds.py -p ${data_eval_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_eval_name_filter}.ds $l -d $dataset_eval_file
echo "P2Rank dataset file created for $dataset_eval_name_filter"
#python3 ${python_scripts_path}create_p2rank_ds.py -p ${data_eval2_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_eval2_name_filter}.ds $l -d $dataset_eval2_file
#echo "P2Rank dataset file created for $dataset_eval2_name_filter"

#todo udelat najednou cely list se vsemi featurami, potom uz jen colit sloupce v trenovani pranku pomoci feat_csv_columns

#create input files for P2Rank Custom feature
python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_train_file -i $data_train_dir -o $data_train_dir/P2Rank/$experiment_name -t $threads -f $features_list
echo "P2Rank custom feature files created for $dataset_train_name_filter"
#python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_test_file -i $data_test_dir -o $data_test_dir/P2Rank/$experiment_name -t $threads -f $features_list
#echo "P2Rank custom feature files created for $dataset_test_name_filter"
python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_eval_file -i $data_eval_dir -o $data_eval_dir/P2Rank/$experiment_name -t $threads -f $features_list
echo "P2Rank custom feature files created for $dataset_eval_name_filter"
#python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_eval2_file -i $data_eval2_dir -o $data_eval2_dir/P2Rank/$experiment_name -t $threads -f $features_list #todo smazat
#echo "P2Rank custom feature files created for $dataset_eval2_name_filter"

#############################################
#   P2Rank MODEL TRAINING AND EVALUATION    #
#############################################

extra_features='(csv)'

dataset_prank_train=${P2Rank_path}/datasets/${dataset_train_name_filter}.ds
dataset_prank_test=${P2Rank_path}/datasets/${dataset_test_name_filter}.ds
dataset_prank_eval=${P2Rank_path}/datasets/${dataset_eval_name_filter}.ds
dataset_prank_eval2=${P2Rank_path}/datasets/${dataset_eval2_name_filter}.ds

old_IFS=$IFS
IFS=','
#train and evaluate new model for each feature separately

IFS=$old_IFS

#copy input files with feature values to a common directory
[ -e ${P2Rank_path}/custom_feature/ ] && rm -r ${P2Rank_path}/custom_feature/
mkdir ${P2Rank_path}/custom_feature/
cp $data_train_dir/P2Rank/$experiment_name/*.csv ${P2Rank_path}/custom_feature/
#	cp $data_test_dir/P2Rank/$experiment_name/*.csv ${P2Rank_path}/custom_feature/
cp $data_eval_dir/P2Rank/$experiment_name/*.csv ${P2Rank_path}/custom_feature/
#	cp $data_eval2_dir/P2Rank/$experiment_name/*.csv ${P2Rank_path}/custom_feature/ || true
echo "Files with feature values copied to ${P2Rank_path}/custom_feature"

#train
bash ${P2Rank_path}/p2rank/prank traineval -t $dataset_prank_train -e $dataset_prank_eval \
	-label ${experiment_name} \
	-rf_features 6 -rf_trees 100 -rf_depth 0 \
	-threads $threads -delete_models 0 -loop ${loops} -seed 42 \
	-classifier FastRandomForest -feature_importances 1 \
	-features $extra_features -atom_table_features '()' \
	-feat_csv_directories ",${P2Rank_path}/custom_feature," \
	-feat_csv_columns "($features_list)" \
	-feat_csv_ignore_missing 1 \
	-average_feat_vectors true -avg_weighted true
echo "P2Rank model training finished."

##evaluate
#bash ${P2Rank_path}/p2rank/prank eval-predict $dataset_prank_eval \
#		-label ${experiment_name}_${feature} \
#		-threads $threads -delete_models 0 -loop ${loops} -seed 42 \
#		-model ${P2Rank_path}/p2rank/test_output/traineval_${dataset_train_name_filter}_${dataset_test_name_filter}_${experiment_name}_${feature}/runs/seed.42/FastRandomForest.model \
#		-features $extra_features \
#		-feat_csv_directories ",${P2Rank_path}/custom_feature," \
#		-feat_csv_columns "($feature)" \
#		-feat_csv_ignore_missing 1 || true
#	echo "P2Rank model evaluation finished."
#
#	#evaluate 2
#	bash ${P2Rank_path}/p2rank/prank eval-predict $dataset_prank_eval2 \
#		-label ${experiment_name}_${feature} \
#		-threads $threads -delete_models 0 -loop ${loops} -seed 42 \
#		-model ${P2Rank_path}/p2rank/test_output/traineval_${dataset_train_name_filter}_${dataset_test_name_filter}_${experiment_name}_${feature}/runs/seed.42/FastRandomForest.model \
#		-features $extra_features \
#		-feat_csv_directories ",${P2Rank_path}/custom_feature," \
#		-feat_csv_columns "($feature)" \
#		-feat_csv_ignore_missing 1 || true
#
#	echo "P2Rank model evaluation 2 finished."
#
