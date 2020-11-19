set -e

#todo log experiment

python_scripts_path=./pythonScripts/
P2Rank_path="/home/katebrich/Documents/diplomka/P2Rank"
main_folder="/home/katebrich/Documents/diplomka/GitHub"

experiment_name="TEST"
threads=4
loops=1

filter="_filter_p2rank"

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

##########################################
#            PREPARE DATA                #
##########################################

#features_list="unp_disulfid,unp_PTM,unp_glycosylation,unp_lipidation,unp_mod_res,unp_variation,unp_topology,unp_sec_str,unp_non_standard,unp_natural_variant,unp_compbias,aa,aa_pairs,hydropathy,polarity,polarity_binary,charged,aromaticity,mol_weight,H_bond_atoms,HSE_up,HSE_down,exposureCN,bfactor,bfactor_CA,pdbekb_sec_str,pdbekb_conservation,dynamine,efoldmine,depth,mobiDB,conservation,phi_angle,psi_angle,cis_peptide"

features_list="pdbekb_conservation,unp_PTM"

echo "Running experiment $experiment_name..." #todo log

#create dataset files for P2Rank
#python3 ${python_scripts_path}create_prank_ds.py -d ${data_train_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_train_name_filter}.ds
#echo "P2Rank dataset file created for $dataset_train_name_filter"
#python3 ${python_scripts_path}create_prank_ds.py -d ${data_test_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_test_name_filter}.ds
#echo "P2Rank dataset file created for $dataset_test_name_filter"
#python3 ${python_scripts_path}create_prank_ds.py -d ${data_eval_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_eval_name_filter}.ds
#echo "P2Rank dataset file created for $dataset_eval_name_filter"

old_IFS=$IFS
IFS=','

#todo udelat najednou cely list se vsemi featurami, potom uz jen colit sloupce v trenovani pranku pomoci feat_csv_columns

##create input files for P2Rank Custom feature
#for feature in $features_list; do
#	python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_train_file -i $data_train_dir -o $data_train_dir/P2Rank/$feature -t $threads -f $feature
#	echo "Feature $feature: P2Rank custom feature files created for $dataset_train_name_filter"
#	python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_test_file -i $data_test_dir -o $data_test_dir/P2Rank/$feature -t $threads -f $feature
#	echo "Feature $feature: P2Rank custom feature files created for $dataset_test_name_filter"
#	python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_eval_file -i $data_eval_dir -o $data_eval_dir/P2Rank/$feature -t $threads -f $feature
#	echo "Feature $feature: P2Rank custom feature files created for $dataset_eval_name_filter"
#done

#############################################
#   P2Rank MODEL TRAINING AND EVALUATION    #
#############################################

extra_features='(chem,volsite,protrusion,bfactor,csv)'

dataset_prank_train=${P2Rank_path}/datasets/${dataset_train_name_filter}.ds
dataset_prank_test=${P2Rank_path}/datasets/${dataset_test_name_filter}.ds
dataset_prank_eval=${P2Rank_path}/datasets/${dataset_eval_name_filter}.ds

#train and evaluate new model for each feature separately
for feature in $features_list; do

	IFS=$old_IFS

	#copy input files with feature values to a common directory
	[ -e ${P2Rank_path}/custom_feature/ ] && rm -r ${P2Rank_path}/custom_feature/
	mkdir ${P2Rank_path}/custom_feature/
	cp $data_train_dir/P2Rank/$feature/*.csv ${P2Rank_path}/custom_feature/
	cp $data_test_dir/P2Rank/$feature/*.csv ${P2Rank_path}/custom_feature/

	echo "Datasets $dataset_train_name_filter and $dataset_test_name_filter: Feature $feature: files with feature values copied to ${P2Rank_path}/custom_feature"

	#todo vyzkouset feat_csv_ignore_missing

	#train
	bash ${P2Rank_path}/p2rank/prank traineval -t $dataset_prank_train -e $dataset_prank_test \
		-label ${experiment_name}_${feature} \
		-rf_features 6 -rf_trees 200 -rf_depth 0 \
		-threads $threads -delete_models 0 -loop ${loops} -seed 42 \
		-classifier FastRandomForest -feature_importances 1 \
		-features $extra_features \
		-csv_file_feature_directories ",${P2Rank_path}/custom_feature," \
		-feat_csv_columns "($feature)" \
		-feat_csv_ignore_missing 1
	#-average_feat_vectors true -avg_weighted true

	echo "P2Rank model training finished."

	#copy input files with feature values to a common directory
	[ -e ${P2Rank_path}/custom_feature/ ] && rm -r ${P2Rank_path}/custom_feature/
	mkdir ${P2Rank_path}/custom_feature/
	cp $data_eval_dir/P2Rank/$feature/*.csv ${P2Rank_path}/custom_feature/

	echo "Dataset $dataset_eval_name_filter: Feature $feature: files with feature values copied to ${P2Rank_path}/custom_feature"

	#evaluate
	bash ${P2Rank_path}/p2rank/prank eval-predict $dataset_prank_eval \
		-label ${experiment_name} \
		-threads $threads -delete_models 0 -loop 1 -seed 42 \
		-model ${P2Rank_path}/p2rank/test_output/traineval_${dataset_train_name_filter}_${dataset_test_name_filter}_${experiment_name}_${feature}/runs/seed.42/FastRandomForest.model \
		-features $extra_features \
		-csv_file_feature_directories ",${P2Rank_path}/custom_feature," \
		-feat_csv_columns $feature \
		-feat_csv_ignore_missing 1

	echo "P2Rank model evaluation finished."

done
