set -e

#todo log experiment

python_scripts_path=./pythonScripts/
P2Rank_path="/home/katebrich/Documents/diplomka/P2Rank"
main_folder="/home/katebrich/Documents/diplomka/GitHub"

experiment_name="11_09_test"

dataset_train_name="chen11"
dataset_train_file="${main_folder}/datasets/${dataset_train_name}.txt"
data_train_dir="${P2Rank_path}/datasets/${experiment_name}/${dataset_train_name}"

dataset_test_name="joined"
dataset_test_file="${main_folder}/datasets/${dataset_test_name}.txt"
data_test_dir="${P2Rank_path}/datasets/${experiment_name}/${dataset_test_name}"

dataset_eval_name="coach420"
dataset_eval_file="${main_folder}/datasets/${dataset_eval_name}.txt"
data_eval_dir="${P2Rank_path}/datasets/${experiment_name}/${dataset_eval_name}"

threads=4

##########################################
#            PREPARE DATA                #
##########################################

tasks="A"
features_list="pdbekb_conservation"

#copy conservation data
mkdir -p $data_train_dir/conservation
cp -r ${main_folder}/conservation_results/$dataset_train_name/*.json $data_train_dir/conservation
mkdir -p $data_test_dir/conservation
cp -r ${main_folder}/conservation_results/$dataset_test_name/*.json $data_test_dir/conservation
mkdir -p $data_eval_dir/conservation
cp -r ${main_folder}/conservation_results/$dataset_eval_name/*.json $data_eval_dir/conservation

#download, compute binding sites, features and statistical analysis
python3 ${python_scripts_path}pipeline.py -d $dataset_train_file -o $data_train_dir -t $tasks -m $threads -f $features_list
python3 ${python_scripts_path}pipeline.py -d $dataset_test_file -o $data_test_dir -t $tasks -m $threads -f $features_list
python3 ${python_scripts_path}pipeline.py -d $dataset_eval_file -o $data_eval_dir -t $tasks -m $threads -f $features_list

#create dataset file for P2Rank
python3 ${python_scripts_path}create_prank_ds.py -d ${data_train_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_train_name}.ds
python3 ${python_scripts_path}create_prank_ds.py -d ${data_test_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_test_name}.ds
python3 ${python_scripts_path}create_prank_ds.py -d ${data_eval_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_eval_name}.ds

IFS=','

#create input files for P2Rank Custom feature
for feature in $features_list; do
	python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_train_file -i $data_train_dir -o $data_train_dir/P2Rank/$feature -t $threads -f $feature
	python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_test_file -i $data_test_dir -o $data_test_dir/P2Rank/$feature -t $threads -f $feature
	python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_eval_file -i $data_eval_dir -o $data_eval_dir/P2Rank/$feature -t $threads -f $feature
done

#############################################
#   P2Rank MODEL TRAINING AND EVALUATION    #
#############################################

extra_features='(chem.volsite.protrusion.bfactor.csv_file_atom_feature)'
dataset_prank_train=${P2Rank_path}/datasets/${dataset_train_name}.ds
dataset_prank_test=${P2Rank_path}/datasets/${dataset_test_name}.ds
dataset_prank_eval=${P2Rank_path}/datasets/${dataset_eval_name}.ds

#train and evaluate new model for each feature separately
for feature in $features_list; do

	#copy input files with feature values to a common directory
	[ -e ${P2Rank_path}/custom_feature/ ] && rm -r ${P2Rank_path}/custom_feature/
	mkdir ${P2Rank_path}/custom_feature/
	cp $data_train_dir/P2Rank/$feature/*.csv ${P2Rank_path}/custom_feature/
	cp $data_test_dir/P2Rank/$feature/*.csv ${P2Rank_path}/custom_feature/

	#train
	bash ${P2Rank_path}/p2rank/prank traineval -t $dataset_prank_train -e $dataset_prank_test \
		-label $experiment_name \
		-rf_features 6 -rf_trees 200 \
		-threads $threads -delete_models 0 -loop 1 -seed 42 \
		-classifier FastRandomForest -feature_importances 1 \
		-extra_features $extra_features \
		-csv_file_feature_directories ",${P2Rank_path}/custom_feature,"

	#copy input files with feature values to a common directory
	[ -e ${P2Rank_path}/custom_feature/ ] && rm -r ${P2Rank_path}/custom_feature/
	mkdir ${P2Rank_path}/custom_feature/
	cp $data_eval_dir/P2Rank/$feature/*.csv ${P2Rank_path}/custom_feature/

	#evaluate
	bash ${P2Rank_path}/p2rank/prank eval-predict $dataset_prank_eval \
		-label $experiment_name \
		-threads $threads -delete_models 0 -loop 1 -seed 42 \
		-model ${P2Rank_path}/p2rank/test_output/traineval_${dataset_train_name}_${dataset_test_name}_${experiment_name}/runs/seed.42/FastRandomForest.model \
		-extra_features $extra_features \
		-csv_file_feature_directories ",${P2Rank_path}/custom_feature,"

done
