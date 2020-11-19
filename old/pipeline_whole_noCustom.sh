set -e

#todo log experiment

python_scripts_path=./pythonScripts/
P2Rank_path="/home/katebrich/Documents/diplomka/P2Rank"
main_folder="/home/katebrich/Documents/diplomka/GitHub"

experiment_name="noFilter_noCustom"
threads=4

dataset_train_name="chen11"
dataset_train_file="${main_folder}/datasets/${dataset_train_name}.txt"
data_train_dir="${P2Rank_path}/datasets/final/${dataset_train_name}"

dataset_test_name="joined"
dataset_test_file="${main_folder}/datasets/${dataset_test_name}.txt"
data_test_dir="${P2Rank_path}/datasets/final/${dataset_test_name}"

dataset_eval_name="coach420"
dataset_eval_file="${main_folder}/datasets/${dataset_eval_name}.txt"
data_eval_dir="${P2Rank_path}/datasets/final/${dataset_eval_name}"

#############################################
#   P2Rank MODEL TRAINING AND EVALUATION    #
#############################################

extra_features='(chem.volsite.protrusion.bfactor)'

dataset_prank_train=${P2Rank_path}/datasets/${dataset_train_name}.ds
dataset_prank_test=${P2Rank_path}/datasets/${dataset_test_name}.ds
dataset_prank_eval=${P2Rank_path}/datasets/${dataset_eval_name}.ds

#create dataset file for P2Rank
python3 ${python_scripts_path}create_prank_ds.py -d ${data_train_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_train_name}.ds
echo "P2Rank dataset file created for $dataset_train_name"
python3 ${python_scripts_path}create_prank_ds.py -d ${data_test_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_test_name}.ds
echo "P2Rank dataset file created for $dataset_test_name"
python3 ${python_scripts_path}create_prank_ds.py -d ${data_eval_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_eval_name}.ds
echo "P2Rank dataset file created for $dataset_eval_name"

#train
bash ${P2Rank_path}/p2rank/prank traineval -t $dataset_prank_train -e $dataset_prank_test \
	-label ${experiment_name} \
	-rf_features 6 -rf_trees 200 \
	-threads $threads -delete_models 0 -loop 1 -seed 42 \
	-classifier FastRandomForest -feature_importances 1 \
	-extra_features $extra_features

echo "P2Rank model training finished."

#evaluate
bash ${P2Rank_path}/p2rank/prank eval-predict $dataset_prank_eval \
	-label ${experiment_name} \
	-threads $threads -delete_models 0 -loop 1 -seed 42 \
	-model ${P2Rank_path}/p2rank/test_output/traineval_${dataset_train_name}_${dataset_test_name}_${experiment_name}/runs/seed.42/FastRandomForest.model \
	-extra_features $extra_features

echo "P2Rank model evaluation finished."
