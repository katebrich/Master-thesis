set -e
python_scripts_path=./pythonScripts/
P2Rank_path="/home/katebrich/Documents/diplomka/P2Rank"
tasks=t,e

train_ds=chen11
test_ds=joined
eval_ds=coach420
train_label=_10_27
test_label=_10_27
eval_label=_10_27
label=__custom_rfFeatures6_extraAsConfig

threads=4

features_list=unp_PTM
features_dir=${P2Rank_path}/datasets/${dataset_name}${label}/p2rank

#extra_features='(chem.protrusion.bfactor.volsite.csv_file_atom_feature)'
extra_features='(chem.volsite.protrusion.bfactor)'

dataset_prank_train=${P2Rank_path}/datasets/${train_ds}.ds
dataset_prank_test=${P2Rank_path}/datasets/${test_ds}.ds
dataset_prank_eval=${P2Rank_path}/datasets/${eval_ds}.ds

#if [[ $tasks == *"t"* ]]; then
#    cp -a $data_dir_train/p2rank/$label/. ${p2rank_location}/custom_feature/${label}
#    cp -a $data_dir_eval/p2rank/$label/. ${p2rank_location}/custom_feature/${label}
#fi

if [[ $tasks == *"t"* ]]; then
    bash ${P2Rank_path}/p2rank/prank traineval -t $dataset_prank_train -e $dataset_prank_test \
        -label $label \
        -rf_features 6 -rf_trees 200 \
        -threads $threads -delete_models 0 -loop 1 -seed 42 \
        -classifier FastRandomForest -feature_importances 1 \
        -extra_features $extra_features #\
    #-csv_file_feature_directories ",${P2Rank_path}/datasets/customFeature,"
    #-csv_file_feature_directories ",${P2Rank_path}/datasets/${train_ds}${train_label}/p2rank,${P2Rank_path}/datasets/${test_ds}${test_label}/p2rank,"
    #-c ${p2rank_location}/p2rank/config/thesis
fi

if [[ $tasks == *"e"* ]]; then
    bash ${P2Rank_path}/p2rank/prank eval-predict $dataset_prank_eval \
        -label $label \
        -threads $threads -delete_models 0 -loop 1 -seed 42 \
        -model ${P2Rank_path}/p2rank/test_output/traineval_${train_ds}_${test_ds}_${label}/runs/seed.42/FastRandomForest.model \
        -extra_features $extra_features
    #-csv_file_feature_directories ",${P2Rank_path}/datasets/customFeature,"
fi
