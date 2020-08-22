set -e
python_scripts_path=./pythonScripts/

#todo parsovani argumentu
dataset_file="/home/katebrich/Documents/diplomka/datasets/chen11.txt" #povinny
out_dir="/home/katebrich/Documents/diplomka/datasets/TEST_pipeline"   #nepovinny. Kdyz neni zadan, vytvori se v umisteni dataset_file
features_dir="${out_dir}/features"                                    #nepovinny. Default podslozka out_dir
filter_ligands=none                                                   #todo none, water, small molecules, given IDs, MOAD,...
features_list=pdbkb_conservation                                      #povinny
threads=4
#log_file="${out_dir}/run.log"
#strict=true #todo
analysis_dir="${out_dir}/analysis"

#if [ -d "$out_dir" ]; then
#    rm -rf "$out_dir"
#fi

#python3 ${python_scripts_path}download_dataset.py -d $dataset_file -o $out_dir -l $filter_ligands -t $threads

python3 ${python_scripts_path}create_mappings_cache.py -d $dataset_file -o $out_dir/mappings -t $threads

#python3 ${python_scripts_path}compute_ligand_binding_sites.py -d $dataset_file -i ${out_dir}/PDB -o ${out_dir}/lbs -t $threads

oldIFS=$IFS
IFS=','

#for feature in $features_list; do
#    python3 ${python_scripts_path}compute_feature.py -f $feature -d $dataset_file -i $out_dir -o $features_dir/$feature -t $threads
#done

#mkdir -p $analysis_dir #todo dovnitr skriptu

#for feature in $features_list; do
#todo check jestli ta slozka s hodnotami featury existuje
#    python3 ${python_scripts_path}run_analysis.py -d $dataset_file -l $out_dir/lbs -v $features_dir/$feature -o $analysis_dir/$feature -f $feature -t $threads
#done

IFS=$oldIFS #todo je to potreba?

#cp ./run.log ${out_dir}/run.log
