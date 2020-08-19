set -e
python_scripts_path=./pythonScripts/

#todo parsovani argumentu
dataset_file="/home/katebrich/Documents/diplomka/datasets/test.txt" #povinny
out_dir="/home/katebrich/Documents/diplomka/datasets/pipeline_test" #nepovinny. Kdyz neni zadan, vytvori se v umisteni dataset_file
features_dir="${out_dir}/features"                                  #nepovinny. Default podslozka out_dir
filter_ligands=none                                                 #todo none, water, small molecules, given IDs, MOAD,...
features_list=pdbkb_conservation                                    #povinny
threads=4
#log_file="${out_dir}/run.log"
#strict=true #todo

if [ -d "$out_dir" ]; then
    rm -rf "$out_dir"
fi

#stahnout PDB a FASTA soubory pro dany dataset z databaze
#bash download_dataset.sh -d $dataset_file -o $out_dir -l $filter_ligands -t $threads -s
python3 ${python_scripts_path}download_dataset.py -d $dataset_file -o $out_dir -l $filter_ligands -t $threads

#vytvorit cache mapovani author residue number : pdbe sequence number
python3 ${python_scripts_path}create_mappings_cache.py -d $dataset_file -o $out_dir/mappings -t $threads

#vypocitat binding sites
#bash compute_binding_sites.sh -d $dataset_file -i ${out_dir}/PDB -o ${out_dir}/lbs -t $threads
python3 ${python_scripts_path}compute_ligand_binding_sites.py -d $dataset_file -i ${out_dir}/PDB -o ${out_dir}/lbs -t $threads

#vypocitat hodnoty zadanych featur
#bash compute_features.sh -f $features_list -d $dataset_file -i $out_dir -o $out_dir/features -t $threads -s #todo udelat volitelne, ze se vypocitaji?

oldIFS=$IFS
IFS=','

for feature in $features_list; do
    python3 ${python_scripts_path}compute_feature.py -f $feature -d $dataset_file -i $out_dir -o $features_dir/$feature -t $threads
done

analysis_dir="${out_dir}/analysis"
#mkdir -p $analysis_dir #todo dovnitr skriptu

#spustit hypothesis test pro kazdou featuru, ulozit vysledky
for feature in $features_list; do
    #todo check jestli ta slozka s hodnotami featury existuje
    bash run_analysis.sh -d $dataset_file -l $out_dir/lbs -v $features_dir/$feature -o $analysis_dir/$feature -f $feature -t $threads -s
done

#restore IFS
IFS=$oldIFS #todo je to potreba?

#copy log to output directory
cp ./run.log ${out_dir}/run.log
