set -e
python_scripts_path=./pythonScripts/

#todo parsovani argumentu
dataset_file="/home/katebrich/Documents/diplomka/datasets/holo4k(mlig).txt"
data_dir="/home/katebrich/Documents/diplomka/datasets/holo4k(mlig)"
filter_ligands=none #todo none, water, small molecules, given IDs, MOAD,...
features_list=pKa_COOH,pKa_NH3,molecular_weight,unp_variants,hydropathy,unp_PTM,unp_glycosylation
threads=4
#strict=true #todo

#stahnout PDB a FASTA soubory pro dany dataset z databaze
bash download_dataset.sh -d $dataset_file -o $data_dir -l $filter_ligands -t $threads

#vytvorit cache mapovani author residue number : pdbe sequence number
python3 ${python_scripts_path}create_mappings_cache.py -d $dataset_file -o $data_dir/mappings -t $threads

#vypocitat binding sites
bash compute_binding_sites.sh -d $dataset_file -i ${data_dir}/PDB -o ${data_dir}/lbs -t $threads

#vypocitat hodnoty zadanych featur
bash compute_features.sh -f $features_list -d $dataset_file -i $data_dir -o $data_dir/features -t $threads -s

mkdir -p $data_dir/analysis

#spustit hypothesis test pro kazdou featuru, ulozit vysledky
oldIFS=$IFS
IFS=','
for feature in $features_list; do
    bash run_analysis.sh -d $dataset_file -l $data_dir/lbs -v $data_dir/features/$feature -o $data_dir/analysis/$feature -f $feature -t $threads -s
done
