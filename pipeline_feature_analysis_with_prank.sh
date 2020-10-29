set -e
python_scripts_path=./pythonScripts/
P2Rank_path=/home/katebrich/Documents/diplomka/P2Rank_with_csv_feature

tasks="f,p" #"r,d,m,l,f,a,p"
#todo parsovani argumentu
dataset_name="joined"
label="_10_29"

dataset_file="/home/katebrich/Documents/diplomka/datasets/${dataset_name}.txt" #povinny
data_dir_name=${dataset_name}${label}
data_dir="${P2Rank_path}/datasets/${data_dir_name}" #nepovinny. Kdyz neni zadan, vytvori se v umisteni dataset_file
features_dir="${data_dir}/features"                 #nepovinny. Default podslozka data_dir
filter_ligands=p2rank                               #todo none, water, small molecules, given IDs, MOAD,...
features_list="unp_PTM,unp_glycosylation,unp_lipidation,unp_mod_res,unp_variation,unp_topology,unp_sec_str,unp_non_standard,unp_natural_variant,unp_compbias,pdbekb_conservation,pdbekb_sec_str,aa,aa_pairs,hydropathy,polarity,polarity_binary,charged,aromaticity,mol_weight,H_bond_atoms,dynamine_website,dynamine_funPDBe,efoldmine_funPDBe,mobiDB,HSE,HSE_down,exposureCN,bfactor,bfactor_Calpha,depth,phi_angle,psi_angle,cis_peptide"
threads=4
#log_file="${data_dir}/run.log"
#strict=true #todo
analysis_dir="${data_dir}/analysis"

rm -f ./run.log

if [[ $tasks == *"r"* ]]; then
    #!!!!!POZOR!!!!!!!!#
    if [ -d "$data_dir" ]; then
        rm -rf "$data_dir"
    fi
fi

if [[ $tasks == *"d"* ]]; then
    python3 ${python_scripts_path}download_dataset.py -d $dataset_file -o $data_dir -l $filter_ligands -t $threads
fi

if [[ $tasks == *"m"* ]]; then
    python3 ${python_scripts_path}create_mappings_cache.py -d $dataset_file -o $data_dir/mappings -t $threads
fi

if [[ $tasks == *"l"* ]]; then
    python3 ${python_scripts_path}compute_ligand_binding_sites.py -d $dataset_file -i ${data_dir} -t $threads
fi

#oldIFS=$IFS #todo smazat?
IFS=','

if [[ $tasks == *"f"* ]]; then
    for feature in $features_list; do
        python3 ${python_scripts_path}compute_feature.py -f $feature -d $dataset_file -i $data_dir -o $features_dir/$feature -t $threads
    done
fi

if [[ $tasks == *"a"* ]]; then
    mkdir -p $analysis_dir #todo dovnitr skriptu
    for feature in $features_list; do
        #todo check jestli ta slozka s hodnotami featury existuje
        python3 ${python_scripts_path}run_analysis.py -d $dataset_file -l $data_dir/lbs -v $features_dir/$feature -o $analysis_dir/$feature -f $feature -t $threads
    done
fi

#IFS=$oldIFS #todo je to potreba?

if [ -f ./run.log ]; then
    cp ./run.log ${data_dir}/run.log #todo zkopirovat i kdyz ten skript rpedtim spadne na chybu a nedobehne to sem
fi

if [[ $tasks == *"p"* ]]; then
    python3 ${python_scripts_path}create_prank_ds.py -d ${data_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_name}.ds

    python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_file -i $data_dir -o $data_dir/p2rank -t $threads -f $features_list
fi
