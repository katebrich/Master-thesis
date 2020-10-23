set -e
python_scripts_path=./pythonScripts/

#todo parsovani argumentu
dataset_name="chen11"
dataset_file="/home/katebrich/Documents/diplomka/datasets/${dataset_name}.txt" #povinny
data_dir="/home/katebrich/Documents/diplomka/datasets/${dataset_name}_10_23"   #nepovinny. Kdyz neni zadan, vytvori se v umisteni dataset_file
features_dir="${data_dir}/features"                                            #nepovinny. Default podslozka data_dir
filter_ligands=p2rank                                                          #todo none, water, small molecules, given IDs, MOAD,...
features_list=dynamine_website                                                 #unp_PTM,unp_glycosylation,unp_lipidation,unp_mod_res,unp_variation,unp_topology,unp_sec_str,unp_non_standard,unp_natural_variant,unp_compbias,pdbekb_conservation,aa,aa_pairs,hydropathy,polarity,polarity_binary,charged,aromaticity,mol_weight,H_bond_atoms,dynamine_website,dynamine_funPDBe,efoldmine_funPDBe,mobiDB,HSE,HSE_down,exposureCN,bfactor,bfactor_Calpha,depth,phi_angle,psi_angle,cis_peptide
threads=4
#log_file="${data_dir}/run.log"
#strict=true #todo
analysis_dir="${data_dir}/analysis"

rm -f ./run.log

#if [ -d "$data_dir" ]; then
#    rm -rf "$data_dir"
#fi

#python3 ${python_scripts_path}download_dataset.py -d $dataset_file -o $data_dir -l $filter_ligands -t $threads

#python3 ${python_scripts_path}create_mappings_cache.py -d $dataset_file -o $data_dir/mappings -t $threads

#python3 ${python_scripts_path}compute_ligand_binding_sites.py -d $dataset_file -i ${data_dir} -t $threads

oldIFS=$IFS
IFS=','

for feature in $features_list; do
    python3 ${python_scripts_path}compute_feature.py -f $feature -d $dataset_file -i $data_dir -o $features_dir/$feature -t $threads
done

mkdir -p $analysis_dir #todo dovnitr skriptu

for feature in $features_list; do
    #todo check jestli ta slozka s hodnotami featury existuje
    python3 ${python_scripts_path}run_analysis.py -d $dataset_file -l $data_dir/lbs -v $features_dir/$feature -o $analysis_dir/$feature -f $feature -t $threads
done

IFS=$oldIFS #todo je to potreba?

cp ./run.log ${data_dir}/run.log #todo zkopirovat i kdyz ten skript rpedtim spadne na chybu a nedobehne to sem
