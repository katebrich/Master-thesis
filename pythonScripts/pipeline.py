import os
import shutil
import logger

from DatasetDownloader import DatasetDownloader
from LigandBindingSitesComputer import LigandBindingSitesComputer
from MappingsComputer import MappingsComputer
from FeaturesComputer import FeaturesComputer
from AnalysisComputer import AnalysisComputer
from Config import Config


P2Rank_path="/home/katebrich/Documents/diplomka/P2Rank" #todo config

config_path="config.json"

dataset_name="test"
label="_11_07"
dataset_file=f"/home/katebrich/Documents/diplomka/GitHub/datasets/{dataset_name}.txt"
data_dir_name= f"{dataset_name}{label}"
output_dir= f"{P2Rank_path}/datasets/{data_dir_name}"

filter_ligands="p2rank"                              #todo none, water, small molecules, given IDs, MOAD,...

config = Config(config_path)
features_list = "unp_PTM"  #config.get_all_feature_names()       #"unp_PTM,unp_glycosylation,unp_lipidation,unp_mod_res,unp_variation,unp_topology,unp_sec_str,unp_non_standard,unp_natural_variant,unp_compbias,pdbekb_conservation,pdbekb_sec_str,aa,aa_pairs,hydropathy,polarity,polarity_binary,charged,aromaticity,mol_weight,H_bond_atoms,dynamine,efoldmine,mobiDB,HSE_up,HSE_down,exposureCN,bfactor,bfactor_CA,depth,phi_angle,psi_angle,cis_peptide"
features_list = features_list.split(',')

threads=4

downloads_dir = output_dir
pdb_dir = f"{output_dir}/PDB"
fasta_dir = f"{output_dir}/FASTA"
lbs_dir = f"{output_dir}/lbs"
mappings_dir = f"{output_dir}/mappings"
features_dir=f"{output_dir}/features"
analysis_dir=f"{output_dir}/analysis"

distance_threshold = 4
SASA_threshold = 0.1

remove = True

#todo config parametrem, predavat rovnou tridam

features_computed = []
analysis_computed = []
processed = False

tasks="L"
tasks=tasks.split(',')

log_file = "run.log"
#remove old log
os.remove(log_file)
logger = logger.get_logger(os.path.basename(__file__))

try:
    while processed == False:
        processed = True
        if ('R' in tasks): #todo jen debug
            if (remove and os.path.exists(output_dir)):
                shutil.rmtree(output_dir)
            tasks.remove('R')
        if ('D' in tasks):
            dd = DatasetDownloader(dataset_file, downloads_dir) #todo filter ligands?
            dd.run(threads)
            tasks.remove('D')
        if ('M' in tasks):
            mc = MappingsComputer(dataset_file, mappings_dir)
            mc.run(threads)
            tasks.remove('M')
        if ('L' in tasks):
            if not os.path.exists(pdb_dir):
                tasks.append('D')
                processed = False
                continue
            if not os.path.exists(mappings_dir):
                tasks.append('M')
                processed = False
                continue
            lc = LigandBindingSitesComputer(dataset_file, lbs_dir, mappings_dir, pdb_dir, distance_threshold, SASA_threshold)
            lc.run(threads)
            tasks.remove('L')
        if ('F' in tasks):
            if not os.path.exists(downloads_dir):
                tasks.append('D')
                processed = False
                continue
            if not os.path.exists(mappings_dir):
                tasks.append('M')
                processed = False
                continue
            fc = FeaturesComputer(dataset_file, downloads_dir, config)
            for feature in features_list:
                if not feature in features_computed:
                    feature_output_dir = os.path.join(features_dir, feature)
                    fc.run(feature, feature_output_dir,  threads)
                    features_computed.append(feature)
            tasks.remove('F')
        if ('A' in tasks):
            if not os.path.exists(lbs_dir):
                tasks.append('L')
                processed = False
                continue
            ac = AnalysisComputer(analysis_dir, lbs_dir, config)
            for feature in features_list:
                if not feature in analysis_computed:
                    feature_dir = os.path.join(features_dir, feature)
                    if not os.path.exists(feature_dir):
                        tasks.append('F')
                        processed = False
                        continue
                    ac.run(feature, feature_dir)
                    analysis_computed.append(feature)
            tasks.remove('A')
except Exception as ex:
    logger.exception(f"{ex}", exc_info=True)
finally:
    #copy log to output_dir
    #todo proc se nekopiruje?
    shutil.copyfile(log_file, os.path.join(output_dir, log_file)) #todo co kdyz nebude zadany output_dir?


'''

if [[ $tasks == *"p"* ]]; then
    python3 ${python_scripts_path}create_prank_ds.py -d ${data_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_name}.ds

    python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_file -i $data_dir -o $data_dir/p2rank -t $threads -f $features_list
fi

'''