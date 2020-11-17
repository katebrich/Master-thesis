import getopt
import os
import shutil
import sys

import Logger

#log_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "run.log")
#remove old log
if os.path.exists(Logger.get_log_path()):
    os.remove(Logger.get_log_path())
logger = Logger.get_logger(os.path.basename(__file__))

from DatasetDownloader import DatasetDownloader
from LigandBindingSitesComputer import LigandBindingSitesComputer
from MappingsComputer import MappingsComputer
from FeaturesComputer import FeaturesComputer
from AnalysisComputer import AnalysisComputer
from Config import Config



#default values
threads=4
config_path=os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.json")
tasks="A"
features_list="x"
distance_threshold = 4
dataset_file=""
output_dir=""


P2Rank_path="/home/katebrich/Documents/diplomka/P2Rank"
dataset_name="test"
experiment_name="debug"
dataset_file=f"/home/katebrich/Documents/diplomka/GitHub/datasets/{dataset_name}.txt"
output_dir= f"{P2Rank_path}/datasets/{experiment_name}/{dataset_name}"
#features_list = "unp_disulfid"  #config.get_all_feature_names()       #"unp_PTM,unp_glycosylation,unp_lipidation,unp_mod_res,unp_variation,unp_topology,unp_sec_str,unp_non_standard,unp_natural_variant,unp_compbias,pdbekb_conservation,pdbekb_sec_str,aa,aa_pairs,hydropathy,polarity,polarity_binary,charged,aromaticity,mol_weight,H_bond_atoms,dynamine,efoldmine,mobiDB,HSE_up,HSE_down,exposureCN,bfactor,bfactor_CA,depth,phi_angle,psi_angle,cis_peptide"
features_list = "unp_disulfid" #"aa,aa_pairs,hydropathy,polarity,polarity_binary,charged,aromaticity,mol_weight,H_bond_atoms,HSE_up,HSE_down,exposureCN,bfactor,bfactor_CA,pdbekb_sec_str,pdbekb_conservation,dynamine,efoldmine,depth,mobiDB,phi_angle,psi_angle,cis_peptide,lbs,aa_ratio,conservation,unp_variation"


#parse arguments: #todo check
try:
    opts, args = getopt.getopt(sys.argv[1:], 'd:o:t:m:f:')
except getopt.GetoptError as err:
    logger.error(err) #unknown option or missing argument
    sys.exit(1)
for opt, arg in opts:
    if opt in ("-d", "--dataset"):
        dataset_file = arg
    elif opt in ("-o", "--output_dir"):
        output_dir = arg
    elif opt in ("-t", "--tasks"):
        tasks = arg #todo check possible values
    elif opt in ("-m", "--threads"): #todo check if threads >= 1, int
        threads = arg
    elif opt in ("-f", "--features"):
        features_list = arg

if (dataset_file == ""):
    logger.error("Dataset must be specified.")
    sys.exit(1)
    #todo print help
if (output_dir == ""):
    logger.error("Output directory must be specified.")
    sys.exit(1)
    #todo print help

config = Config(config_path)
if (features_list == "" or features_list == "x"): #todo
    features_list = config.get_all_feature_names()
else:
    features_list = features_list.split(',')

tasks=tasks.split(',') #todo check spravny format

downloads_dir = output_dir
pdb_dir = f"{output_dir}/PDB"
fasta_dir = f"{output_dir}/FASTA"
lbs_dir = f"{output_dir}/lbs"
mappings_dir = f"{output_dir}/mappings"
features_dir=f"{output_dir}/features"
analysis_dir=f"{output_dir}/analysis" #todo

#todo config parametrem, predavat rovnou tridam

features_computed = []
analysis_computed = []
processed = False
try:
    for t in tasks:
        if t != 'R' and t != 'D' and t != 'M' and t != 'L' and t != 'F' and t != 'A':
            raise ValueError(f"Unrecognized task '{t}'.")
            #todo print help

    #todo log command line parameters

    # features check
    for feature in features_list:
        if not config.is_feature_defined(feature):
            raise ValueError(f"Feature {feature} not defined in {config_path}")

    while processed == False:
        processed = True
        if ('R' in tasks): #todo jen debug
            if (os.path.exists(output_dir)):
                shutil.rmtree(output_dir)
            tasks.remove('R')
        if ('D' in tasks):
            dd = DatasetDownloader(dataset_file, downloads_dir)
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
            lc = LigandBindingSitesComputer(dataset_file, lbs_dir, mappings_dir, pdb_dir, distance_threshold)
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
            if not os.path.exists(features_dir):
                tasks.append('F')
                processed = False
                continue
            if processed == False:
                continue
            ac = AnalysisComputer(analysis_dir, lbs_dir, config)
            for feature in features_list:
                if not feature in analysis_computed:
                    feature_dir = os.path.join(features_dir, feature)
                    ac.run(feature, feature_dir)
                    analysis_computed.append(feature)
            tasks.remove('A')
            ac.write_summary()
except Exception as ex:
    logger.exception(f"{ex}", exc_info=True)
finally:
    #copy log to output_dir
    log_path = os.path.join(output_dir, "run.log")
    if os.path.exists(log_path):
        #append to existing log
        with open(log_path, "a") as f:
            with open(Logger.get_log_path(), "r") as log:
                f.write('\n')
                f.write("------------------------------------------------------------------------------------------")
                f.write('\n')
                f.write('\n')
                f.write(log.read())
    else:
        #create new log file in output directory
        shutil.copyfile(Logger.get_log_path(), log_path)


#todo prank
'''

if [[ $tasks == *"p"* ]]; then
    python3 ${python_scripts_path}create_prank_ds.py -d ${data_dir}/PDB -o ${P2Rank_path}/datasets/${dataset_name}.ds

    python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_file -i $data_dir -o $data_dir/p2rank -t $threads -f $features_list
fi

'''