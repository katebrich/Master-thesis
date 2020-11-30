import getopt
import os
import shutil
import sys
import Logger

from AnalysisPipeline.DatasetDownloader import DatasetDownloader
from AnalysisPipeline.LigandBindingSitesComputer import LigandBindingSitesComputer
from AnalysisPipeline.MappingsComputer import MappingsComputer
from AnalysisPipeline.FeaturesComputer import FeaturesComputer
from AnalysisPipeline.AnalysisComputer import AnalysisComputer
from Config import Config

#remove old log
if os.path.exists(Logger.get_log_path()):
    os.remove(Logger.get_log_path())
logger = Logger.get_logger(os.path.basename(__file__))

#default values
threads=1
config_path=os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.json")
tasks="A"
features_list="."
lbs_distance_threshold = 4 # in Angstroms
dataset_file=""
output_dir=""
sample_size = 0 # 0 = do not sample, take all available data
iterations = 1
balance_binding_ratio = False # if True, take the same sample size for binding/nonbinding residues, regardless of their ratio in the original dataset
draw_plots = True
alpha = 0.05

#todo smazat
'''
P2Rank_path="/home/katebrich/Documents/diplomka/P2Rank"
dataset_name="mix_filter_MOAD"
tasks="A"
tasks = tasks.split(',')
dataset_file=f"/home/katebrich/Documents/diplomka/GitHub/datasets/{dataset_name}.txt"
output_dir= f"{P2Rank_path}/datasets/{dataset_name}"
features_list = "x" #"aa,hydropathy,polarity,polarity_binary,charged,aromaticity,mol_weight,H_bond_atoms,HSE_up,HSE_down,exposureCN,bfactor,bfactor_CA,pdbekb_sec_str,pdbekb_conservation,dynamine,efoldmine,depth,mobiDB,phi_angle,psi_angle,cis_peptide,lbs,aa_ratio,conservation,unp_variation"
sample_size = 500
iterations = 1000
balance_binding_ratio = True
'''

#parse arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], 'd:o:t:m:f:s:i:a:p:b:c:l:')
except getopt.GetoptError as err:
    #todo print help
    logger.error(err) #unknown option or missing argument
    sys.exit(1)

try:
    for opt, arg in opts:
        if opt in ("-d", "--dataset"):
            dataset_file = arg
        elif opt in ("-o", "--output_dir"):
            output_dir = arg
        elif opt in ("-t", "--tasks"):
            tasks = arg.upper().split(',')
        elif opt in ("-m", "--threads"):
            threads = int(arg)
        elif opt in ("-f", "--features"):
            features_list = arg
        elif opt in ("-s", "--sample_size"):
            sample_size = int(arg)
        elif opt in ("-i", "--iterations"):
            iterations = int(arg)
        elif opt in ("-c", "--config_path"): #todo check if it works
            config_path = arg
        elif opt in ("-l", "--lbs_distance_threshold"):
            lbs_distance_threshold = float(arg) #todo check if decimal works
        elif opt in ("-b", "--balance_binding_ratio"):
            if arg == '0' or arg == "false" or arg == "False":
                balance_binding_ratio = False
            elif arg == '1' or arg == "true" or arg == "True":
                balance_binding_ratio = True
            else:
                raise ValueError(f"Option '-b' ('--balance_binding_ratio') has invalid argument: '{arg}'. Possible values: 0, 1, true, false, True, False.")
        elif opt in ("-p", "--draw_plots"):
            if arg == '0' or arg == "false" or arg == "False":
                draw_plots = False
            elif arg == '1' or arg == "true" or arg == "True":
                draw_plots = True
            else:
                raise ValueError(
                    f"Option '-p' ('--draw_plots') has invalid argument: '{arg}'. Possible values: 0, 1, true, false, True, False.")
        elif opt in ("-a", "--alpha"):
            alpha = int(arg)

    if (dataset_file == ""):
        raise ValueError("Argument '-d' ('--dataset') is compulsory.")
    if (output_dir == ""):
        raise ValueError("Argument '-o' ('--output_dir') is compulsory.")

    #check valid tasks values
    for t in tasks:
        if t != 'D' and t != 'M' and t != 'L' and t != 'F' and t != 'A':
            raise ValueError(f"Unrecognized task '{t}'.")

    #init config
    config = Config(config_path)

    #parse features
    if (features_list == "" or features_list == "."):
        features_list = config.get_all_feature_names() #if no specific features, take all from the config
    else:
        features_list = features_list.split(',')
        # features check
        for feature in features_list:
            if not config.is_feature_defined(feature):
                raise ValueError(f"Feature {feature} not defined in {config_path}")

    #set directories structure
    downloads_dir = output_dir
    pdb_dir = f"{output_dir}/PDB"
    fasta_dir = f"{output_dir}/FASTA"
    lbs_dir = f"{output_dir}/lbs"
    mappings_dir = f"{output_dir}/mappings"
    features_dir=f"{output_dir}/features"
    analysis_dir=f"{output_dir}/analysis_{sample_size}"

    features_computed = []
    analysis_computed = []
    processed = False

    while processed == False:
        processed = True
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
            lc = LigandBindingSitesComputer(dataset_file, lbs_dir, mappings_dir, pdb_dir, lbs_distance_threshold)
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
            ac = AnalysisComputer(analysis_dir, lbs_dir, features_dir, features_list, config)
            ac.run(sample_size, iterations, balance_binding_ratio, draw_plots, alpha, threads)
            ac.write_summary()
            tasks.remove('A')
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
