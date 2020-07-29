# Downloads all structures defined in the dataset_file.
# Creates a directory with output_dir path with all the downloaded data.
# If the output directory is not defined, it is created in the dataset_file location. It has a name as the dataset_file plus unique uuid.
# If the output directory is defined and it does not exist, it is created.
# If the output directory is defined, it exists, it is nonempty and the -f (force) option is not given, the program ends with an error.
# Ligands filtering...

python_scripts_path=./pythonScripts/

OPTIND=1 # reset
# init parameters with defaults
dataset_file=""
filter_ligands=none #todo none, water, small molecules, given IDs, MOAD,...
output_dir=""
force=false

show_help() {
    echo TODO PRINT HELP
}

### parse arguments:
while getopts ":h?d:o:l:f" opt; do
    case "$opt" in
        h)
            show_help
            exit
            ;;
        \?)
            echo "ERROR: Unknown option $OPTARG"
            show_help
            exit
            ;;
        :)
            echo "ERROR: Missing argument for option $OPTARG"
            show_help
            exit
            ;;
        d) dataset_file=$OPTARG ;;
        o) output_dir=$OPTARG ;;
        l) filter_ligands=$OPTARG ;;
        f) force=true ;;
    esac
done
shift $((OPTIND - 1))

### check arguments:
#if dataset file missing, exit
[ -z "$dataset_file" ] && echo "ERROR: Dataset file path must be entered." && show_help && exit

if [ -z "$output_dir" ]; then # if output dir argument missing, create it in dataset file location with unique name
    output_dir=${dataset_file%.*}_$(uuidgen)
    mkdir "$output_dir" && echo "INFO: Output directory ${output_dir} created."
else                                # if the output directory argument was given
    if [ ! -d "$output_dir" ]; then # if the directory does not exist yet, it is created.
        mkdir "$output_dir" && echo "INFO: Output directory ${output_dir} created."
    elif [ -n "$(ls -A ${output_dir})" ]; then # it exists and it is not empty
        if [ "$force" = true ]; then
            rm -rf "$output_dir"
            mkdir "$output_dir" && echo "INFO: Output directory ${output_dir} created."
        else # If it is nonempty and the -f (force) option is not given, the program ends.
            echo "ERROR: Given output directory is not empty and option -f was not given."
            show_help
            exit
        fi
    fi
fi

#todo check ligands option

echo "INFO: Downloading dataset ${dataset_file} to ${output_dir} started..."

python3 ${python_scripts_path}download_dataset.py -d $dataset_file -o $output_dir -l $ligands
