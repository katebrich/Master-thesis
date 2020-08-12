python_scripts_path=./pythonScripts/

OPTIND=1 # reset
# init parameters with defaults
dataset_file=""
input_dir=""
output_dir=""
force=false

show_help() {
    echo TODO PRINT HELP
}

### parse arguments:
while getopts ":h?d:i:o:f" opt; do
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
        i) input_dir=$OPTARG ;;
        o) output_dir=$OPTARG ;;
        f) force=true ;;
    esac
done
shift $((OPTIND - 1))

### check arguments:
[ -z "$dataset_file" ] && echo "ERROR: Dataset file path must be entered." && show_help && exit

[ -z "$input_dir" ] && echo "ERROR: Input directory must be entered." && show_help && exit
[ ! -d "$input_dir" ] && echo "ERROR: Input directory ${input_dir} does not exist." && show_help && exit

if [ -z "$output_dir" ]; then      # if output dir argument missing, create it in input dir location
    output_dir=${input_dir}/../lbs #todo check
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

echo "INFO: Computing ligand binding sites started..."
python3 ${python_scripts_path}compute_ligand_binding_sites.py -d $dataset_file -i $input_dir -o $output_dir
echo "INFO: Computing ligand binding sites finished." #todo co kdyz tam je chyba? nepsat finished
