python_scripts_path=./pythonScripts/

OPTIND=1 # reset
# init parameters with defaults
dataset_file=""
lbs_dir=""
feature_vals_dir=""
output_dir=""
strict=false
threads=1
feature_name=""

show_help() {
    echo TODO PRINT HELP
}

### parse arguments:
while getopts ":h?d:l:o:t:v:f:s" opt; do
    case "$opt" in
        h)
            show_help
            exit
            ;;
        \?)
            echo "ERROR: Unknown option $OPTARG"
            show_help
            exit 1
            ;;
        :)
            echo "ERROR: Missing argument for option $OPTARG"
            show_help
            exit 1
            ;;
        d) dataset_file=$OPTARG ;;
        l) lbs_dir=$OPTARG ;;
        v) feature_vals_dir=$OPTARG ;;
        o) output_dir=$OPTARG ;;
        f) feature_name=$OPTARG ;;
        t) threads=$OPTARG ;;
        s) strict=true ;;
    esac
done
shift $((OPTIND - 1))

### check arguments:
[ -z "$dataset_file" ] && echo "ERROR: Dataset file path must be entered." && show_help && exit

[ -z "$feature_name" ] && echo "ERROR: Feature name must be entered." && show_help && exit

[ -z "$lbs_dir" ] && echo "ERROR: Directory with ligand binding sites data must be entered." && show_help && exit
[ ! -d "$lbs_dir" ] && echo "ERROR: Directory with ligand binding sites data ${lbs_dir} does not exist." && show_help && exit

[ -z "$feature_vals_dir" ] && echo "ERROR: Directory with feature values must be entered." && show_help && exit
[ ! -d "$feature_vals_dir" ] && echo "ERROR: Directory with feature values ${feature_vals_dir} does not exist." && show_help && exit

if [ -z "$output_dir" ]; then # if output dir argument missing, create it in feature_vals_dir location
    #todo chyba..nebo label...?
    output_dir=${feature_vals_dir}/../analysis #todo check
    mkdir "$output_dir" && echo "INFO: Output directory ${output_dir} created."
else                                # if the output directory argument was given
    if [ ! -d "$output_dir" ]; then # if the directory does not exist yet, it is created.
        mkdir "$output_dir" && echo "INFO: Output directory ${output_dir} created."
    elif [ -n "$(ls -A ${output_dir})" ]; then # it exists and it is not empty
        if [ "$strict" = true ]; then
            rm -rf "$output_dir"
            mkdir "$output_dir" && echo "INFO: Output directory ${output_dir} was emptied."
        else # If it is nonempty and the -s (strict) option is not given, the program ends.
            echo "ERROR: Given output directory ${output_dir} is not empty and option -s was not given."
            show_help
            exit 1
        fi
    fi
fi
#todo check threads argument

echo "INFO: Computing feature ${feature_name} started..."
python3 ${python_scripts_path}run_analysis.py -d $dataset_file -l $lbs_dir -v $feature_vals_dir -o $output_dir -f $feature_name -t $threads
echo "INFO: Computing feature ${feature_name} finished."
