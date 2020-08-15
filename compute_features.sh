set -e

python_scripts_path=./pythonScripts/

OPTIND=1 # reset
# init parameters with defaults
dataset_file=""
input_dir=""
output_dir=""
strict=false
threads=1
features_list=""

show_help() {
    echo TODO PRINT HELP
}

### parse arguments:
while getopts ":h?d:i:o:t:f:s" opt; do
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
        i) input_dir=$OPTARG ;;
        o) output_dir=$OPTARG ;;
        t) threads=$OPTARG ;;
        f) features_list=$OPTARG ;;
        s) strict=true ;;
    esac
done
shift $((OPTIND - 1))

### check arguments:
[ -z "$dataset_file" ] && echo "ERROR: Dataset file path must be entered." && show_help && exit

[ -z "$input_dir" ] && echo "ERROR: Input directory must be entered." && show_help && exit
[ ! -d "$input_dir" ] && echo "ERROR: Input directory ${input_dir} does not exist." && show_help && exit

if [ -z "$output_dir" ]; then           # if output dir argument missing, create it in input dir location
    output_dir=${input_dir}/../features #todo check
    mkdir -p "$output_dir" && echo "INFO: Output directory ${output_dir} created."
fi
#todo check threads argument
#todo check features argument

oldIFS=$IFS
IFS=','
for feature in $features_list; do
    feature_output_dir=${output_dir}/${feature}
    if [ ! -d "$feature_output_dir" ]; then # if the directory does not exist yet, it is created.
        mkdir -p "$feature_output_dir" && echo "INFO: Output directory ${feature_output_dir} created."
    elif [ -n "$(ls -A ${feature_output_dir})" ]; then # it exists and it is not empty
        if [ "$strict" = true ]; then
            rm -rf "$feature_output_dir"
            mkdir -p "$feature_output_dir" && echo "INFO: Output directory ${feature_output_dir} was emptied."
        else # If it is nonempty and the -s (strict) option is not given, the program ends.
            echo "ERROR: The output directory ${feature_output_dir} is not empty and option -s was not given."
            show_help
            exit 1
        fi
    fi

    echo "INFO: Computing feature ${feature} started..."
    python3 ${python_scripts_path}compute_feature.py -f $feature -d $dataset_file -i $input_dir -o $feature_output_dir -t $threads
    echo "INFO: Computing feature ${feature} finished."
done

#restore IFS
IFS=$oldIFS #todo je to potreba?
