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

if [ -z "$output_dir" ]; then      # if output dir argument missing, create it in input dir location
    output_dir=${input_dir}/../lbs #todo check
    mkdir -p "$output_dir" && echo "INFO: Output directory ${output_dir} created."
else                                # if the output directory argument was given
    if [ ! -d "$output_dir" ]; then # if the directory does not exist yet, it is created.
        mkdir -p "$output_dir" && echo "INFO: Output directory ${output_dir} created."
    elif [ -n "$(ls -A ${output_dir})" ]; then # it exists and it is not empty
        if [ "$strict" = true ]; then
            rm -rf "$output_dir"
            mkdir -p "$output_dir" && echo "INFO: Output directory ${output_dir} was emptied."
        else # If it is nonempty and the -s (strict) option is not given, the program ends.
            echo "ERROR: Given output directory ${output_dir} is not empty and option -s was not given."
            show_help
            exit 1
        fi
    fi
fi
#todo check threads argument
#todo check features
#todo jestli existuje vsechny potrebne slozky s features?

echo "INFO: Creating prank custom feature files..."
python3 ${python_scripts_path}create_p2rank_custom_feature_files.py -d $dataset_file -i $input_dir -o $output_dir -t $threads -f $features_list
echo "INFO: Creating prank custom feature files finished."
