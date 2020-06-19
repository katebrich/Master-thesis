output_dir="PTM_data"
rm -r $output_dir 
mkdir $output_dir

bash get_list_proteins.sh -d coach420 | xargs -I% python pythonScripts/get_PTM.py -p % | awk -v out_dir="PTM_data" -f to_p2rank_file.awk 
