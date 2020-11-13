lbs_dir=""
feature_dir=""
output_dir=""
feature_type="binary"
feature_name="" #default - name of feature_dir folder
threads=4

#todo vytvorit config

python3 pipeline.py -o $output_dir -l $lbs_dir -f $feature_dir -v $feature_type -t $threads
