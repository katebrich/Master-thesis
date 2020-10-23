#!/bin/bash

docker run -v /home/katebrich/Documents/diplomka/data/conservation:/data/conservation -u $(id -u ${USER}):$(id -g ${USER}) -it conservation

cd /data/conservation

for f in ./FASTA/*.fasta; do
    file_name="${f##*/}"        #with extension
    file_name="${file_name%.*}" #without extension
    file_name=${file_name}.json
    if [ ! -e ./results/$file_name ]; then
        echo python3 /opt/conservation/calculate_conservation.py --input $f --output ./temp/$file_name
        echo mv ./temp/$file_name ./results/$file_name
    fi
done
