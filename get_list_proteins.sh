#!/bin/bash

usage() { echo "Usage: $0 [-d <string>]" 1>&2; exit 1; }

while getopts ":d:" o; do
    case "${o}" in
        d)
            d=${OPTARG} 
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

case "${d}" in
        coach420)
    	    cd /home/katebrich/Documents/diplomka/P2Rank/datasets/coach420 ;
	    for f in *.pdb ; do echo ${f: : -4} ; done 
            ;;
        *)
            echo unknown dataset
            ;;
esac

