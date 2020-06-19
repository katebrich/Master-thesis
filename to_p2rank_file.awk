BEGIN { RS = "" ; FS = "\n" }
{
	out_file=(out_dir "/" $1 $2 ".pdb.csv")
	print "\"chain\",\"ins. code\",\"seq. code\",\"PTM\"" > out_file

	split($7, vals, ",");
	i=1
        for (res=$5; res<=$6; res++)
        {
		printf "\"%s\",,\"%s\",\"%s\"\n", $2, res, vals[i] > out_file
		i++
        }
	
	
}
