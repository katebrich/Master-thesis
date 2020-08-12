oldIFS=$IFS
IFS=','
for feature in $features; do
    echo "INFO: Computing feature ${feature} started..."
    python3 ${python_scripts_path}compute_ligand_binding_sites.py #todo
    echo "INFO: Computing feature ${feature} finished."
done

#restore IFS
IFS=$oldIFS #todo je to potreba?
