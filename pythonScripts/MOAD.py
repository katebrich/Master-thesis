import csv
import os
import uuid

import pandas as pd
from helper import restAPI_get

class MOAD:
    MOAD_file_path = "moad.csv" #todo f"temp_moad_{uuid.uuid1()}"
    moad_dict = {}

    def __init__(self):
        self.__download_MOAD_file() #todo

    def __del__(self):
        pass
        #if os.path.exists(self.MOAD_file_path):
        #    os.remove(self.MOAD_file_path)

    def __download_MOAD_file(self):
        url = "http://bindingmoad.org/files/csv/every.csv"
        #response = restAPI_get(url)
        #with open(self.MOAD_file_path, 'wb') as file:
        #    file.write(response)
        with open(self.MOAD_file_path, 'r') as file:
            csv_reader = csv.reader(file, delimiter=',')
            pdb_id = ""
            ligands = []
            for row in csv_reader:
                if (row[2] != ""): #beginning of new structure
                    if (pdb_id != ""): #save previous structure and clear variables
                        self.moad_dict[pdb_id] = ligands
                    ligands = []
                    pdb_id = row[2]
                elif (row[3] != ""): #ligand; append to current structure
                    ligands.append((row[3], row[4])) # (ligandCode:chain:residueNumber, validity)
                else: #save structure and clear variables
                    if (pdb_id != ""):
                        self.moad_dict[pdb_id] = ligands
                    ligands = []
                    pdb_id = ""


    def get_relevant_ligands(self, pdb_id, chain_id):
        try:
            ligands = self.moad_dict[pdb_id.upper()]
        except KeyError: #Excluded from MOAD - no valid ligands or not an x-ray structure or low resolution etc.
            return None
        result = []
        for l in ligands:
            validity = l[1]
            if (validity == "valid"):
                res_code, chain, res_num = l[0].split(':')
                if (chain == chain_id):
                    result.append((res_code, res_num))
        return result
