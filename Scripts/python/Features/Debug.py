import os

import numpy as np


class LBS():
    def get_values(self, data_dir, pdb_id, chain_id):
        # get ligand binding sites values
        file = os.path.join(data_dir, "lbs", f"{pdb_id}{chain_id}.txt")
        lbs = tuple(map(tuple, np.genfromtxt(file, delimiter=' ', dtype=None)))
        return lbs