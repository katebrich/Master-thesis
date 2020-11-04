import json
import os

class Config():
    parsed_config = ""

    def __init__(self, config_path):
        with open(config_path, 'r') as config:
            self.parsed_config = json.load(config)

    def get_feature_function(self, feature_name):
        return self.parsed_config["features"][feature_name]["import_path"]

    def get_all_feature_names(self):
        return list((self.parsed_config["features"]).keys()) #todo

    def get_feature_type(self, feature_name):
        return self.parsed_config["features"][feature_name]["type"]
