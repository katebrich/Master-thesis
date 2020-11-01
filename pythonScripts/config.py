import json
import os

config_path=os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.json")

with open(config_path, 'r') as config:
    parsed_config = json.load(config)

def get_feature_path(feature_name):
    return parsed_config["features"][feature_name]["import_path"]
