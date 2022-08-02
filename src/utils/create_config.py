import sys
import os
import glob
import re
import yaml
from collections import defaultdict
from collections import OrderedDict
from ruamel.yaml import YAML

run_id = sys.argv[1]
folder = sys.argv[2]
files=glob.glob(folder + "/*R1*")

yaml_dict2 = defaultdict(lambda: defaultdict(dict))

for R1 in files:
    # get fastq
    R2 = re.sub("R1","R2", R1)
    UMI = re.sub("_R1", ".UMI", R1).replace(".gz", "")
    file = os.path.basename(R1)
    S_ID = file.replace("_R1.fastq.gz","")
    # create dictionary
    yaml_dict2['Run']=run_id
    yaml_dict2['Samples'][S_ID]['fwd'] = R1
    yaml_dict2['Samples'][S_ID]['rev'] = R2

# sort diction
yaml_dict2 = OrderedDict(yaml_dict2)

# create and save yaml file
#with open("/mnt/nfs/shared/000000-Shared/groups/MedGen/TP53_smk/src/config/config_" + run_id + ".yaml", 'w') as outfile:
#    yaml.dump(yaml_dict2, outfile, default_flow_style=False)
with open("./config/config_" + run_id + ".yaml", 'w') as outfile:
    yaml.dump(yaml_dict2, outfile, default_flow_style=False)
