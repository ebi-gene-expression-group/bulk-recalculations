#!/usr/bin/env python3

import sys
import yaml

if len(sys.argv) == 3:
    input_sp = str(sys.argv[1])
    YAML_file = str(sys.argv[2])
else:
    print('Incorrect number of arguments passed')
    sys.exit(0)

def check_sp(input):
    with open(YAML_file , 'r') as stream:
        try:
            yaml_file=yaml.safe_load(stream)
            if yaml_file != None and input_sp in yaml_file['species']:
                # if input exists in dictionary then return corrected species name
                return ''.join(yaml_file['species'][input_sp])
            else:
                # otherwise return unmodified
                return input_sp
        except yaml.YAMLError as exc:
            print("There was a problem reading the YAML file")
            print(exc)

output_sp = check_sp(input_sp)

print(output_sp)

