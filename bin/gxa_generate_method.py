#!/usr/bin/env python3
# Generate a simple tsv formatted table from pipeline_info/nf_core_rnaseq_software_mqc_versions.yml of nf-core/rnaseq workflow

import yaml
import sys

if len(sys.argv) < 2:
    sys.stderr.write("Usage: yaml_to_tsv.py nf_core_rnaseq_software_mqc_versions.yml\n")
    sys.exit(1)

with open(sys.argv[1]) as f:
    data = yaml.safe_load(f)

# Output header
print("Tool\tVersion")

for module, tools in data.items():
    # Skip workflow block
    if module == "Workflow":
        continue

    # tools may be dict or other
    if isinstance(tools, dict):
        for tool, ver in tools.items():
            print(f"{tool}\t{ver}")
