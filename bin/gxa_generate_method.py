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

seen = set()

for module, tools in data.items():
    # tools may be dict or other
    if isinstance(tools, dict):
        for tool, ver in tools.items():
            # normalize values to avoid duplicates caused by formatting differences
            tool_str = str(tool).strip() if tool is not None else ""
            ver_str = str(ver).strip() if ver is not None else ""
            line = f"{tool_str}\t{ver_str}"

            # skip completely empty entries
            if not tool_str and not ver_str:
                continue

            # ensure we don't print duplicate Tool\tVersion lines
            if line in seen:
                continue
            seen.add(line)

            print(line)
