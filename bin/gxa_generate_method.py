#!/usr/bin/env python3
# Generate a simple tsv formatted table from pipeline_info/nf_core_rnaseq_software_mqc_versions.yml of nf-core/rnaseq workflow

import yaml
import sys
import re

if len(sys.argv) < 2:
    sys.stderr.write("Usage: gxa_generate_method.py nf_core_rnaseq_software_mqc_versions.yml\n")
    sys.exit(1)

PIPELINE_ORDER = (
    "nf-core/rnaseq",
    "nextflow",
    "fastqc",
    "umi_tools",
    "fastp",
    "trim_galore",
    "cutadapt",
    "bbduk",
    "bbsplit",
    "sortmerna",
    "star",
    "hisat2",
    "bowtie2",
    "samtools",
    "picard",
    "rsem",
    "salmon",
    "kallisto",
    "subread",
    "featurecounts",
    "stringtie",
    "gffread",
    "bedtools",
    "qualimap",
    "rseqc",
    "preseq",
    "dupradar",
    "deeptools",
    "multiqc",
)

BASIC_TOOLS = {
    "awk",
    "bash",
    "bioconductor",
    "bzip2",
    "conda",
    "coreutils",
    "curl",
    "data",
    "description",
    "git",
    "grep",
    "gzip",
    "id",
    "java",
    "jq",
    "make",
    "markdown",
    "nextflowconfig",
    "numpy",
    "openjdk",
    "pandas",
    "parentid",
    "parentname",
    "perl",
    "pigz",
    "plottype",
    "python",
    "python3",
    "pythonabi",
    "pyyaml",
    "r",
    "rbase",
    "rscript",
    "sectionhref",
    "sectionname",
    "sed",
    "softwareversions",
    "tar",
    "wget",
    "xz",
    "yaml",
    "zlib",
}

BASIC_PREFIXES = (
    "bioconductor-",
    "lib",
    "perl-",
    "py-",
    "python-",
    "r-",
)

ORDER_INDEX = {
    re.sub(r"[^a-z0-9]+", "", tool.lower()): index
    for index, tool in enumerate(PIPELINE_ORDER)
}


def normalize_tool(tool):
    return re.sub(r"[^a-z0-9]+", "", str(tool).strip().lower())


def is_scalar(value):
    return value is None or isinstance(value, (str, int, float, bool))


def should_skip_tool(tool):
    tool_str = str(tool).strip()
    tool_lower = tool_str.lower()
    normalized = normalize_tool(tool_str)

    if not normalized or normalized in BASIC_TOOLS:
        return True

    return any(tool_lower.startswith(prefix) for prefix in BASIC_PREFIXES)


def iter_tool_versions(node):
    if isinstance(node, dict):
        for key, value in node.items():
            if is_scalar(value):
                yield key, value
            else:
                yield from iter_tool_versions(value)
    elif isinstance(node, list):
        for item in node:
            yield from iter_tool_versions(item)


with open(sys.argv[1]) as f:
    data = yaml.safe_load(f) or {}

# Output header
print("Tool\tVersion")

seen = set()
records = []

for first_seen, (tool, version) in enumerate(iter_tool_versions(data)):
    tool_str = str(tool).strip() if tool is not None else ""
    version_str = str(version).strip() if version is not None else ""
    normalized = normalize_tool(tool_str)

    if should_skip_tool(tool_str) or not version_str:
        continue

    if normalized in seen:
        continue
    seen.add(normalized)

    records.append(
        {
            "tool": tool_str,
            "version": version_str,
            "order": ORDER_INDEX.get(normalized, len(ORDER_INDEX)),
            "first_seen": first_seen,
        }
    )

for record in sorted(records, key=lambda item: (item["order"], item["first_seen"])):
    print(f"{record['tool']}\t{record['version']}")
