from pathlib import Path
from argparse import Namespace

# Store all programs and their requirements used by pipelines.
requirements = {
        'macs2': '2.1.1.20160309',
        'samtools': '1.9',
        'bcl2fastq': '2.20.0.422',
        "gcc": "9.2.0",
        "bedtools": "2.30.0-gcc-9.2.0",
        "fastqc": "0.11.9",
        'java': "jdk-17.0.7",
        'trimmomatic': '0.36',
        "picard": "2.8.0",
        "sambamba": "0.7.1"
    }

# Create namespace containing names for program and their respective path on O2
o2_paths = {}

for prog in requirements:
    version = requirements[prog]
    o2_paths[prog] = Path(f'/n/app/{prog}/{version}/bin/{prog}')

progs = Namespace(**o2_paths)

trimmomatic_jar_path = "/n/app/trimmomatic/0.36/bin/trimmomatic-0.36.jar"
