from pathlib import Path


class Aligner:

    requires = {
        "gcc": "6.2.0",
        "samtools": "0.1.19",
        "bowtie2": "2.2.9",
        "bedtools": "2.26.0",
        "fastqc": "0.11.5",
        "java": "jdk-1.8u112",
        "trimmomatic": "0.36",
        "picard": "2.8.0",
        "sambamba": "0.7.1"
    }

    def atac_alignment_script(self, r1_fastq: Path, r2_fastq: Path) -> str:
        """
        Gets an alignment script to perform
        :param fastq_path:
        :param assay:
        :return:
        """

        s = f"""
        
        mkdir
        java -jar {trimmomatic_jar_path} PE {r1_fastq} {r2_fastq} ${file/R1/P1} ${file/R1/U1} ${file/R1/P2} ${file/R1/U2} ILLUMINACLIP:/n/app/trimmomatic/0.36/bin/adapters/TruSeq3-PE.fa:2:30:10 CROP:30
        
        
        """




