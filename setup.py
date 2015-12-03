from setuptools import setup, find_packages
from glob import glob

import bio_bits

setup(
    name = bio_bits.__projectname__,
    version = bio_bits.__release__,
    packages = find_packages(),
    author = bio_bits.__authors__,
    author_email = bio_bits.__authoremails__,
    description = bio_bits.__description__,
    license = "GPLv2",
    keywords = "biopython split fasta concat",
    entry_points = {
        'console_scripts': [
            'rename_fasta = bio_bits.rename_fasta:main',
            'vcfcat = bio_bits.vcfcat_main:main',
            'amos2fastq  = bio_bits.amos2fastq_main:main',
            'group_references = bio_bits.group_references:main',
            'beast_checkpoint = bio_bits.beast_checkpoint:main',
            'beast_wrapper = bio_bits.beast_wrapper:beast_wrapper',
            'beast_est_time = bio_bits.beast_wrapper:beast_est_time',
            'ctleptop = bio_bits.ctleptop:main',
            'parallel_blast = bio_bits.parallel_blast:main',
            'version = bio_bits.version:main',
            'degen = bio_bits.degen:main'
            #'sequence_concat = bio_bits.sequence_concat:main',
            #'sequence_files_concat = bio_bits.sequence_files_concat:main',
            #'sequence_split = bio_bits_old.sequence_split:main',
            #'cat_sequences = bio_bits.cat_sequences:main',
            #'phyml_seqrename = bio_bits.phyml_seqrename:main',
            #'raxmlrunner = bio_bits.raxmlrunner:main',
            #'phymlrunner = bio_bits.phymlrunner:main',
            #'seaview_phyml_renamer = bio_bits.seaview_phyml_renamer:main',
        ],
    },
)
