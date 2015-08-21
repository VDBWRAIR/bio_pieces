from setuptools import setup, find_packages
from glob import glob

import bio_pieces

setup(
    name = bio_pieces.__projectname__,
    version = bio_pieces.__release__,
    packages = find_packages(),
    author = bio_pieces.__authors__,
    author_email = bio_pieces.__authoremails__,
    description = bio_pieces.__description__,
    license = "GPLv2",
    keywords = "biopython split fasta concat",
    entry_points = {
        'console_scripts': [
            'rename_fasta = bio_pieces.rename_fasta:main',
            'vcfcat = bio_pieces.vcfcat_main:main',
            'amos2fastq  = bio_pieces.amos2fastq_main:main',
            'group_references = bio_pieces.group_references:main',
            'beast_checkpoint = bio_pieces.beast_checkpoint:main',
            'beast_wrapper = bio_pieces.beast_wrapper:beast_wrapper',
            'beast_est_time = bio_pieces.beast_wrapper:beast_est_time',
            'parallel_blast = bio_pieces.parallel_blast:main',
            'version = bio_pieces.version:main',
            'degen = bio_pieces.degen:main'
            #'sequence_concat = bio_pieces.sequence_concat:main',
            #'sequence_files_concat = bio_pieces.sequence_files_concat:main',
            #'sequence_split = bio_pieces_old.sequence_split:main',
            #'cat_sequences = bio_pieces.cat_sequences:main',
            #'phyml_seqrename = bio_pieces.phyml_seqrename:main',
            #'raxmlrunner = bio_pieces.raxmlrunner:main',
            #'phymlrunner = bio_pieces.phymlrunner:main',
            #'seaview_phyml_renamer = bio_pieces.seaview_phyml_renamer:main',
        ],
    },
)
