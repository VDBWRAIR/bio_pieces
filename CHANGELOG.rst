=========
CHANGELOG
=========

Version 1.2.0
-------------

* Renamed project to bio_bits to fix naming issue with other project
* GPL License added
* dgen_regions script added
* parallel_blast added

Version 1.1.0
-------------

* Renamed parse_contigs to group_references to better name functionality
* group_references now supports bam files

Version 1.0.0
-------------

* Version bump. Starting here we will employ semantic versioning
* Added version script to get version from project

Version 0.1.0
-------------

* Started project over to setup for Continuous Integration testing
* Added rename_fasta that can rename fasta sequence identifiers based
  on a input rename file
* Added travis, coveralls, readthedocs
* Added amos file parser that is specific to Ray assembler amos format
* Added format functionality for amos classes such that it is easy to
  convert to different formats
* Added amos2fastq to pull sequences out of AMOS files organized by their contigs.
* Added vcfcat.py, a commandline app for filtering and comparing vcf files.  
* Completed documentation for vcfcat
* Added beast_checkpoint script and documentation
* Added beast_wrapper script that prints estimated time column in beast output
* Added beast_est_time script that allows you to easily get estimated time left
  from already running beast run
