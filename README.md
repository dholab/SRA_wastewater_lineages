# Summary
An automated pipeline to download and aggregate Sars-cov-2 wastewater sequencing data from National Center for Biotechnology Information's (NCBI) Sequence Read Archive (SRA) and discover variants of concern within the data and lineages
- This pipeline leverages the University of Wisconsin-Madison's distributed computing platform, CHTC (Center for High Throughput Computing), to quickly retrieve and analyze results
- As of December 20, 2022, more than 38,000 lineages have been downloaded and analyzed with this pipeline.
# Prerequisites:
## Python Prerequisite libraries:
- pandas>=1.1.5
  - https://pandas.pydata.org/
  - Version used: 1.1.5
## additional installed executables
These should in the $PATH environment variable directory(s)
- minimap2>=2.17
  - https://github.com/lh3/minimap2
  - Package used: minimap2-2.17.tar.bz2
- samtools>=1.9
  - http://www.htslib.org/
  - Package used: samtools-1.9.tar.bz2
- histlab>=1.9
  - http://www.htslib.org/
- bbmap>=39.01
  - https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/
  - Package used: BBMap_39.01.tar.gz
## Workflow manager
- Remote pipeline management system for CHTC and remote servers.
  - https://github.com/DABAKER165/chtc_serpent_v2
  - git checkout urL: https://github.com/DABAKER165/chtc_serpent_v2.git

# SRA CRYPTIC CHTC driven Auto Downloader and Aggregator
There are 5 components :
- CHTC serpent (for managing jobs) (should be dual authentication compatible)
- SRA_Query (to get era number)
  - This is the SRA Poller.
  - It can be set up to run every x hours and get new data
  - Unfortunately you must download all to get the last files, so you have a 1GB file to get the last 100 samples due to the non-standard database set up of SRA.
- SRA Stager
  - Based on the query, figures out what already ran, needs to be run and the projected size of the files
  - Creates batch jobs, based on file size that will optimize node usage, with upto 30 files running across 25 nodes at a time.
  - It limits the SRA files to 5 GB for multiple files per node submissions, or 1 file if the file size > 5GB
- SRA Downloader that runs in CHTC,
  - Downloads the SRA files
- Variant/CRAM/Fasta deduplication (CHTC)
  - then runs a deduplication script, an alignment, and a variant extractor on the same CHTC node as when it downloads it.
  - Exports the results in compressed formats (tsv.gz, CRAM)
  - Tars the results (as the files are compressed.
- Aggregator unpacks tar/compressed files then identifies VOC's (variant of concerns) based on hard coded inputs of the files
  - This packs the files into different folders (so we don't hit too many files per folder).
  - It then runs a variant of concern extractor.

# key files:
## SRA QUERY
### launch_sra_query.sh & sra_query.py
- This queries from the SRA NCBI server the SRA accession numbers with wastewater (spelled correctly)
- This is ran once like follows:
```
python3 ~/github/sra_cryptic_lineages/sra_poller/sra_query.py \
--out_dir=/Volumes/T7/sra_cryptic_loop/sra_query
```
###sra_query folder:
- Done_SRAs.csv.gz
  - This was what was done by Marc Johnson and Gregory Devon of Missouri
- downloaded_list.txt
  - this is the individual SRA numbers that have been downloaded
- downloaded_group_list.txt
  - This is what group of files has been downloaded
- SRA_meta.tsv
  - this is the file created by the sra_query that has key info (size location) of the SRA download.
## config files
Th is is what serpent uses to run it.

### config_example.json
	- Rename this file something you can find (SRA_121213.json)
	- Change the file path in this file to where your config file is for cryptic lineages
	- move this file to your serpent config folder i.e. /Users/username/chtc_serpent/config/SRA_221213.json
```
{
  "all": {
    "local": {
      "default_config_dir": "/Users/username/github/sra_cryptic_lineages/config"
    }
  }
}
```

### chtc.json
- change any paths with "##" to where your paths are and remove the ##
### default.json
- change any paths with "##" to where your paths are and remove the ##

## Post Aggergation and extraction:
### extract_results_aggregate.py
- Run this after your chtc results are finished.
- This untars the results files
- it is triggered automatically as part of the CHTC serpent workflow

###  extract_multivariants.py
- run after extract_results_aggregate.py finishes
- it can be run as long as the files are completely created (not partial)
- check the command output to make sure the workflow is not running (is sleeping) and stopped (ctrl-c) prior to running
```
python3 ~/github/sra_cryptic_lineages/modules/extract_multivariants.py \
--untar_dir=/Volumes/T7/serpent/sra_cryptic_lineages/SRA_221213/out_files/extract_results \
--multivariant_py_path=~/github/sra_cryptic_lineages/static_files/multivariant.py
```

# Acknowledgments:
## Software:
minimap2:
- Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37:4572-4574. doi:10.1093/bioinformatics/btab705

samtools:
- Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]

BBMap:
- Bushnell, Brian. BBMap: A Fast, Accurate, Splice-Aware Aligner. United States: N. p., 2014.

### Python packages
Pandas: 
```
@software{reback2020pandas,
    author       = {The pandas development team},
    title        = {pandas-dev/pandas: Pandas},
    month        = feb,
    year         = 2020,
    publisher    = {Zenodo},
    version      = {latest},
    doi          = {10.5281/zenodo.3509134},
    url          = {https://doi.org/10.5281/zenodo.3509134}
}
@InProceedings{ mckinney-proc-scipy-2010,
  author    = { {W}es {M}c{K}inney },
  title     = { {D}ata {S}tructures for {S}tatistical {C}omputing in {P}ython },
  booktitle = { {P}roceedings of the 9th {P}ython in {S}cience {C}onference },
  pages     = { 56 - 61 },
  year      = { 2010 },
  editor    = { {S}t\'efan van der {W}alt and {J}arrod {M}illman },
  doi       = { 10.25080/Majora-92bf1922-00a }
}
```
ncbi/sra-tools:
- https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
- SRA Toolkit Development Team

## Development team:
Baker, David A <sup>1*</sup>, Gregory, Devon A <sup>2*</sup> Johnson, Marc <sup>2</sup>, O'Connor David H. <sup>1</sup>
- <sup>*</sup> denotes most significant contributions
- <sup>1</sup> University of Wisconsin-Madison, School of Public Health and Medicine
- <sup>2</sup> University of Missouri, School of Medicine
