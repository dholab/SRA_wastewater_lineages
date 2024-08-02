# SRA Wastewater Lineage Finder Refactor

### Refactor Roadmap

- [x] SRA accession fetching
- [x] SRA accession dereplicating and mapping
- [ ] Introduce more sophisticated data structures like dataclasses to reduce file I/O overall.
- [ ] Pure Python command line interface with type-checking and a PyPI-compatible executable entrypoint
- [ ] SAM refiner/haplotyping identification library/libraries
- [ ] Variant extractor library/libraries
- [ ] DuckDB database setup
- [ ] Potential distributed backend via Dask
- [ ] Project and environment management with Poetry or similar PyPI-compatible system
- [ ] Docker rebuild and push to Docker Hub

---

## Deprecated from prior to refactor

### Summary

1. This workflow will query SRA for the latest samples, then align and determine variants phased to a unique reads.
2. Unlike a typical VCF, this reports the unique reads in the sample. If a sample has multiple variants on the same read, it will be treated a collection of "phased variants" rather than individual variants.
3. It has been automated and has successfully processed ~140,000 SRA samples for wastewater.

### Workflow overview.

1. sra_query_daily.py to get the latest SRA entries and compare to already completed files
   - This ensures if the script is still running and you kick off another query, it will not double download files.
2. stage_sra_loop_daily.py
   - This adds a list of files needing to be queried, in batches to CHTC, based on the projected file sizes
   - This will send up to a set number of files per node (30 is default), and also if the total input file size per SRA query is under 5GB.
   - If the entry is >5GB it will only query up to one file per node.
   - This is because 1.) Large files can time out or take very long to process, beyond time limits of the nodes.
   - By limiting file sizes we mitigate this source of error from happening and affecting multiple SRA entries
3. sra_cryptic_loop.sh
   - This runs in CHTC and will run a single node workflow of the steps to download, align, and aggregate the unique sequence/ phased variants
   - This runs a hardcoded pipeline of the following scripts: derep.py, multivariant.poy, SAM_Refiner.py, SRA_fetch.py and Variant_extractor.py
   - This loops all SRA entries that were sent to the node.
4. Extract results_aggregate.py
   - Files are packaged as a tar or compressed tar (tar.gz)
   - These files need to be sorted by type and
5. Aggregate Results
   - This then takes the extracted results and condenses them into fewer files that are easier to back up and save.
6. Agggregate Unique Sequences
   - This will aggregate the unique sequences into a single file if it is above the 1.4% allele fraction threshold

### Running the workflow

- If you are not useing a orchestrator to run you pipeline you will only need to run the the modules/sra_cryptic_loop.sh script

```bash
sra_name=SRA12345678
/path/to/repository/modules/sra_cryptic_loop.sh -s ${sra_name}

```

- You can also use an orchestrator such as CHTC Serpent
  - https://github.com/DABAKER165/chtc_serpent_v2
- To run the orchestrator you will need access to a HPCcondor HPC platform.

### References

- BBMap
  - BBMap – Bushnell B. – sourceforge.net/projects/bbmap/
  - Deployed Version 39.01
- Minimap2
  - Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191
  - Deployed version minimap2-2.17
- Samtools
  - Twelve years of SAMtools and BCFtools
    Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li
    GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008
  - Deployed version 1.9

# Contributors

- David A. Baker
  - Automation
    - sra_query_daily.py
    - sra_cryptic_loop.sh
    - state_sra_loop_daily.py poller
    - chtc_serpent for HTC computing.
  - Aggregation
    - aggregate_results
    - extract_results_aggregate.py
- Devon Gregory
  - Download from SRA and Create unique sequencing files.
    - derep.py, multivariant.poy, SAM_Refiner.py, SRA_fetch.py and Variant_extractor.py were wrtitten by
