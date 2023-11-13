# Apollo-variant-typing
Interpretation of variants identified in fungal genomes.

## Pipeline information
* **Author(s):**            Boas van der Putten, Roxanne Wolthuis
* **Organization:**         Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
* **Department:**           Infektieziekteonderzoek, Diagnostiek en Laboratorium Surveillance (IDS), Bacteriologie (BPD)
* **Start date:**           07 - 04 - 2023
* **Commissioned by:**      Thijs Bosch, Julia (Jianhua) Zhang

## About this project
Apollo-variant-typing uses output generated with apollo-mapping to type *Candida auris* and *Aspergillus fumigatus* genomes.

*Candida auris* typing currently includes:
- Extraction of AMR mutations, based on reference data present in `files/cauris/resistance_list.csv`.
- Prediction of *C. auris* clade using [`auriclass`](https://github.com/rivm-bioinformatics/auriclass).

*Aspergillus fumigatus* typing will be implemented in the near future.

The apollo-variant-typing pipeline is created with the [juno-template](https://github.com/RIVM-bioinformatics/juno-template) and [juno-library](https://github.com/RIVM-bioinformatics/juno-library).

The pipeline uses the following tools:
1. [bcftools](https://samtools.github.io/bcftools/bcftools.html): predict the effects of genomic variants using the `csq` subcommand.
2. [GATK](https://gatk.broadinstitute.org/hc/en-us): comprehensive toolkit for genomic variant analysis. Currently, only the VariantsToTable tool is used to convert VCF to tab-separated values.
3. [Picard](https://broadinstitute.github.io/picard/): comprehensive toolkit for NGS data handling. Currently, only the `SamToFastq` tool is used to extract FastQ data from a BAM file. 
4. [AuriClass](https://github.com/rivm-bioinformatics/auriclass): to predict *Candida auris* clade from FastQ data.

## Prerequisities
* Linux environment
* mamba
* Python 3.11

## Installation
1. Clone the repository and enter the cloned directory
```
git clone https://github.com/RIVM-bioinformatics/apollo-variant-typing.git
cd apollo-variant-typing
```
2. Create & activate mamba environment.
```
mamba env update -f envs/apollo_variant_typing.yaml
```

3. Create & activate apollo environment.
```
conda activate apollo_variant_typing
```

4. Example of run:
```
python apollo_variant_typing.py -i [input] -o [output] -s [species]
```

## Parameters & Usage
```
usage: apollo_variant_typing.py [-h] -i DIR [-o DIR] [-w DIR] [-ex FILE] [-p PATH] [-l] [-tl INT] [-u] [-n] [-q QUEUE] [--no-containers] [--snakemake-args [SNAKEMAKE_ARGS ...]] [-m FILE]
                                [-s GENUS SPECIES] [-d DIR] [--presets-path PATH]

Apollo-variant-typing for interpretation of variants identified in fungal genomes.

options:
  -h, --help            show this help message and exit
  -i DIR, --input DIR   Relative or absolute path to the input directory. It must contain all the raw reads (fastq) files for all samples to be processed (not in subfolders).
  -o DIR, --output DIR  Relative or absolute path to the output directory. If none is given, an 'output' directory will be created in the current directory.
  -w DIR, --workdir DIR
                        Relative or absolute path to the working directory. If none is given, the current directory is used.
  -ex FILE, --exclusionfile FILE
                        Path to the file that contains samplenames to be excluded.
  -p PATH, --prefix PATH
                        Conda or singularity prefix. Basically a path to the place where you want to store the conda environments or the singularity images.
  -l, --local           If this flag is present, the pipeline will be run locally (not attempting to send the jobs to an HPC cluster**). The default is to assume that you are working on a cluster.
                        **Note that currently only LSF clusters are supported.
  -tl INT, --time-limit INT
                        Time limit per job in minutes (passed as -W argument to bsub). Jobs will be killed if not finished in this time.
  -u, --unlock          Unlock output directory (passed to snakemake).
  -n, --dryrun          Dry run printing steps to be taken in the pipeline without actually running it (passed to snakemake).
  -q QUEUE, --queue QUEUE
                        Name of the queue that the job will be submitted to if working on a cluster.
  --no-containers       Use conda environments instead of containers.
  --snakemake-args [SNAKEMAKE_ARGS ...]
                        Extra arguments to be passed to snakemake API (https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html).
  -m FILE, --metadata FILE
                        Relative or absolute path to the metadata csv file. If provided, it must contain at least one column named 'sample' with the name of the sample (same than file name but
                        removing the suffix _R1.fastq.gz), a column called 'genus' and a column called 'species'. The genus and species provided will be used to choose the serotyper and the MLST
                        schema(s).If a metadata file is provided, it will overwrite the --species argument for the samples present in the metadata file.
  -s GENUS SPECIES, --species GENUS SPECIES
                        Species name (any species in the metadata file will overwrite this argument). It should be given as two words (e.g. --species Candida auris)
  -d DIR, --db_dir DIR  Relative or absolute path to the directory that contains the databases for all the tools used in this pipeline or where they should be downloaded. Default is:
                        /mnt/db/apollo/variant-typing
  --presets-path PATH   Relative or absolute path to custom presets.yaml to use. If none is provided, the default (config/presets.yaml) is used.
```

## Explanation of the output
* **cauris_typing** (if *C. auris* was analysed): Files containing *C. auris*-specific typing results, such as AMR mutation reports and clade predictions.
* **audit_trail**: Logs of conda, git and the pipeline, a sample sheet, the used parameters and a snakemake report.
* **prepared_files**: Generated indices of files, necessary for typing analyses.
* **log**: Log with output and error file from the cluster for each Snakemake rule/step that is performed.


## Issues
* The default confifuration of this pipeline only works on the RIVM cluster. Paths to reference data can be specified on the command line. Cluster integration is currently only implemented for IBM LSF (`bsub`). To run without submission of jobs to a cluster, specify `--local` on the command line.

## License
This pipeline is licensed with a AGPL3 license. Detailed information can be found inside the 'LICENSE' file in this repository.

## Contact
* **Contact person:**       IDS-Bioinformatics
* **Email:**                ids-bioinformatics@rivm.nl  

## Contribution guidelines
Apollo pipelines use a [feature branch workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow). To work on features, create a branch from the `main` branch to make changes to. This branch can be merged to the main branch via a pull request. Hotfixes for bugs can be committed to the `main` branch.

Please adhere to the [conventional commits](https://www.conventionalcommits.org/) specification for commit messages. These commit messages can be picked up by [release please](https://github.com/googleapis/release-please) to create meaningful release messages.
