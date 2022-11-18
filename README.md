# Artemis.jl benchmark

This is a benchmark of speed for finding off-target sites for a given gRNA.
Currently this pipeline builds and searches for off-target sites within the whole genome for four software packages: 

* ARTEMIS.jl
* CRISPRITZ
* Cas-offinder
* FlashFry

# Installation 

To execute the pipeline first build an environment (as per snakemake recommendations this step uses mamba):

`mamba env create --name benchmark-pipeline --file environment.yaml`

Activate the environment:

`conda activate benchmark-pipeline`

To run the pipeline end-to-end use the command:

`snakemake --cores 1`

To run specific rules consult the rule you wish to execute and copy one of it's output files and run:

`snakemake --cores 1 path/to/output.file`

# Requirements

The data/ folder should contain three files before running the pipeline. Please make sure to download hg38v34.fa and place it here. We provide files curated_guides_wo_PAM.txt (list of guides to find off-targets for) and 20bp-NGG-SpCas9.txt (a template for guide format required by CRISPRITZ), while the rest of remaining required input files are generated. 

Due to the time it takes to execute several of the steps of the pipeline, by default it only runs search for distance = 1. This is indicated in the config.yaml file. However, if you wish to run the pipeline for a different, or several, distance(s) simply modify the config file accordingly. Include several distances by standard python array notation. Due to the particularities of the software, running the pipeline without cas-offinder is also default behavior as indicated by the variable "cas" in config.yaml. Simply change this to True if you know that cas-offinder is compatible to your platform. 

Please note that particularly the indexing step is a time-costly process and so running the pipeline from end-to-end will take some time. Furthermore, output files for indexing and searching for greater distances grow to a considerable size, so please make sure to run the pipelines in a directory with ample free storage space. 


# Possible issues

* CRISPRITz only recognizes .fa files and not .fasta files as input.

* Cas-offinder requires OpenCL, which is highly system and hardware specific, though it might already be installed and thus Cas-offinder ready to run. To check this the following command can be run to build and check Cas-offinder. If the output recognizes devices under "Available devices list" Cas-offinder should be good to go. 

* ARTEMIS.jl generates a large number of files in its build step, depending on the prefix size for some databases. This can sometimes cause the computer to run out of memory. You could try to make sure many files are allowed for creation, on linux its setting `ulimit -n`.
On Ubuntu to fix it I have set `DefaultLimitNOFILE=524288:524288` to /etc/systemd/system.conf and /etc/systemd/user.conf.