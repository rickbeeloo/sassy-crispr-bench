# Artemis.jl benchmark

This is a benchmark of speed for finding off-target sites for a given gRNA.
Currently this pipeline builds and searches for off-target sites within the whole genome for four software packages: 

* ARTEMIS.jl
* CRISPRITZ
* Cas-offinder

# Results

Results of this pipeline are available for inspection in the `summary.txt` file. There is no need to rerun all the benchmark.

# Installation 

Make sure to install mamba form mambaforge:

```
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```

To execute the pipeline first build an environment (as per snakemake recommendations this step uses mamba):

`mamba env create --name artemis-benchmark --file environment.yaml`

Activate the environment:

`conda activate artemis-benchmark`

To run the pipeline end-to-end use the command:

`snakemake --cores 1`

To run specific rules consult the rule you wish to execute and copy one of it's output files and run:

`snakemake --cores 1 path/to/output.file`

# Requirements

`config.yaml` contains default settings for the pipeline to run over. This pipeline should download all the files and work as is with simple `snakemake --cores 1`. It might take a loong time to perform all calculations, be aware. You need large disc space 
of at least 60GB to perform all calculations.

# Possible issues

* CRISPRITz only recognizes .fa files and not .fasta files as input.

* Cas-offinder requires OpenCL, which is highly system and hardware specific, though it might already be installed and thus Cas-offinder ready to run. To check this the following command can be run to build and check Cas-offinder. If the output recognizes devices under "Available devices list" Cas-offinder should be good to go. 

* ARTEMIS.jl generates a large number of files in its build step, depending on the prefix size for some databases. This can sometimes cause the computer to run out of memory. You could try to make sure many files are allowed for creation, on linux its setting `ulimit -n`.
On Ubuntu to fix it I have set `DefaultLimitNOFILE=524288:524288` to /etc/systemd/system.conf and /etc/systemd/user.conf.