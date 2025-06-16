# Sassy crispr benchmark

This file was modified from the original [Chopoff.jl benchmark](https://git.app.uib.no/valenlab/chopoff-benchmark) to include Sassy, credits to making their benchmarks available. Modifications:
* Added `Sassy crispr`
* Added timings for each Chopoff edit distance build (instead of re-using max edit database)
* Added Rust to environment.yaml
* Removed Cas-offinder and CRISPRITz

# Running

To execute the pipeline first build an environment (as per snakemake recommendations this step uses mamba):

`mamba env create --name sassy-benchmark --file environment.yaml` 

Activate the environment:

`conda activate sassy-benchmark`

To run the pipeline end-to-end use the command:

`snakemake --cores 1`

To run specific rules consult the rule you wish to execute and copy one of it's output files and run:

`snakemake --cores 1 path/to/output.file`


# Results
See `summary.txt` file

