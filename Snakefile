import sys

configfile: "config.yaml"

rule all:
    input:
	    [
			expand("out_dir/artemis_out/results/{db_kind}_{prefix}_{dist}_time.csv", db_kind=config["artemis"], prefix=config["prefix"], dist=config["dist"])
		]
		
		

rule fa_index:
    input:
        "data/hg38v34.fa"
    output:
        "data/hg38v34.fa.fai"
    shell:
        "samtools faidx {input}"


rule download_genome:
    output:
        "data/hg38v34.fa"
    shell: 
        """
        wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz
        gunzip GRCh38.primary_assembly.genome.fa.gz
        mv GRCh38.primary_assembly.genome.fa {output}
        """


## ARTEMIS.jl

rule clone_and_build_artemis:
	output:
		"soft/ARTEMIS.jl/build/bin/ARTEMIS"
	shell:
		"""
        rm -rf soft/ARTEMIS.jl
		git clone https://github.com/JokingHero/ARTEMIS.jl soft/ARTEMIS.jl
		cd soft/ARTEMIS.jl
        ./build_standalone.sh
        """


rule artemis_build_trio:
	input:
		soft="soft/ARTEMIS.jl/build/bin/ARTEMIS",
		idx="data/hg38v34.fa.fai",
		genome="data/hg38v34.fa"
	output:
		db="out_dir/artemis_out/db/{db_kind}_{prefix}_{dist}/{db_kind}.bin"
	shell:
		"soft/ARTEMIS.jl/build/bin/ARTEMIS build "
        "--name {wildcards.db_kind}_{wildcards.prefix}_{wildcards.dist}_Cas9_hg38v34 "
        "--genome {input.genome} "
        "-o out_dir/artemis_out/db/{wildcards.db_kind}_{wildcards.prefix}_{wildcards.dist}/ "
        "--distance {wildcards.dist} "
        "--motif Cas9 {wildcards.db_kind} --prefix_length {wildcards.prefix}"


rule artemis_run_trio:
	input:
	    soft="soft/ARTEMIS.jl/build/bin/ARTEMIS",
		db=rules.artemis_build_trio.output.db,
		guides="data/curated_guides_wo_PAM.txt"
	output:
		res="out_dir/artemis_out/results/{db_kind}_{prefix}_{dist}.csv",
		time="out_dir/artemis_out/results/{db_kind}_{prefix}_{dist}_time.csv"
	shell:
		"{{ time {input.soft} "
		"search out_dir/artemis_out/db/{wildcards.db_kind}_{wildcards.prefix}_{wildcards.dist}/ '{wildcards.db_kind}' {input.guides} {output.res} "
		"--distance {wildcards.dist}; }} 2> {output.time}"