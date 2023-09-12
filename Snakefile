import sys

configfile: "config.yaml"


rule all:
    input:
        [
            expand("out_dir/artemis_out/results/{db_kind}_{prefix}_{dist}_time.csv", 
                db_kind=config["artemis"], prefix=config["prefix"], dist=config["dist"]),
            expand("out_dir/crispritz_out/results/crispritz_{dist}_time.csv", 
                dist=config["dist"]), 
            expand("out_dir/cas-offinder_out/results/casoffinder_{dist}_time.txt",
                dist=config["dist"]),
            #"out_dir/artemis_out/results/esMax_9_4_time.csv", 
            #"out_dir/artemis_out/results/esMin_9_4_time.csv", 
            #"out_dir/artemis_out/results/es1_9_1_time.csv", 
        ]


rule dag:
    output:
        "dag.pdf"
    shell:
        "snakemake --dry-run --cores 1 --dag | dot -Tpdf > dag.pdf"


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


# this should only run for the largest distance once - other distances can then reuse
rule artemis_build_trio:
    input:
        soft="soft/ARTEMIS.jl/build/bin/ARTEMIS",
        idx="data/hg38v34.fa.fai",
        genome="data/hg38v34.fa"
    output:
        db=str("out_dir/artemis_out/db/{db_kind}_{prefix}_" f"{config['max_dist']}" "/{db_kind}.bin")
    shell:
        "export JULIA_NUM_THREADS={config[threads_build]}; "
        "soft/ARTEMIS.jl/build/bin/ARTEMIS build "
        "--name {wildcards.db_kind}_{wildcards.prefix}_{config[max_dist]}_Cas9_hg38v34 "
        "--genome {input.genome} "
        "-o out_dir/artemis_out/db/{wildcards.db_kind}_{wildcards.prefix}_{config[max_dist]}/ "
        "--distance {config[max_dist]} "
        "--motif Cas9 {wildcards.db_kind} --prefix_length {wildcards.prefix}"


rule artemis_run_trio:
    input:
        soft="soft/ARTEMIS.jl/build/bin/ARTEMIS",
        db=str("out_dir/artemis_out/db/{db_kind}_{prefix}_" f"{config['max_dist']}" "/{db_kind}.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/artemis_out/results/{db_kind}_{prefix}_{dist}.csv",
        time="out_dir/artemis_out/results/{db_kind}_{prefix}_{dist}_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'artemis {wildcards.db_kind} {wildcards.dist} %e %U %S' {input.soft} "
        "search "
        "--database out_dir/artemis_out/db/{wildcards.db_kind}_{wildcards.prefix}_{config[max_dist]}/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance {wildcards.dist} "
        "{wildcards.db_kind}; }} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


rule early_stopping_max:
    input:
        soft="soft/ARTEMIS.jl/build/bin/ARTEMIS",
        db=str("out_dir/artemis_out/db/linearDB_9_4/linearDB.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/artemis_out/results/esMax_9_4.csv",
        time="out_dir/artemis_out/results/esMax_9_4_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'artemis esMax 4 %e %U %S' {input.soft} "
        "search "
        "--database out_dir/artemis_out/db/linearDB_9_4/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance 4 "
        "linearDB "
        "--early_stopping 1000000 1000000 1000000 1000000 1000000; "
        "}} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"
 
 
rule early_stopping_min:
    input:
        soft="soft/ARTEMIS.jl/build/bin/ARTEMIS",
        db=str("out_dir/artemis_out/db/linearDB_9_4/linearDB.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/artemis_out/results/esMin_9_4.csv",
        time="out_dir/artemis_out/results/esMin_9_4_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'artemis esMin 4 %e %U %S' {input.soft} "
        "search "
        "--database out_dir/artemis_out/db/linearDB_9_4/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance 4 "
        "linearDB "
        "--early_stopping 1 1 10 50 200; "
        "}} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


rule early_stopping_1:
    input:
        soft="soft/ARTEMIS.jl/build/bin/ARTEMIS",
        db=str("out_dir/artemis_out/db/linearDB_9_4/linearDB.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/artemis_out/results/es1_9_4.csv",
        time="out_dir/artemis_out/results/es1_9_4_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'artemis es1 1 %e %U %S' {input.soft} "
        "search "
        "--database out_dir/artemis_out/db/linearDB_8_3/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance 1 "
        "linearDB "
        "--early_stopping 1 1; "
        "}} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


## CRISPRITz
# installed through conda environment
rule split_genome:
    input:
        genome="data/hg38v34.fa"
    output:
        "data/chrom_split/stdin.part_chr1.fa"
    shell:
        """
        seqkit split -O data/chrom_split -i < {input.genome}
        rename 's/.fasta$/.fa/' data/chrom_split/*.fasta
        """


# also build once for largest distance
rule crispritz_index:
    input:
        "data/chrom_split/stdin.part_chr1.fa",
        pam="data/20bp-NGG-SpCas9.txt"
    output:
        "genome_library/NGG_{max_dist}_hg38v34_{max_dist}_ref/NGG_chr1 1_1.bin"
    shell:
        "crispritz.py index-genome hg38v34_{config[max_dist]}_ref data/chrom_split/ {input.pam}  -bMax {config[max_dist]} -th {config[threads_run]}"


rule crispritz_search:
    input:
        index=expand("genome_library/NGG_{max_dist}_hg38v34_{max_dist}_ref/NGG_chr1 1_1.bin", max_dist=config["max_dist"]),
        pam="data/20bp-NGG-SpCas9.txt",
        guides="data/curated_guides_wo_PAM.txt"
    output:
        data="out_dir/crispritz_out/results/crispritz_{dist}.targets.txt",
        time="out_dir/crispritz_out/results/crispritz_{dist}_time.csv"
    shell:
        "mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time -f 'CRISPRitz search {wildcards.dist} %e %U %S' crispritz.py "
        "search "
        "genome_library/NGG_{config[max_dist]}_hg38v34_{config[max_dist]}_ref/ {input.pam} {input.guides} "
        "out_dir/crispritz_out/results/crispritz_{wildcards.dist} "
        "-index hg38v34_{config[max_dist]}_ref -mm {wildcards.dist} -bMax {wildcards.dist} "
        "-bDNA {wildcards.dist} -bRNA {wildcards.dist} -th {config[threads_run]} -r; }} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


# Cas-OFFinder
# we don't install through conda as there exists newer version on github 
# its unstable, but native bulges support should give it a chance in speed competition
# install 

rule download_casoff:
    output:
        "soft/build/cas-offinder"
    shell:
        """
        wget https://github.com/snugel/cas-offinder/releases/download/3.0.0b3/cas-offinder_linux_x86_64.zip -O soft/cas-offinder_linux_x86-64.zip
        unzip -q soft/cas-offinder_linux_x86-64 -d ./soft
        rm soft/cas-offinder_linux_x86-64.zip
        sudo chmod 764 {output}
        """


rule input_prep_casoff:
    input:
        "data/curated_guides_wo_PAM.txt"
    output:
        dist="data/cas_guides_input_dist{dist}.txt"
    shell:
        """
        touch {output.dist}
        echo data/hg38v34.fa >> {output.dist}
        echo NNNNNNNNNNNNNNNNNNNNNGG {wildcards.dist} {wildcards.dist} >> {output.dist}
        cat {input} | while read line; do echo ${{line}}'NNN {wildcards.dist}'; done >> {output.dist}
        """


rule run_casoff:
    input:
        soft="soft/build/cas-offinder",
        guides="data/cas_guides_input_dist{dist}.txt"
    output:
        offt="out_dir/cas-offinder_out/results/casoffinder_{dist}.txt",
        time="out_dir/cas-offinder_out/results/casoffinder_{dist}_time.txt"
    shell:
        "mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time -f 'cas-offinder GPU {wildcards.dist} %e %U %S' ./soft/build/cas-offinder "
        "{input.guides} G {output.offt}; }} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"
