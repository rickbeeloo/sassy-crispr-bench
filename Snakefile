import sys

configfile: "config.yaml"

wildcard_constraints:
    db_kind="linearDB|treeDB|motifDB",
    prefix="5|6|7|8|9|10",
    dist="0|1|2|3|4|5|6",
    restrict_to_len="14|15|16|17|18|19|20"

rule all:
    input:
        [  
            expand("out_dir/chopoff_out/results/prefixHashDB_{restrict_to_len}_{dist}_time.csv", 
                restrict_to_len=config["restrict_to_len"], dist=config["dist"]),
            #"out_dir/chopoff_out/results/linearHashDB_with_es_max_8_16_3_time.csv",
            #"out_dir/chopoff_out/results/linearHashDB_with_es_min_8_16_3_time.csv"
            #expand("out_dir/chopoff_out/results/bffDB_{restrict_to_len}_{dist}_time.csv", restrict_to_len=config["restrict_to_len"], dist=config["dist"])
            #expand("out_dir/chopoff_out/results/{db_kind}_{prefix}_{dist}_time.csv", 
            #    db_kind=config["chopoff"], prefix=config["prefix"], dist=config["dist"]),
            #expand("out_dir/crispritz_out/results/crispritz_{dist}_time.csv", 
            #    dist=config["dist"]), 
            #expand("out_dir/cas-offinder_out/results/casoffinder_{dist}_time.txt",
            #    dist=config["dist"]),
            #"out_dir/chopoff_out/results/esMax_9_4_time.csv", 
            #"out_dir/chopoff_out/results/esMin_9_4_time.csv", 
            #"out_dir/chopoff_out/results/esOne_9_4_time.csv",
            #"out_dir/chopoff_out/results/dictDB_time.csv",
            #"out_dir/chopoff_out/results/hashDB_right_time.csv",
            #"out_dir/chopoff_out/results/hashDB_left_time.csv"
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


## CHOPOFF.jl

rule clone_and_build_chopoff:
    output:
        "soft/CHOPOFF.jl/build/bin/CHOPOFF"
    shell:
        """
        rm -rf soft/CHOPOFF.jl
        git clone https://github.com/JokingHero/CHOPOFF.jl soft/CHOPOFF.jl
        cd soft/CHOPOFF.jl
        ./build_standalone.sh
        """


# this should only run for the largest distance once - other distances can then reuse
rule chopoff_build_trio:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        idx="data/hg38v34.fa.fai",
        genome="data/hg38v34.fa"
    output:
        db=str("out_dir/chopoff_out/db/{db_kind}_{prefix}_" f"{config['max_dist']}" "/{db_kind}.bin")
    shell:
        "export JULIA_NUM_THREADS={config[threads_build]}; "
        "soft/CHOPOFF.jl/build/bin/CHOPOFF build "
        "--name {wildcards.db_kind}_{wildcards.prefix}_{config[max_dist]}_Cas9_hg38v34 "
        "--genome {input.genome} "
        "-o out_dir/chopoff_out/db/{wildcards.db_kind}_{wildcards.prefix}_{config[max_dist]}/ "
        "--distance {config[max_dist]} "
        "--motif Cas9 {wildcards.db_kind} --prefix_length {wildcards.prefix}"


rule chopoff_run_trio:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=str("out_dir/chopoff_out/db/{db_kind}_{prefix}_" f"{config['max_dist']}" "/{db_kind}.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/chopoff_out/results/{db_kind}_{prefix}_{dist}.csv",
        time="out_dir/chopoff_out/results/{db_kind}_{prefix}_{dist}_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff {wildcards.db_kind} {wildcards.dist} %e %U %S' {input.soft} "
        "search "
        "--database out_dir/chopoff_out/db/{wildcards.db_kind}_{wildcards.prefix}_{config[max_dist]}/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance {wildcards.dist} "
        "{wildcards.db_kind}; }} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


rule chopoff_build_prefixHashDB:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        idx="data/hg38v34.fa.fai",
        genome="data/hg38v34.fa"
    output:
        db=str("out_dir/chopoff_out/db/prefixHashDB_{restrict_to_len}_" f"{config['max_dist']}" "/prefixHashDB.bin")
    shell:
        "export JULIA_NUM_THREADS={config[threads_build]}; "
        "soft/CHOPOFF.jl/build/bin/CHOPOFF build "
        "--name prefixHashDB_{wildcards.restrict_to_len}_{config[max_dist]}_Cas9_hg38v34 "
        "--genome {input.genome} "
        "-o out_dir/chopoff_out/db/prefixHashDB_{wildcards.restrict_to_len}_{config[max_dist]}/ "
        "--distance {config[max_dist]} "
        "--motif Cas9 prefixHashDB --hash_length {wildcards.restrict_to_len}"


rule chopoff_run_prefixHashDB:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=str("out_dir/chopoff_out/db/prefixHashDB_{restrict_to_len}_" f"{config['max_dist']}" "/prefixHashDB.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/chopoff_out/results/prefixHashDB_{restrict_to_len}_{dist}.csv",
        time="out_dir/chopoff_out/results/prefixHashDB_{restrict_to_len}_{dist}_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff prefixHashDB_{wildcards.restrict_to_len} {wildcards.dist} %e %U %S' {input.soft} "
        "search "
        "--database out_dir/chopoff_out/db/prefixHashDB_{wildcards.restrict_to_len}_{config[max_dist]}/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance {wildcards.dist} "
        "prefixHashDB; }} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


rule chopoff_build_linearHashDB:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        idx="data/hg38v34.fa.fai",
        genome="data/hg38v34.fa"
    output:
        db=str("out_dir/chopoff_out/db/linearHashDB_{prefix}_{restrict_to_len}_" f"{config['max_dist']}" "/linearHashDB.bin")
    shell:
        "export JULIA_NUM_THREADS={config[threads_build]}; "
        "soft/CHOPOFF.jl/build/bin/CHOPOFF build "
        "--name linearHashDB_{wildcards.prefix}_{wildcards.restrict_to_len}_{config[max_dist]}_Cas9_hg38v34 "
        "--genome {input.genome} "
        "-o out_dir/chopoff_out/db/linearHashDB_{wildcards.prefix}_{wildcards.restrict_to_len}_{config[max_dist]}/ "
        "--distance {config[max_dist]} "
        "--motif Cas9 linearHashDB --prefix_length {wildcards.prefix} --hash_length {wildcards.restrict_to_len}"


rule chopoff_run_linearHashDB:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=str("out_dir/chopoff_out/db/linearHashDB_{prefix}_{restrict_to_len}_" f"{config['max_dist']}" "/linearHashDB.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/chopoff_out/results/linearHashDB_{prefix}_{restrict_to_len}_{dist}.csv",
        time="out_dir/chopoff_out/results/linearHashDB_{prefix}_{restrict_to_len}_{dist}_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff linearHashDB_{wildcards.restrict_to_len} {wildcards.dist} %e %U %S' {input.soft} "
        "search "
        "--database out_dir/chopoff_out/db/linearHashDB_{wildcards.prefix}_{wildcards.restrict_to_len}_{config[max_dist]}/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance {wildcards.dist} "
        "linearHashDB; }} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


rule chopoff_run_linearHashDB_with_es_max:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=str("out_dir/chopoff_out/db/linearHashDB_8_16_4/linearHashDB.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/chopoff_out/results/linearHashDB_with_es_max_8_16_3.csv",
        time="out_dir/chopoff_out/results/linearHashDB_with_es_max_8_16_3_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff linearHashDB_with_es_max_8_16 3 %e %U %S' {input.soft} "
        "search "
        "--database out_dir/chopoff_out/db/linearHashDB_8_16_4/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance 3 "
        "linearHashDB "
        "--early_stopping 10000 10000 10000 10000; "
        " }} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


rule chopoff_run_linearHashDB_with_es_min:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=str("out_dir/chopoff_out/db/linearHashDB_8_16_4/linearHashDB.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/chopoff_out/results/linearHashDB_with_es_min_8_16_3.csv",
        time="out_dir/chopoff_out/results/linearHashDB_with_es_min_8_16_3_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff linearHashDB_with_es_min_8_16 3 %e %U %S' {input.soft} "
        "search "
        "--database out_dir/chopoff_out/db/linearHashDB_8_16_4/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance 3 "
        "linearHashDB "
        "--early_stopping 1 1 1 1; "
        " }} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


rule early_stopping_max:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=str("out_dir/chopoff_out/db/linearDB_9_4/linearDB.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/chopoff_out/results/esMax_9_4.csv",
        time="out_dir/chopoff_out/results/esMax_9_4_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff esMax 4 %e %U %S' {input.soft} "
        "search "
        "--database out_dir/chopoff_out/db/linearDB_9_4/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance 4 "
        "linearDB "
        "--early_stopping {wildcards.es}; "
        "}} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"
 
 
rule early_stopping_min:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=str("out_dir/chopoff_out/db/linearDB_9_4/linearDB.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/chopoff_out/results/esMin_9_4.csv",
        time="out_dir/chopoff_out/results/esMin_9_4_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff esMin 4 %e %U %S' {input.soft} "
        "search "
        "--database out_dir/chopoff_out/db/linearDB_9_4/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance 4 "
        "linearDB "
        "--early_stopping {wildcards.es}; "
        "}} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


rule early_stopping_one:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=str("out_dir/chopoff_out/db/linearDB_9_4/linearDB.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/chopoff_out/results/esOne_9_4.csv",
        time="out_dir/chopoff_out/results/esOne_9_4_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff esOne 1 %e %U %S' {input.soft} "
        "search "
        "--database out_dir/chopoff_out/db/linearDB_9_4/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance 1 "
        "linearDB "
        "--early_stopping 1 1; "
        "}} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


rule chopoff_build_dictDB:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        idx="data/hg38v34.fa.fai",
        genome="data/hg38v34.fa"
    output:
        db="out_dir/chopoff_out/db/dictDB/dictDB.bin"
    shell:
        "export JULIA_NUM_THREADS={config[threads_build]}; mkdir -p $(dirname {output.db}); "
        "soft/CHOPOFF.jl/build/bin/CHOPOFF build "
        "--name dictDB_default_Cas9_hg38v34 "
        "--genome {input.genome} "
        "--output {output.db} "
        "--distance 1 "
        "--motif Cas9 dictDB"


rule chopoff_build_hashDB:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        idx="data/hg38v34.fa.fai",
        genome="data/hg38v34.fa"
    output:
        db="out_dir/chopoff_out/db/hashDB/hashDB.bin"
    shell:
        "export JULIA_NUM_THREADS={config[threads_build]}; mkdir -p $(dirname {output.db}); "
        "soft/CHOPOFF.jl/build/bin/CHOPOFF build "
        "--name hashDB_default_Cas9_hg38v34 "
        "--genome {input.genome} "
        "--output {output.db} "
        "--distance 1 "
        "--motif Cas9 hashDB"


rule chopoff_run_dictDB:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=str("out_dir/chopoff_out/db/dictDB/dictDB.bin"),
        guides="data/all_genes_only_guides.txt"
    output:
        res="out_dir/chopoff_out/results/dictDB.csv",
        time="out_dir/chopoff_out/results/dictDB_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff dictDB 1 %e %U %S' {input.soft} "
        "estimate "
        "--database {input.db} "
        "--guides {input.guides} "
        "--output {output.res}; "
        "}} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


rule chopoff_run_hashDB_right:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=str("out_dir/chopoff_out/db/hashDB/hashDB.bin"),
        guides="data/all_genes_only_guides.txt"
    output:
        res="out_dir/chopoff_out/results/hashDB_right.csv",
        time="out_dir/chopoff_out/results/hashDB_right_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff hashDB_right 1 %e %U %S' {input.soft} "
        "estimate "
        "--database {input.db} "
        "--guides {input.guides} "
        "--output {output.res} "
        "--right; "
        "}} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


rule chopoff_run_hashDB_left:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=str("out_dir/chopoff_out/db/hashDB/hashDB.bin"),
        guides="data/all_genes_only_guides.txt"
    output:
        res="out_dir/chopoff_out/results/hashDB_left.csv",
        time="out_dir/chopoff_out/results/hashDB_left_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff hashDB_left 1 %e %U %S' {input.soft} "
        "estimate "
        "--database {input.db} "
        "--guides {input.guides} "
        "--output {output.res}; "
        "}} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


#BINARY FUSE FILTER + FM-index
rule chopoff_build_bff:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        idx="data/hg38v34.fa.fai",
        genome="data/hg38v34.fa"
    output:
        db=str("out_dir/chopoff_out/db/bffDB_{restrict_to_len}_{dist}/BinaryFuseFilterDB.bin")
    shell:
        "export JULIA_NUM_THREADS={config[threads_build]}; "
        "soft/CHOPOFF.jl/build/bin/CHOPOFF build "
        "--name bffDB_{wildcards.restrict_to_len}_{wildcards.dist}_Cas9_hg38v34 "
        "--genome {input.genome} "
        "-o out_dir/chopoff_out/db/bffDB_{wildcards.restrict_to_len}_{wildcards.dist}/ "
        "--distance {wildcards.dist} "
        "--motif Cas9 bffDB --restrict_to_len {wildcards.restrict_to_len}"


rule chopoff_build_fmi:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        idx="data/hg38v34.fa.fai",
        genome="data/hg38v34.fa"
    output:
        db=str("out_dir/chopoff_out/db/fmi/genomeInfo.bin")
    shell:
        "export JULIA_NUM_THREADS={config[threads_build]}; "
        "soft/CHOPOFF.jl/build/bin/CHOPOFF build "
        "--name fmi_Cas9_hg38v34 "
        "--genome {input.genome} "
        "-o out_dir/chopoff_out/db/fmi/ "
        "--distance 3 "
        "--motif Cas9 fmi"


rule chopoff_run_bff:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=[str("out_dir/chopoff_out/db/bffDB_{restrict_to_len}_{dist}/BinaryFuseFilterDB.bin"), "out_dir/chopoff_out/db/fmi/genomeInfo.bin"],
        genome="data/hg38v34.fa",
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/chopoff_out/results/bffDB_{restrict_to_len}_{dist}.csv",
        time="out_dir/chopoff_out/results/bffDB_{restrict_to_len}_{dist}_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff bffDB_{wildcards.restrict_to_len} {wildcards.dist} %e %U %S' {input.soft} "
        "search "
        "--database out_dir/chopoff_out/db/bffDB_{wildcards.restrict_to_len}_{wildcards.dist}/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance {wildcards.dist} bffDB "
        "--fmiDB out_dir/chopoff_out/db/fmi/ --genome {input.genome}; }} "
        "2> {output.time};"
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
