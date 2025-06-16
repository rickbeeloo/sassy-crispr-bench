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
            expand("out_dir/chopoff_out/results/{db_kind}_{prefix}_{dist}_time.csv", 
                db_kind=config["chopoff"], prefix=config["prefix"], dist=config["dist"]),
            expand("out_dir/swoffinder_out/results/swoffinder_{dist}_time.txt",
                dist=config["dist"]),
            expand("out_dir/sassy_out/results/sassy_{dist}_time.txt",
                dist=config["dist"]),
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


## Sassy
rule clone_and_build_sassy:
    output:
        "soft/sassy/sassy"
    shell:
        """
        rm -rf soft/sassy
        git clone https://github.com/RagnarGrootKoerkamp/sassy.git soft/sassy
        cd soft/sassy
        cargo build --release
        """
rule run_sassy:
    input:
        soft="soft/sassy/sassy",
        guides="data/curated_guides_wo_PAM.txt",
        genome="data/hg38v34.fa"
    output:
        res="out_dir/sassy_out/results/sassy_{dist}.txt",
        time="out_dir/sassy_out/results/sassy_{dist}_time.txt"
    threads: config["threads_run"]
    shell:
        "mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time -f 'sassy {wildcards.dist} %e %U %S' "
        "{input.soft} crispr "
        "-g {input.guides} "
        "-k {wildcards.dist} "
        "-t {input.genome} "
        "-o {output.res} "
        "-j {threads} "
        "-n 0.0"
        "; }} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"


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
        db=str("out_dir/chopoff_out/db/{db_kind}_{prefix}_{dist}/{db_kind}.bin")
    shell:
        "export JULIA_NUM_THREADS={config[threads_build]}; "
        "soft/CHOPOFF.jl/build/bin/CHOPOFF build "
        "--name {wildcards.db_kind}_{wildcards.prefix}_{wildcards.dist}_Cas9_hg38v34 "
        "--genome {input.genome} "
        "-o out_dir/chopoff_out/db/{wildcards.db_kind}_{wildcards.prefix}_{wildcards.dist}/ "
        "--distance {wildcards.dist} "
        "--motif Cas9 {wildcards.db_kind} --prefix_length {wildcards.prefix}"


rule chopoff_run_trio:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=str("out_dir/chopoff_out/db/{db_kind}_{prefix}_{dist}/{db_kind}.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/chopoff_out/results/{db_kind}_{prefix}_{dist}.csv",
        time="out_dir/chopoff_out/results/{db_kind}_{prefix}_{dist}_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff {wildcards.db_kind} {wildcards.dist} %e %U %S' {input.soft} "
        "search "
        "--database out_dir/chopoff_out/db/{wildcards.db_kind}_{wildcards.prefix}_{wildcards.dist}/ "
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
        db=str("out_dir/chopoff_out/db/prefixHashDB_{restrict_to_len}_{dist}/prefixHashDB.bin")
    shell:
        "export JULIA_NUM_THREADS={config[threads_build]}; "
        "soft/CHOPOFF.jl/build/bin/CHOPOFF build "
        "--name prefixHashDB_{wildcards.restrict_to_len}_{wildcards.dist}_Cas9_hg38v34 "
        "--genome {input.genome} "
        "-o out_dir/chopoff_out/db/prefixHashDB_{wildcards.restrict_to_len}_{wildcards.dist}/ "
        "--distance {wildcards.dist} "
        "--motif Cas9 prefixHashDB --hash_length {wildcards.restrict_to_len}"


rule chopoff_run_prefixHashDB:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        db=str("out_dir/chopoff_out/db/prefixHashDB_{restrict_to_len}_{dist}/prefixHashDB.bin"),
        guides="data/curated_guides_wo_PAM.txt"
    output:
        res="out_dir/chopoff_out/results/prefixHashDB_{restrict_to_len}_{dist}.csv",
        time="out_dir/chopoff_out/results/prefixHashDB_{restrict_to_len}_{dist}_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_run]}; mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time  -f 'chopoff prefixHashDB_{wildcards.restrict_to_len} {wildcards.dist} %e %U %S' {input.soft} "
        "search "
        "--database out_dir/chopoff_out/db/prefixHashDB_{wildcards.restrict_to_len}_{wildcards.dist}/ "
        "--guides {input.guides} "
        "--output {output.res} "
        "--distance {wildcards.dist} "
        "prefixHashDB; }} 2> {output.time};"
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


# SWOFfinder
rule clone_and_build_swoffinder:
    output:
        "soft/SWOffinder/bin/SmithWatermanOffTarget/SmithWatermanOffTargetSearchAlign.class"
    shell:
        """
        rm -rf soft/SWOffinder
        git clone https://github.com/OrensteinLab/SWOffinder.git soft/SWOffinder
        cd soft/SWOffinder
        javac -d bin SmithWatermanOffTarget/*.java
        """

rule run_swoff:
    input:
        soft="soft/SWOffinder/bin/SmithWatermanOffTarget/SmithWatermanOffTargetSearchAlign.class",
        guides="data/curated_guides_with_PAM.txt",
        genome="data/hg38v34.fa"
    output:
        time="out_dir/swoffinder_out/results/swoffinder_{dist}_time.txt"
    shell:
        "mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time -f 'swoffinder search {wildcards.dist} %e %U %S' "
        "java -cp soft/SWOffinder/bin SmithWatermanOffTarget.SmithWatermanOffTargetSearchAlign "
        "{input.genome} {input.guides} out_dir/swoffinder_out/results/swoffinder_{wildcards.dist} {wildcards.dist} {wildcards.dist} {wildcards.dist} "
        "{wildcards.dist} {config[threads_run]} false 23 NGG false"
        "; }} 2> {output.time};"
        "tail -1 {output.time} >> summary.txt;"

