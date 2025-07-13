import sys

configfile: "config.yaml"

wildcard_constraints:
    dist="0|1|2|3|4|5|6",
    restrict_to_len="14|15|16|17|18|19|20"

rule all:
    input:
        # Sassy timing files for all dist
        expand("out_dir/sassy_out/results/sassy_{dist}_time.txt",
               dist=config["dist"]),
        
        # Chopoff prefixHashDB build time files for all restrict_to_len and dist
        expand("out_dir/chopoff_out/results/prefixHashDB_build_{restrict_to_len}_{dist}_time.csv", 
               restrict_to_len=config["restrict_to_len"], dist=config["dist"]),
        
        # Chopoff database files themselves to ensure building happens
        expand("out_dir/chopoff_out/db/prefixHashDB_{restrict_to_len}_{dist}/prefixHashDB.bin",
               restrict_to_len=config["restrict_to_len"], dist=config["dist"]),
        
        # Chopoff search time files
        expand("out_dir/chopoff_out/results/prefixHashDB_{restrict_to_len}_{dist}_time.csv", 
               restrict_to_len=config["restrict_to_len"], dist=config["dist"]),
        
        # SWOffinder timing files for all dist
        expand("out_dir/swoffinder_out/results/swoffinder_{dist}_time.txt",
               dist=config["dist"]),

# rule all:
#     input:
#         expand("out_dir/chopoff_out/results/prefixHashDB_build_{restrict_to_len}_{dist}_time.csv", 
#                restrict_to_len=config["restrict_to_len"], dist=config["dist"])


rule dag:
    output:
        "dag.pdf"
    shell:
        "snakemake --dry-run --cores 1 --dag | dot -Tpdf > dag.pdf"


rule fa_index:
    input:
        "data/chm13v2.0.fa"
    output:
        "data/chm13v2.0.fa.fai"
    shell:
        "samtools faidx {input}"


rule download_genome:
    output:
        "data/chm13v2.0.fa"
    shell: 
        """
        wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
        gunzip chm13v2.0.fa.gz
        mv chm13v2.0.fa {output}
        """

# Sassy
rule clone_and_build_sassy:
    output:
        "soft/sassy/sassy"
    shell:
        """
        rm -rf soft/sassy
        git clone https://github.com/RagnarGrootKoerkamp/sassy.git soft/sassy
        cd soft/sassy
        cargo build --release
        cp target/release/sassy ./sassy
        """

rule run_sassy:
    input:
        soft="soft/sassy/sassy",
        guides="data/curated_guides_with_PAM.txt",
        genome="data/chm13v2.0.fa"
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

rule chopoff_build_prefixHashDB:
    input:
        soft="soft/CHOPOFF.jl/build/bin/CHOPOFF",
        idx="data/chm13v2.0.fa.fai",
        genome="data/chm13v2.0.fa"
    output:
        db=str("out_dir/chopoff_out/db/prefixHashDB_{restrict_to_len}_{dist}/prefixHashDB.bin"),
        time="out_dir/chopoff_out/results/prefixHashDB_build_{restrict_to_len}_{dist}_time.csv"
    shell:
        "export JULIA_NUM_THREADS={config[threads_build]}; "
        "mkdir -p $(dirname {output.time}); touch {output.time}; "
        "{{ /usr/bin/time -f 'chopoff_build prefixHashDB_{wildcards.restrict_to_len} {wildcards.dist} %e %U %S' "
        "soft/CHOPOFF.jl/build/bin/CHOPOFF build "
        "--name prefixHashDB_{wildcards.restrict_to_len}_{wildcards.dist}_Cas9_chm13v2.0 "
        "--genome {input.genome} "
        "-o out_dir/chopoff_out/db/prefixHashDB_{wildcards.restrict_to_len}_{wildcards.dist}/ "
        "--distance {wildcards.dist} "
        "--motif Cas9 prefixHashDB --hash_length {wildcards.restrict_to_len}"
        "; }} 2> {output.time}; "
        "tail -1 {output.time} >> summary.txt;"


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


# installed through conda environment
rule split_genome:
    input:
        genome="data/chm13v2.0.fa"
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
        genome="data/chm13v2.0.fa"
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

