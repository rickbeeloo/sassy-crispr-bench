import sys

configfile: "config.yaml"

def get_input(wildcards):
	input_list=[]
	if config["cas"]:
		input_list = (expand("out_dir/cas-offinder_out/out_dist{dist}_time.txt", dist=config["dist"]))
	return input_list

# Separate build and benchmark steps!!!

rule all:
	input:
		get_input,
		"genome_library/NGG_3_hg38v34_ref/NGG_chr10 10_1.bin",
		expand("out_dir/crispritz_out/crispritz_dist{dist}_time.txt", dist=config["dist"]),
		"soft/FlashFry-assembly-1.12.jar",
		"out_dir/flash/hg38_database",
		"out_dir/flash/hg38_database.header",
		"data/curated_guides_w_PAM.fasta",
		expand("out_dir/flash/flash_out.mm_{dist}.time.txt", dist=config["dist"]),
		"soft/build/cas-offinder",
		expand("data/cas_guides_input_dist{dist}.txt", dist=config["dist"]),
		#expand("out_dir/cas-offinder_out/out_dist{dist}_time.txt", dist=config["dist"]),
		#"data/cas_guides_input_dist1.txt",
		#"data/cas_guides_input_dist2.txt",
		#"data/cas_guides_input_dist3.txt",
		#"out_dir/cas-offinder_out/out_dist1.txt",
		#"out_dir/cas-offinder_out/out_dist2.txt",
		#"out_dir/cas-offinder_out/out_dist3.txt",
		#"out_dir/cas-offinder_out/out_dist2_time.txt",
		#"out_dir/cas-offinder_out/out_dist3_time.txt",
		"soft/CRISPRofftargetHunter.jl/build_standalone.sh",
		"soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter",
		expand("out_dir/coh_out/results/{coh}_dist_{dist}_time.txt", coh=config["coh"], dist=config["dist"])

envvars:
	"GIT_TOKEN"

#################################################
#
#	CRISPRITZ
#
#################################################


rule split_gen:
	input:
		genome="data/hg38v34.fa",
		pam="data/20bp-NGG-SpCas9.txt"
	output:
		split=directory("data/chrom_split"),
		index=directory("genome_library"),
		file="genome_library/NGG_4_hg38v34_ref/NGG_chr10 10_1.bin"
	run:
		shell("seqkit split -O {output.split} -i < {input.genome}")
		shell("rename 's/.fasta$/.fa/' {output.split}/*.fasta")
		shell("crispritz.py index-genome hg38v34_ref {output.split} {input.pam}  -bMax 4")

rule create_pam:
	output:
		"data/20bp-NGG-SpCas9.txt"
	shell:
		"echo NNNNNNNNNNNNNNNNNNNNNGG 3 >> {output}"

rule crisp_search_dist1:
        input:
                rules.split_gen.output.file,
                pam="data/20bp-NGG-SpCas9.txt",
                guides="data/curated_guides_wo_PAM.txt"
        output:
                "out_dir/crispritz_out/test_guide_dist1.out.hg38.extended_profile.xls",
                "out_dir/crispritz_out/test_guide_dist1.out.hg38.profile_complete.xls",
                "out_dir/crispritz_out/test_guide_dist1.out.hg38.profile_dna.xls",
                "out_dir/crispritz_out/test_guide_dist1.out.hg38.profile_rna.xls",
                "out_dir/crispritz_out/test_guide_dist1.out.hg38.profile.xls",
                "out_dir/crispritz_out/test_guide_dist1.out.hg38.targets.txt",
                time="out_dir/crispritz_out/crispritz_dist1_time.txt"
        shell:
                "{{ time crispritz.py search genome_library/NGG_3_hg38v34_ref/ {input.pam} {input.guides} "
                "out_dir/crispritz_out/test_guide_dist1.out.hg38 -index -mm 1 -bMax 1 -bDNA 1 -bRNA 1 -t -score hg38v34_ref; }} 2> {output.time}"

rule crisp_search_dist2:
        input:
                pam="data/20bp-NGG-SpCas9.txt",
                guides="data/curated_guides_wo_PAM.txt"
        output:
                "out_dir/crispritz_out/test_guide_dist2.out.hg38.extended_profile.xls",
                "out_dir/crispritz_out/test_guide_dist2.out.hg38.profile_complete.xls",
                "out_dir/crispritz_out/test_guide_dist2.out.hg38.profile_dna.xls",
                "out_dir/crispritz_out/test_guide_dist2.out.hg38.profile_rna.xls",
                "out_dir/crispritz_out/test_guide_dist2.out.hg38.profile.xls",
                "out_dir/crispritz_out/test_guide_dist2.out.hg38.targets.txt",
                time="out_dir/crispritz_out/crispritz_dist2_time.txt"
        shell:
                "{{ time crispritz.py search genome_library/NGG_3_hg38v34_ref/ {input.pam} {input.guides} "
                "out_dir/crispritz_out/test_guide_dist2.out.hg38 -index -mm 2 -bMax 2 -bDNA 2 -bRNA 2 -t -score hg38v34_ref; }} 2> {output.time}"

rule crisp_search_dist3:
        input:
                pam="data/20bp-NGG-SpCas9.txt",
                guides="data/curated_guides_wo_PAM.txt"
        output:
                #expand("out_dir/crispritz_out/test_guide.out.hg38.{type}.{ext}", type=["extended_profile","profile_complete","profile_dna","profile_rna","profile","targets"], ext=["xls","txt"])
                "out_dir/crispritz_out/test_guide_dist3.out.hg38.extended_profile.xls",
                "out_dir/crispritz_out/test_guide_dist3.out.hg38.profile_complete.xls",
                "out_dir/crispritz_out/test_guide_dist3.out.hg38.profile_dna.xls",
                "out_dir/crispritz_out/test_guide_dist3.out.hg38.profile_rna.xls",
                "out_dir/crispritz_out/test_guide_dist3.out.hg38.profile.xls",
                "out_dir/crispritz_out/test_guide_dist3.out.hg38.targets.txt",
                time="out_dir/crispritz_out/crispritz_dist3_time.txt"
        shell:
                "{{ time crispritz.py search genome_library/NGG_3_hg38v34_ref/ {input.pam} {input.guides} "
                "out_dir/crispritz_out/test_guide_dist3.out.hg38 -index -mm 3 -bMax 3 -bDNA 3 -bRNA 3 -t -score hg38v34_ref; }} 2> {output.time}"

rule crisp_search_dist4:
        input:
                pam="data/20bp-NGG-SpCas9.txt",
                guides="data/curated_guides_wo_PAM.txt"
        output:
                #expand("out_dir/crispritz_out/test_guide.out.hg38.{type}.{ext}", type=["extended_profile","profile_complete","profile_dna","profile_rna","profile","targets"], ext=["xls","txt"])
                "out_dir/crispritz_out/test_guide_dist4.out.hg38.extended_profile.xls",
                "out_dir/crispritz_out/test_guide_dist4.out.hg38.profile_complete.xls",
                "out_dir/crispritz_out/test_guide_dist4.out.hg38.profile_dna.xls",
                "out_dir/crispritz_out/test_guide_dist4.out.hg38.profile_rna.xls",
                "out_dir/crispritz_out/test_guide_dist4.out.hg38.profile.xls",
                "out_dir/crispritz_out/test_guide_dist4.out.hg38.targets.txt",
                time="out_dir/crispritz_out/crispritz_dist4_time.txt"
        shell:
                "{{ time crispritz.py search genome_library/NGG_3_hg38v34_ref/ {input.pam} {input.guides} "
                "out_dir/crispritz_out/test_guide.out.hg38 -index -mm 4 -bMax 3 -bDNA 4 -bRNA 4 -t -score hg38v34_ref; }} 2> {output.time}"
	
#################################################
#
#	FlashFry
#
#################################################


rule download_flashfry:
       output:
               "soft/FlashFry-assembly-1.12.jar"
       shell:
               "wget https://github.com/mckennalab/FlashFry/releases/download/1.12/FlashFry-assembly-1.12.jar -O {output}"

rule flash_index:
       input:
	       rules.download_flashfry.output,
               genome="data/hg38v34.fa"
       output:
               db="out_dir/flash/hg38_database",
               dbh="out_dir/flash/hg38_database.header",
	       temp=directory("out_dir/temp")
       run:
               shell("mkdir out_dir/temp")
               shell("java -Xmx6000m -jar soft/FlashFry-assembly-1.12.jar index --tmpLocation ./out_dir/temp --database {output.db} --reference {input.genome} --enzyme spcas9ngg")


rule flash_prep:
       input:
               "data/curated_guides_wo_PAM.txt"
       output:
               "data/curated_guides_w_PAM.fasta"
       run:
               shell("count=0; cat {input} | while read line; do count=$((count+1)); printf '>'$count'\n'${{line}}'GGG\n'; done >> {output};")

rule flash_discover:
       input:
               db="out_dir/flash/hg38_database",
               guide="data/curated_guides_w_PAM.fasta"
       output:
               result="out_dir/flash/flash_out.mm_1.off_targets",
	       time="out_dir/flash/flash_out.mm_1.time.txt"
       shell:
               "{{ time java -Xmx4g -jar soft/FlashFry-assembly-1.12.jar discover --database {input.db} --fasta {input.guide} --output {output.result} --maxMismatch 1 ; }} 2> {output.time}"

rule flash_discover_mm2:
       input:
               db="out_dir/flash/hg38_database",
               guide="data/curated_guides_w_PAM.fasta"
       output:
               result="out_dir/flash/flash_out.mm_2.off_targets",
               time="out_dir/flash/flash_out.mm_2.time.txt"
       shell:
               "{{ time java -Xmx4g -jar soft/FlashFry-assembly-1.12.jar discover --database {input.db} --fasta {input.guide} --output {output.result} --maxMismatch 2 ; }} 2> {output.time}"

rule flash_discover_mm3:
       input:
               db="out_dir/flash/hg38_database",
               guide="data/curated_guides_w_PAM.fasta"
       output:
               result="out_dir/flash/flash_out.mm_3.off_targets",
               time="out_dir/flash/flash_out.mm_3.time.txt"
       shell:
               "{{ time java -Xmx4g -jar soft/FlashFry-assembly-1.12.jar discover --database {input.db} --fasta {input.guide} --output {output.result} --maxMismatch 3 ; }} 2> {output.time}"

rule flash_discover_mm4:
       input:
               db="out_dir/flash/hg38_database",
               guide="data/curated_guides_w_PAM.fasta"
       output:
               result="out_dir/flash/flash_out.mm_4.off_targets",
               time="out_dir/flash/flash_out.mm_4.time.txt"
       shell:
               "{{ time java -Xmx4g -jar soft/FlashFry-assembly-1.12.jar discover --database {input.db} --fasta {input.guide} --output {output.result} --maxMismatch 4 ; }} 2> {output.time}"

rule flash_score:
	input:
		targets="out_dir/flash/flash_out.off_targets",
		db="out_dir/flash/hg38_database"
	output:
		"out_dir/flash/flash_out.off_targets.scored"
	shell:
		"java -Xmx4g -jar soft/FlashFry-assembly-1.12.jar score --input {input.targets} --output {output} --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013 --database {input.db}"

#################################################
#
#       cas-offinder
#
#################################################

rule download_casoff:
	output:
		"soft/build/cas-offinder"
	run:
		shell("wget https://github.com/snugel/cas-offinder/releases/download/3.0.0b3/cas-offinder_linux_x86_64.zip -O soft/cas-offinder_linux_x86-64.zip")
		shell("unzip -q soft/cas-offinder_linux_x86-64 -d ./soft")
		shell("rm soft/cas-offinder_linux_x86-64.zip")
		shell("sudo chmod 764 {output}")

rule input_prep_cas:
	input:
		"data/curated_guides_wo_PAM.txt"
	output:
		dist1="data/cas_guides_input_dist1.txt",
		dist2="data/cas_guides_input_dist2.txt",
		dist3="data/cas_guides_input_dist3.txt",
		dist4="data/cas_guides_input_dist4.txt"
	run:
		shell("touch {output.dist1}; echo data/hg38v34.fa >> {output.dist1}; cut -d ' ' -f 1 data/20bp-NGG-SpCas9.txt >> {output.dist1}; sed -i '$ s/$/ 1 1/' {output.dist1}; cat {input} | while read line; do echo ${{line}}'NNN 1'; done >> {output.dist1};")
		shell("touch {output.dist2}; echo data/hg38v34.fa >> {output.dist2}; cut -d ' ' -f 1 data/20bp-NGG-SpCas9.txt >> {output.dist2}; sed -i '$ s/$/ 2 2/' {output.dist2}; cat {input} | while read line; do echo ${{line}}'NNN 2'; done >> {output.dist2};")
		shell("touch {output.dist3}; echo data/hg38v34.fa >> {output.dist3}; cut -d ' ' -f 1 data/20bp-NGG-SpCas9.txt >> {output.dist3}; sed -i '$ s/$/ 3 3/' {output.dist3}; cat {input} | while read line; do echo ${{line}}'NNN 3'; done >> {output.dist3};")
		shell("touch {output.dist4}; echo data/hg38v34.fa >> {output.dist4}; cut -d ' ' -f 1 data/20bp-NGG-SpCas9.txt >> {output.dist4}; sed -i '$ s/$/ 4 4/' {output.dist4}; cat {input} | while read line; do echo ${{line}}'NNN 4'; done >> {output.dist4};")

rule run_cas_dist1:
	input:
		"data/cas_guides_input_dist1.txt"
	output:
		"out_dir/cas-offinder_out/out_dist1.txt",
		time="out_dir/cas-offinder_out/out_dist1_time.txt"
	run:
		shell("{{ ./soft/build/cas-offinder {input} G {output}; }} 2> {output.time}")

rule run_cas_dist2:
        input:
                "data/cas_guides_input_dist2.txt"
        output:
                "out_dir/cas-offinder_out/out_dist2.txt",
		time="out_dir/cas-offinder_out/out_dist2_time.txt"
        run:
                shell("{{ time ./soft/build/cas-offinder {input} G {output}; }} 2> {output.time}")

rule run_cas_dist3:
        input:
                "data/cas_guides_input_dist1.txt"
        output:
                "out_dir/cas-offinder_out/out_dist3.txt",
		time="out_dir/cas-offinder_out/out_dist3_time.txt"
        run:
                shell("{{ time ./soft/build/cas-offinder {input} G {output}; }} 2> {output.time}")

rule run_cas_dist4:
        input:
                "data/cas_guides_input_dist4.txt"
        output:
                "out_dir/cas-offinder_out/out_dist4.txt",
                time="out_dir/cas-offinder_out/out_dist4_time.txt"
        run:
                shell("{{ time ./soft/build/cas-offinder {input} G {output}; }} 2> {output.time}")

#################################################
#
#       CRISPRofftargetHunter - Initial build
#
#       All modules of COH implementation below,
#	Motif, Tree, and Linear, depend on these 
#	steps, and so I symbolically separate 
#	them from the rest of the COH rules. 
#
#	We also resolve the issues with cloning
#	the COH repo by pre-emptively clearing 
#	the location to which we download the 
#	repo.  
#
#################################################

rule clone_n_build_coh:
	output:
		"soft/CRISPRofftargetHunter.jl/build_standalone.sh",
		bin="soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter"
	params:
		x=os.environ["GIT_TOKEN"]
	run:
		shell("[ -d 'soft/CRISPRofftargetHunter.jl' ] && echo 'Directory soft/CRISPRofftargetHunter.jl exists.'")
		shell("rm -rf soft/CRISPRofftargetHunter.jl")
		shell("git clone https://{params.x}@github.com/JokingHero/CRISPRofftargetHunter.jl soft/CRISPRofftargetHunter.jl")
		shell("cd soft/CRISPRofftargetHunter.jl; ./build_standalone.sh")

rule fa_index:
	input:
		"data/hg38v34.fa"
	output:
		"data/hg38v34.fa.fai"
	shell:
		"samtools faidx {input}"

#################################################
#
#       CRISPRofftargetHunter - MotifDB
#
#       For some reason experiencing issues with
#       this when trying to clone and execute
#       build_standalone.sh. Assumed working 
#       from a funtctioning computer.
#
#################################################

rule build_coh_db:
	input:
		soft="soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter",
		genome="data/hg38v34.fa",
		idx="data/hg38v34.fa.fai"
	output:
		"out_dir/coh_out/db/motifDB_2/AAAAAAA.bin"
	shell:
		"./{input.soft} build --name test_motif --genome {input.genome} --output out_dir/coh_out/db/motifDB_2/ --motif Cas9 motifDB --prefix_length 7"

rule coh_search_dist_1:
	input:
		rules.build_coh_db.output,
		guides="data/curated_guides_wo_PAM.txt"
	output:
		result="out_dir/coh_out/results/mot_dist_1_out.csv",
		detail="out_dir/coh_out/results/mot_dist_1_detail.csv",
		time="out_dir/coh_out/results/mot_dist_1_time.txt"
	shell:
		"{{ time ./soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter "
		"search out_dir/coh_out/db/motifDB_2/ 'motifDB' {input.guides} {output.result} "
		"--distance 1 --detail {output.detail}; }} 2> {output.time}"

rule coh_search_dist_2:
        input:
                "data/curated_guides_wo_PAM.txt"
        output:
                result="out_dir/coh_out/results/mot_dist_2_out.csv",
                detail="out_dir/coh_out/results/mot_dist_2_detail.csv",
                time="out_dir/coh_out/results/mot_dist_2_time.txt"
        shell:
                "{{ time ./soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter "
		"search out_dir/coh_out/db/motifDB_2/ 'motifDB' {input} {output.result} "
		"--distance 2 --detail {output.detail}; }} 2> {output.time}"

rule coh_search_dist_3:
        input:
                "data/curated_guides_wo_PAM.txt"
        output:
                result="out_dir/coh_out/results/mot_dist_3_out.csv",
                detail="out_dir/coh_out/results/mot_dist_3_detail.csv",
                time="out_dir/coh_out/results/mot_dist_3_time.txt"
        shell:
                "{{ time ./soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter "
                "search out_dir/coh_out/db/motifDB_2/ 'motifDB' {input} {output.result} "
                "--distance 3 --detail {output.detail}; }} 2> {output.time}"

rule coh_search_dist_4:
    input:
        "data/curated_guides_wo_PAM.txt"
    output:
        result="out_dir/coh_out/results/mot_dist_4_out.csv",
        detail="out_dir/coh_out/results/mot_dist_4_detail.csv",
        time="out_dir/coh_out/results/mot_dist_4_time.txt"
    shell:
        "{{ time ./soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter "
    	"search out_dir/coh_out/db/motifDB_2/ 'motifDB' {input} {output.result} "
    	"--distance 4 --detail {output.detail}; }} 2> {output.time}"

#################################################
#
#       CRISPRofftargetHunter - TreeDB
#
#       TODO
#
#################################################

rule build_tree_db:
	input:
		soft="soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter",
		idx="data/hg38v34.fa.fai",
		genome="data/hg38v34.fa"
	output:
		"out_dir/coh_out/db/treeDB/AAAAAAA.bin"
	shell:
		"soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter build --name tree --genome {input.genome} --output out_dir/coh_out/db/treeDB/ --motif Cas9 treeDB --prefix_length 7"

rule coh_search_tree_dist_1:
	input:
		rules.build_tree_db.output,
		guides="data/curated_guides_wo_PAM.txt"
	output:
		result="out_dir/coh_out/results/tree_dist_1_out.csv",
		detail="out_dir/coh_out/results/tree_dist_1_detail.csv",
		time="out_dir/coh_out/results/tree_dist_1_time.txt"
	shell:
		"{{ time ./soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter "
		"search out_dir/coh_out/db/treeDB/ 'treeDB' {input.guides} {output.result} "
		"--distance 1 --detail {output.detail}; }} 2> {output.time}"

rule coh_search_tree_dist_2:
	input:
		rules.build_tree_db.output,
		guides="data/curated_guides_wo_PAM.txt"
	output:
		result="out_dir/coh_out/results/tree_dist_2_out.csv",
		detail="out_dir/coh_out/results/tree_dist_2_detail.csv",
		time="out_dir/coh_out/results/tree_dist_2_time.txt"
	shell:
		"{{ time ./soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter "
		"search out_dir/coh_out/db/treeDB/ 'treeDB' {input.guides} {output.result} "
		"--distance 2 --detail {output.detail}; }} 2> {output.time}"

rule coh_search_tree_dist_3:
	input:
		rules.build_tree_db.output,
		guides="data/curated_guides_wo_PAM.txt"
	output:
		result="out_dir/coh_out/results/tree_dist_3_out.csv",
		detail="out_dir/coh_out/results/tree_dist_3_detail.csv",
		time="out_dir/coh_out/results/tree_dist_3_time.txt"
	shell:
		"{{ time ./soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter "
		"search out_dir/coh_out/db/treeDB/ 'treeDB' {input.guides} {output.result} "
		"--distance 3 --detail {output.detail}; }} 2> {output.time}"

rule coh_search_tree_dist_4:
	input:
		rules.build_tree_db.output,
		guides="data/curated_guides_wo_PAM.txt"
	output:
		result="out_dir/coh_out/results/tree_dist_4_out.csv",
		detail="out_dir/coh_out/results/tree_dist_4_detail.csv",
		time="out_dir/coh_out/results/tree_dist_4_time.txt"
	shell:
		"{{ time ./soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter "
		"search out_dir/coh_out/db/treeDB/ 'treeDB' {input.guides} {output.result} "
		"--distance 4 --detail {output.detail}; }} 2> {output.time}"

#################################################
#
#       CRISPRofftargetHunter - LinearDB
#
#       TODO
#
#################################################

rule build_linear_db:
	input:
		soft="soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter",
		idx="data/hg38v34.fa.fai",
		genome="data/hg38v34.fa"
	output:
		"out_dir/coh_out/db/linearDB/AAAAAAA.bin"
	shell:
		"soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter build --name tree --genome {input.genome} --output out_dir/coh_out/db/linearDB/ --motif Cas9 linearDB --prefix_length 7"

rule coh_search_lin_dist_1:
	input:
		rules.build_linear_db.output,
		guides="data/curated_guides_wo_PAM.txt"
	output:
		result="out_dir/coh_out/results/linear_dist_1_out.csv",
		detail="out_dir/coh_out/results/linear_dist_1_detail.csv",
		time="out_dir/coh_out/results/linear_dist_1_time.txt"
	shell:
		"{{ time ./soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter "
		"search out_dir/coh_out/db/linearDB/ 'linearDB' {input.guides} {output.result} "
		"--distance 1 --detail {output.detail}; }} 2> {output.time}"

rule coh_search_lin_dist_2:
	input:
		rules.build_linear_db.output,
		guides="data/curated_guides_wo_PAM.txt"
	output:
		result="out_dir/coh_out/results/linear_dist_2_out.csv",
		detail="out_dir/coh_out/results/linear_dist_2_detail.csv",
		time="out_dir/coh_out/results/linear_dist_2_time.txt"
	shell:
		"{{ time ./soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter "
		"search out_dir/coh_out/db/linearDB/ 'linearDB' {input.guides} {output.result} "
		"--distance 2 --detail {output.detail}; }} 2> {output.time}"

rule coh_search_lin_dist_3:
	input:
		rules.build_linear_db.output,
		guides="data/curated_guides_wo_PAM.txt"
	output:
		result="out_dir/coh_out/results/linear_dist_3_out.csv",
		detail="out_dir/coh_out/results/linear_dist_3_detail.csv",
		time="out_dir/coh_out/results/linear_dist_3_time.txt"
	shell:
		"{{ time ./soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter "
		"search out_dir/coh_out/db/linearDB/ 'linearDB' {input.guides} {output.result} "
		"--distance 3 --detail {output.detail}; }} 2> {output.time}"

rule coh_search_lin_dist_4:
	input:
		rules.build_linear_db.output,
		guides="data/curated_guides_wo_PAM.txt"
	output:
		result="out_dir/coh_out/results/linear_dist_4_out.csv",
		detail="out_dir/coh_out/results/linear_dist_4_detail.csv",
		time="out_dir/coh_out/results/linear_dist_4_time.txt"
	shell:
		"{{ time ./soft/CRISPRofftargetHunter.jl/build/bin/CRISPRofftargetHunter "
		"search out_dir/coh_out/db/linearDB/ 'linearDB' {input.guides} {output.result} "
		"--distance 4 --detail {output.detail}; }} 2> {output.time}"
