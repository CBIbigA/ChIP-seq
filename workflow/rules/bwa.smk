rule bwa_index:
	input:
		config["genome"]["dir"]+config["genome"]["name"]+config["genome"]["ext"],
	output:
		idx=multiext(config["genome"]["dir"]+config["genome"]["name"]+config["genome"]["ext"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
	log:
		"logs/bwa_index/genome.log",
	params:
		algorithm="bwtsw"
	conda: "../envs/bwa_mem.yaml"
	shell:
		"bwa index -a {params.construction_algorithm} {input} 2> {log}"

rule bwa_mem:
	input:
		fastq=config["general"]["sample_dir"]+"/{prefix}"+config["raw_spec"]["ext"],
		genome=config["genome"]["dir"]+config["genome"]["name"]+config["genome"]["ext"],
		rules.bwa_index.output
	output:
		temp(config["general"]["experiment_name"]+"/mapping/sam/{prefix}.sam")
	params:
		custom=config["bwa"]["custom"]
	threads: config["general"]["threads"]
	benchmark :
		config["general"]["experiment_name"]+"/benchmarks/bwa_mem/{prefix}.txt"
	priority: 50
	message: "##RUNNING : bwa mem for {input.fastq}"
	conda: "../envs/bwa_mem.yaml"
	log:
		"logs/bwa_mem/{prefix}.log"
	shell:
		"bwa mem {params.custom} "
		"-t {threads} "
		"{input.genome} {input.fastq} "
		"> {output} 2> {log}"
