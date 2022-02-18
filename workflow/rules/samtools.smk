rule samtools_faidx:
	input:
		genome=config["genome"]["dir"]+config["genome"]["name"]+config["genome"]["ext"]
	output:
		config["genome"]["dir"]+config["genome"]["name"]+config["genome"]["ext"]+".fai"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools faidx -o {output} {input}"

rule sam_to_bam:
	input:
		sam=config["general"]["experiment_name"]+"/mapping/sam/{prefix}.sam",
		genome=config["genome"]["dir"]+config["genome"]["name"]+config["genome"]["ext"]+".fai"
	output:
		config["general"]["experiment_name"]+"/mapping/bam/raw/{prefix}.bam"
	params:
		quality=config["samtools"]["quality"],
		custom=config["samtools"]["custom"]
	benchmark :
		config["general"]["experiment_name"]+"/benchmarks/sam_to_bam/{prefix}.txt"
	priority: 50
	threads: config["general"]["threads"]
	conda: "../envs/samtools.yaml"
	message: "##RUNNING : samtools view for {input.sam}"
	shell:
		"samtools view "
		"{params.custom} -@ {threads} "
		"-b -S "
		"-q {params.quality} "
		"-t {input.genome} "
		"-o {output} "
		"{input.sam}"

rule samtools_sort:
	input:
		bam=config["general"]["experiment_name"]+"/mapping/bam/raw/{prefix}.bam"
	output:
		config["general"]["experiment_name"]+"/mapping/bam/sorted/{prefix}.sorted.bam"
	benchmark :
		config["general"]["experiment_name"]+"/benchmarks/samtools_sort/{prefix}.txt"
	priority: 50
	threads: config["general"]["threads"]
	message: "##RUNNING : samtools sort {input.bam}"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools sort -@ {threads} "
		"-o {output} "
		"{input.bam}"


rule samtools_sortn:
	input:
		bam=config["general"]["experiment_name"]+"/mapping/bam/sorted/{prefix}.sorted.bam"
	output:
		temp(config["general"]["experiment_name"]+"/mapping/bam/sorted/{prefix}.nsorted.bam")
	benchmark :
		config["general"]["experiment_name"]+"/benchmarks/samtools_sortn/{prefix}.txt"
	priority: 50
	message: "##RUNNING : samtools sort -n {input.bam}"
	conda: "../envs/samtools.yaml"
	threads: config["general"]["threads"]
	shell:
		"samtools sort -@ {threads} -n -o {output}"
		"{input.bam} "

rule samtools_fixmate:
	input:
		bam=config["general"]["experiment_name"]+"/mapping/bam/sorted/{prefix}.sorted.bam"
	output:
		temp(config["general"]["experiment_name"]+"/mapping/bam/sorted/{prefix}.fixmate.bam")
	benchmark :
		config["general"]["experiment_name"]+"/benchmarks/samtools_fixmate/{prefix}.txt"
	priority: 50
	message: "##RUNNING : samtools fixmate -m {input.bam}"
	threads: config["general"]["threads"]
	conda: "../envs/samtools.yaml"
	shell:
		"samtools fixmate -m -@ {threads}"
		"{input.bam} "
		"{output}"

rule samtools_markdups:
	input:
		bam=config["general"]["experiment_name"]+"/mapping/bam/sorted/{prefix}.fixmate.bam"
	output:
		config["general"]["experiment_name"]+"/mapping/bam/rmdups/{prefix}.rmdups.bam"
	benchmark :
		config["general"]["experiment_name"]+"/benchmarks/samtools_markdups/{prefix}.txt"
	priority: 50
	message: "##RUNNING : samtools markdup {input.bam}"
	threads: config["general"]["threads"]
	conda: "../envs/samtools.yaml"
	shell:
		"samtools markdup -@ {threads} -r "
		"{input.bam} "
		"{output}"

rule samtools_index:
	input:
		bam=config["general"]["experiment_name"]+"/mapping/bam/{samtype}/{prefix}.{samtype}.bam"
	output:
		config["general"]["experiment_name"]+"/mapping/bam/{samtype}/{prefix}.{samtype}.bai"
	benchmark :
		config["general"]["experiment_name"]+"/benchmarks/samtools_index/{prefix}.{samtype}.txt"
	priority: 50
	threads: config["general"]["threads"]
	conda: "../envs/samtools.yaml"
	message: "##RUNNING : samtools index {input}"
	shell:
		"samtools index -@ {threads} {input} {output}"


rule samtools_stats:
	input:
		bam=config["general"]["experiment_name"]+"/mapping/bam/{samtype}/{prefix}.{samtype}.bam",
		bai=config["general"]["experiment_name"]+"/mapping/bam/{samtype}/{prefix}.{samtype}.bai"
	output:
		config["general"]["experiment_name"]+"/QC/STATS/stats/{prefix}.{samtype}.stats"
	benchmark :
		config["general"]["experiment_name"]+"/benchmarks/samtools_stats/{prefix}.{samtype}.txt"
	priority: 50
	message: "##RUNNING : samtools stats {input}"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools stats {input.bam} > {output}"

rule samtools_idxstats:
	input:
		bam=config["general"]["experiment_name"]+"/mapping/bam/{samtype}/{prefix}.{samtype}.bam",
		bai=config["general"]["experiment_name"]+"/mapping/bam/{samtype}/{prefix}.{samtype}.bai"
	output:
		config["general"]["experiment_name"]+"/QC/STATS/idxstats/{prefix}.{samtype}.idxstats"
	benchmark :
		config["general"]["experiment_name"]+"/benchmarks/samtools_idxstats/{prefix}.{samtype}.txt"
	priority: 50
	message: "##RUNNING : samtools idxstats {input}"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools idxstats {input.bam} > {output}"

rule samtools_flagstat:
	input:
		bam=config["general"]["experiment_name"]+"/mapping/bam/{samtype}/{prefix}.{samtype}.bam",
		bai=config["general"]["experiment_name"]+"/mapping/bam/{samtype}/{prefix}.{samtype}.bai"
	output:
		config["general"]["experiment_name"]+"/QC/STATS/flagstat/{prefix}.{samtype}.flagstat"
	benchmark :
		config["general"]["experiment_name"]+"/benchmarks/samtools_flagstat/{prefix}.{samtype}.txt"
	priority: 50
	message: "##RUNNING : samtools flagstat {input}"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools flagstat {input.bam} > {output}"

