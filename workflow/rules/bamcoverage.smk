rule bamCoverage:
	input:
		bam=config["general"]["experiment_name"]+"/mapping/bam/{samtype}/{prefix}.{samtype}.bam",
		bai=config["general"]["experiment_name"]+"/mapping/bam/{samtype}/{prefix}.{samtype}.bai"
	threads : config["general"]["threads"]
	output:
		protected(config["general"]["experiment_name"]+"/mapping/BIGWIG/{prefix}_{samtype}.{normalization}.bw")
	params:
		normalization="{normalization}",
		binsize=config["bamcoverage"]["binsize"]
	benchmark :
		config["general"]["experiment_name"]+"/benchmarks/bamCoverage/{prefix}.{normalization}.{samtype}.txt"
	priority: 50
	message : "##RUNNING : Rscript to make BIGWIG with {input}"
	conda: "../envs/deeptools.yaml"
	shell:
		"bamCoverage -b {input.bam} -o {output} -of bigwig -bs {params.binsize} -p {threads} --exactScaling --normalizeUsing {params.normalization} "