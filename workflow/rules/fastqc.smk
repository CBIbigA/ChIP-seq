rule fastqc_raw:
	input:
		fastq=config["general"]["sample_dir"]+"/{prefix}"+config["raw_spec"]["ext"]
	output:
		config["general"]["experiment_name"]+"/QC/RAW/{prefix}_fastqc.html"
	params:
		dir=config["general"]["experiment_name"]+"/QC/RAW/",
		before=config["general"]["experiment_name"]+"/QC/RAW/{prefix}_fastqc.zip",
		after=config["general"]["experiment_name"]+"/QC/RAW/{prefix}_raw_fastqc.zip"
	threads: config["general"]["threads"]
	benchmark :
		config["general"]["experiment_name"]+"/benchmarks/fastqc_raw/{prefix}.txt"
	priority: 100
	conda: "../envs/fastqc.yaml"
	message : "##RUNNING : fastqc for {input.fastq}"
	shell: 
		"fastqc -q -t {threads} --outdir {params.dir} {input.fastq}"

rule rename_fastqc_raw:
	input:
		sample=expand(config["general"]["experiment_name"]+"/QC/RAW/{sample}_fastqc.html",sample=config["general"]["samples"])
	output : 
		sample=expand(config["general"]["experiment_name"]+"/QC/RAW/{sample}_raw_fastqc.zip",sample=config["general"]["samples"])
	params:
		dirraw = os.getcwd()+"/"+config["general"]["experiment_name"]+"/QC/RAW"
	run:
		shell("rename 's/_fastqc.zip/_raw_fastqc.zip/' {params.dirraw}/*_fastqc.zip")


rule fastqc_bam:
	input:
		bam=config["general"]["experiment_name"]+"/mapping/bam/{samtype}/{prefix}.{samtype}.bam"
	output:
		config["general"]["experiment_name"]+"/QC/bam/{prefix}.{samtype}_fastqc.html"
	params:
		dir=config["general"]["experiment_name"]+"/QC/bam/"
	threads: config["general"]["threads"]
	benchmark :
		config["general"]["experiment_name"]+"/benchmarks/fastqc_bam/{prefix}.{samtype}.txt"
	priority: 10
	message : "##RUNNING : fastqc for {input.bam}"
	conda: "../envs/fastqc.yaml"
	shell: "fastqc -q -t {threads} --outdir {params.dir} {input.bam}"


rule multi_qc:
	input:
		expand(config["general"]["experiment_name"]+"/QC/RAW/{sample}_raw_fastqc.zip",sample=config["general"]["samples"]),
		expand(config["general"]["experiment_name"]+"/QC/bam/{sample}.{samtype}_fastqc.html",samtype=samtype,sample=config["general"]["samples"]),
		expand(config["general"]["experiment_name"]+"/QC/STATS/{type}/{sample}.{samtype}.{type}",samtype=samtype,type = ["stats","idxstats","flagstat"],sample=config["general"]["samples"])
	output : 
		config["general"]["experiment_name"]+"/QC/MULTIQC/"+config["general"]["experiment_name"]+"_multiqc_report.html",
	params:
		title = config["general"]["experiment_name"],
		conf = config["fastqc"]["multi_qc_path"],
		output = config["general"]["experiment_name"]+"/QC/MULTIQC",
		filename = config["general"]["experiment_name"]+"_multiqc_report.html"
	priority: 50
	message : "##RUNNING : MultiQC"
	conda: "../envs/fastqc.yaml"
	shell:
		"export LC_ALL=C.UTF-8 && "
		"export LANG=C.UTF-8 && "
		"multiqc {params.title} --config {params.conf} --title {params.title} -o {params.output} --filename {params.filename}"

