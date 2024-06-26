process COLLECT_BASECALLED {
	tag "COLLECT_BASECALLED on $sample.name using $task.cpus CPUs and $task.memory memory"
	label "small_process"

	input:
	val(sample)

	output:
	tuple val(sample), path("*.fastq.gz")

	script:
	"""
	echo COLLECT_BASECALLED $sample.name
	cp  /mnt/share/710000-CEITEC/713000-cmm/713003-pospisilova/base/sequencing_results/primary_data/*${sample.run}/raw_fastq/${sample.name}* ./
	"""
} 


process TRIMMING {
	tag "trimming on $sample.name using $task.cpus CPUs and $task.memory memory"
	label "small_process"
	
	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("*.fastq.gz")

	script:
	""" 
	cutadapt -m 50 -o ${sample.name}.trimmed.R1.fastq.gz -p ${sample.name}.trimmed.R2.fastq.gz $reads
	"""
}

process FIRST_ALIGN_BAM {
	tag "first align on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
	label "medium_cpu"
	label "large_mem"

	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("${sample.name}.sorted.bam")
	tuple val(sample), path("${sample.name}.sorted.bai")

	script:
	rg = "\"@RG\\tID:${sample.name}\\tSM:${sample.name}\\tLB:${sample.name}\\tPL:ILLUMINA\""
	"""
	bwa mem -R ${rg} -t 4 ${params.refindex} $reads \
	| samtools view -Sb -o - - | samtools sort -o ${sample.name}.sorted.bam
	samtools index ${sample.name}.sorted.bam ${sample.name}.sorted.bai	
	"""
}

process FIRST_QC {
	tag "first QC on $sample.name using $task.cpus CPUs and $task.memory memory"
	label "smallest_process"
	container 'registry.gitlab.ics.muni.cz:443/450402/btk_k8s:16'	
	
	input:
	tuple val(sample), path(bam)

	output:
	path "*"

	script:
	"""
	samtools flagstat $bam > ${sample.name}.flagstat
	samtools stats $bam > ${sample.name}.samstats
	picard BedToIntervalList I=${params.covbed} O=${sample.name}.interval_list SD=${params.ref}.dict
	picard CollectHsMetrics I=$bam BAIT_INTERVALS=${sample.name}.interval_list TARGET_INTERVALS=${sample.name}.interval_list R=${params.ref}.fasta O=${sample.name}.aln_metrics
	"""
}

process MARK_DUPLICATES {
	tag "Mark duplicates on $sample.name using $task.cpus CPUs and $task.memory memory"
	label "small_process"

	input:
	tuple val(sample), path(bam)

	output:
	path "*.txt"
	tuple val(sample), path("*first.md.bam")
	path "*.bai"

	script:
	"""
	picard MarkDuplicates I=$bam M=${sample.name}.MD.metrics.txt O=${sample.name}.first.md.bam
	samtools index ${sample.name}.first.md.bam
	"""
}

process MULTIQC {
	tag "MultiQC using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/multiqc_reports/", mode:'copy'
	label "smallest_process"
	container "staphb/multiqc:1.19"
	// container 'registry.gitlab.ics.muni.cz:443/450402/tp53_nf:5'

	input:
	path '*'

	output:
	path '*.html'

	script:
	"""
	export LC_ALL=C.UTF-8
	export LANG=C.UTF-8
	multiqc . -n MultiQC-"`date +"%d-%m-%Y"`".html
	"""
}

process MUTECT2 {
	tag "MUTECT2 on $sample.name using $task.cpus CPUs and $task.memory memory"
	
	input:
	tuple val(sample), path(bam)
	
	output:
	tuple val(sample), path ('*.vcf')
	path '*'

	script:
	"""
	gatk Mutect2 --reference ${params.ref}.fasta --input ${bam} --annotation StrandArtifact --min-base-quality-score 10 --intervals $params.covbed -bamout ${sample.name}.bamout.bam --output ${sample.name}.mutect.vcf
	"""
}

process FILTER_MUTECT {
	tag "filter mutect on $sample.name using $task.cpus CPUs and $task.memory memory"
	label "smallest_process"

	input:
	tuple val(sample), path(vcf_input)
	
	output:
	tuple val(sample), path ('*.vcf')

	script:
	"""
	gatk FilterMutectCalls -V $vcf_input -O ${sample.name}.mutect.filt.vcf
	"""

}

process NORMALIZE_MUTECT {
	tag "normalize filtered mutect on $sample.name using $task.cpus CPUs $task.memory"
	label "smallest_process"

	input:
	tuple val(sample), path(vcf_input)
	
	output:
	tuple val(sample), path ('*.vcf')

	script:
	"""
	bcftools norm -m-both $vcf_input > ${sample.name}.mutect.filt.norm.vcf
	"""
}

process ANNOTATE_MUTECT {
	tag "annotate mutect on $sample.name using $task.cpus CPUs $task.memory"
	container "ensemblorg/ensembl-vep:release_108.0"
	publishDir "${params.outDirectory}/${sample.run}/variants/", mode:'copy'
	label "smallest_process"

	input:
	tuple val(sample), path(vcf_input)
	
	output:
	tuple val(sample), path('*.vcf')

	script:
	"""
	vep -i $vcf_input --cache --cache_version 108 --dir_cache $params.vep \
	--fasta ${params.ref}.fasta --merged --mane_select --offline --vcf --everything -o ${sample.name}.mutect2.filt.norm.vep.vcf
	"""	
}

process JOIN_VARS_TO_FILE {
	tag "JOIN_VARS_TO_FILE  using $task.cpus CPUs $task.memory"
	publishDir "${params.variantsDBdir}/", mode:'copy'
	label "smallest_process"

	input:
	path "VcfsToMerge"
	
	output:
	path "DB.bed"

	script:
	"""
	for vcf_file in $VcfsToMerge; do bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE]\\n' "\$vcf_file" >> temp.bed; done
	awk '!seen[\$0]++' temp.bed ${params.variantsDBdir}/DB.bed >> DB.bed
	# make sure there are no duplicates
	"""	
}

process CREATE_FULL_TABLE {
	tag "creating full table on $sample.name using $task.cpus CPUs $task.memory"
		publishDir "${params.outDirectory}/${sample.run}/create_full_table/", mode:'copy'

	label "smallest_process"

	input:
	tuple val(sample), path(vcf_input)
	
	output:
	tuple val(sample), path("${sample.name}.mutect2.filt.norm.vep.full.csv")

	script:
	"""
	python $params.vcftbl simple --build GRCh38 -i $vcf_input -t ${sample.name} -o ${sample.name}.mutect2.filt.norm.vep.full.csv
	"""	
}

process ALTER_FULL_TABLE {
	tag "ALTER_FULL_TABLE on $sample.name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/variants/", mode:'copy'
	label "smallest_process"

	input:
	tuple val(sample), path(variantsTableCsv), path(joinedTsv)//, val(ntotal)
	
	output:
	tuple val(sample), path("${sample.name}.final.csv")

	script:
	"""
	echo alter full table on $sample.name
	python ${params.mergetables} --table $variantsTableCsv --varlist $joinedTsv --outname ${sample.name}.final.csv
	"""	
}

process COVERAGE {
	tag "calculating coverage on $sample.name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/coverage/", mode:'copy'

	input:
	tuple val(sample), path(bam)
	
	output:
	tuple val(sample), path('*.PBcov.txt')

	script:
	"""
	bedtools coverage -abam $params.covbed -b $bam -d > ${sample.name}.PBcov.txt
	"""	
}



process COVERAGE_R {
	tag "R coverage on $sample.name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/coverage/", mode:'copy'
	label "smallest_process"

	input:
	tuple val(sample), path(pbcov)

	script:
	"""
	Rscript --vanilla $params.coverstat $pbcov ${params.outDirectory}/${sample.run}/coverage/${sample.name}.perexon_stat.txt
	"""
}

workflow {

runlist = channel.fromList(params.samples)
rawfastq = COLLECT_BASECALLED(runlist)

trimmed	= TRIMMING(rawfastq)
sortedbam	= FIRST_ALIGN_BAM(trimmed)
qc_files	= FIRST_QC(sortedbam[0])
qcdup_file	= MARK_DUPLICATES(sortedbam[0])
MULTIQC(qc_files.collect())
raw_vcf	= MUTECT2(qcdup_file[1]) //markdup.bam 
filtered	= FILTER_MUTECT(raw_vcf[0])
normalized	= NORMALIZE_MUTECT(filtered)
annotated	= ANNOTATE_MUTECT(normalized)
full_table	= CREATE_FULL_TABLE(annotated)
Vcf_paths = normalized.map({it -> [it[1]]})
Vcf_paths_collected = Vcf_paths.collect()
joined_vars = JOIN_VARS_TO_FILE(Vcf_paths_collected)
ALTER_FULL_TABLE(full_table.combine(joined_vars))
pbcov = COVERAGE(sortedbam[0])
COVERAGE_R(pbcov)
}
