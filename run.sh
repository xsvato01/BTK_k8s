nextflow kuberun xsvato01/BTK_k8s -r main -pod-image 'cerit.io/nextflow/nextflow:22.06.1'\
	-resume -c zaloha_nextflow.config --datain /mnt/shared/MedGen/sequencing_results/primary_data/230106_TP53_20230106/raw_fastq/
