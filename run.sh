nextflow kuberun xsvato01/BTK_k8s -r main -pod-image 'cerit.io/nextflow/nextflow:22.06.1'\
	-resume -c zaloha_nextflow.config --datain /mnt/shared/MedGen/sequencing_results/primary_data/220708_TP53_20220708/raw_fastq
