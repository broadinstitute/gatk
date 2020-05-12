qsub -cwd -V -b y -pe smp_pe 16 \
    ./gatk --java-options "-Dlog4j.configurationFile=./log4j2-test.xml" \
           EstimateDragstrParameters \
	   --reference /Users/valentin/Analysis/dragstr/download/DRAGstr/references/hg38_alt_aware.fa \
	   -I ../out/NA12878_illumina_full_str_bam/NA12878_illumina_full_str_bam.bam \
	   --output ../out/NA12878_illumina_full_str_bam/NA12878_illumina_full_str_bam.gatk-dragstr-model.txt \
	   --sampling-min-mq 60 \
	   --sampling-loci-path /Users/valentin/Analysis/dragstr/download/DRAGstr/references/hg38_alt_aware.dragstr-sites.2.zip \
	   --parallel --threads 16 --verbosity DEBUG
