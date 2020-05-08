workflow create_jg_tsv {
	String gvcf
	File interval_list
	String chr
	String basename = basename(gvcf, ".vcf.gz")


	call make_tsv {
	  input:
	    gvcf = gvcf,
	    basename = basename,
	    chr = chr,
	    interval_list = interval_list
	}
}

task make_tsv {
	String gvcf
	String basename
	String chr
	File interval_list

	command {
    set -e
    gatk --java-options "-Xms2500m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      BlahVariantWalker \
      -V ${gvcf} \
      -L ${interval_list} \
      --vet-out-path gs://broad-dsp-spec-ops/scratch/jsoto/dalio_exome_bq_uploads/${chr}/vet/${basename}_vet.tsv \
      --pet-out-path gs://broad-dsp-spec-ops/scratch/jsoto/dalio_exome_bq_uploads/${chr}/pet/${basename}_pet.tsv
  }

  runtime {
    docker: "jsotobroad/gatk:jg_1.1"
    preemptible: 0
    memory: "3.5 GB"
    cpu: "1"
    disks: "local-disk 20 HDD"
  }

}
