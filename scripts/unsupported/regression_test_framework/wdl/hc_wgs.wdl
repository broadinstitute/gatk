import "hc_subworkflow.wdl" as  hc
import "genotypeconcordance.wdl" as gc

workflow HaplotypeCallerComparison {

    String picard_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.1-1504795437" # for picard, need to use the normal genomes in the cloud docker image

    File wgs_calling_interval_list = "gs://broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list"
  	Int haplotype_scatter_count = 50
  	Int break_bands_at_multiples_of = 1000000

  	Array[String] bams = [
    "gs://fc-483ee185-781e-4bb3-857e-8e06ae03ea3a/1e0003b4-e480-4118-b36b-69e31b1ee194/CramToBamWorkflow/9dbe7281-bbef-4a7f-a7d6-25e4e8d3b225/call-CramToBam/attempt-2/G94982.bam",
    "gs://broad-dsde-methods/cromwell-execution-26/CramToBamWorkflow/0acf2c44-b267-4249-82be-8d98029c237d/call-CramToBam/shard-1/G96830.bam",
    "gs://fc-483ee185-781e-4bb3-857e-8e06ae03ea3a/1e0003b4-e480-4118-b36b-69e31b1ee194/CramToBamWorkflow/6a4d50f2-e8d3-48c8-b3b0-a693f3319485/call-CramToBam/G96831.bam",
    "gs://broad-dsde-methods/cromwell-execution-26/CramToBamWorkflow/0acf2c44-b267-4249-82be-8d98029c237d/call-CramToBam/shard-0/G96832.bam"
  	]
  	Array[String] samples = ["G94982","G96830", "G96831", "G96832"]
  	String gc_sample = "NA12878"

  	File nist_truth = "gs://broad-dsde-methods/skwalker/NIST_truth.vcf.gz"
  	File nist_truth_index = "gs://broad-dsde-methods/skwalker/NIST_truth.vcf.gz.tbi"
  	Array[File] nist_intervals = [ "gs://broad-dsde-methods/skwalker/NIST_highconfidenceregions.interval_list", wgs_calling_interval_list]

  	Int disk_size = 400

  	File picard_jar = "gs://broad-dsde-methods/skwalker/picard_2.12.1.jar"
  	File gatk_jar = "gs://broad-dsde-methods/skwalker/GATK_3.8.jar"

    File ref_dict = "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    File ref_fasta = "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File ref_fasta_index = "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

	call ScatterIntervalList {
	    input:
	      interval_list = wgs_calling_interval_list,
	      scatter_count = haplotype_scatter_count,
	      break_bands_at_multiples_of = break_bands_at_multiples_of,
        picard_docker = picard_docker
  	}

  	scatter (i in range(length(bams))) {

  		call hc.HaplotypeCallerSubWorkflow {
  			input:
	  			intervals = ScatterIntervalList.out, 
	  			num_intervals = ScatterIntervalList.interval_count,
	  			bam = bams[i],
	  			base_file_name = samples[i],
	  			disk_size = disk_size,
          picard_docker = picard_docker
  		}

  		call gc.GenotypeConcordance as GenotypeConcordance3vNIST {
  			input:
	  			call_vcf = HaplotypeCallerSubWorkflow.gatk3_vcf,
	  			call_index = HaplotypeCallerSubWorkflow.gatk3_vcf_index,
	  			call_sample = "NA12878",
	  			truth_vcf = nist_truth,
	  			truth_index = nist_truth_index,
	  			truth_sample = "HG001",
	  			output_name = samples[i],
	  			intervals = nist_intervals,
	  			jar = picard_jar,
	  			disk_size = disk_size
  		}

  		call gc.GenotypeConcordance as GenotypeConcordance4vNIST {
  			input:
	  			call_vcf = HaplotypeCallerSubWorkflow.gatk4_vcf,
	  			call_index = HaplotypeCallerSubWorkflow.gatk4_vcf_index,
	  			call_sample = "NA12878",
	  			truth_vcf = nist_truth,
	  			truth_index = nist_truth_index,
	  			truth_sample = "HG001",
	  			output_name = samples[i],
	  			intervals =  nist_intervals,
	  			jar = picard_jar,
	  			disk_size = disk_size
  		}

  		call gc.GenotypeConcordance as GenotypeConcordance3v4 {
  			input:
	  			call_vcf = HaplotypeCallerSubWorkflow.gatk4_vcf,
	  			call_index = HaplotypeCallerSubWorkflow.gatk4_vcf_index,
	  			call_sample = "NA12878",
	  			truth_vcf = HaplotypeCallerSubWorkflow.gatk3_vcf,
	  			truth_index = HaplotypeCallerSubWorkflow.gatk3_vcf_index,
	  			truth_sample = "NA12878",
	  			output_name = samples[i],
	  			intervals = [wgs_calling_interval_list],
	  			jar = picard_jar,
	  			disk_size = disk_size

  		}

  		call VariantsToTable as VariantsToTable3 {
  			input:
  				vcf = HaplotypeCallerSubWorkflow.gatk3_vcf,
  				vcf_index =  HaplotypeCallerSubWorkflow.gatk3_vcf_index,
          jar = gatk_jar,
          output_name = samples[i] + ".table",
          ref = ref_fasta,
          ref_index = ref_fasta_index,
          ref_dict = ref_dict,
        	disk_size = disk_size
  		}

  		call VariantsToTable as VariantsToTable4 {
  			input:
          vcf = HaplotypeCallerSubWorkflow.gatk4_vcf,
    			vcf_index =  HaplotypeCallerSubWorkflow.gatk4_vcf_index,
          jar = gatk_jar,
          output_name = samples[i] + ".table",
          ref = ref_fasta,
          ref_index = ref_fasta_index,
          ref_dict = ref_dict,
        	disk_size = disk_size
  		}

  		call GCVariantsToTable as GCVariantsToTable3 {
  			input:
  				vcf = GenotypeConcordance3vNIST.output_vcf,
  				vcf_index =  GenotypeConcordance3vNIST.output_vcf_index,
          jar = gatk_jar,
          output_name = samples[i] + ".table",
          ref = ref_fasta,
          ref_index = ref_fasta_index,
          ref_dict = ref_dict,
        	disk_size = disk_size
  		}

  		call GCVariantsToTable as GCVariantsToTable4 {
  			input:
  				vcf = GenotypeConcordance4vNIST.output_vcf,
  				vcf_index =  GenotypeConcordance4vNIST.output_vcf_index,
          jar = gatk_jar,
          output_name = samples[i] + ".table",
          ref = ref_fasta,
          ref_index = ref_fasta_index,
          ref_dict = ref_dict,
        	disk_size = disk_size
  		}

  		call GCVariantsToTable as GCVariantsToTable3v4 {
  			input:
  				vcf = GenotypeConcordance3v4.output_vcf,
  				vcf_index =  GenotypeConcordance3v4.output_vcf_index,
          jar = gatk_jar,
          output_name = samples[i] + ".table",
          ref = ref_fasta,
          ref_index = ref_fasta_index,
          ref_dict = ref_dict,
        	disk_size = disk_size
  		}
  	}
}


# This task calls picard's IntervalListTools to scatter the input interval list into scatter_count sub interval lists
# Note that the number of sub interval lists may not be exactly equal to scatter_count.  There may be slightly more or less.
# Thus we have the block of python to count the number of generated sub interval lists.
task ScatterIntervalList {
  File interval_list
  Int scatter_count
  Int break_bands_at_multiples_of
  String picard_docker

  command <<<
    set -e
    mkdir out
    java -Xms1g -jar /usr/gitc/picard.jar \
      IntervalListTools \
      SCATTER_COUNT=${scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=${break_bands_at_multiples_of} \
      INPUT=${interval_list} \
      OUTPUT=out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals))
    CODE
  >>>
  output {
    Array[File] out = glob("out/*/*.interval_list")
    Int interval_count = read_int(stdout())
  }
  runtime {
    docker: picard_docker
    memory: "2 GB"
  }
}

task VariantsToTable {
    File vcf
    File vcf_index
    File jar
    File ref
    File ref_index
    File ref_dict
    String output_name
    Int disk_size

    command {
        java -jar ${jar} -T VariantsToTable -V ${vcf} -R ${ref} -L chr20:10000000-25000000 -o ${output_name} -SMA -F CHROM -F POS \
        -F REF -F ALT -F QUAL -F AC -F AF -F AN -F BaseQRankSum -F ClippingRankSum -F DB -F DP -F DS \
        -F END -F ExcessHet -F FS -F MLEAC -F HaplotypeScore -F InbreedingCoeff -F MLEAF -F MQ \
        -F MQRankSum -F ReadPosRankSum -F SOR -F QD -GF GT -GF AD -GF DP -GF GQ -GF PL -GF PGT -GF PID \
        -GF RGQ -GF SB
    }

    runtime {
         memory: "3 GB"
         disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_table = "${output_name}"
    }
}

task GCVariantsToTable {
    File vcf
    File vcf_index
    File jar
    File ref
    File ref_index
    File ref_dict
    String output_name
    Int disk_size

    command {
        java -jar ${jar} -T VariantsToTable -V ${vcf} -R ${ref} -L chr20:10000000-25000000 -o ${output_name} -SMA -F CHROM -F POS \
        -F REF -F ALT -F QUAL -F FILTER -F CONC_ST
    }

    runtime {
         memory: "3 GB"
         cpu: "1"
         disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_table = "${output_name}"
    }
}
