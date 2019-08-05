import argparse
import re
import gzip

def main(vcf, out, dist):
	with gzip.open(vcf,"r") as inp, open(out, "w") as out:

		out_header = "sample\tchrom\tpos\tref\talt\tfilt\tgt\tpgt\tdp\tgq\tpid\tvar_type\tsample1\tchrom1\tpos1\tref1\talt1\tfilt1\tgt1\tpgt1\tdp1\tgq1\tpid1\tvar_type1\tdistance\tsame_pid\trelative_phase\n"
		out.write(out_header)
		
		evaluated_chrom = "0"
		
		for variant_line in inp:
			if variant_line.startswith("##"):
				continue
			if variant_line.startswith("#CHROM"):
				header = variant_line.strip().split("\t")
				samples = header[9:len(header)]
				
				#Doing this to be able to print variants where the next variant is less than x bp away
				last_variant_info_dict = {key: list() for key in samples}
				continue
			
			vcf_elems = variant_line.strip().split("\t")
			chrom, pos , rsid,ref, alt, qual, filt, info, format_field = vcf_elems[0:9]		

			genotypes = vcf_elems[9:len(vcf_elems)]
			genotypes = [zip(format_field.split(':'), x.split(':')) for x in genotypes]
			

			#Only evaluating biallelic variants
			if len(alt.split(",")) > 1: 
				continue

			if len(alt) != len(ref):
				var_type = "indel"
			if len(alt) == len(ref):
				var_type = "SNV"

			if chrom != evaluated_chrom:
				print "starting chromosome " + chrom
				evaluated_chrom = chrom 
				#Reinitialize dictionary for every chromosome, likely not super necessary since VCF is in order 
				last_variant_dict = dict.fromkeys(samples,0)
			
			for ind in range(len(genotypes)):
				genotype_dic  = dict(genotypes[ind])
				sample =  samples[ind]

				pgt = genotype_dic.get("PGT",None)
				pid = genotype_dic.get("PID",None)
				if pid is None and pgt is not None:
					print "\t".join([sample,chrom,pos,ref,alt,filt, gt]) + "\n"
					raise ValueError("PGT without PID")
					break
				elif pid is not None and pgt is None:
					print "\t".join([sample,chrom,pos,ref,alt,filt, gt]) + "\n"
					raise ValueError("PID without PGT")
					break

			 	gt = genotype_dic .get("GT",None)
			 	dp = genotype_dic .get("DP",None)
			 	gq = genotype_dic .get("GQ",None)


				if gt == "0/1":
					last_variant_for_sample = last_variant_dict[sample] 
	 		 		current_variant = pos
	 		 		distance = int(current_variant)-int(last_variant_for_sample)
					
	 		 		line_to_print = [sample,chrom,pos,ref,alt,filt, gt,pgt, str(dp),str(gq), pid, var_type]

			 		if distance <= int(dist):
			 			last_variant_for_sample_info = last_variant_info_dict[sample]

			 			last_pgt = last_variant_for_sample_info[7]
		 		 		last_pid = last_variant_for_sample_info[10]

		 		 		same_pid = "NA"
		 		 		relative_phase = "NA"
		 		 		if pid is not None and pid != "." and last_pid is not None and last_pid != "." and last_pid != "None":
		 		 			if pid == last_pid:
		 		 				same_pid = "True"
		 		 				if pgt == last_pgt:
		 		 					relative_phase = "cis"
		 		 				else:
		 		 					relative_phase = "trans"
		 		 			else:
		 		 				same_pid = "False"
		 		 		else:
		 		 			if pid is None:
		 		 				line_to_print[10] = "None"
		 		 				line_to_print[7] = "None"
		 		 			if last_pid is None:
		 		 				last_variant_for_sample_info[10] = "None"
		 		 				last_variant_for_sample_info[7] = "None"

		 		 		pair_info = [str(distance), same_pid, relative_phase]

						final_print = last_variant_for_sample_info + line_to_print + pair_info
						out.write("\t".join(final_print) + "\n")
			 		
			 		last_variant_dict[sample] = current_variant
			 		last_variant_info_dict[sample] = line_to_print



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='''MNPs''')
	parser.add_argument("-vcf",help="vcf")
	parser.add_argument("-out", help = "output file (.tsv or .txt)")
	parser.add_argument("-dist",help="max distance")
	args = parser.parse_args()
	main(args.vcf, args.out, args.dist)






