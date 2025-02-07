version 1.0

struct Reference {
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File wgs_calling_interval_list
    File exome_calling_interval_list
    File dbsnp_vcf
    File dbsnp_vcf_index

    # Resource files - used only for VQSR and VETS:
    # Probably should be done in another struct with more of a free form list and associated parameters to support VQSR/VETS
    File axiomPoly_resource_vcf
    File axiomPoly_resource_vcf_index
    File hapmap_resource_vcf
    File hapmap_resource_vcf_index
    File mills_resource_vcf
    File mills_resource_vcf_index
    File omni_resource_vcf
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf
    File one_thousand_genomes_resource_vcf_index
}