version 1.0

struct Reference {
    String reference_fasta
    String reference_fasta_index
    String reference_dict
    String wgs_calling_interval_list
    String exome_calling_interval_list
    String dbsnp_vcf
    String dbsnp_vcf_index

    # Resource files - used only for VQSR and VETS:
    # Probably should be done in another struct with more of a free form list and associated parameters to support VQSR/VETS
    String axiomPoly_resource_vcf
    String axiomPoly_resource_vcf_index
    String hapmap_resource_vcf
    String hapmap_resource_vcf_index
    String mills_resource_vcf
    String mills_resource_vcf_index
    String omni_resource_vcf
    String omni_resource_vcf_index
    String one_thousand_genomes_resource_vcf
    String one_thousand_genomes_resource_vcf_index
}
