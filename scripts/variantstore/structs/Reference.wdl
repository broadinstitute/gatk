version 1.0

struct Reference {
    String reference_version              # 38 or CUSTOM - a magic value sent to several java programs to handle contig naming
    String reference_fasta          # Path to the reference FASTA file
    String reference_fasta_index    # Path to the reference FASTA index file
    String reference_dict           # Path to the reference dictionary file

    String wgs_calling_interval_list
    String exome_calling_interval_list
}
