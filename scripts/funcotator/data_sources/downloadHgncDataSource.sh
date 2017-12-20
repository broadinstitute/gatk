#!/usr/bin/env bash

# Downloads the HGNC data source from the HGNC website.

curl 'https://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_locus_type&col=gd_locus_group&col=gd_prev_sym&col=gd_prev_name&col=gd_aliases&col=gd_name_aliases&col=gd_pub_chrom_map&col=gd_date_mod&col=gd_date_sym_change&col=gd_date_name_change&col=gd_pub_acc_ids&col=gd_enz_ids&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=gd_pubmed_ids&col=gd_pub_refseq_ids&col=family.id&col=family.name&col=gd_ccds_ids&col=gd_vega_ids&col=md_eg_id&col=md_mim_id&col=md_refseq_id&col=md_prot_id&col=md_ensembl_id&col=md_ucsc_id&status=Approved&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit' > hgnc_download_$(date +%b%d%Y).tsv

