# -*- coding: utf-8 -*-
import uuid
import time
import datetime

import csv
import json
import argparse

JOB_IDS = set()

#
# CONSTANTS
#
VAT_TABLE_PREFIX = "vat_"
SAMPLES_PER_PARTITION = 4000

## TODO in the future do I want to map the original json key to a function!?!?!?!?

vat_nirvana_positions_dictionary = {
  "position": "position", # required  TODO pull this out! check how nirvana handles "when the vcf position is not just the preceeding base for the variant"
}

vat_nirvana_variants_dictionary = {
  "vid": "vid", # required
  "contig": "chromosome", # required
  "ref_allele": "refAllele", # required
  "alt_allele": "altAllele", # required
  "variant_type": "variantType", # required
  "genomic_location": "hgvsg", # required
  "dbsnp_rsid": "dbsnp",  # nullable  -- TODO is this always a single val array?
}

vat_nirvana_transcripts_dictionary = {
  "transcript": "transcript", # nullable
  "gene_symbol": "hgnc", # nullable
  "transcript_source": "source", # nullable
  "aa_change": "hgvsp", # nullable
  "consequence": "consequence", # nullable -- TODO check on this one. May want this to be "[]"
  "dna_change": "hgvsc", # nullable
  "exon_number": "exons", # nullable
  "intron_number": "introns", # nullable
  "splice_distance": "hgvsc", # nullable -- TODO Additional processing:  Extract the splice distance from the hgvsc field
  "entrez_gene_id": "geneId", # nullable
  # "hgnc_gene_id": "hgncid", # nullable -- TODO ignore for now Lees notes dont match up here: genes.hgncid?
  "is_canonical_transcript": "isCanonical" # nullable -- (and lets make the nulls false)
}

vat_nirvana_gvs_alleles_dictionary = {
  "gvs_all_ac": "x", # required
  "gvs_all_an": "x", # required
  "gvs_all_af": "x" # required
}

vat_nirvana_revel_dictionary = {
  "revel": "score" # nullable
}

vat_nirvana_splice_ai_dictionary = {
  "splice_ai_acceptor_gain_score": "acceptorGainScore", # nullable
  "splice_ai_acceptor_gain_distance": "acceptorGainDistance", # nullable
  "splice_ai_acceptor_loss_score": "acceptorLossScore", # nullable
  "splice_ai_acceptor_loss_distance": "acceptorLossDistance", # nullable
  "splice_ai_donor_gain_score": "donorGainScore", # nullable
  "splice_ai_donor_gain_distance": "donorGainDistance", # nullable
  "splice_ai_donor_loss_score": "donorLossScore", # nullable
  "splice_ai_donor_loss_distance": "donorLossDistance" # nullable
}

vat_nirvana_clinvar_dictionary = {
  "clinvar_classification": "significance", # nullable
  "clinvar_last_updated": "lastUpdatedDate", # nullable
  "clinvar_phenotype": "phenotypes" # nullable -- currently here: "phenotypes":["not specified"] <-- lets talk arrays!
}

vat_nirvana_gnomad_dictionary = {
  "gnomad_all_af": "allAf", # nullable
  "gnomad_all_ac": "allAc", # nullable
  "gnomad_all_an": "allAn", # nullable
  "gnomad_max_af": "afrAf", # nullable THIS NEEDS MORE THAN JUST A MAPPING
  "gnomad_max_ac": "afrAc", # nullable THIS NEEDS MORE THAN JUST A MAPPING
  "gnomad_max_an": "afrAn", # nullable THIS NEEDS MORE THAN JUST A MAPPING
  "gnomad_max_subpop": "x" #TODO need to choose the correct mapping. Lee says maybe drop for now
}

vat_nirvana_omim_dictionary = { # TODO or should this be vat_nirvana_genes_dictionary ?
  #"hgnc_gene_id": "hgncid", # nullable genes.hgncid NOT HERE?!??!?
  # "gene_omim_id": "mimNumber", # nullable genes.omim.mimNumber HARD CODED FOR NOW
  "omim_phenotypes_id": "mimNumber", # nullable <-- lets talk arrays! genes.omim.phenotypes.mimNumber
  "omim_phenotypes_name": "phenotype" #TODO where is the omim stuff?????? <-- lets talk arrays! genes.omim.phenotypes.phenotype
}

def make_annotated_json_row(row_position, variant_line, transcript_line): # would it be better to not pass the transcript_line since its dupe data?
    row = {}
    row["position"] = row_position # this is a required field -- do we want validation? (what about validation for all the variants_fieldnames?)

    for vat_variants_fieldname in vat_nirvana_variants_dictionary.keys():  # like "contig"
      nirvana_variants_fieldname = vat_nirvana_variants_dictionary.get(vat_variants_fieldname)
      variant_fieldvalue = variant_line.get(nirvana_variants_fieldname)
      row[vat_variants_fieldname] = variant_fieldvalue

    if transcript_line != None:
      for vat_transcripts_fieldname in vat_nirvana_transcripts_dictionary.keys():  # like "transcript"
        nirvana_transcripts_fieldname = vat_nirvana_transcripts_dictionary.get(vat_transcripts_fieldname)
        transcript_fieldvalue = transcript_line.get(nirvana_transcripts_fieldname)
        if nirvana_transcripts_fieldname == "isCanonical" and transcript_fieldvalue != True: # oooof this is ugly
          transcript_fieldvalue = False
        row[vat_transcripts_fieldname] = transcript_fieldvalue

    for vat_gvs_alleles_fieldname in vat_nirvana_gvs_alleles_dictionary.keys():  # like "gvs_all_ac"
      nirvana_gvs_alleles_fieldname = vat_nirvana_gvs_alleles_dictionary.get(vat_gvs_alleles_fieldname)
      gvs_alleles_fieldvalue = variant_line.get(nirvana_gvs_alleles_fieldname)
      row[vat_gvs_alleles_fieldname] = gvs_alleles_fieldvalue

    if variant_line.get("spliceAI") != None:
      splice_ai_line = variant_line["spliceAI"][0] # TODO I am making the huge assumption that we are only grabbing 1
      for vat_splice_ai_fieldname in vat_nirvana_splice_ai_dictionary.keys():  # like "splice_ai_acceptor_gain_score"
        nirvana_splice_ai_fieldname = vat_nirvana_splice_ai_dictionary.get(vat_splice_ai_fieldname)
        splice_ai_fieldvalue = splice_ai_line.get(nirvana_splice_ai_fieldname)
        row[vat_splice_ai_fieldname] = splice_ai_fieldvalue

    if variant_line.get("gnomad") != None:
      for vat_gnomad_fieldname in vat_nirvana_gnomad_dictionary.keys():  # like "gnomad_all_af"
        nirvana_gnomad_fieldname = vat_nirvana_gnomad_dictionary.get(vat_gnomad_fieldname)
        gnomad_fieldvalue = variant_line.get("gnomad").get(nirvana_gnomad_fieldname) # not a list like the others
        row[vat_gnomad_fieldname] = gnomad_fieldvalue

    if variant_line.get("clinvar") != None:
      clinvar_lines = variant_line["clinvar"] # TODO I am making the huge assumption that this is correctly pulling the RCV one
      for clinvar_correct_line in clinvar_lines:
        # get the clinvar line with the id that starts with RCV
        if clinvar_correct_line.get("id")[:2] == "RCV": # TODO, does this need to be 3, am I being a dummy here?
              for vat_clinvar_fieldname in vat_nirvana_clinvar_dictionary.keys():  # like "clinvar_classification"
                nirvana_clinvar_fieldname = vat_nirvana_clinvar_dictionary.get(vat_clinvar_fieldname)
                clinvar_fieldvalue = clinvar_correct_line.get(nirvana_clinvar_fieldname)
                row[vat_clinvar_fieldname] = clinvar_fieldvalue

    if variant_line.get("revel") != None:
      row["revel"] = variant_line.get("revel").get("score")

    return row

def make_annotation_jsons(annotated_json, output_json, output_genes_json):
  output_file=open(output_json, 'w')
  output_genes_file=open(output_genes_json, 'w')
  json_data = open(annotated_json, 'r')
  data = json.load(json_data)
  annotated_position_lines = data['positions']
  annotated_gene_lines = data['genes'] # Is there any reason to split up the genes and positions json creation into two helper scripts?
  for annotated_json_line in annotated_position_lines:
    position=annotated_json_line.get('position')  # this is a required field -- do we want validation?
    variants=annotated_json_line.get('variants')
    # row for each variant - transcript
    # so let's start with each variant
    for variant in variants:
      # remember that we want one for each variant-transcript and variant-null for variants without transcripts
      if variant.get("transcripts") == None:
        # then we make a special row
        row = make_annotated_json_row(position, variant, None)
        back_json=json.dumps(row)
        output_file.write(back_json)
        output_file.write("\n")
      else:
        transcript_lines = variant.get("transcripts")
        for transcript in transcript_lines:
          row = make_annotated_json_row(position, variant, transcript)
          back_json=json.dumps(row)
          output_file.write(back_json)
          output_file.write("\n")
  output_file.close()
  omim_fieldnames = list(vat_nirvana_omim_dictionary.keys())
  for gene_line in annotated_gene_lines:
    if gene_line.get("omim") != None:
      row = {}
      row["gene_symbol"] = gene_line.get("name") # TODO throw an error if it's None? We dont care if there's no omim tho, right? right! cuz nothin would be there to add in the join
      omim_line = gene_line["omim"][0] # TODO I am making the huge assumption that we are only grabbing 1
      row["gene_omim_id"] = omim_line.get("mimNumber")
      if omim_line.get("phenotypes") != None:
        phenotypes_line = omim_line["phenotypes"][0] # TODO I am making the huge assumption that we are only grabbing 1
        for vat_omim_fieldname in omim_fieldnames:  # like "mimNumber", "phenotype"
          nirvana_omim_fieldname = vat_nirvana_omim_dictionary.get(vat_omim_fieldname)
          omim_fieldvalue = phenotypes_line.get(nirvana_omim_fieldname)
          row[vat_omim_fieldname] = omim_fieldvalue
      back_json=json.dumps(row)
      output_genes_file.write(back_json)
      output_genes_file.write("\n")
  output_genes_file.close()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(allow_abbrev=False, description='Create BQ load friendly jsons for VAT creation')
  parser.add_argument('--annotated_json',type=str, help='nirvana created annotation json', required=True)
  parser.add_argument('--output_vt_json',type=str, help='name of the vt json', required=True)
  parser.add_argument('--output_genes_json',type=str, help='name of the genes json', required=True)

  # Execute the parse_args() method
  args = parser.parse_args()

  make_annotation_jsons(args.annotated_json,
                        args.output_vt_json,
                        args.output_genes_json)

#make_annotation_jsons("ralpha1.json", "aou_alpha1_shard_annotations_bq_load.json", "aou_alpha1_shard_genes_bq_load.json")
# TODO add this to the ah_var_store docker!!!
# TODO do I want to gsutil cp up the files into a bucket here?