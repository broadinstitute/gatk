# -*- coding: utf-8 -*-
import uuid
import time
from datetime import datetime

import csv
import json
import ijson
import gzip
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
  "consequence": "consequence", # nullable -- this is now an array
  "dna_change": "hgvsc", # nullable
  "exon_number": "exons", # nullable
  "intron_number": "introns", # nullable
  # "splice_distance": "hgvsc", # nullable -- Lee has pushed this out of p0 for now
  "gene_id": "geneId", # nullable
  # "entrez_gene_id": "geneId", # nullable
  # "hgnc_gene_id": "hgncid", # nullable --  Lee has pushed this out of p0 for now
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
  # "hgnc_gene_id": "hgncid", # nullable genes.hgncid NOT HERE -- Lee has pushed this out of p0 for now
  # "gene_omim_id": "mimNumber", # nullable genes.omim.mimNumber HARD CODED FOR NOW
  "omim_phenotypes_id": "mimNumber", # nullable <-- lets talk arrays! genes.omim.phenotypes.mimNumber
  "omim_phenotypes_name": "phenotype" # nullable
}

significance_ordering = [ # TODO I have lowercased this---check that this is okay. Nirvana seems to have diff cases :P
  "benign",
  "likely benign",
  "uncertain significance",
  "aikely pathogenic",
  "aathogenic",
  "drug response",
  "association",
  "risk factor",
  "protective",
  "affects",
  "conflicting data from submitters",
  "other",
  "not provided",
  "'-'"
 ]


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
        row[vat_transcripts_fieldname] = transcript_fieldvalue

      if (test_clinvar.jsonvariant_line.get("spliceAI") != None) and (transcript_line.get("hgnc") != None):
        splice_ai_list = variant_line["spliceAI"]
        for splice_ai_obj in splice_ai_list:
          # get the splice AI value that matches to the transcript_line transcripts.hgnc to "spliceAI.hgnc"
          if splice_ai_obj.get("hgnc") == transcript_line["hgnc"]:
            for vat_splice_ai_fieldname in vat_nirvana_splice_ai_dictionary.keys():  # like "splice_ai_acceptor_gain_score"
              nirvana_splice_ai_fieldname = vat_nirvana_splice_ai_dictionary.get(vat_splice_ai_fieldname)
              splice_ai_fieldvalue = splice_ai_obj.get(nirvana_splice_ai_fieldname)
              row[vat_splice_ai_fieldname] = splice_ai_fieldvalue

    for vat_gvs_alleles_fieldname in vat_nirvana_gvs_alleles_dictionary.keys():  # like "gvs_all_ac"
      nirvana_gvs_alleles_fieldname = vat_nirvana_gvs_alleles_dictionary.get(vat_gvs_alleles_fieldname)
      gvs_alleles_fieldvalue = variant_line.get(nirvana_gvs_alleles_fieldname)
      row[vat_gvs_alleles_fieldname] = gvs_alleles_fieldvalue

    if variant_line.get("gnomad") != None:
      for vat_gnomad_fieldname in vat_nirvana_gnomad_dictionary.keys():  # like "gnomad_all_af"
        nirvana_gnomad_fieldname = vat_nirvana_gnomad_dictionary.get(vat_gnomad_fieldname)
        gnomad_fieldvalue = variant_line.get("gnomad").get(nirvana_gnomad_fieldname) # not a list like the others
        row[vat_gnomad_fieldname] = gnomad_fieldvalue

    if variant_line.get("clinvar") != None:
      clinvar_lines = variant_line["clinvar"]
      significance_values = [] # ordered by Benign, Likely Benign, Uncertain significance, Likely pathogenic, Pathogenic # https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
      updated_dates = [] # grab the most recent
      phenotypes = [] # ordered alphabetically
      clinvar_ids = [] # For easy validation downstream
      for clinvar_RCV_line in clinvar_lines:
        # get the clinvar lines with the id that starts with RCV
        if clinvar_RCV_line.get("id")[:3] == "RCV":
          clinvar_ids.append(clinvar_RCV_line.get("id"))
          significance_values.extend([x.lower() for x in clinvar_RCV_line.get("significance")])
          updated_dates.append(clinvar_RCV_line.get("lastUpdatedDate"))
          phenotypes.extend(clinvar_RCV_line.get("phenotypes"))
      # We want to collect all the significance values and order them by the significance_ordering list
      # So I will loop through the significance_ordering values and check for matching values in the significance_values list and put them in a new list
      ordered_significance_values = [] # how do we deal with cases issues?
      for value in significance_ordering:
        if value in significance_values:
          ordered_significance_values.append(value) # this adds the id to the end of the list
      values_not_accounted_for = list(set(ordered_significance_values).difference(significance_values))
      ordered_significance_values.extend(values_not_accounted_for) # add any values that aren't in significance_ordering to the end
      row["clinvar_id"] = clinvar_ids # array
      row["clinvar_classification"] = ordered_significance_values # special sorted array
      updated_dates.sort(key=lambda date: datetime.strptime(date, "%Y-%m-%d")) # note: method is in-place, and returns None
      row["clinvar_last_updated"] = updated_dates[-1] # most recent date
      row["clinvar_phenotype"] = sorted(phenotypes) # union of all phenotypes

    if variant_line.get("revel") != None:
      row["revel"] = variant_line.get("revel").get("score")

    return row


def make_positions_json(annotated_json, output_json):
  output_file=gzip.open(output_json, 'w')

  if annotated_json.endswith(".gz"):
      json_data = gzip.open(annotated_json, 'rb')
  else:
      json_data = open(annotated_json, 'rb')

  positions = ijson.items(json_data, 'positions.item', use_float=True)

  #TODO: filter down to Ensembl transcripts only.
  # so we can no longer just check if the transcripts exist or not---need to see if they exist and then if they have transcript_source of "Ensembl"
  # TODO is it better to do the transcript_source filtering here in the python, or in BQ? MAYBE BQ is cleaner?

  for p in positions:
    position=p['position']  # this is a required field -- do we want validation?
    variants=p['variants']
    # row for each variant - transcript
    # so let's start with each variant
    for variant in variants:
      # remember that we want one for each variant-transcript and variant-null for variants without transcripts
      if variant.get("transcripts") == None:
        # then we make a special row
        row = make_annotated_json_row(position, variant, None)
        json_str = json.dumps(row) + "\n"
        json_bytes = json_str.encode('utf-8')
        output_file.write(json_bytes)
      else:
        transcript_lines = variant.get("transcripts")
        # Collect all the transcript sources and check for if they contain Ensembl <-- this might be a good place for optimization
        sources = [transcript.get('source') for transcript in transcript_lines]
        if "Ensembl" in sources:
          for transcript in transcript_lines:
            if transcript.get('source') == "Ensembl":
              row = make_annotated_json_row(position, variant, transcript)
              json_str = json.dumps(row) + "\n"
              json_bytes = json_str.encode('utf-8')
              output_file.write(json_bytes)
        else:
          # if there are transcripts, but they are not Ensembl, we now only want one row in the VAT, not one row per transcript
          row = make_annotated_json_row(position, variant, None)
          json_str = json.dumps(row) + "\n"
          json_bytes = json_str.encode('utf-8')
          output_file.write(json_bytes)
  output_file.close()

def make_genes_json(annotated_json, output_genes_json):
  output_genes_file=gzip.open(output_genes_json, 'w')

  if annotated_json.endswith(".gz"):
      json_data = gzip.open(annotated_json, 'rb')
  else:
      json_data = open(annotated_json, 'rb')

  genes = ijson.items(json_data, 'genes.item', use_float=True)

  omim_fieldnames = list(vat_nirvana_omim_dictionary.keys())
  for gene_line in genes:
    if gene_line.get("omim") != None:
      row = {}
      row["gene_symbol"] = gene_line.get("name")
      omim_line = gene_line["omim"][0] # TODO I am making the huge assumption that we are only grabbing 1
      if len(gene_line.get("omim")) > 1:
        print(gene_line.get("name"),len(gene_line.get("omim")))
        raise ValueError("An assumption about the possible count of omim values is incorrect.")
      row["gene_omim_id"] = omim_line.get("mimNumber")
      if omim_line.get("phenotypes") != None:
        phenotypes = omim_line["phenotypes"]
        for vat_omim_fieldname in omim_fieldnames:  # like "mimNumber", "phenotype"
          # both of these are arrays -- TODO would it be clearer if this was just hard coded?
          omim_field_array=[]
          for phenotype in phenotypes:
            nirvana_omim_fieldname = vat_nirvana_omim_dictionary.get(vat_omim_fieldname)
            omim_fieldvalue = phenotype.get(nirvana_omim_fieldname, "")
            omim_field_array.append(omim_fieldvalue)
          row[vat_omim_fieldname] = omim_field_array
      json_str = json.dumps(row) + "\n"
      json_bytes = json_str.encode('utf-8')
      output_genes_file.write(json_bytes)
  output_genes_file.close()

def make_annotation_jsons(annotated_json, output_json, output_genes_json):
  make_positions_json(annotated_json, output_json)
  # we've already read the whole file once so we have to open it again
  # TODO: cleanup closing of file handles
  make_genes_json(annotated_json, output_genes_json)
  # TODO: should we be taking on the ".json.gz" to all the file names?


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

# TODO do I want to gsutil cp up the files into a bucket here?