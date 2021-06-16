# -*- coding: utf-8 -*-
import uuid
import time
import datetime

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig
from google.oauth2 import service_account

import csv
import json


JOB_IDS = set()

#
# CONSTANTS
#
VAT_TABLE_PREFIX = "vat_"
SAMPLES_PER_PARTITION = 4000

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
  "gnomad_max_af": "afrAf", # nullable
  "gnomad_max_ac": "afrAc", # nullable
  "gnomad_max_an": "afrAn", # nullable
  "gnomad_max_subpop": "x" #TODO need to choose the correct mapping. Lee says maybe drop for now
}

vat_nirvana_omim_dictionary = { # TODO or should this be vat_nirvana_genes_dictionary ?
  "gene_omim_id": "mimNumber", # nullable
  "omim_phenotypes_id": "phenotype", # nullable <-- lets talk arrays!
  "omim_phenotypes_name": "x" #TODO where is the omim stuff?????? <-- lets talk arrays!
}

def make_fieldnames():
    positions_fieldnames = list(vat_nirvana_positions_dictionary.keys())
    variants_fieldnames = list(vat_nirvana_variants_dictionary.keys())
    transcripts_fieldnames = list(vat_nirvana_transcripts_dictionary.keys())
    gvs_alleles_fieldnames = list(vat_nirvana_gvs_alleles_dictionary.keys())
    revel_fieldnames = list(vat_nirvana_revel_dictionary)
    splice_ai_fieldnames = list(vat_nirvana_splice_ai_dictionary.keys())
    clinvar_fieldnames = list(vat_nirvana_clinvar_dictionary.keys())
    gnomad_fieldnames = list(vat_nirvana_gnomad_dictionary.keys())
    omim_fieldnames = list(vat_nirvana_omim_dictionary.keys())

    return positions_fieldnames, variants_fieldnames, transcripts_fieldnames, gvs_alleles_fieldnames, revel_fieldnames, splice_ai_fieldnames, clinvar_fieldnames, gnomad_fieldnames, omim_fieldnames



def make_annotated_row(row_position, variant_line, transcript_line): # would it be better to not pass the transcript_line since its dupe data?
    positions_fieldnames, variants_fieldnames, transcripts_fieldnames, gvs_alleles_fieldnames, revel_fieldnames, splice_ai_fieldnames, clinvar_fieldnames, gnomad_fieldnames, omim_fieldnames = make_fieldnames()
    row = {}
    row["position"] = row_position # this is a required field -- do we want validation? (what about validation for all the variants_fieldnames?)

    for vat_variants_fieldname in variants_fieldnames:  # like "contig"
      nirvana_variants_fieldname = vat_nirvana_variants_dictionary.get(vat_variants_fieldname)
      variant_fieldvalue = variant_line.get(nirvana_variants_fieldname)
      row[vat_variants_fieldname] = variant_fieldvalue

    if transcript_line != None:
      for vat_transcripts_fieldname in transcripts_fieldnames:  # like "transcript"
        nirvana_transcripts_fieldname = vat_nirvana_transcripts_dictionary.get(vat_transcripts_fieldname)
        transcript_fieldvalue = transcript_line.get(nirvana_transcripts_fieldname)
        if nirvana_transcripts_fieldname = "isCanonical" & transcript_fieldvalue != True: # oooof this is ugly
          transcript_fieldvalue = False
        row[vat_transcripts_fieldname] = transcript_fieldvalue

    for vat_gvs_alleles_fieldname in gvs_alleles_fieldnames:  # like "gvs_all_ac"
      nirvana_gvs_alleles_fieldname = vat_nirvana_gvs_alleles_dictionary.get(vat_gvs_alleles_fieldname)
      gvs_alleles_fieldvalue = variant_line.get(nirvana_gvs_alleles_fieldname)
      row[vat_gvs_alleles_fieldname] = gvs_alleles_fieldvalue

    if variant_line.get("gnomad") != None:
      for vat_gnomad_fieldname in gnomad_fieldnames:  # like "gnomad_all_af"
        nirvana_gnomad_fieldname = vat_nirvana_gnomad_dictionary.get(vat_gnomad_fieldname)
        gnomad_fieldvalue = variant_line.get("gnomad").get(nirvana_gnomad_fieldname) # not a list like the others
        row[vat_gnomad_fieldname] = gnomad_fieldvalue

    if variant_line.get("spliceAI") != None:
      splice_ai_line = variant_line["spliceAI"][0] # TODO I am making the huge assumption that we are only grabbing 1
      for vat_splice_ai_fieldname in splice_ai_fieldnames:  # like "splice_ai_acceptor_gain_score"
        nirvana_splice_ai_fieldname = vat_nirvana_splice_ai_dictionary.get(vat_splice_ai_fieldname)
        splice_ai_fieldvalue = splice_ai_line.get(nirvana_splice_ai_fieldname)
        row[vat_splice_ai_fieldname] = splice_ai_fieldvalue

    if variant_line.get("clinvar") != None:
      clinvar_line = variant_line["clinvar"][0] # TODO I am making the huge assumption that we are only grabbing 1
      for vat_clinvar_fieldname in clinvar_fieldnames:  # like "clinvar_classification"
        nirvana_clinvar_fieldname = vat_nirvana_clinvar_dictionary.get(vat_clinvar_fieldname)
        clinvar_fieldvalue = clinvar_line.get(nirvana_clinvar_fieldname)
        row[vat_clinvar_fieldname] = clinvar_fieldvalue

    if variant_line.get("genes") != None:
      omim_line = variant_line["genes"][0] # TODO I am making the huge assumption that we are only grabbing 1
      for vat_omim_fieldname in omim_fieldnames:  # like "clinvar_classification"
        nirvana_omim_fieldname = vat_nirvana_omim_dictionary.get(vat_omim_fieldname)
        omim_fieldvalue = omim_line.get(nirvana_omim_fieldname)
        row[vat_omim_fieldname] = omim_fieldvalue

    if variant_line.get("revel") != None:
      row["revel"] = variant_line.get("revel").get("score")

    return row

def make_annotation_table(annotated_json, output_csv):
    with open(output_csv, 'w', newline='') as csvfile:
        positions_fieldnames, variants_fieldnames, transcripts_fieldnames, gvs_alleles_fieldnames, revel_fieldnames, splice_ai_fieldnames, clinvar_fieldnames, gnomad_fieldnames, omim_fieldnames = make_fieldnames()
        fieldnames = positions_fieldnames
        fieldnames.extend(variants_fieldnames)
        fieldnames.extend(transcripts_fieldnames)
        fieldnames.extend(gvs_alleles_fieldnames)
        fieldnames.extend(revel_fieldnames)
        fieldnames.extend(splice_ai_fieldnames)
        fieldnames.extend(clinvar_fieldnames)
        fieldnames.extend(gnomad_fieldnames)
        fieldnames.extend(omim_fieldnames)
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        # in the future can we simplify this pre-processing?
        json_data = open(annotated_json)
        data = json.load(json_data)
        annotated_json_lines = data["positions"]
        for annotated_json_line in annotated_json_lines:
          row_position = annotated_json_line.get("position") # this is a required field -- do we want validation?
          # want to write a csv row for each variant - transcript
          # so let's start with each variant
          variant_lines = annotated_json_line.get("variants")  # this is a required field -- do we want validation?
          for variant_line in variant_lines:
            if variant_line.get("transcripts") == None:
              # then we make a special row
              print("can I handle this?")
              row = make_annotated_row(row_position, variant_line, None)
              print("could I handle this?")
              print(row)
              writer.writerow(row)
            else:
              transcript_lines = variant_line.get("transcripts")
              for transcript_line in transcript_lines:
                row = make_annotated_row(row_position, variant_line, transcript_line)
                print(row)
                writer.writerow(row)

"""               row = {}
              row["position"] = annotated_json_line.get("position") # this is a required field -- do we want validation? (what about validation for all the variants_fieldnames?)

              for vat_variants_fieldname in variants_fieldnames:  # like "contig"
                nirvana_variants_fieldname = vat_nirvana_variants_dictionary.get(vat_variants_fieldname)
                variant_fieldvalue = variant_line.get(nirvana_variants_fieldname)
                row[vat_variants_fieldname] = variant_fieldvalue

              for vat_transcripts_fieldname in transcripts_fieldnames:  # like "transcript"
                nirvana_transcripts_fieldname = vat_nirvana_transcripts_dictionary.get(vat_transcripts_fieldname)
                transcript_fieldvalue = transcript_line.get(nirvana_transcripts_fieldname)
                row[vat_transcripts_fieldname] = transcript_fieldvalue

              for vat_gvs_alleles_fieldname in gvs_alleles_fieldnames:  # like "gvs_all_ac"
                nirvana_gvs_alleles_fieldname = vat_nirvana_gvs_alleles_dictionary.get(vat_gvs_alleles_fieldname)
                gvs_alleles_fieldvalue = variant_line.get(nirvana_gvs_alleles_fieldname)
                row[vat_gvs_alleles_fieldname] = gvs_alleles_fieldvalue

              if variant_line.get("gnomad") != None:
                for vat_gnomad_fieldname in gnomad_fieldnames:  # like "gnomad_all_af"
                   nirvana_gnomad_fieldname = vat_nirvana_gnomad_dictionary.get(vat_gnomad_fieldname)
                   gnomad_fieldvalue = variant_line.get("gnomad").get(nirvana_gnomad_fieldname) # not a list like the others
                   row[vat_gnomad_fieldname] = gnomad_fieldvalue

              if variant_line.get("spliceAI") != None:
                splice_ai_line = variant_line["spliceAI"][0] # TODO I am making the huge assumption that we are only grabbing 1
                for vat_splice_ai_fieldname in splice_ai_fieldnames:  # like "splice_ai_acceptor_gain_score"
                  nirvana_splice_ai_fieldname = vat_nirvana_splice_ai_dictionary.get(vat_splice_ai_fieldname)
                  splice_ai_fieldvalue = splice_ai_line.get(nirvana_splice_ai_fieldname)
                  row[vat_splice_ai_fieldname] = splice_ai_fieldvalue

              if variant_line.get("clinvar") != None:
                clinvar_line = variant_line["clinvar"][0] # TODO I am making the huge assumption that we are only grabbing 1
                for vat_clinvar_fieldname in clinvar_fieldnames:  # like "clinvar_classification"
                  nirvana_clinvar_fieldname = vat_nirvana_clinvar_dictionary.get(vat_clinvar_fieldname)
                  clinvar_fieldvalue = clinvar_line.get(nirvana_clinvar_fieldname)
                  row[vat_clinvar_fieldname] = clinvar_fieldvalue

              if variant_line.get("genes") != None:
                omim_line = variant_line["genes"][0] # TODO I am making the huge assumption that we are only grabbing 1
                for vat_omim_fieldname in omim_fieldnames:  # like "clinvar_classification"
                  nirvana_omim_fieldname = vat_nirvana_omim_dictionary.get(vat_omim_fieldname)
                  omim_fieldvalue = omim_line.get(nirvana_omim_fieldname)
                  row[vat_omim_fieldname] = omim_fieldvalue

              if variant_line.get("revel") != None:
                row["revel"] = variant_line.get("revel").get("score")


                print(row)
                writer.writerow(row)"""


""" if __name__ == '__main__':
  parser = argparse.ArgumentParser(allow_abbrev=False, description='Extract a cohort from BigQuery Variant Store ')
  parser.add_argument('--fq_petvet_dataset',type=str, help='project.dataset location of pet/vet data', required=True)


  # Execute the parse_args() method
  args = parser.parse_args()

  make_annotation_table(annotated_json) """

make_annotation_table("hello_did_I_annotate.json", "vat_annotations.csv")
