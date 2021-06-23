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

## TODO in the future I want to map the original json key to a function!

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
        if nirvana_transcripts_fieldname == "isCanonical" and transcript_fieldvalue != True: # oooof this is ugly
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
      clinvar_lines = variant_line["clinvar"] # TODO I am making the huge assumption that this is correctly pulling the RCV one
      for clinvar_correct_line in clinvar_lines:
        # get the clinvar line with the id that starts with RCV
        if clinvar_correct_line.get("id")[:2] == "RCV":
              for vat_clinvar_fieldname in clinvar_fieldnames:  # like "clinvar_classification"
                nirvana_clinvar_fieldname = vat_nirvana_clinvar_dictionary.get(vat_clinvar_fieldname)
                clinvar_fieldvalue = clinvar_correct_line.get(nirvana_clinvar_fieldname)
                row[vat_clinvar_fieldname] = clinvar_fieldvalue

    if variant_line.get("revel") != None:
      row["revel"] = variant_line.get("revel").get("score")

    return row

def make_annotation_tables(annotated_json):
    variant_transcript_table_output_csv = "aou_alpha1_shard_annotations_vat.csv" # TODO get this from the annotated_json name
    genes_output_csv = "aou_alpha1_shard_annotations_vat_genes.csv" # TODO get this from the annotated_json name
  #  make_variant_transcript_table(annotated_json, variant_transcript_table_output_csv)
    make_genes_table(annotated_json, genes_output_csv)

def make_genes_table(annotated_json, output_csv): # why is this taking so long? I think we will want the reading to be done as a single script to make both tables
    with open(output_csv, 'w', newline='') as csvfile:
      omim_fieldnames = list(vat_nirvana_omim_dictionary.keys())
      fieldnames = omim_fieldnames
      fieldnames.extend(["gene_symbol","gene_omim_id"])
      writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
      writer.writeheader()
      json_data = open(annotated_json)
      data = json.load(json_data)
      gene_lines = data["genes"] # we want to write this out elsewhere?
      for gene_line in gene_lines:
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
          writer.writerow(row)

def make_variant_transcript_table(annotated_json, output_csv):
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
        # fieldnames.extend(omim_fieldnames)
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
              row = make_annotated_row(row_position, variant_line, None)
              writer.writerow(row)
            else:
              transcript_lines = variant_line.get("transcripts")
              for transcript_line in transcript_lines:
                row = make_annotated_row(row_position, variant_line, transcript_line)
                writer.writerow(row)


def make_annotation_json(annotated_json, output_json):
  output_file=open(output_json, 'w')
  json_data = open(annotated_json, 'r')
  data = json.load(json_data)
  annotated_json_lines = data['positions']
  for annotated_json_line in annotated_json_lines:
    position=annotated_json_line.get('position')
    variants=annotated_json_line.get('variants')
    for variant in variants:
      make_variant_line={} # remember that we want one for each variant-transcript and variant-null for variants without transcripts
      make_variant_line['position']=position
      for value in vat_nirvana_variants_dictionary.values():
        make_variant_line[value]=variant.get(value)
      if variant.get("revel") != None: # TODO check with Lee to see if this spacing is safe -- if there are no transcripts will there always never be a revel, spliceAT etc?
        make_variant_line["revel"] = variant.get("revel").get("score")
      if variant.get("spliceAI") != None:
        spliceAI_first_item = variant.get("spliceAI")[0]
        make_variant_line["spliceAI"]={}
        for value in vat_nirvana_splice_ai_dictionary.values():
          make_variant_line["spliceAI"][value]=spliceAI_first_item.get(value) # TODO I am making the huge assumption that we are only grabbing the first
      if variant.get("clinvar") != None:
        clinvar_items = variant.get("clinvar")
        make_variant_line["clinvar"]={}
        for clinvar_item in clinvar_items:
          if clinvar_item.get("id")[:3] == "RCV": # would it be better to save this processing for in BQ?
            for value in vat_nirvana_clinvar_dictionary.values():
              make_variant_line["clinvar"][value]=clinvar_item.get(value)
      if variant.get("gnomad") != None:
        make_variant_line["gnomad"]={}
        for value in vat_nirvana_gnomad_dictionary.values():
          make_variant_line["gnomad"][value]=variant.get("gnomad").get(value)
      if variant.get("transcripts") == None:
        back_json=json.dumps(make_variant_line)
        output_file.write(back_json)
        output_file.write("\n")
      else:
        transcript_lines = variant.get("transcripts")
        for transcript_line in transcript_lines:
          make_transcript_line=make_variant_line # remember that we want one for each variant-transcript and variant-null for variants without transcripts
          for value in vat_nirvana_transcripts_dictionary.values():
            make_transcript_line[value]=transcript_line.get(value)
          back_json=json.dumps(make_transcript_line)
          output_file.write(back_json)
          output_file.write("\n")
  output_file.close()

def make_annotation_json_better(annotated_json, output_json):
  output_file=open(output_json, 'w')
  json_data = open(annotated_json, 'r')
  data = json.load(json_data)
  annotated_position_json_lines = data['positions']
  annotated_genes_json_lines = data['genes']
  for annotated_json_line in annotated_position_json_lines:
      position=annotated_json_line.get('position')
      variants=annotated_json_line.get('variants')

  row_position = annotated_json_line.get("position") # this is a required field -- do we want validation?
            # want to write a csv row for each variant - transcript
            # so let's start with each variant
            variant_lines = annotated_json_line.get("variants")  # this is a required field -- do we want validation?
            for variant_line in variant_lines:
              if variant_line.get("transcripts") == None:
                # then we make a special row
                row = make_annotated_row(row_position, variant_line, None)
                writer.writerow(row)
              else:
                transcript_lines = variant_line.get("transcripts")
                for transcript_line in transcript_lines:
                  row = make_annotated_row(row_position, variant_line, transcript_line)
                  writer.writerow(row)


  for annotated_json_line in annotated_position_json_lines:
    position=annotated_json_line.get('position')
    variants=annotated_json_line.get('variants')
    for variant in variants:
      make_variant_line={} # remember that we want one for each variant-transcript and variant-null for variants without transcripts
      make_variant_line['position']=position
      for variant_entry in vat_nirvana_variants_dictionary:
        #this needs help
        make_variant_line[variant_entry.key()]=variant.get(variant_entry.value())
      if variant.get("revel") != None: # TODO check with Lee to see if this spacing is safe -- if there are no transcripts will there always never be a revel, spliceAT etc?
        make_variant_line["revel"] = variant.get("revel").get("score")
      if variant.get("spliceAI") != None:
        spliceAI_first_item = variant.get("spliceAI")[0]
        make_variant_line["spliceAI"]={}
        for value in vat_nirvana_splice_ai_dictionary.values():
          make_variant_line["spliceAI"][value]=spliceAI_first_item.get(value) # TODO I am making the huge assumption that we are only grabbing the first
      if variant.get("clinvar") != None:
        clinvar_items = variant.get("clinvar")
        make_variant_line["clinvar"]={}
        for clinvar_item in clinvar_items:
          if clinvar_item.get("id")[:3] == "RCV": # would it be better to save this processing for in BQ?
            for value in vat_nirvana_clinvar_dictionary.values():
              make_variant_line["clinvar"][value]=clinvar_item.get(value)
      if variant.get("gnomad") != None:
        make_variant_line["gnomad"]={}
        for value in vat_nirvana_gnomad_dictionary.values():
          make_variant_line["gnomad"][value]=variant.get("gnomad").get(value)
      if variant.get("transcripts") == None:
        back_json=json.dumps(make_variant_line)
        output_file.write(back_json)
        output_file.write("\n")
      else:
        transcript_lines = variant.get("transcripts")
        for transcript_line in transcript_lines:
          make_transcript_line=make_variant_line # remember that we want one for each variant-transcript and variant-null for variants without transcripts
          for value in vat_nirvana_transcripts_dictionary.values():
            make_transcript_line[value]=transcript_line.get(value)
          back_json=json.dumps(make_transcript_line)
          output_file.write(back_json)
          output_file.write("\n")
  output_file.close()

""" if __name__ == '__main__':
  parser = argparse.ArgumentParser(allow_abbrev=False, description='Extract a cohort from BigQuery Variant Store ')
  parser.add_argument('--fq_petvet_dataset',type=str, help='project.dataset location of pet/vet data', required=True)


  # Execute the parse_args() method
  args = parser.parse_args()

  make_annotation_tables(annotated_json) """

make_annotation_tables("ralpha1.json")
#make_annotation_json_better("ralpha1.json", "aou_alpha1_shard_annotations_bq_load.json")
#make_annotation_json("ralpha1.json", "aou_alpha1_shard_annotations_bq_load.json")
