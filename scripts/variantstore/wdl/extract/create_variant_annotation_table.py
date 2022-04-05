# -*- coding: utf-8 -*-
import uuid
import time
from datetime import datetime

import csv
import json
import ijson
import gzip
import argparse

vat_nirvana_positions_dictionary = {
  "position": "position", # required
  "ref_allele": "refAllele", # required
  "alt_allele": "altAlleles", # required
}

vat_nirvana_variants_dictionary = {
  "vid": "vid", # required
  "contig": "chromosome", # required
  # "ref_allele": "refAllele", # required # do we want to grab both to compare?
  # "alt_allele": "altAllele", # required
  "variant_type": "variantType", # required
  "genomic_location": "hgvsg", # required
  "dbsnp_rsid": "dbsnp",  # nullable
}

vat_nirvana_transcripts_dictionary = {
  "transcript": "transcript", # nullable
  "gene_symbol": "hgnc", # nullable
  "transcript_source": "source", # nullable
  "aa_change": "hgvsp", # nullable
  "consequence": "consequence", # nullable
  "dna_change_in_transcript": "hgvsc", # nullable
  "exon_number": "exons", # nullable
  "intron_number": "introns", # nullable
  "gene_id": "geneId", # nullable
  "is_canonical_transcript": "isCanonical" # nullable
}

vat_nirvana_gvs_alleles_dictionary = {
  "gvs_all_an": "AN", # required
  "gvs_all_ac": "AC", # required
  "gvs_all_af": "AF" # required
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
  "clinvar_phenotype": "phenotypes" # nullable
}

vat_nirvana_gnomad_dictionary = {
  "gnomad_all_af": "allAf", # nullable
  "gnomad_all_ac": "allAc", # nullable
  "gnomad_all_an": "allAn", # nullable
  "gnomad_failed_filter": "failedFilter" #nullable
}

vat_nirvana_omim_dictionary = {
  "omim_phenotypes_id": "mimNumber", # nullable
  "omim_phenotypes_name": "phenotype" # nullable
}

significance_ordering = [
  "benign",
  "likely benign",
  "uncertain significance",
  "likely pathogenic",
  "pathogenic",
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

gnomad_ordering = [
 "afr",
 "amr",
 "asj",
 "eas",
 "fin",
 "nfr",
 "oth",
 "sas"
]

gvs_subpopulations = [
  "afr",
  "amr",
  "eas",
  "eur",
  "mid",
  "oth",
  "sas"
]

def check_filtering(variant):
    # skip any row (with a warning) if no gvsAnnotations exist
    if variant.get("gvsAnnotations") == None: # <-- enum since we need this to be in tandem with the custom annotations header / template
      print("WARNING: There has been an error in creating custom annotations for AC/AF/AN", variant.get("vid"))
      return False
    # skip any row (with a warning) if the AC value is 0
    elif variant["gvsAnnotations"].get("AC") == 0:
      print("WARNING: Its AC is 0 so we are dropping this variant", variant.get("vid"))
      return False
    # skip any row (with a warning) if AC, AN or AF is missing
    elif variant["gvsAnnotations"].get("AC") == None:
      print("WARNING: There has been an error-- there is no AC value---should AN be 0 for this variant?", variant.get("vid"))
      return False
    elif variant["gvsAnnotations"].get("AN") == None:
      print("WARNING: There has been an error-- there is an AC value---but no AN value", variant.get("vid"))
      return False
    elif variant["gvsAnnotations"].get("AF") == None:
      print("WARNING: There has been an error-- there is an AC value---but no AF value", variant.get("vid"))
      return False
    else:
      return True

def get_gnomad_subpop(gnomad_obj):
    row = {}
    max_ac = None
    max_an = None
    max_af = None
    max_subpop = ""
    for gnomad_subpop in gnomad_ordering: # since we cycle through this in order, if there is a tie, we just ignore it because we wouldn't chose it anyway based on order
      subpop_af_key = "".join([gnomad_subpop, "Af"])
      subpop_ac_key = "".join([gnomad_subpop, "Ac"])
      subpop_an_key = "".join([gnomad_subpop, "An"])
      subpop_af_val = gnomad_obj.get(subpop_af_key) # note that these can be null if there is no value in the annotations. They will be null in the VAT
      subpop_ac_val = gnomad_obj.get(subpop_ac_key)
      subpop_an_val = gnomad_obj.get(subpop_an_key)
      # here we set the subpopulation ac/an/af values
      row["_".join(["gnomad", gnomad_subpop, "an"])] = subpop_an_val
      row["_".join(["gnomad", gnomad_subpop, "ac"])] = subpop_ac_val
      row["_".join(["gnomad", gnomad_subpop, "af"])] = subpop_af_val
      if subpop_af_val != None and (max_af == None or subpop_af_val > max_af): # this will also set the first max_af value
        max_subpop = gnomad_subpop
        max_ac = subpop_ac_val
        max_an = subpop_an_val
        max_af = subpop_af_val
    # here we set the MAX subpopulation ac/an/af values
    row["gnomad_max_subpop"] = max_subpop
    row["gnomad_max_ac"] = max_ac
    row["gnomad_max_an"] = max_an
    row["gnomad_max_af"] = max_af
    return row

def get_subpopulation_calculations(subpop_annotations):
    row = {}
    max_ac = None
    max_an = None
    max_af = None
    max_sc = None
    max_subpop = ""
    for gvs_subpop in gvs_subpopulations:
      subpop_ac_val = subpop_annotations.get("_".join(["AC", gvs_subpop]))
      subpop_an_val = subpop_annotations.get("_".join(["AN", gvs_subpop]))
      subpop_af_val = subpop_annotations.get("_".join(["AF", gvs_subpop]))
       # note the assumption is made that AC_Hom must be even because by it's nature it means there are two, but there could be an error
      subpop_sc_val = int(subpop_annotations.get("_".join(["AC_Hom", gvs_subpop])) / 2 ) + subpop_annotations.get("_".join(["AC_Het", gvs_subpop])) + subpop_annotations.get("_".join(["AC_Hemi", gvs_subpop]))
      # here we set the subpopulation ac/an/af values
      row["_".join(["gvs", gvs_subpop, "ac"])] = subpop_ac_val
      row["_".join(["gvs", gvs_subpop, "an"])] = subpop_an_val
      row["_".join(["gvs", gvs_subpop, "af"])] = subpop_af_val
      row["_".join(["gvs", gvs_subpop, "sc"])] = subpop_sc_val
      if subpop_af_val != None and (max_af == None or subpop_af_val > max_af): # this will also set the first max_af value
        max_subpop = gvs_subpop
        max_ac = subpop_ac_val
        max_an = subpop_an_val
        max_af = subpop_af_val
        max_sc = subpop_sc_val
    # here we set the MAX subpopulation ac/an/af values
    row["gvs_max_subpop"] = max_subpop
    row["gvs_max_ac"] = max_ac
    row["gvs_max_an"] = max_an
    row["gvs_max_af"] = max_af
    row["gvs_max_sc"] = max_sc
    return row

def make_annotated_json_row(row_position, row_ref, row_alt, variant_line, transcript_line):
    row = {}
    row["position"] = row_position # this is a required field
    row["ref_allele"] = row_ref # this is a required field
    row["alt_allele"] = row_alt # this is a required field

    for vat_variants_fieldname in vat_nirvana_variants_dictionary.keys():  # like "contig"
      nirvana_variants_fieldname = vat_nirvana_variants_dictionary.get(vat_variants_fieldname)
      variant_fieldvalue = variant_line.get(nirvana_variants_fieldname)
      row[vat_variants_fieldname] = variant_fieldvalue

    if transcript_line != None:
      for vat_transcripts_fieldname in vat_nirvana_transcripts_dictionary.keys():  # like "transcript"
        nirvana_transcripts_fieldname = vat_nirvana_transcripts_dictionary.get(vat_transcripts_fieldname)
        transcript_fieldvalue = transcript_line.get(nirvana_transcripts_fieldname)
        row[vat_transcripts_fieldname] = transcript_fieldvalue

      if (variant_line.get("spliceAI") != None) and (transcript_line.get("hgnc") != None):
        splice_ai_list = variant_line["spliceAI"]
        for splice_ai_obj in splice_ai_list:
          # get the splice AI value that matches to the transcript_line transcripts.hgnc to "spliceAI.hgnc"
          if splice_ai_obj.get("hgnc") == transcript_line["hgnc"]:
            for vat_splice_ai_fieldname in vat_nirvana_splice_ai_dictionary.keys():  # like "splice_ai_acceptor_gain_score"
              nirvana_splice_ai_fieldname = vat_nirvana_splice_ai_dictionary.get(vat_splice_ai_fieldname)
              splice_ai_fieldvalue = splice_ai_obj.get(nirvana_splice_ai_fieldname)
              row[vat_splice_ai_fieldname] = splice_ai_fieldvalue

    if variant_line.get("gnomad") != None:
      for vat_gnomad_fieldname in vat_nirvana_gnomad_dictionary.keys():  # like "gnomad_all_af"
        nirvana_gnomad_fieldname = vat_nirvana_gnomad_dictionary.get(vat_gnomad_fieldname)
        gnomad_fieldvalue = variant_line.get("gnomad").get(nirvana_gnomad_fieldname)
        row[vat_gnomad_fieldname] = gnomad_fieldvalue
      gnomad_row = get_gnomad_subpop(variant_line["gnomad"])
      row.update(gnomad_row)

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
      ordered_significance_values = []
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

    gvs_annotations = variant_line["gvsAnnotations"]

    for vat_gvs_alleles_fieldname in vat_nirvana_gvs_alleles_dictionary.keys():  # like "gvs_all_ac"
      nirvana_gvs_alleles_fieldname = vat_nirvana_gvs_alleles_dictionary.get(vat_gvs_alleles_fieldname)
      gvs_alleles_fieldvalue = gvs_annotations.get(nirvana_gvs_alleles_fieldname)
      row[vat_gvs_alleles_fieldname] = gvs_alleles_fieldvalue
      row["gvs_all_sc"] = int(gvs_annotations.get("AC_Hom") / 2 ) + gvs_annotations.get('AC_Het')


    subpopulation_info = get_subpopulation_calculations(gvs_annotations)
    row.update(subpopulation_info)

    return row


def make_positions_json(annotated_json, output_json):
  output_file=gzip.open(output_json, 'w')

  if annotated_json.endswith(".gz"):
      json_data = gzip.open(annotated_json, 'rb')
  else:
      json_data = open(annotated_json, 'rb')

  positions = ijson.items(json_data, 'positions.item', use_float=True)

  for p in positions:
    position=p['position']
    refAllele=p['refAllele'] # ref_allele
    altAlleles=p['altAlleles'] # this is an array that we need to split into each alt_allele
    variants=p['variants']
    if '*' in altAlleles:
      altAlleles.remove('*')
    ## there should be the same number of altAlleles as variants, right?
    assert len(altAlleles)==len(variants)
    # row for each variant - transcript
    # so let's start with each variant
    for index, variant in enumerate(variants):
      # check if it's a row we care about
      if check_filtering(variant) is False:
        continue
      # remember that we want one for each variant-transcript and variant-null for variants without transcripts
      if variant.get("transcripts") == None:
        # then we make a special row
        row = make_annotated_json_row(position, refAllele, altAlleles[index], variant, None)
        json_str = json.dumps(row) + "\n"
        json_bytes = json_str.encode('utf-8')
        output_file.write(json_bytes)
      else:
        transcript_lines = variant.get("transcripts")
        # Collect all the transcript sources and check for if they contain Ensembl
        sources = [transcript.get('source') for transcript in transcript_lines]
        if "Ensembl" in sources:
          for transcript in transcript_lines:
            if transcript.get('source') == "Ensembl":
              row = make_annotated_json_row(position, refAllele, altAlleles[index], variant, transcript)
              json_str = json.dumps(row) + "\n"
              json_bytes = json_str.encode('utf-8')
              output_file.write(json_bytes)
        else:
          # if there are transcripts, but they are not Ensembl, we now only want one row in the VAT, not one row per transcript
          row = make_annotated_json_row(position, refAllele, altAlleles[index], variant, None)
          json_str = json.dumps(row) + "\n"
          json_bytes = json_str.encode('utf-8')
          output_file.write(json_bytes)
  output_file.close()
  json_data.close()

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
      omim_line = gene_line["omim"][0]
      if len(gene_line.get("omim")) > 1:
        print("WARNING: An assumption about the possible count of omim values is incorrect.", gene_line.get("name"),len(gene_line.get("omim")))
      row["gene_omim_id"] = omim_line.get("mimNumber")
      if omim_line.get("phenotypes") != None:
        phenotypes = omim_line["phenotypes"]
        for vat_omim_fieldname in omim_fieldnames:  # like "mimNumber", "phenotype"
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
  json_data.close()


def make_annotation_jsons(annotated_json, output_json, output_genes_json):
  make_positions_json(annotated_json, output_json)
  # we've already read the whole file once so we have to open it again
  make_genes_json(annotated_json, output_genes_json)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(allow_abbrev=False, description='Create BQ load friendly jsons for VAT creation')
  parser.add_argument('--annotated_json',type=str, help='nirvana created annotation json', required=True)
  parser.add_argument('--output_vt_json',type=str, help='name of the vt json', required=True)
  parser.add_argument('--output_genes_json',type=str, help='name of the genes json', required=True)

  args = parser.parse_args()

  make_annotation_jsons(args.annotated_json,
                        args.output_vt_json,
                        args.output_genes_json)
