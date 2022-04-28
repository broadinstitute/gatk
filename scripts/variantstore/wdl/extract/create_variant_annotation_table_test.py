# -*- coding: utf-8 -*-
import uuid
import time
from datetime import datetime

import csv
import json
#import ijson
import gzip
import argparse
import unittest

from create_variant_annotation_table import make_annotated_json_row

class TestMakeAnnotatedJsonRow(unittest.TestCase):
    def test_make_annotated_json_row_success(self):
        actual = make_annotated_json_row(
          row_position=117559590,
          row_ref="ATCT",
          row_alt="A",
          variant_line=variantLine,
          transcript_line=None)
        expected = {"tank_a": ["shark", "tuna"]}
        self.maxDiff=None
        self.assertEqual(actual, expected)

variantLine = {
          "vid": "7-117559590-ATCT-A",
          "chromosome": "chr7",
          "begin": 117559591,
          "end": 117559593,
          "refAllele": "TCT",
          "altAllele": "-",
          "variantType": "deletion",
          "hgvsg": "NC_000007.14:g.117559592_117559594del",
          "regulatoryRegions": [
            {
              "id": "ENSR00000217281",
              "type": "open_chromatin_region",
              "consequence": [
                "regulatory_region_variant"
              ]
            }
          ],
          "clinvar": [
            {
              "id": "VCV000007105.11",
              "reviewStatus": "practice guideline",
              "significance": [
                "pathogenic"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            },
            {
              "id": "VCV000007105.11",
              "reviewStatus": "practice guideline",
              "significance": [
                "pathogenic"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            },
            {
              "id": "RCV000007524.9",
              "variationId": 7105,
              "reviewStatus": "no assertion criteria provided",
              "alleleOrigins": [
                "germline"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "phenotypes": [
                "Bronchiectasis with or without elevated sweat chloride 1, modifier of"
              ],
              "medGenIds": [
                "CN258830"
              ],
              "omimIds": [
                "602421.0001"
              ],
              "significance": [
                "risk factor"
              ],
              "lastUpdatedDate": "2020-03-01",
              "pubMedIds": [
                "1370875",
                "1377276",
                "1380673",
                "1384321",
                "1536179",
                "1673094",
                "1756602",
                "1997384",
                "2210767",
                "2220803",
                "2236053",
                "2300168",
                "7537148",
                "9439669",
                "15367919",
                "17692578",
                "18507830",
                "20595578",
                "20705837",
                "22981120",
                "25981758",
                "26618866"
              ],
              "isAlleleSpecific": True
            },
            {
              "id": "VCV000007105.11",
              "reviewStatus": "practice guideline",
              "significance": [
                "pathogenic"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            },
            {
              "id": "RCV000417138.1",
              "variationId": 7105,
              "reviewStatus": "reviewed by expert panel",
              "alleleOrigins": [
                "germline"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "phenotypes": [
                "ivacaftor / lumacaftor response - Efficacy"
              ],
              "medGenIds": [
                "CN240582"
              ],
              "omimIds": [
                "602421.0001"
              ],
              "significance": [
                "drug response"
              ],
              "lastUpdatedDate": "2020-03-01",
              "pubMedIds": [
                "21825083",
                "21976485",
                "22992668",
                "24973281",
                "25981758",
                "27214033",
                "27298017",
                "27334259",
                "27805836",
                "27898234",
                "28325531",
                "28606620",
                "29126871",
                "29327948",
                "29451946"
              ],
              "isAlleleSpecific": True
            },
            {
              "id": "VCV000007105.11",
              "reviewStatus": "practice guideline",
              "significance": [
                "pathogenic"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            },
            {
              "id": "RCV000058929.15",
              "variationId": 7105,
              "reviewStatus": "criteria provided, multiple submitters, no conflicts",
              "alleleOrigins": [
                "germline"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "phenotypes": [
                "not provided"
              ],
              "medGenIds": [
                "CN517202"
              ],
              "omimIds": [
                "602421.0001"
              ],
              "significance": [
                "pathogenic"
              ],
              "lastUpdatedDate": "2020-03-01",
              "pubMedIds": [
                "25741868",
                "26467025"
              ],
              "isAlleleSpecific": True
            },
            {
              "id": "VCV000634837.1",
              "reviewStatus": "reviewed by expert panel",
              "significance": [
                "pathogenic"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            },
            {
              "id": "RCV000785641.1",
              "variationId": 634837,
              "reviewStatus": "reviewed by expert panel",
              "alleleOrigins": [
                "germline"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "phenotypes": [
                "Cystic fibrosis"
              ],
              "medGenIds": [
                "C0010674"
              ],
              "omimIds": [
                "219700",
                "602421.0001"
              ],
              "orphanetIds": [
                "586"
              ],
              "significance": [
                "pathogenic"
              ],
              "lastUpdatedDate": "2020-03-01",
              "pubMedIds": [
                "23974870"
              ],
              "isAlleleSpecific": True
            },
            {
              "id": "VCV000007105.11",
              "reviewStatus": "practice guideline",
              "significance": [
                "pathogenic"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            },
            {
              "id": "RCV000007523.21",
              "variationId": 7105,
              "reviewStatus": "practice guideline",
              "alleleOrigins": [
                "germline",
                "maternal",
                "unknown",
                "paternal",
                "inherited"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "phenotypes": [
                "Cystic fibrosis"
              ],
              "medGenIds": [
                "C0010674"
              ],
              "omimIds": [
                "219700",
                "602421.0001"
              ],
              "orphanetIds": [
                "586"
              ],
              "significance": [
                "pathogenic"
              ],
              "lastUpdatedDate": "2020-03-01",
              "pubMedIds": [
                "1370875",
                "1377276",
                "1380673",
                "1381146",
                "1384321",
                "1536179",
                "1673094",
                "1756602",
                "1997384",
                "2210767",
                "2220803",
                "2236053",
                "2300168",
                "2475911",
                "2570460",
                "7537148",
                "9135274",
                "9439669",
                "10782933",
                "11280952",
                "15141088",
                "15246977",
                "15367919",
                "15371902",
                "17206681",
                "17692578",
                "18507830",
                "20595578",
                "20705837",
                "22981120",
                "23974870",
                "24033266",
                "24559724",
                "25741868",
                "25981758",
                "26618866",
                "28492532"
              ],
              "isAlleleSpecific": True
            },
            {
              "id": "RCV001000022.1",
              "variationId": 7105,
              "reviewStatus": "criteria provided, single submitter",
              "alleleOrigins": [
                "germline"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "phenotypes": [
                "not specified"
              ],
              "medGenIds": [
                "CN169374"
              ],
              "omimIds": [
                "602421.0001"
              ],
              "significance": [
                "pathogenic"
              ],
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            },
            {
              "id": "VCV000007105.11",
              "reviewStatus": "practice guideline",
              "significance": [
                "pathogenic"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            },
            {
              "id": "VCV000007105.11",
              "reviewStatus": "practice guideline",
              "significance": [
                "pathogenic"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            },
            {
              "id": "RCV000119038.4",
              "variationId": 7105,
              "reviewStatus": "no assertion criteria provided",
              "alleleOrigins": [
                "germline"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "phenotypes": [
                "Hereditary pancreatitis"
              ],
              "medGenIds": [
                "C0238339"
              ],
              "omimIds": [
                "167800",
                "602421.0001"
              ],
              "orphanetIds": [
                "676"
              ],
              "significance": [
                "pathogenic"
              ],
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            },
            {
              "id": "RCV000626693.1",
              "variationId": 7105,
              "reviewStatus": "criteria provided, single submitter",
              "alleleOrigins": [
                "unknown"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "phenotypes": [
                "Recurrent pancreatitis"
              ],
              "medGenIds": [
                "C4551632"
              ],
              "omimIds": [
                "602421.0001"
              ],
              "significance": [
                "pathogenic"
              ],
              "lastUpdatedDate": "2020-03-01",
              "pubMedIds": [
                "25741868"
              ],
              "isAlleleSpecific": True
            },
            {
              "id": "RCV000626692.1",
              "variationId": 7105,
              "reviewStatus": "criteria provided, single submitter",
              "alleleOrigins": [
                "unknown"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "phenotypes": [
                "Duodenal stenosis"
              ],
              "medGenIds": [
                "C0238093"
              ],
              "omimIds": [
                "602421.0001"
              ],
              "significance": [
                "likely pathogenic"
              ],
              "lastUpdatedDate": "2020-03-01",
              "pubMedIds": [
                "25741868"
              ],
              "isAlleleSpecific": True
            },
            {
              "id": "VCV000007105.11",
              "reviewStatus": "practice guideline",
              "significance": [
                "pathogenic"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            },
            {
              "id": "VCV000007105.11",
              "reviewStatus": "practice guideline",
              "significance": [
                "pathogenic"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            },
            {
              "id": "RCV001004459.1",
              "variationId": 7105,
              "reviewStatus": "criteria provided, single submitter",
              "alleleOrigins": [
                "germline"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "phenotypes": [
                "Cystic fibrosis",
                "Congenital bilateral aplasia of vas deferens from CFTR mutation"
              ],
              "medGenIds": [
                "C0010674",
                "C0403814"
              ],
              "omimIds": [
                "219700",
                "277180",
                "602421.0001"
              ],
              "orphanetIds": [
                "586",
                "48"
              ],
              "significance": [
                "pathogenic"
              ],
              "lastUpdatedDate": "2020-03-01",
              "pubMedIds": [
                "25741868"
              ],
              "isAlleleSpecific": True
            },
            {
              "id": "VCV000007105.11",
              "reviewStatus": "practice guideline",
              "significance": [
                "pathogenic"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            },
            {
              "id": "RCV000624683.1",
              "variationId": 7105,
              "reviewStatus": "criteria provided, single submitter",
              "alleleOrigins": [
                "germline"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "phenotypes": [
                "Inborn genetic diseases"
              ],
              "medGenIds": [
                "C0950123"
              ],
              "omimIds": [
                "602421.0001"
              ],
              "significance": [
                "pathogenic"
              ],
              "lastUpdatedDate": "2020-03-01",
              "pubMedIds": [
                "2378364",
                "11737931",
                "12833420",
                "12940920",
                "14685937",
                "15221447",
                "15482777",
                "15619636",
                "18463704",
                "18639722",
                "21152102",
                "21658649",
                "23457292",
                "23951356",
                "23974870",
                "26631874",
                "28129809",
                "28129811"
              ],
              "isAlleleSpecific": True
            },
            {
              "id": "RCV000211188.1",
              "variationId": 7105,
              "reviewStatus": "reviewed by expert panel",
              "alleleOrigins": [
                "germline"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "phenotypes": [
                "ivacaftor response - Efficacy"
              ],
              "medGenIds": [
                "CN236562"
              ],
              "omimIds": [
                "602421.0001"
              ],
              "significance": [
                "drug response"
              ],
              "lastUpdatedDate": "2020-03-01",
              "pubMedIds": [
                "19846789",
                "21602569",
                "22293084",
                "22383668",
                "22942289",
                "22992668",
                "23891399",
                "25148434",
                "26968770",
                "29099333"
              ],
              "isAlleleSpecific": True
            },
            {
              "id": "VCV000007105.11",
              "reviewStatus": "practice guideline",
              "significance": [
                "pathogenic"
              ],
              "refAllele": "TCT",
              "altAllele": "-",
              "lastUpdatedDate": "2020-03-01",
              "isAlleleSpecific": True
            }
          ],
          "dbsnp": [
            "rs113993960"
          ],
          "gnomad": {
            "coverage": 36,
            "allAf": 0.007172,
            "allAn": 282630,
            "allAc": 2027,
            "allHc": 1,
            "afrAf": 0.002604,
            "afrAn": 24958,
            "afrAc": 65,
            "afrHc": 0,
            "amrAf": 0.003811,
            "amrAn": 35426,
            "amrAc": 135,
            "amrHc": 0,
            "easAf": 0,
            "easAn": 19948,
            "easAc": 0,
            "easHc": 0,
            "finAf": 0.002433,
            "finAn": 25074,
            "finAc": 61,
            "finHc": 0,
            "nfeAf": 0.012384,
            "nfeAn": 129034,
            "nfeAc": 1598,
            "nfeHc": 1,
            "asjAf": 0.005594,
            "asjAn": 10368,
            "asjAc": 58,
            "asjHc": 0,
            "sasAf": 0.001993,
            "sasAn": 30608,
            "sasAc": 61,
            "sasHc": 0,
            "othAf": 0.006792,
            "othAn": 7214,
            "othAc": 49,
            "othHc": 0,
            "maleAf": 0.007099,
            "maleAn": 153256,
            "maleAc": 1088,
            "maleHc": 1,
            "femaleAf": 0.007258,
            "femaleAn": 129374,
            "femaleAc": 939,
            "femaleHc": 0,
            "controlsAllAf": 0.006306,
            "controlsAllAn": 120200,
            "controlsAllAc": 758
          },
          "oneKg": {
            "allAf": 0.003994,
            "afrAf": 0,
            "amrAf": 0.010086,
            "easAf": 0,
            "eurAf": 0.008946,
            "sasAf": 0.00409,
            "allAn": 5008,
            "afrAn": 1322,
            "amrAn": 694,
            "easAn": 1008,
            "eurAn": 1006,
            "sasAn": 978,
            "allAc": 20,
            "afrAc": 0,
            "amrAc": 7,
            "easAc": 0,
            "eurAc": 9,
            "sasAc": 4
          },
          "topmed": {
            "allAf": 0.008521,
            "allAn": 125568,
            "allAc": 1070,
            "allHc": 1
          },
          "transcripts": [
            {
              "transcript": "ENST00000003084.10",
              "source": "Ensembl",
              "bioType": "protein_coding",
              "codons": "aTCTtt/att",
              "aminoAcids": "F/-",
              "cdnaPos": "1652-1654",
              "cdsPos": "1520-1522",
              "exons": "11/27",
              "proteinPos": "508",
              "geneId": "ENSG00000001626",
              "hgnc": "CFTR",
              "consequence": [
                "inframe_deletion"
              ],
              "hgvsc": "ENST00000003084.10:c.1521_1523del",
              "hgvsp": "ENSP00000003084.6:p.(Phe508del)",
              "isCanonical": True,
              "proteinId": "ENSP00000003084.6"
            },
            {
              "transcript": "NM_000492.3",
              "source": "RefSeq",
              "bioType": "protein_coding",
              "codons": "aTCTtt/att",
              "aminoAcids": "F/-",
              "cdnaPos": "1652-1654",
              "cdsPos": "1520-1522",
              "exons": "11/27",
              "proteinPos": "508",
              "geneId": "1080",
              "hgnc": "CFTR",
              "consequence": [
                "inframe_deletion"
              ],
              "hgvsc": "NM_000492.3:c.1521_1523del",
              "hgvsp": "NP_000483.3:p.(Phe508del)",
              "isCanonical": True,
              "proteinId": "NP_000483.3"
            },
            {
              "transcript": "ENST00000426809.5",
              "source": "Ensembl",
              "bioType": "protein_coding",
              "codons": "aTCTtt/att",
              "aminoAcids": "F/-",
              "cdnaPos": "1430-1432",
              "cdsPos": "1430-1432",
              "exons": "10/26",
              "proteinPos": "478",
              "geneId": "ENSG00000001626",
              "hgnc": "CFTR",
              "consequence": [
                "inframe_deletion"
              ],
              "hgvsc": "ENST00000426809.5:c.1431_1433del",
              "hgvsp": "ENSP00000389119.1:p.(Phe478del)",
              "proteinId": "ENSP00000389119.1"
            },
            {
              "transcript": "ENST00000472848.1",
              "source": "Ensembl",
              "bioType": "processed_transcript",
              "geneId": "ENSG00000001626",
              "hgnc": "CFTR",
              "consequence": [
                "upstream_gene_variant"
              ]
            },
            {
              "transcript": "ENST00000441019.1",
              "source": "Ensembl",
              "bioType": "antisense_RNA",
              "geneId": "ENSG00000232661",
              "hgnc": "CFTR-AS1",
              "consequence": [
                "downstream_gene_variant"
              ],
              "isCanonical": True
            }
          ],
          "gvsAnnotations":{
            "AC_Hom": 100, "AC_Het": 100, "AC_Hemi": 100,
            "AC_afr":10, "AN_afr": 10, "AF_afr": 10, "AC_Hom_afr": 10, "AC_Het_afr": 10, "AC_Hemi_afr":10,
            "AC_amr":10, "AN_amr": 10, "AF_amr": 10, "AC_Hom_amr": 10, "AC_Het_amr": 10, "AC_Hemi_amr":10,
            "AC_eas":10, "AN_eas": 10, "AF_eas": 10, "AC_Hom_eas": 10, "AC_Het_eas": 10, "AC_Hemi_eas":10,
            "AC_eur":10, "AN_eur": 10, "AF_eur": 10, "AC_Hom_eur": 10, "AC_Het_eur": 10, "AC_Hemi_eur":10,
            "AC_mid":10, "AN_mid": 10, "AF_mid": 10, "AC_Hom_mid": 10, "AC_Het_mid": 10, "AC_Hemi_mid":10,
            "AC_oth":10, "AN_oth": 10, "AF_oth": 10, "AC_Hom_oth": 10, "AC_Het_oth": 10, "AC_Hemi_oth":10,
            "AC_sas":10, "AN_sas": 10, "AF_sas": 10, "AC_Hom_sas": 10, "AC_Het_sas": 10, "AC_Hemi_sas":10}
        }



