import json
import ijson
import gzip
import argparse
import logging
import sys

vat_nirvana_omim_dictionary = {
    "omim_phenotypes_id": "mimNumber",      # nullable
    "omim_phenotypes_name": "phenotype"     # nullable
}


def make_genes_json(annotated_json, output_genes_json):
    output_genes_file=gzip.open(output_genes_json, 'w')

    if annotated_json.endswith(".gz"):
        json_data = gzip.open(annotated_json, 'rb')
    else:
        json_data = open(annotated_json, 'rb')

    logging.info(f"Loading the genes json data")
    genes = ijson.items(json_data, 'item', use_float=True)
    logging.info(f"Done loading the genes json data")

    gene_count = 0
    for gene_line in genes:
        gene_count += 1
        logging.debug(f"gene_line: {gene_line}")
        if gene_line.get("omim") != None:
            row = {}
            row["gene_symbol"] = gene_line.get("name")
            omim_line = gene_line["omim"][0]
            if len(gene_line.get("omim")) > 1:
                logging.warning("WARNING: An assumption about the possible count of omim values is incorrect.", gene_line.get("name"),len(gene_line.get("omim")))
            row["gene_omim_id"] = omim_line.get("mimNumber")
            if omim_line.get("phenotypes") != None:
                phenotypes = omim_line["phenotypes"]
                for vat_omim_fieldname in vat_nirvana_omim_dictionary.keys():  # like "mimNumber", "phenotype"
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

    if gene_count == 0:
        logging.warning(f"WARNING: Found no items in annotated json file: {annotated_json}")
        sys.exit(0)

def make_annotation_json(annotated_json, output_genes_json, logging):
    logging.info("Making the genes json")
    make_genes_json(annotated_json, output_genes_json)
    logging.info("Done")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Create BQ load friendly json for VAT genes table creation')
    parser.add_argument('--annotated_json', type=str, help='nirvana created annotation json', required=True)
    parser.add_argument('--output_genes_json', type=str, help='name of the genes json', required=True)

    args = parser.parse_args()

    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')

    make_annotation_json(args.annotated_json,
                         args.output_genes_json,
                         logging)
