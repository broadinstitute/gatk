package org.broadinstitute.hellbender.tools.variantdb.arrays;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;

import java.util.ArrayList;
import java.util.Arrays;


/**
 * Expected headers for the  uncompressed array table
 *     sample_id, // req
 *     GT
 *     DS
 *     GP1
 *     GP2
 *     GP3
 */

public enum ImputedFieldEnum {
    sample_id {
        public String getColumnValue(VariantContext variant, String sampleId) {
            return sampleId;
        }
    },

    GT {
        public String getColumnValue(final VariantContext variant, String sampleId) {
            IndexedAlleleList<Allele> alleleList = new IndexedAlleleList<>(variant.getAlleles());
            ArrayList<Integer> allele_indices = new ArrayList<Integer>();

            for (Allele allele : variant.getGenotype(0).getAlleles()) {
                allele_indices.add(alleleList.indexOfAllele(allele));
            }
            if (allele_indices.size() != 2){
                throw new IllegalArgumentException("GT doesnt have two alleles");
            }
            String separator = variant.getGenotype(0).isPhased() ? VCFConstants.PHASED : VCFConstants.UNPHASED;
            return StringUtils.join(allele_indices, separator);
        }
    },

    DS {
        public String getColumnValue(final VariantContext variant, String sampleId) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("DS"));
        }
    },
    GP1 {
        public String getColumnValue(final VariantContext variant, String sampleId) {
            return getIthElementFromAttributeList(variant.getGenotype(0).getExtendedAttribute("GP").toString(), 1);
        }
    },
    GP2 {
        public String getColumnValue(final VariantContext variant, String sampleId) {
            return getIthElementFromAttributeList(variant.getGenotype(0).getExtendedAttribute("GP").toString(), 2);
        }
    },
    GP3 {
        public String getColumnValue(final VariantContext variant, String sampleId) {
            return getIthElementFromAttributeList(variant.getGenotype(0).getExtendedAttribute("GP").toString(), 3);
        }
    };

    private static String getIthElementFromAttributeList(String GPStringValue, int i) {
        String response = null;
        if (GPStringValue != null && !GPStringValue.isEmpty()) {
            String[] vals = GPStringValue.split(",");
            if (vals.length >= i) {
                response = vals[i-1];
            }
        }
        return response;
    }

    public String getColumnValue(final VariantContext variant, String sampleId) {
        throw new IllegalArgumentException("Not implemented");
    }

    public static ImputedFieldEnum[] getUncompressedRawArrayFieldEnums() {
        return new ImputedFieldEnum[] { sample_id, GT, DS, GP1, GP2, GP3 };
    }
//    public static ImputedArrayFieldEnum[] getCompressedRawArrayFieldEnums() {
//        return new ImputedArrayFieldEnum[] { basic_array_data, raw_array_data };
//    }
}
