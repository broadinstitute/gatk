package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Expected headers for the Variant Table (VET)
 *     location, // req
 *     sample, // req
 *     reference_bases, // req
 *     alternate_bases_alt, // req
 *     alternate_bases_AS_RAW_MQ, // req
 *     alternate_bases_AS_RAW_MQRankSum,
 *     alternate_bases_AS_QUALapprox, // req
 *     alternate_bases_AS_RAW_ReadPosRankSum,
 *     alternate_bases_AS_SB_TABLE, // req
 *     alternate_bases_AS_VarDP, // req
 *     call_genotype, // req
 *     call_AD,
 *     call_DP, // Laura says consider removing for now-- so similar to AS_VarDP
 *     call_GQ, // req
 *     call_PGT,
 *     call_PID,
 *     call_PL // req
 *
 */
public enum ExomeFieldEnum {
    // This where the validation step (required vs not) lives  -- fail if there is missing data for a required field
    // and just leave it empty if not required

    location, // Required-- encoded chromosome and position
    sample, // Required-- sample Id for sample

    ref { // Required
        public String getColumnValue(final VariantContext variant) {
            final String referenceBase = variant.getReference().getBaseString();
            if (referenceBase == null) {
                throw new IllegalArgumentException("Cannot be missing required value for reference_bases"); // TODO, should this be UserException too?
            }
            return referenceBase;
        }
    },

    alt { // remove "<NON_REF>"
        //TODO what if this field is null and if <NON_REF> is not there--throw an error
        public String getColumnValue(final VariantContext variant) {
            List<String> outList = new ArrayList<>();
            for(Allele a : variant.getAlternateAlleles()) {
                if (!a.isNonRefAllele()) { // TODO unit test this
                    outList.add(a.getDisplayString());
                }
            }
            return String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, outList);
        }
    },

    AS_RAW_MQ {
        // Required
        // can strip off the first one?
        // TODO sci notation?
        public String getColumnValue(final VariantContext variant) {
            String out = getAttribute(variant, GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, null);
            if (out == null) {
                throw new UserException("Cannot be missing required value for alternate_bases.AS_RAW_MQ");
            }
            if (out.endsWith("|0.00")) {
                out = out.substring(0, out.length() - 5);
                String[] outValues = out.split("\\|");
                out = Arrays
                        .stream(outValues)
                        .map(val -> val.endsWith(".00") ? val.substring(0, val.length() - 3) : val)
                        .collect(Collectors.joining(VCFConstants.PHASED));
            } else {
                throw new UserException("Expected AS_RAW_MQ value to end in |0.00");
            }
            return out;
        }
    },

    AS_RAW_MQRankSum { // TODO -- maybe rely on 1/1 for call_GT, also get rid of the | at the beginning
        public String getColumnValue(final VariantContext variant) {
            String out =  getAttribute(variant, GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY, "");
            if (out.contentEquals("||") || out.contentEquals("|||") ) {
                out = " "; //  TODO is this better than null?
                return out;
            }
            if (out.startsWith("|")) {
                out = out.substring(1);
            } else {
                throw new UserException("Expected AS_RAW_MQRankSum value to begin with a |");
            }
            if (out.endsWith("|NaN")) {
                out = out.substring(0, out.length() - 4);
            } else {
                throw new UserException("Expected AS_RAW_MQRankSum value to be ||, ||| or to end in |NaN");
            }
            return out;

        }
    },

    AS_QUALapprox { // Required
        public String getColumnValue(final VariantContext variant) {
            //TODO find a constant for "AS_QUALapprox"
            String out = getAttribute(variant, "AS_QUALapprox", null);
            if (out == null) {
                throw new UserException("Cannot be missing required value for alternate_bases.AS_QUALapprox");
            }
            if (out.startsWith("|")) {
                out = out.substring(1);
            } else {
                throw new UserException("Expected AS_RAW_MQRankSum value to begin with a |");
            }
            // check to see if there are two or three values, make sure last is smallest, throw out last
            List<String> outList = Arrays.asList(out.split("\\|"));
            if (outList.size() == 2 | outList.size() == 3) { // check length of array -- needs to be 2 or 3
                if (outList.lastIndexOf(Collections.min(outList)) == outList.size() - 1) { // this should be the smallest value
                    out = StringUtils.join(outList.subList(0, outList.size() - 1), VCFConstants.PHASED);
                } else {
                    throw new UserException(String.format("Expected the final value of AS_QUALapprox to be the smallest at %d", variant.getStart()));
                }
            } else {
                throw new UserException("Expected AS_QUALapprox to have two or three values");
            }
            return out;
        }
    },

    AS_RAW_ReadPosRankSum {  // TODO -- maybe rely on 1/1 for call_GT
        public String getColumnValue(final VariantContext variant) {
            String out =  getAttribute(variant, GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY, "");
            if (out.contentEquals("||") || out.contentEquals("|||") ) {
                out = " "; // TODO is this better than null?
                return out;
            }
            if (out.startsWith("|")) {
                out = out.substring(1);
            } else {
                throw new UserException("Expected AS_RAW_ReadPosRankSum value to begin with a |");
            }
            if (out.endsWith("|NaN")) {
                out = out.substring(0, out.length() - 4);
            } else {
                throw new UserException("Expected AS_RAW_ReadPosRankSum value to be ||, ||| or to end in |NaN");
            }
            return out;
        }
    },

    AS_SB_TABLE { // Required
        public String getColumnValue(final VariantContext variant) {
            String out = getAttribute(variant, GATKVCFConstants.AS_SB_TABLE_KEY, null);
            if (out == null) {
                throw new UserException("Cannot be missing required value for alternate_bases.AS_SB_TABLE");
            }
            if (out.endsWith("|0,0")) {
                out = out.substring(0, out.length() - 4);
            } else {
                throw new UserException("Expected AS_SB_TABLE value to end in |0,0");
            }
            return out;
        }
    },

    AS_VarDP { // Required
        public String getColumnValue(final VariantContext variant) {
            //TODO find a constant for "AS_VarDP"
            String out = getAttribute(variant, "AS_VarDP", null);
            if (out == null) {
                throw new UserException("Cannot be missing required value for alternate_bases.AS_VarDP");
            }
            if (out.endsWith("|0")) {
                out = out.substring(0, out.length() - 2);
            } else {
                throw new UserException("Expected AS_VarDP value to end in |0");
            }
            return out;
        }
    },

    call_GT {
        public String getColumnValue(final VariantContext variant) {
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

    call_AD {
        public String getColumnValue(final VariantContext variant) {
            String out = variant.getGenotype(0).hasAD() ? Arrays.stream(variant.getGenotype(0).getAD())
                    .mapToObj(String::valueOf)
                    .collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)) : "";
            if (out.endsWith(",0")) {
                out = out.substring(0, out.length() - 2);
            } else {
                throw new UserException("Expected call_AD to have a final value of 0");
            }
            return out;
        }
    },

    // call_DP { // TODO we can drop whole column since it is similar to AS_VarDP
    // TODO come up with a check-- looks like we can drop this whole one!!!!!
    // public String getColumnValue(final VariantContext variant) {
    // return variant.getGenotype(0).hasDP() ? String.valueOf(variant.getGenotype(0).getDP()): "";
    // }
    // },

    call_GQ { // Required
        public String getColumnValue(final VariantContext variant) {
            if (!variant.getGenotype(0).hasGQ()) {
                throw new UserException("Cannot be missing required value for call.GQ");
            }
            return  String.valueOf(variant.getGenotype(0).getGQ());
        }
    },

    call_PGT {
        public String getColumnValue(final VariantContext variant) {
            return variant.getGenotype(0).hasAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) ? String.valueOf(variant.getGenotype(0).getAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY)) : "";
        }
    },

    call_PID {
        public String getColumnValue(final VariantContext variant) {
            return variant.getGenotype(0).hasAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY) ? String.valueOf(variant.getGenotype(0).getAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY)) : "";
        }
    },

    call_PL {
        public String getColumnValue(final VariantContext variant) {
            return variant.getGenotype(0).hasPL() ? Arrays.stream(variant.getGenotype(0).getPL())
                    .mapToObj(String::valueOf)
                    .collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)) : "";
        }
    };

    public String getColumnValue(final VariantContext variant) {
        throw new IllegalArgumentException("Not implemented");
    }

    private static String getAttribute(VariantContext vc, String key, String defaultValue){
        Object attr = vc.getAttribute(key);
        if ( attr == null ) return defaultValue;
        if ( attr instanceof String ) return (String)attr;
        if ( attr instanceof List) return StringUtils.join((List)attr, VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
        return String.valueOf(attr); // throws an exception if this isn't a string
    }
}

