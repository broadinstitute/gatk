package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.CommonCode;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Expected headers for the Variant Table (VET)
 *     sample_id, // req
 *     location, // req
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
public enum VetFieldEnum {
    // This where the validation step (required vs not) lives  -- fail if there is missing data for a required field
    // and just leave it empty if not required

    sample_id, // Required-- sample Id for sample
    location, // Required-- encoded chromosome and position

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
            String outNotAlleleSpecific = getAttribute(variant, GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, null);
            String outNotAlleleSpecificAndOld = getAttribute(variant, GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED, null);
            if (out == null && outNotAlleleSpecific == null && outNotAlleleSpecificAndOld == null) {
                throw new UserException("Cannot be missing required value for alternate_bases.AS_RAW_MQ, RAW_MQandDP or RAW_MQ.");
            }
            if (out != null) {
                if(!out.endsWith("|0.00")) {
                    logger.warn("Expected AS_RAW_MQ value to end in |0.00. value is: " + out + " for variant " + variant.toString());
                }
                out = out.substring(0, out.lastIndexOf("|"));
                String[] outValues = out.split("\\|");
                out = Arrays
                        .stream(outValues)
                        .map(val -> val.endsWith(".00") ? val.substring(0, val.length() - 3) : val)
                        .collect(Collectors.joining(VCFConstants.PHASED));
                return out;
            // If we have gvcfs that are not allele specific from GATK4 we'll get RAW_MQandDP.
            // We can drop DP here and use AS_VarDP when finalizing RMS Mapping Quality
            } else {
                String outValue;
                if (outNotAlleleSpecific != null) {
                    String[] outValues = outNotAlleleSpecific.split(",");
                    if (outValues.length !=2) {
                        throw new UserException("Expected RAW_MQandDP to be two values separated by a comma.");
                    }
                    // First value is MQ the second is DP. Use the only MQ value we have for all alleles since we're faking allele specific annotations.
                    outValue = outValues[0];
                } else {
                    outValue = outNotAlleleSpecificAndOld;
                }
                // Spread MQ accross multiple alleles.
                //TODO: check how hail does this and change this to do the same thing
                double mq = Double.parseDouble(outValue);
                if (variant.getAlleles().size() == 3) {
                    outNotAlleleSpecific = (int) mq / 2 + VCFConstants.PHASED + (int) mq / 2;
                } else if (variant.getAlleles().size() == 4) {
                    outNotAlleleSpecific = (int) mq / 3 + VCFConstants.PHASED + (int) mq / 3 + VCFConstants.PHASED + (int) mq / 3;
                } else {
                    throw new UserException("Expected diploid sample to either have 3 alleles (ref, alt, non-ref) or 4 alleles (ref, alt 1, alt 2, non-ref)");
                }
                return outNotAlleleSpecific;
            }
        }
    },

    AS_RAW_MQRankSum { // TODO -- maybe rely on 1/1 for call_GT, also get rid of the | at the beginning
        public String getColumnValue(final VariantContext variant) {
            String out =  getAttribute(variant, GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY, null);
            if ( out == null || out.contentEquals("||") || out.contentEquals("|||") ) {
                out = "";
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
                // for now remove the last value even if not NaN
                out = out.substring(0, out.lastIndexOf("|"));
                //throw new UserException("Expected AS_RAW_MQRankSum value to be ||, ||| or to end in |NaN");
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
            String out =  getAttribute(variant, GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY, null);
            if (out == null || out.contentEquals("||") || out.contentEquals("|||") ) {
                out = "";
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
                // for now remove the last value even if not NaN
                out = out.substring(0, out.lastIndexOf("|"));
                //throw new UserException("Expected AS_RAW_ReadPosRankSum value to be ||, ||| or to end in |NaN");
            }
            return out;
        }
    },

    AS_SB_TABLE { // Required
        public String getColumnValue(final VariantContext variant) {
            String out = getAttribute(variant, GATKVCFConstants.AS_SB_TABLE_KEY, null);
            String outNotAlleleSpecific = variant.getGenotype(0).getExtendedAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, null).toString();
            if (out == null) {
                String[] outValues = outNotAlleleSpecific.split(",");
                if (variant.getAlleles().size() == 3) {
                    outNotAlleleSpecific = outValues[0] + "," + outValues[1] + "|" + outValues[2] + "," + outValues[3];
                } else if (variant.getAlleles().size() == 4) {
                    int sbPosSpread = Integer.parseInt(outValues[2]) / 2;
                    int sbNegSpread = Integer.parseInt(outValues[3]) / 2;
                    outNotAlleleSpecific = outValues[0] + "," + outValues[1] + "|" + sbPosSpread + "," + sbNegSpread + "|" + sbPosSpread + "," + sbNegSpread;
                } else {
                    throw new UserException("Expected diploid sample to either have 3 alleles (ref, alt, non-ref) or 4 alleles (ref, alt 1, alt 2, non-ref)");
                }
                return outNotAlleleSpecific;
            }
            if (out.endsWith("|0,0")) {
                out = out.substring(0, out.length() - 4);
            } else {
                // for now remove the last value even if not NaN
                out = out.substring(0, out.lastIndexOf("|"));
                //throw new UserException("Expected AS_SB_TABLE value to end in |0,0");
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
                // for now remove the last value even if not NaN
                out = out.substring(0, out.lastIndexOf("|"));
                //throw new UserException("Expected AS_VarDP value to end in |0");
            }
            return out;
        }
    },

    call_GT {
        public String getColumnValue(final VariantContext variant) {
            // TODO how is missing handled?
            return CommonCode.getGTString(variant);
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
                // for now remove the last value even if not NaN
                out = out.substring(0, out.lastIndexOf("|"));
                //throw new UserException("Expected call_AD to have a final value of 0");
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

    static final Logger logger = LogManager.getLogger(CreateVariantIngestFiles.class);

    private static String getAttribute(VariantContext vc, String key, String defaultValue){
        Object attr = vc.getAttribute(key);
        if ( attr == null ) return defaultValue;
        if ( attr instanceof String ) return (String)attr;
        if ( attr instanceof List) return StringUtils.join((List)attr, VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
        return String.valueOf(attr); // throws an exception if this isn't a string
    }
}

