package org.broadinstitute.hellbender.tools.variantdb.ingest.arrays;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import java.util.Set;

/**
 * Expected headers for the Variant Table (VET)
 *     position, // req
 *     sample, // req
 *     ref, // req
 *     alt, // req
 *     call_GT, // req
 *     call_NORMX, // req
 *     call_NORMY, // req
 *     call_BAF // b allele fraction --> AD proxy
 *     call_LRR // Log R ratio --> intensity value instead of DP
 *     call_IGC // IGC < --- instead of the GQ
 *
 */

public enum RawArrayFieldEnum {
    // This where the validation step (required vs not) lives  -- fail if there is missing data for a required field
    // and just leave it empty if not required

    sample_id, // Required-- sample Id for sample


    rsid { // Required
        public String getColumnValue(final VariantContext variant) {
            final String siteId = variant.getID();
            if (siteId == null) {
                throw new IllegalArgumentException("Cannot be missing required value for site_name"); // TODO, should this be UserException too?
            }
            return siteId;
        }
    },

    filter {
        public String getColumnValue(final VariantContext variant) {
            Set<String> outList = variant.getFilters();
            return outList.isEmpty() ? "null" : String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, outList);
        }
    },

    call_GT_encoded {
        public String getColumnValue(final VariantContext variant) {
            Genotype g = variant.getGenotype(0);
            RawArrayTsvCreator.GT_encoding gt = RawArrayTsvCreator.GT_encoding.MISSING;
            if (g.isHomRef()) {
                gt = RawArrayTsvCreator.GT_encoding.HOM_REF;
            } else if (g.isHomVar()) {
                gt = RawArrayTsvCreator.GT_encoding.HOM_VAR;
            } else if (g.isHetNonRef()) {
                gt = RawArrayTsvCreator.GT_encoding.HET_NON_REF;
            } else if (g.isHet()) {
                gt = RawArrayTsvCreator.GT_encoding.HET;
            }
            return gt.getValue();
        }
    },

    call_NORMX { // variant.getGenotypes().get(0).getExtendedAttributes().get("NORMX")
        public String getColumnValue(final VariantContext variant) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("NORMX"));
        }
    },
    call_NORMY { // variant.getGenotypes().get(0).getExtendedAttributes().get("NORMY")
        public String getColumnValue(final VariantContext variant) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("NORMY"));
        }
    },
    call_BAF { // variant.getGenotypes().get(0).getExtendedAttributes().get("BAF")
        public String getColumnValue(final VariantContext variant) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("BAF"));
        }
    },
    call_LRR { // variant.getGenotypes().get(0).getExtendedAttributes().get("LRR")
        public String getColumnValue(final VariantContext variant) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("LRR"));
        }
    },
    call_IGC { // Required instead of GQ variant.getGenotypes().get(0).getExtendedAttributes().get("IGC")
        public String getColumnValue(final VariantContext variant) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("IGC"));
        }
    };

    public String getColumnValue(final VariantContext variant) {
        throw new IllegalArgumentException("Not implemented");
    }

//    private static String getAttribute(VariantContext vc, String key, String defaultValue){
//        Object attr = vc.getAttribute(key);
//        if ( attr == null ) return defaultValue;
//        if ( attr instanceof String ) return (String)attr;
//        if ( attr instanceof List) return StringUtils.join((List)attr, VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
//        return String.valueOf(attr); // throws an exception if this isn't a string
//    }
}
