package org.broadinstitute.hellbender.tools.variantdb.ingest.arrays;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;

import java.util.Set;

/**
 * Expected headers for the Variant Table (VET)
 *     sample, // req
 *     probe_id, // req
 *     filter,
 *     call_GT_encoded,
 *     call_NORMX, // intensity
 *     call_NORMY, // intensity
 *     call_BAF // b allele fraction --> AD proxy
 *     call_LRR // Log R ratio --> intensity value instead of DP
 *
 */

public enum RawArrayFieldEnum {
    // This where the validation step (required vs not) lives  -- fail if there is missing data for a required field
    // and just leave it empty if not required

    sample_id, // Required-- sample Id for sample


    probe_id { // Required
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo) {
            return String.valueOf(probeInfo.probeId);
        }
    },

    filter {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo) {
            Set<String> outList = variant.getFilters();
            return outList.isEmpty() ? "null" : String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, outList);
        }
    },

    call_GT_encoded {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo) {
            Genotype g = variant.getGenotype(0);
            RawArrayTsvCreator.GT_encoding gt = RawArrayTsvCreator.GT_encoding.MISSING;
            if (g.isHomRef() || g.isHomVar()) {
                Allele allele = g.getAllele(0);
                if (allele.basesMatch(probeInfo.alleleA)) {
                    gt = RawArrayTsvCreator.GT_encoding.AA;
                } else if (allele.basesMatch(probeInfo.alleleB)) {
                    gt = RawArrayTsvCreator.GT_encoding.BB;
                } else {
                    throw new IllegalStateException("allele: " + allele + " must match either A: " + probeInfo.alleleA + " or B: " + probeInfo.alleleB);
                }
            } else if (g.isHet()) {
                gt = RawArrayTsvCreator.GT_encoding.AB;
            }
            return gt.getValue();
        }
    },

    call_NORMX {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("NORMX"));
        }
    },
    call_NORMY {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("NORMY"));
        }
    },
    call_BAF {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("BAF"));
        }
    },
    call_LRR {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("LRR"));
        }
    };

    public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo) {
        throw new IllegalArgumentException("Not implemented");
    }
}
