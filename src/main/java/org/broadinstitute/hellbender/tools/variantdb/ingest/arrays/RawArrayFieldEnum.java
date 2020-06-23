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
 *     GT_encoded,
 *     NORMX, // intensity
 *     NORMY, // intensity
 *     BAF // b allele fraction --> AD proxy
 *     LRR // Log R ratio --> intensity value instead of DP
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

    GT_encoded {
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

    NORMX {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("NORMX"));
        }
    },
    NORMY {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("NORMY"));
        }
    },
    BAF {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("BAF"));
        }
    },
    LRR {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("LRR"));
        }
    };

    public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo) {
        throw new IllegalArgumentException("Not implemented");
    }
}
