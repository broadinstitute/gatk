package org.broadinstitute.hellbender.tools.variantdb.arrays;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;


/**
 * Expected headers for the  uncompressed array table
 *     sample, // req
 *     probe_id, // req
 *     GT_encoded,
 *     NORMX, // intensity
 *     NORMY, // intensity
 *     BAF // b allele fraction --> AD proxy
 *     LRR // Log R ratio --> intensity value instead of DP
 * 
 * Headers for the compressed array table
 *     basic_array_data
 *     raw_array_data
 */

public enum RawArrayFieldEnum {
    sample_id {
        public String getColumnValue(VariantContext variant, ProbeInfo probeInfo, String sampleId) {
            return sampleId;
        }
    },

    // This where the validation step (required vs not) lives  -- fail if there is missing data for a required field
    // and just leave it empty if not required
    basic_array_data {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo, String sampleId) {
            String gt = GT_encoded.getColumnValue(variant, probeInfo, sampleId);
            BasicArrayData.ArrayGenotype agt;
            if (".".equals(gt)) {
                agt = BasicArrayData.ArrayGenotype.NO_CALL;
            } else {
                agt = BasicArrayData.ArrayGenotype.valueOf(gt);
            }
            BasicArrayData d = new BasicArrayData(Integer.parseInt(sampleId), (int) probeInfo.probeId, agt);
            return String.valueOf(d.encode());
        }
    },

    raw_array_data {
        private Float convert(String s) {
            if (s == null || "".equals(s) || "null".equals(s) ) {
                return null;
            } else {
                return Float.parseFloat(s);                
            }
        }

        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo, String sampleId) {
            String normx = NORMX.getColumnValue(variant, probeInfo, sampleId);
            String normy = NORMY.getColumnValue(variant, probeInfo, sampleId);
            String baf = BAF.getColumnValue(variant, probeInfo, sampleId);
            String lrr = LRR.getColumnValue(variant, probeInfo, sampleId);

            RawArrayData d = new RawArrayData(convert(normx),
                                              convert(normy),
                                              convert(baf),
                                              convert(lrr)
                                              );
            return String.valueOf(d.encode());
        }
    },

    probe_id { // Required
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo, String sampleId) {
            return String.valueOf(probeInfo.probeId);
        }
    },

    GT_encoded {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo, String sampleId) {
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
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo, String sampleId) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("NORMX"));
        }
    },
    NORMY {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo, String sampleId) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("NORMY"));
        }
    },
    BAF {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo, String sampleId) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("BAF"));
        }
    },
    LRR {
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo, String sampleId) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("LRR"));
        }
    };

    public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo, String sampleId) {
        throw new IllegalArgumentException("Not implemented");
    }

    public static RawArrayFieldEnum[] getUncompressedRawArrayFieldEnums() {
        return new RawArrayFieldEnum[] { sample_id, probe_id, GT_encoded, NORMX, NORMY, BAF, LRR };
    }
    public static RawArrayFieldEnum[] getCompressedRawArrayFieldEnums() {
        return new RawArrayFieldEnum[] { basic_array_data, raw_array_data };
    }
}
