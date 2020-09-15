package org.broadinstitute.hellbender.tools.variantdb.arrays;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.variantdb.CommonCode;
import org.broadinstitute.hellbender.tools.variantdb.arrays.tables.ProbeInfo;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


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

    probe_id { // Required
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo, String sampleId) {
            return String.valueOf(probeInfo.probeId);
        }
    },

    GT_encoded { // Required
        public String getColumnValue(final VariantContext variant, ProbeInfo probeInfo, String sampleId) {
            List<Integer> alleleIndexes = CommonCode.getGTAlleleIndexes(variant);

            RawArrayTsvCreator.GT_encoding gt = RawArrayTsvCreator.GT_encoding.MISSING;
            if (alleleIndexes.size() == 2) {
                Set<Integer> uniqueAlleleIndexes = new HashSet<>(alleleIndexes);

                if (uniqueAlleleIndexes.size() == 1) {
                    // we know it's HOM something
                    if (uniqueAlleleIndexes.contains(0)) {
                        gt = RawArrayTsvCreator.GT_encoding.HOM_REF;
                    } else if (uniqueAlleleIndexes.contains(1)) {
                        gt = RawArrayTsvCreator.GT_encoding.HOM_VAR;
                    } else if (uniqueAlleleIndexes.contains(2)) {
                        gt = RawArrayTsvCreator.GT_encoding.HOM_ALT2;
                    }
                } else {
                    // we know its het
                    if (uniqueAlleleIndexes.containsAll(new HashSet<>(Arrays.asList(0, 1)))) {
                        gt = RawArrayTsvCreator.GT_encoding.HET0_1;
                    } else if (uniqueAlleleIndexes.containsAll(new HashSet<>(Arrays.asList(1, 2))))
                        gt = RawArrayTsvCreator.GT_encoding.HET1_2;
                }
            }
            return gt == RawArrayTsvCreator.value_to_drop ? "null" : gt.getValue();
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
}
