package org.broadinstitute.hellbender.tools.variantdb.arrays;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.variantdb.CommonCode;
import org.broadinstitute.hellbender.tools.variantdb.arrays.tables.ProbeInfo;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


/**
 * Expected headers for the  uncompressed array table
 *     sample, // required
 *     probe_id, // required
 *     GT_encoded, // required
 *     NORMX, // intensity
 *     NORMY, // intensity
 *     BAF // b allele fraction --> AD proxy
 *     LRR // Log R ratio --> intensity value instead of DP
 * 
 */

public enum RawArrayFieldEnum {
    // fail if there is missing data for a required field
    // and return the string "null" if there is missing data for an optional field
    // (in the bq import, it will convert the "null" to an actual null value in the database

    sample_id,
    probe_id,
    GT_encoded { // Required
        public String getColumnValue(final VariantContext variant) {
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
            } else {
                logger.warn("Found " + alleleIndexes.size() + " alleles instead of 2. Not processing variant \t" + variant);
            }
            return gt == RawArrayTsvCreator.value_to_drop ? "null" : gt.getValue();
        }
    },

    NORMX { // Required
        public String getColumnValue(final VariantContext variant) {
            Object value = variant.getGenotype(0).getExtendedAttribute("NORMX");
            if (value == null) {
                throw new IllegalStateException("Missing required value NORMX for variant:  \t" + variant);
            }
            return  String.valueOf(value);
        }
    },
    NORMY { // Required
        public String getColumnValue(final VariantContext variant) {
            Object value = variant.getGenotype(0).getExtendedAttribute("NORMY");
            if (value == null) {
                throw new IllegalStateException("Missing required value NORMY for variant:  \t" + variant);
            }
            return  String.valueOf(value);
        }
    },
    BAF {
        public String getColumnValue(final VariantContext variant) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("BAF"));
        }
    },
    LRR {
        public String getColumnValue(final VariantContext variant) {
            return  String.valueOf(variant.getGenotype(0).getExtendedAttribute("LRR"));
        }
    };

    private static final Logger logger = LogManager.getLogger(RawArrayFieldEnum.class);

    public String getColumnValue(final VariantContext variant) {
        throw new IllegalArgumentException("Not implemented");
    }

    public static RawArrayFieldEnum[] getUncompressedRawArrayFieldEnums() {
        return new RawArrayFieldEnum[] { sample_id, probe_id, GT_encoded, NORMX, NORMY, BAF, LRR };
    }
}
