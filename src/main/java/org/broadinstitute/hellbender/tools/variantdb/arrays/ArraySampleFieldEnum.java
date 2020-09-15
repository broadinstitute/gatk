package org.broadinstitute.hellbender.tools.variantdb.arrays;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.variantdb.CommonCode;
import org.broadinstitute.hellbender.tools.variantdb.arrays.tables.ProbeInfo;

import java.util.*;



public enum ArraySampleFieldEnum {
    sample_id,
    sample_name,


    // This where the validation step (required vs not) lives  -- fail if there is missing data for a required field
    // and just leave it empty if not required

    NUM_ASSAYS {
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },
    NUM_NON_FILTERED_ASSAYS{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },
    NUM_FILTERED_ASSAYS{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },
    NUM_ZEROED_OUT_ASSAYS{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },
    NUM_SNPS{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },
    NUM_INDELS{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },
    NUM_CALLS{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },
    NUM_AUTOCALL_CALLS{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },
    NUM_NO_CALLS{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },
    NUM_IN_DB_SNP{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },
    NOVEL_SNPS{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },

    PCT_DBSNP{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },
    CALL_RATE{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },
    AUTOCALL_CALL_RATE{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    },
    NUM_SINGLETONS{
        public String getColumnValue(final Map<String, String> metricsMap) {
            return metricsMap.get(this.name());
        }
    };
    
    public String getColumnValue(final Map<String, String> metricsMap) {
        throw new IllegalArgumentException("Not implemented");
    }

}
