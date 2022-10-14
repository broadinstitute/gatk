package org.broadinstitute.hellbender.tools.gvs.filtering;

import java.io.File;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.io.Resource;

public class ExtractFeaturesBQ {
    private final static String FEATURE_EXTRACT_QUERY_RESOURCE =
        "org/broadinstitute/hellbender/tools/gvs/filtering/feature_extract.sql";

    private final static String FEATURE_EXTRACT_USER_DEFINED_FUNCTIONS =
        "org/broadinstitute/hellbender/tools/gvs/filtering/udf_freq_table.sql";

    private final static String FEATURE_EXTRACT_USER_DEFINED_FUNCTION_MEDIAN =
            "org/broadinstitute/hellbender/tools/gvs/filtering/udf_median.sql";

    private final static String VQSR_TRAINING_SITES_TABLE =
        "broad-dsp-spec-ops.joint_genotyping_ref.vqsr_training_sites_*";

    public static String getVQSRFeatureExtractQueryString(final TableReference altAllele, final TableReference sampleList,
                                                          final Long minLocation, final Long maxLocation, final int hqGenotypeGQThreshold, final int hqGenotypeDepthThreshold, final double hqGenotypeABThreshold) {
        String locationStanza =
            ((minLocation != null)?"AND location >= " + minLocation + " \n":"") +
            ((maxLocation != null)?"AND location < " + maxLocation + " \n":"");

        try {
            File file =  Resource.getResourceContentsAsFile(FEATURE_EXTRACT_QUERY_RESOURCE);
            String query = FileUtils.readFileToString(file, "UTF-8");

            return query
                .replaceAll("@locationStanza", locationStanza)
                .replaceAll("@sample", sampleList.getFQTableName())
                .replaceAll("@altAllele", altAllele.getFQTableName())
                .replaceAll("@hqGenotypeGQThreshold", Double.toString(hqGenotypeGQThreshold))
                .replaceAll("@hqGenotypeDepthThreshold", Double.toString(hqGenotypeDepthThreshold))
                .replaceAll("@hqGenotypeABThreshold", Double.toString(hqGenotypeABThreshold));

        } catch (Exception ioe) {
            throw new GATKException("Unable to read query file from resources", ioe);
        }
    }

    public static String getVQSRFeatureExtractUserDefinedFunctionsString() {
        try {
            File file = Resource.getResourceContentsAsFile(FEATURE_EXTRACT_USER_DEFINED_FUNCTIONS);
            // also grab our definition of median
            File fileMedian = Resource.getResourceContentsAsFile(FEATURE_EXTRACT_USER_DEFINED_FUNCTION_MEDIAN);
            return FileUtils.readFileToString(file, "UTF-8")
                    + System.lineSeparator()
                    + FileUtils.readFileToString(fileMedian, "UTF-8");
        } catch (Exception ioe) {
            throw new GATKException("Unable to read udf file from resources", ioe);
        }
    }
}
