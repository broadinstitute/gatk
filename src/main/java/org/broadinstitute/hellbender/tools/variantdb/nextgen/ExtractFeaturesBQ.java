package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import java.io.File;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.io.Resource;

public class ExtractFeaturesBQ {
    private final static String FEATURE_EXTRACT_QUERY_RESOURCE = 
        "org/broadinstitute/hellbender/tools/variantdb/nextgen/feature_extract.sql";

    private final static String VQSR_TRAINING_SITES_TABLE = 
        "broad-dsp-spec-ops.joint_genotyping_ref.vqsr_training_sites_*";

    public static String getVQSRFeatureExtractQueryString(final TableReference altAllele, final TableReference sampleList,
                                                          final Long minLocation, final Long maxLocation, final boolean trainingSitesOnly) {

        String trainingSitesStanza =
            !trainingSitesOnly?"":
                "AND location IN (SELECT location FROM `" + VQSR_TRAINING_SITES_TABLE + "`)\n";

        String locationStanza = 
            ((minLocation != null)?"AND location >= " + minLocation + " \n":"") +
            ((maxLocation != null)?"AND location < " + maxLocation + " \n":"");

        try {            
            File file =  Resource.getResourceContentsAsFile(FEATURE_EXTRACT_QUERY_RESOURCE);
            String query = FileUtils.readFileToString(file, "UTF-8");

            return query
                .replaceAll("@locationStanza", locationStanza)
                .replaceAll("@trainingSitesStanza", trainingSitesStanza)
                .replaceAll("@sample", sampleList.getFQTableName())
                .replaceAll("@altAllele", altAllele.getFQTableName());

        } catch (Exception ioe) {
            throw new GATKException("Unable to read query file from resources", ioe);
        }
    }
}
