package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.cloud.bigquery.FieldValueList;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryResultAndStatistics;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;

public class SampleInfo {
    private final String projectID;
    private final String datasetName;
    private final String sampleInfoTableName;


    public SampleInfo(String projectID, String datasetName, String sampleInfoTableName) {
        this.projectID = projectID;
        this.datasetName = datasetName;
        this.sampleInfoTableName = sampleInfoTableName;
    }

    public void setSampleInfoIsLoaded(long sampleId) {
        final String query = "UPDATE `" + projectID + "." + datasetName + "." + sampleInfoTableName + "` " +
                " SET " + SchemaUtils.SAMPLE_INFO_IS_LOADED_FIELD_NAME + " = TRUE " +
                " WHERE " + SchemaUtils.SAMPLE_ID_FIELD_NAME + " = " + sampleId;

        BigQueryUtils.executeQuery(projectID, query, false, null);
    }

    public boolean isLoaded(long sampleId) {
        final String query = "SELECT " + SchemaUtils.SAMPLE_INFO_IS_LOADED_FIELD_NAME +
                " FROM `" + projectID + "." + datasetName + "." + sampleInfoTableName + "` " +
                " WHERE " + SchemaUtils.SAMPLE_ID_FIELD_NAME + " = " + sampleId;

        boolean isLoaded = false;
        BigQueryResultAndStatistics results = BigQueryUtils.executeQuery(projectID, query, false, null);
        if (results.result.getTotalRows() > 1) {
            // TODO - is it overkill to throw an exception here?
            throw new GATKException("Found more than 1 row in '" + projectID + "." + datasetName + "." + sampleInfoTableName + "' for sampleId: " + sampleId);
        }
        for ( final FieldValueList row : results.result.iterateAll() ) {
            isLoaded = !row.get(0).isNull() && row.get(0).getBooleanValue();
        }

        return isLoaded;
    }

    public void writeSampleInfo(long sampleId, String sampleName, boolean is_control) {
        final String query = String.format("INSERT INTO `%s.%s.%s` (%s, %s, %s, %s, %s) VALUES (\"%s\", %d, null, %s, null) ", projectID, datasetName, sampleInfoTableName, SchemaUtils.SAMPLE_NAME_FIELD_NAME, SchemaUtils.SAMPLE_ID_FIELD_NAME, SchemaUtils.SAMPLE_INFO_IS_LOADED_FIELD_NAME, SchemaUtils.SAMPLE_INFO_IS_CONTROL_FIELD_NAME, SchemaUtils.SAMPLE_INFO_WITHDRAWN_FIELD_NAME, sampleName, sampleId, is_control);
        BigQueryUtils.executeQuery(projectID, query, false, null);
    }
}
