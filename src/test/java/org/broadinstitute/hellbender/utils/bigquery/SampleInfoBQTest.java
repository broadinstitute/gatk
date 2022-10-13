package org.broadinstitute.hellbender.utils.bigquery;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.gvs.ingest.SampleInfo;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.UUID;

public class SampleInfoBQTest extends GATKBaseTest {

    private static final String BIGQUERY_TEST_PROJECT = "broad-dsde-dev";
    private static final String BIGQUERY_TEST_DATASET = "gatk_bigquery_test_dataset";

    private static final String uuid =  UUID.randomUUID().toString().replace('-', '_');
    private static final String TEMP_TABLE_NAME = "sample_info_test_" + uuid;

    private static final String FQ_TEMP_TABLE_NAME = String.format("%s.%s.%s",
            BIGQUERY_TEST_PROJECT, BIGQUERY_TEST_DATASET, TEMP_TABLE_NAME);

    @BeforeTest(groups = {"cloud"})
    public void beforeTest() {
        BigQueryUtils.executeQuery(
                "CREATE TABLE " + FQ_TEMP_TABLE_NAME + " (" +
                        SchemaUtils.SAMPLE_NAME_FIELD_NAME + " STRING, " +
                        SchemaUtils.SAMPLE_ID_FIELD_NAME + " INT64, " +
                        SchemaUtils.SAMPLE_INFO_IS_LOADED_FIELD_NAME + " BOOLEAN, " +
                        SchemaUtils.SAMPLE_INFO_IS_CONTROL_FIELD_NAME + " BOOLEAN, " +
                        SchemaUtils.SAMPLE_INFO_WITHDRAWN_FIELD_NAME + " TIMESTAMP " +
                        ")", new HashMap<>());
    }

    @AfterTest(groups = {"cloud"})
    public void afterTest() {
        BigQueryUtils.executeQuery("DROP TABLE " + FQ_TEMP_TABLE_NAME, new HashMap<>());
    }

    @Test(groups = {"cloud"})
    public void testUpdateStatus() {
        SampleInfo sampleInfo = new SampleInfo(BIGQUERY_TEST_PROJECT, BIGQUERY_TEST_DATASET, TEMP_TABLE_NAME);

        long sampleId = 1;
        sampleInfo.writeSampleInfo(sampleId, "NA12878", null, true, null);

        Assert.assertFalse(sampleInfo.isLoaded(sampleId));

        sampleInfo.setSampleInfoIsLoaded(sampleId);

        Assert.assertTrue(sampleInfo.isLoaded(sampleId));
    }
}
