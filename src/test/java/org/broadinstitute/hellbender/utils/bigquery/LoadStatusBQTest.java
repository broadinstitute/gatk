package org.broadinstitute.hellbender.utils.bigquery;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.gvs.ingest.LoadStatus;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.UUID;

public class LoadStatusBQTest extends GATKBaseTest {

    private static final String BIGQUERY_TEST_PROJECT = "broad-dsde-dev";
    private static final String BIGQUERY_TEST_DATASET = "gatk_bigquery_test_dataset";

    private static final String uuid =  UUID.randomUUID().toString().replace('-', '_');
    private static final String TEMP_TABLE_NAME = "sample_load_status_test_" + uuid;

    private static final String FQ_TEMP_TABLE_NAME = String.format("%s.%s.%s",
            BIGQUERY_TEST_PROJECT, BIGQUERY_TEST_DATASET, TEMP_TABLE_NAME);

    @BeforeTest(groups = {"cloud"})
    public void beforeTest() {
        BigQueryUtils.executeQuery(
                "CREATE TABLE " + FQ_TEMP_TABLE_NAME + " (" +
                        SchemaUtils.SAMPLE_ID_FIELD_NAME + " INT64, " +
                        SchemaUtils.LOAD_STATUS_FIELD_NAME + " STRING, " +
                        SchemaUtils.LOAD_STATUS_EVENT_TIMESTAMP_NAME + " TIMESTAMP " +

                        ")", new HashMap<>());
    }

    @AfterTest(groups = {"cloud"})
    public void afterTest() {
        BigQueryUtils.executeQuery("DROP TABLE " + FQ_TEMP_TABLE_NAME, new HashMap<>());
    }

    @Test(groups = {"cloud"})
    public void testUpdateStatus() {
        LoadStatus loadStatus = new LoadStatus(BIGQUERY_TEST_PROJECT, BIGQUERY_TEST_DATASET, TEMP_TABLE_NAME);

        Assert.assertEquals(loadStatus.getSampleLoadState(1), LoadStatus.LoadState.NONE);

        loadStatus.writeLoadStatusStarted(1);
        Assert.assertEquals(loadStatus.getSampleLoadState(1), LoadStatus.LoadState.PARTIAL);

        loadStatus.writeLoadStatusFinished(1);
        Assert.assertEquals(loadStatus.getSampleLoadState(1), LoadStatus.LoadState.COMPLETE);

        Assert.assertEquals(loadStatus.getSampleLoadState(2), LoadStatus.LoadState.NONE);

    }
}
