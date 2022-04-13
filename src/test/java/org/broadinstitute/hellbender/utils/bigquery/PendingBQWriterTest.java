package org.broadinstitute.hellbender.utils.bigquery;

import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import com.google.protobuf.Descriptors;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.json.JSONObject;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;
import java.util.concurrent.ExecutionException;

public class PendingBQWriterTest extends GATKBaseTest {
    private static final String BIGQUERY_TEST_PROJECT = "broad-dsde-dev";
    private static final String BIGQUERY_TEST_DATASET = "gatk_bigquery_test_dataset";

    private static final String uuid =  UUID.randomUUID().toString().replace('-', '_');
    private static final String TEMP_TABLE_NAME = "PendingTestTable_" + uuid;
    private static final String FQ_TEMP_TABLE_NAME = String.format("%s.%s.%s",
            BIGQUERY_TEST_PROJECT, BIGQUERY_TEST_DATASET, TEMP_TABLE_NAME);
    private static final String SAMPLE_ID_KEY = "sample_id";
    private static final String SAMPLE_NAME_KEY = "sample_name";


    @BeforeTest(groups = {"cloud"})
    public void beforeTest() {
        Map<String, String> labels = new HashMap<String, String>();
        labels.put("gatktestwriteapi", "testbefore" + uuid);

        BigQueryUtils.executeQuery("CREATE TABLE " + FQ_TEMP_TABLE_NAME + " (" + SAMPLE_ID_KEY + " INT64, " + SAMPLE_NAME_KEY + " STRING)", labels);
    }

    @AfterTest(groups = {"cloud"})
    public void afterTest() {
        Map<String, String> labels = new HashMap<String, String>();
        labels.put("gatktestwriteapi", "testafter" + uuid);

        BigQueryUtils.executeQuery("DROP TABLE " + FQ_TEMP_TABLE_NAME, labels);
    }

    @Test(groups = {"cloud"})
    public void testWriteToTable() {
        PendingBQWriter writer = new PendingBQWriter(BIGQUERY_TEST_PROJECT, BIGQUERY_TEST_DATASET, TEMP_TABLE_NAME);
        JSONObject jsonObject = new JSONObject();
        jsonObject.put(SAMPLE_ID_KEY, 1);
        jsonObject.put(SAMPLE_NAME_KEY,"Sample_1");

        try {
            writer.addJsonRow(jsonObject);
            writer.flushBuffer();
            writer.commitWriteStreams();
        } catch (Descriptors.DescriptorValidationException | ExecutionException | InterruptedException | IOException e) {
            e.printStackTrace();
        }

        Map<String, String> labels = new HashMap<String, String>();
        labels.put("gatktestwriteapi", "testwrite" + uuid);

        TableResult result = BigQueryUtils.executeQuery("SELECT * FROM " + FQ_TEMP_TABLE_NAME, labels);
        checkQueryResults(result, jsonObject);
    }

    private static void checkQueryResults(final TableResult result, final JSONObject expectedValues) {
        int rowCount = 0;
        for ( final FieldValueList row : result.iterateAll() ) {
            rowCount++;

            final long sample_id = row.get(0).getLongValue();
            final String sample_name = row.get(1).getStringValue();

            Assert.assertEquals(sample_id, expectedValues.getInt(SAMPLE_ID_KEY));
            Assert.assertEquals(sample_name, expectedValues.getString(SAMPLE_NAME_KEY));
        }
    }
}
