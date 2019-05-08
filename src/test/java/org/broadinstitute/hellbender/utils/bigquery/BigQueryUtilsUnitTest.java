package org.broadinstitute.hellbender.utils.bigquery;

import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.Map;

/**
 * A class to test the functionality of {@link BigQueryUtils}.
 */
public class BigQueryUtilsUnitTest extends GATKBaseTest {

    private static final String BIGQUERY_TEST_PROJECT = "broad-dsde-dev";
    private static final String BIGQUERY_TEST_DATASET = "gatk_bigquery_test_dataset";
    private static final String BIGQUERY_TEST_TABLE = "gatk_bigquery_test_table1";
    private static final String BIGQUERY_FULLY_QUALIFIED_TABLE = String.format("%s.%s.%s",
            BIGQUERY_TEST_PROJECT, BIGQUERY_TEST_DATASET, BIGQUERY_TEST_TABLE);

    @Test(groups = {"cloud"})
    public void testExecuteQueryAllRecords() {
        final Map<String, String> expectedNamesAndAges = new HashMap<>();
        expectedNamesAndAges.put("Fred", "35");
        expectedNamesAndAges.put("Matilda", "28");
        expectedNamesAndAges.put("Harry", null);

        final String QUERY = String.format("SELECT * FROM `%s`", BIGQUERY_FULLY_QUALIFIED_TABLE);

        final TableResult result = BigQueryUtils.executeQuery(QUERY);

        checkQueryResults(result, expectedNamesAndAges, QUERY);
    }

    @Test(groups = {"cloud"})
    public void testExecuteQueryWithWhereClause() {
        final Map<String, String> expectedNamesAndAges = new HashMap<>();
        expectedNamesAndAges.put("Fred", "35");

        final String QUERY = String.format("SELECT * FROM `%s` WHERE name = 'Fred'", BIGQUERY_FULLY_QUALIFIED_TABLE);

        final TableResult result = BigQueryUtils.executeQuery(QUERY);

        checkQueryResults(result, expectedNamesAndAges, QUERY);
    }

    private static void checkQueryResults(final TableResult result, final Map<String, String> expectedNamesAndAges, final String query) {
        int rowCount = 0;
        for ( final FieldValueList row : result.iterateAll() ) {
            rowCount++;

            final String name = row.get(0).getStringValue();
            final String age = row.get(1).isNull() ? null : row.get(1).getStringValue();

            Assert.assertTrue(expectedNamesAndAges.containsKey(name), "Unexpected name " + name + " returned from query " + query);
            Assert.assertEquals(expectedNamesAndAges.get(name), age, "Wrong age " + age + " returned for name " + name + " in query " + query);
        }

        Assert.assertEquals(rowCount, expectedNamesAndAges.size(), "Wrong number of rows returned from query: " + query);
    }
}
