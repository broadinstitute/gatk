package org.broadinstitute.hellbender.utils.bigquery;

import com.google.cloud.bigquery.FieldValueList;
import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import java.util.*;

/**
 * A class to test the functionality of {@link BigQueryUtils}.
 */
public class BigQueryUtilsUnitTest extends GATKBaseTest {

    private static final String BIGQUERY_TEST_PROJECT = "broad-dsde-dev";
    private static final String BIGQUERY_TEST_DATASET = "gatk_bigquery_test_dataset";
    private static final String BIGQUERY_TEST_TABLE = "gatk_bigquery_test_table1";
    private static final String BIGQUERY_FULLY_QUALIFIED_TABLE = String.format("%s.%s.%s",
            BIGQUERY_TEST_PROJECT, BIGQUERY_TEST_DATASET, BIGQUERY_TEST_TABLE);

    private static final UUID runUuid = UUID.randomUUID();

    @Test(groups = {"cloud"})
    public void testExecuteQueryAllRecords() {
        final String query = String.format("SELECT * FROM `%s`", BIGQUERY_FULLY_QUALIFIED_TABLE);
        Map<String, String> labels = new HashMap<>();
        labels.put("gatktestquery", "getallrecords" + runUuid);
        final BigQueryResultAndStatistics resultAndStatistics = BigQueryUtils.executeQuery(query, labels);

        checkQueryResults(resultAndStatistics, getAllExpectedNamesAndAges(), query);
    }

    @Test(groups = {"cloud"})
    public void testExecuteQueryWithWhereClause() {
        final Map<String, String> expectedNamesAndAges = new HashMap<>();
        expectedNamesAndAges.put("Fred", "35");

        final String query = String.format("SELECT * FROM `%s` WHERE name = 'Fred'", BIGQUERY_FULLY_QUALIFIED_TABLE);
        Map<String, String> labels = new HashMap<String, String>();
        labels.put("gatktestquery", "testwhereclause" + runUuid);
        logger.info("testwhereclause" + runUuid);
        final BigQueryResultAndStatistics resultAndStatistics = BigQueryUtils.executeQuery(query, labels);

        checkQueryResults(resultAndStatistics, expectedNamesAndAges, query);
    }

    @Test(groups = {"cloud"})
    public void testExecuteQueryInBatchMode() {
        final String query = String.format("SELECT * FROM `%s`", BIGQUERY_FULLY_QUALIFIED_TABLE);
        Map<String, String> labels = new HashMap<>();
        labels.put("gatktestquery", "testbatchmode" + runUuid);
        final BigQueryResultAndStatistics resultAndStatistics = BigQueryUtils.executeQuery(query, true, labels);

        checkQueryResults(resultAndStatistics, getAllExpectedNamesAndAges(), query);
    }

    @Test(groups = {"cloud"})
    public void testExecuteQueryInBatchModeWithProject() {
        final String query = String.format("SELECT * FROM `%s`", BIGQUERY_FULLY_QUALIFIED_TABLE);
        Map<String, String> labels = new HashMap<>();
        labels.put("gatktestquery", "testbatchmode" + runUuid);
        final BigQueryResultAndStatistics resultAndStatistics = BigQueryUtils.executeQuery(BIGQUERY_TEST_PROJECT, query, true, labels);

        checkQueryResults(resultAndStatistics, getAllExpectedNamesAndAges(), query);
    }

    @Test(groups = {"cloud"})
    public void testSpecifiedExecuteQuery() {
        final String query = String.format("SELECT * FROM `%s`", BIGQUERY_TEST_TABLE);
        Map<String, String> labels = new HashMap<>();
        labels.put("gatktestquery", "testspecifiedexecutequery" + runUuid);
        final BigQueryResultAndStatistics resultAndStatistics = BigQueryUtils.executeQuery(BigQueryUtils.getBigQueryEndPoint(), BIGQUERY_TEST_PROJECT, BIGQUERY_TEST_DATASET, query, labels);

        checkQueryResults(resultAndStatistics, getAllExpectedNamesAndAges(), query);
    }

    @Test(groups = {"cloud"})
    public void testQueryWithStorageAPI() {
        final Map<String, String> expectedNamesAndAges = getAllExpectedNamesAndAges();

        final String query = String.format("SELECT * FROM `%s`", BIGQUERY_FULLY_QUALIFIED_TABLE);

        final List<String> fieldsToRetrieve = new LinkedList<>();
        final String noUDFs = null;
        fieldsToRetrieve.add("name");
        Map<String, String> labels = new HashMap<>();
        labels.put("gatktestquery", "teststorageapi" + runUuid);
        final StorageAPIAvroReader result = BigQueryUtils.executeQueryWithStorageAPI(query, fieldsToRetrieve, BIGQUERY_TEST_PROJECT, BIGQUERY_TEST_DATASET, noUDFs, false, labels).storageAPIAvroReader;

        int rowCount = 0;
        final Set<String> retrievedNames = new HashSet<>();

        while( result.hasNext() ) {
            GenericRecord row = result.next();
            String name = row.get("name").toString();
            retrievedNames.add(name);
            rowCount++;

            Assert.assertTrue(expectedNamesAndAges.containsKey(name), "Unexpected name " + name + " returned from query " + query);
            Assert.expectThrows(org.apache.avro.AvroRuntimeException.class, () -> row.get("age"));
        }

        final Set<String> expectedNames = new HashSet<>();
        expectedNames.add("Fred");
        expectedNames.add("Matilda");
        expectedNames.add("Harry");

        Assert.assertEquals(retrievedNames, expectedNames, "Retrieved names did not match.");
        Assert.assertEquals(rowCount, expectedNamesAndAges.size(), "Size of returned results does not match expected.");
    }

    @Test(groups = {"cloud"})
    public void testQueryWithEmptyDatasetStorageAPI() {
        final String query = String.format("SELECT * FROM `%s` WHERE false", BIGQUERY_FULLY_QUALIFIED_TABLE);

        final List<String> fieldsToRetrieve = new LinkedList<>();
        final String noUDFs = null;
        fieldsToRetrieve.add("name");
        Map<String, String> labels = new HashMap<>();
        labels.put("gatktestquery", "testapiwithemptydata" + runUuid);
        final StorageAPIAvroReader result = BigQueryUtils.executeQueryWithStorageAPI(query, fieldsToRetrieve, BIGQUERY_TEST_PROJECT, BIGQUERY_TEST_DATASET, noUDFs, false, labels).storageAPIAvroReader;

        int rowCount = 0;
        final Set<String> retrievedNames = new HashSet<>();

        while( result.hasNext() ) {
            Assert.fail("No Result expected");
        }

        Assert.assertTrue(retrievedNames.isEmpty(), "No Result expected");
    }

    @Test(groups = {"cloud"})
    public void testQueryWithNullLabel() {
        final String query = String.format("SELECT * FROM `%s`", BIGQUERY_FULLY_QUALIFIED_TABLE);
        Map<String, String> labels = null;
        final BigQueryResultAndStatistics resultAndStatistics = BigQueryUtils.executeQuery(query, labels);

        checkQueryResults(resultAndStatistics, getAllExpectedNamesAndAges(), query);
    }

    private Map<String, String> getAllExpectedNamesAndAges() {
        final Map<String, String> expectedNamesAndAges = new HashMap<>();
        expectedNamesAndAges.put("Fred", "35");
        expectedNamesAndAges.put("Matilda", "28");
        expectedNamesAndAges.put("Harry", null);
        return expectedNamesAndAges;
    }

    private static void checkQueryResults(final BigQueryResultAndStatistics resultAndStatistics, final Map<String, String> expectedNamesAndAges, final String query) {
        int rowCount = 0;
        for ( final FieldValueList row : resultAndStatistics.result.iterateAll() ) {
            rowCount++;

            final String name = row.get(0).getStringValue();
            final String age = row.get(1).isNull() ? null : row.get(1).getStringValue();

            Assert.assertTrue(expectedNamesAndAges.containsKey(name), "Unexpected name " + name + " returned from query " + query);
            Assert.assertEquals(expectedNamesAndAges.get(name), age, "Wrong age " + age + " returned for name " + name + " in query " + query);
        }

        Assert.assertEquals(rowCount, expectedNamesAndAges.size(), "Wrong number of rows returned from query: " + query);
    }
}
