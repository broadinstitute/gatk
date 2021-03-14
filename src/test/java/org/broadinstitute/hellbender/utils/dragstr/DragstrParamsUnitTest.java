package org.broadinstitute.hellbender.utils.dragstr;

import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

/**
 * Unit test for class {@link DragstrParams}.
 */
public class DragstrParamsUnitTest extends BaseTest {

    public static final String TEST_PARAMS_FILE = "org/broadinstitute/hellbender/utils/dragstr/dna_nexus_novaseq_plus0_0_params.txt";

    @Test
    public void testToString() throws IOException {
        final String inputParamsFile = ("src/test/resources/" +  TEST_PARAMS_FILE);
        final DragstrParams subject = DragstrParamUtils.parse(new GATKPath(inputParamsFile));
        Assert.assertNotNull(subject.toString());
        final File pathFile = new File(subject.toString());
        Assert.assertEquals(pathFile.getCanonicalFile(), new File(inputParamsFile).getCanonicalFile());
    }

    @Test
    public void testReadWriteAndQuery() throws IOException {
        final String inputParamsFile = ("src/test/resources/" +  TEST_PARAMS_FILE);
        final File outputParamsFile = File.createTempFile("dragst-params-test", ".txt");
        outputParamsFile.deleteOnExit();
        final DragstrParams subject = DragstrParamUtils.parse(new GATKPath(inputParamsFile));
        DragstrParamUtils.print(subject, new GATKPath(outputParamsFile.toPath().toString()));
        try (final BufferedReader r1 = new BufferedReader(new FileReader(inputParamsFile));
             final BufferedReader r2 = new BufferedReader(new FileReader(outputParamsFile))) {
            String l1, l2;
            Query query = null;
            int row = 0;
            while ((l1 = readRelevantLine(r1)) != null & (l2 = readRelevantLine(r2)) != null) {

                // We test the read-and-write checking that the content of the new copied params file is the same as the original (except from some irrelevant white-spaces.)
                Assert.assertEquals(l1.trim(), l2.trim());

                // To test the query here we implement a simple and independent parsing of the table rows and check
                // that the values read match the ones return by the subject.

                // When we encounter a table header we update the query lambda appropriately and reset the table row index.
                if (l1.contains("GOP")) {
                    query = subject::gop; row = 0;
                } else if (l1.contains("GCP")) {
                    query = subject::gcp; row = 0;
                } else if (l1.contains("API")) {
                    query = subject::api; row = 0;
                } else if (query != null) {
                    // when we read in a table row, we check that the subject returns the same values for that row/period and repeat length/column.
                    row++;
                    // The way the tables are formatted split return empty string (first or last) element so we filter them out using a stream.
                    final String[] l2parts = Arrays.stream(l2.split("\\s+")).filter(s -> !s.trim().isEmpty()).toArray(String[]::new);
                    Assert.assertEquals(l2parts.length, subject.maximumRepeats(), l2);
                    for (int i = 0; i < l2parts.length; i++) {
                        final double expected = Double.parseDouble(l2parts[i]);
                        final double actual = query.get(row, i + 1);
                        Assert.assertEquals(actual, expected, 0.001);
                    }
                }
            }
            Assert.assertNull(l1);
            Assert.assertNull(l2);
        } finally {
            outputParamsFile.delete();
        }
    }

    private String readRelevantLine(final BufferedReader reader) throws IOException {
        String line;
        while ((line = reader.readLine()) != null) {
            if (!line.startsWith("#") && !line.trim().isEmpty()) {
                return line;
            }
        }
        return null;
    }
    
    @Test(dependsOnMethods = "testReadWriteAndQuery")
    public void additionalSelfConsistencyTests() {
        final String filePath = ("src/test/resources/" +  TEST_PARAMS_FILE);
        final DragstrParams params = DragstrParamUtils.parse(new GATKPath(filePath));
        Assert.assertEquals(params.maximumPeriod(), 8); // we know this for a fact.
        Assert.assertEquals(params.maximumRepeats(), 20); // we know this for a fact.
        Assert.assertEquals(params.maximumLengthInBasePairs(), params.maximumRepeats() * params.maximumPeriod());
        for (int i = 1; i <= params.maximumPeriod(); i++) {
            // For each table, we check that repeat lengths larger than the maximum record return the same value as the maximum.
            // We test maximum + 1 and maximum + N where N >>> 1
            Assert.assertEquals(params.gop(i, params.maximumRepeats()), params.gop(i, params.maximumRepeats() + 1));
            Assert.assertEquals(params.gop(i, params.maximumRepeats()), params.gop(i, params.maximumRepeats() * 31 + 11)); // any large repeat-length would do.
            Assert.assertEquals(params.gcp(i, params.maximumRepeats()), params.gcp(i, params.maximumRepeats() + 1));
            Assert.assertEquals(params.gcp(i, params.maximumRepeats()), params.gcp(i, params.maximumRepeats() * 131 + 911));  // any large repeat-length would do.
            Assert.assertEquals(params.api(i, params.maximumRepeats()), params.api(i, params.maximumRepeats() + 1));
            Assert.assertEquals(params.api(i, params.maximumRepeats()), params.api(i, params.maximumRepeats() * 311 + 1193));  // any large repeat-length would do.
        }
    }

    @FunctionalInterface
    interface Query {
        double get(final int period, final int repeat);
    }
}
