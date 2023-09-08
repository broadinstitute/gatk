package org.broadinstitute.hellbender.tools.genomicsdb;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.LinkedHashMap;
import java.util.Map;
import java.net.URI;
import java.net.URISyntaxException;

public class SampleNameMapUnitTest extends GATKBaseTest {

    private static final String ORDERED_SAMPLE_MAP =    "Sample1\tfile1\n" +
                                                        "Sample2\tfile2\n" +
                                                        "Sample3\tfile3";

    private static final String UNORDERED_SAMPLE_MAP =  "Sample3\tfile3\n" +
                                                        "Sample2\tfile2\n" +
                                                        "Sample1\tfile1\n";

    @DataProvider
    public Object[][] getBadSampleNameMapFiles(){
        return new Object[][]{
                {"Sample1\tsamplePath\n"
                +"Sample1\tsamplePath"},        // duplicate sample name
                {""},                           // empty file
                {"Sample1\t"},                  // 1 column
                {"Sample1"},                    // 1 column no delimiter
                {"\tfile"},                     // empty first token
                {" \tfile"},                    // first token only whitespace
                {"Sample1\tfile1\t"},           // extra tab
                {"Sample1\nfile"},              // newline instead of tab
                {"\t"},                         // only tab
                {"Sample1 file1"},              // 1 column, internal whitespace
                {" Sample1\tfile1"},            // preceding whitespace
                {"Sample1 \tfile1"},            // trailing whitespace
                {"Sample1\tfile1\t"},                      // empty index
                {"Sample1\tfile1\t "},                     // all-whitespace index
                {"Sample1\tfile1\tindex1\textraColumn"},   // 4 columns
                {"Sample1\tfile1\tindex1\t"}               // 4 columns, blank 4th column
        };
    }

    @Test(dataProvider = "getBadSampleNameMapFiles", expectedExceptions = UserException.BadInput.class)
    public void testBadInputFiles(final String text){
        final File sampleFile = IOUtils.writeTempFile(text, "badSampleMapping", ".txt");
        final SampleNameMap sampleMap = new SampleNameMap(sampleFile.toPath());
    }

    @DataProvider
    public Object[][] getGoodSampleNameMapFiles(){
        return new Object[][]{
                // Note: none of these files are real, these are just valid files syntactically

                // normal sample names, no explicit indices
                {"Sample1\tsamplePath1\n" +
                 "Sample2\tsamplePath2",
                        new String[][] {
                                {"Sample1", "samplePath1"},
                                {"Sample2", "samplePath2"}}},

                // normal sample names, explicit indices for all files
                {"Sample1\tsamplePath1\tindexPath1\n" +
                 "Sample2\tsamplePath2\tindexPath2",
                        new String[][] {
                                {"Sample1", "samplePath1", "indexPath1"},
                                {"Sample2", "samplePath2", "indexPath2"}}},

                // normal sample names, explicit indices for some files but not others
                {"Sample1\tsamplePath1\n" +
                 "Sample2\tsamplePath2\tindexPath2",
                        new String[][] {
                                {"Sample1", "samplePath1"},
                                {"Sample2", "samplePath2", "indexPath2"}}},

                // sample names with internal whitespace
                {"Sample1 001\tFile",
                        new String[][] {
                                {"Sample1 001", "File"}}
                },
                
                // leading and trailing whitespace second column
                {"name name\t file1 ",
                        new String[][] {
                                {"name name", "file1"}}
                },

                // leading and trailing whitespace third column
                {"name name\tfile1\t index1 ",
                        new String[][] {
                                {"name name", "file1", "index1"}}
                },

                // leading and trailing whitespace second and third columns
                {"name name\t file1 \t index1 ",
                        new String[][] {
                                {"name name", "file1", "index1"}}
                },
        };
    }

    @Test(dataProvider = "getGoodSampleNameMapFiles")
    public void testValidSampleFiles(final String text, final String[][] expectedEntries){
        final File sampleFile = IOUtils.writeTempFile(text, "goodSampleMapping", ".txt");

        final SampleNameMap sampleMap = new SampleNameMap(sampleFile.toPath());
        final SortedMap<String, URI> outputMap = sampleMap.getSampleNameToVcfPath();

        Assert.assertEquals(outputMap.size(),expectedEntries.length,
                "Wrong number of entries in the Map returned by getSampleNameToVcfPath()");
        Assert.assertEquals(sampleMap.getNumSamples(), expectedEntries.length,
                "Wrong number of samples returned by getNumSamples()");
        boolean expectedIndicesFound = false;

        for ( final String[] expected : expectedEntries ) {
            Assert.assertTrue(outputMap.containsKey(expected[0]));

            Assert.assertEquals(outputMap.get(expected[0]).toString(),expected[1]);
            Assert.assertEquals(sampleMap.getVCFForSample(expected[0]).toString(), expected[1],
                    "Wrong VCF returned by getVCFForSample() for sample " + expected[0]);
            Assert.assertEquals(sampleMap.getVCFForSampleAsPath(expected[0]).toString(), expected[1],
                    "Wrong VCF returned by getVCFForSampleAsPath() for sample " + expected[0]);

            if ( expected.length == 3 ) {
                expectedIndicesFound = true;

                Assert.assertNotNull(sampleMap.getVCFIndexForSample(expected[0]),
                        "No index returned by getVCFIndexForSample() for sample " + expected[0]);
                Assert.assertNotNull(sampleMap.getVCFForSampleAsPath(expected[0]),
                        "No index returned by getVCFForSampleAsPath() for sample " + expected[0]);

                Assert.assertEquals(sampleMap.getVCFIndexForSample(expected[0]).toString(), expected[2],
                        "Wrong index returned by getVCFIndexForSample() for sample " + expected[0]);
                Assert.assertEquals(sampleMap.getVCFIndexForSampleAsPath(expected[0]).toString(), expected[2],
                        "Wrong index returned by getVCFIndexForSampleAsPath() for sample " + expected[0]);
            } else {
                Assert.assertNull(sampleMap.getVCFIndexForSample(expected[0]),
                        "Index unexpectedly returned by getVCFIndexForSample() for sample " + expected[0]);
                Assert.assertNull(sampleMap.getVCFIndexForSampleAsPath(expected[0]),
                        "Index unexpectedly returned by getVCFIndexForSampleAsPath() for sample " + expected[0]);
            }
        }

        Assert.assertEquals(sampleMap.indicesSpecified(), expectedIndicesFound,
                "Wrong value returned by indicesSpecified()");

    }

    // Test to ensure that the "unsorted" map used in subsequent tests is actually unsorted,
    // to guard against future modifications
    @Test
    public void testUnorderedSampleMapIsActuallyUnordered() throws IOException {
        final File sampleFile = IOUtils.writeTempFile(UNORDERED_SAMPLE_MAP, "badSampleMapping", ".txt");
        final List<String> expectedSampleOrdering = Arrays.asList("Sample3", "Sample2", "Sample1");

        try ( final XReadLines lineReader = new XReadLines(sampleFile) ) {
            int lineNumber = 0;
            for ( final String line : lineReader ) {
                final String sampleFromFile = line.split("\\t", -1)[0];
                Assert.assertEquals(sampleFromFile, expectedSampleOrdering.get(lineNumber));
                ++lineNumber;
            }
        }
    }

    @DataProvider
    public Object[][] getSampleMapsForOrderingTest(){
        final Map<String, URI> expectedMap = new LinkedHashMap<>();
        try {
            expectedMap.put("Sample1", new URI("file1"));
            expectedMap.put("Sample2", new URI("file2"));
            expectedMap.put("Sample3", new URI("file3"));
        }
        catch(URISyntaxException e) {
            throw new RuntimeException("Malformed URI "+e.toString());
        }

        final List<String> expectedSampleOrdering = Arrays.asList("Sample1", "Sample2", "Sample3");

        return new Object[][]{
                {ORDERED_SAMPLE_MAP, expectedMap, expectedSampleOrdering},
                {UNORDERED_SAMPLE_MAP, expectedMap, expectedSampleOrdering}
        };
    }

    @Test(dataProvider = "getSampleMapsForOrderingTest")
    public void testSampleOrdering(final String sampleMapText, final Map<String, URI> expectedMap, final List<String> expectedSampleOrdering){
        final File sampleFile = IOUtils.writeTempFile(sampleMapText, "sampleMapping", ".txt");

        final SampleNameMap sampleMap = new SampleNameMap(sampleFile.toPath());
        final SortedMap<String, URI> actualMap = sampleMap.getSampleNameToVcfPath();
        
        Assert.assertEquals(actualMap, expectedMap);
        Assert.assertEquals(sampleMap.getNumSamples(), expectedSampleOrdering.size(), "Wrong number of samples returned by getNumSamples()");
        Assert.assertEquals(sampleMap.getSampleNamesInSortedOrder().size(), expectedSampleOrdering.size(), "Wrong number of samples returned by getSampleNamesInOrder()");

        final Iterator<String> actualSamplesFromMap = actualMap.keySet().iterator();
        final Iterator<String> actualSamplesFromGetter = sampleMap.getSampleNamesInSortedOrder().iterator();

        for ( final String expectedSample : expectedSampleOrdering ) {
            Assert.assertEquals(actualSamplesFromMap.next(), expectedSample,
                    "Wrong sample found in Map returned by getSampleNameToVcfPath()");
            Assert.assertEquals(actualSamplesFromGetter.next(), expectedSample,
                    "Wrong sample found in Set returned by getSampleNamesInOrder()");
        }
        Assert.assertFalse(actualSamplesFromMap.hasNext());
        Assert.assertFalse(actualSamplesFromGetter.hasNext());
    }

    @Test
    public void testIncrementalAddition() throws URISyntaxException {
        // Use the no-arg constructor to start with an empty SampleNameMap, then
        // add samples incrementally:
        final SampleNameMap sampleMap = new SampleNameMap();
        Assert.assertEquals(sampleMap.getNumSamples(), 0);
        Assert.assertTrue(sampleMap.getSampleNamesInSortedOrder().isEmpty());
        Assert.assertTrue(sampleMap.getSampleNameToVcfPath().isEmpty());
        Assert.assertFalse(sampleMap.indicesSpecified());

        sampleMap.addSample("Sample3", new URI("file3"));
        Assert.assertEquals(sampleMap.getNumSamples(), 1);
        Assert.assertEquals(sampleMap.getSampleNamesInSortedOrder().size(), 1);
        Assert.assertEquals(sampleMap.getSampleNamesInSortedOrder().get(0), "Sample3");
        Assert.assertEquals(sampleMap.getSampleNameToVcfPath().size(), 1);
        Assert.assertEquals(new ArrayList<>(sampleMap.getSampleNameToVcfPath().keySet()), Arrays.asList("Sample3"));
        Assert.assertEquals(sampleMap.getVCFForSample("Sample3").toString(), "file3");
        Assert.assertEquals(sampleMap.getVCFForSampleAsPath("Sample3").toString(), "file3");
        Assert.assertFalse(sampleMap.indicesSpecified());

        sampleMap.addSample("Sample1", new URI("file1"), new URI("index1"));
        Assert.assertEquals(sampleMap.getNumSamples(), 2);
        Assert.assertEquals(sampleMap.getSampleNamesInSortedOrder().size(), 2);
        Assert.assertEquals(sampleMap.getSampleNamesInSortedOrder().get(0), "Sample1");
        Assert.assertEquals(sampleMap.getSampleNamesInSortedOrder().get(1), "Sample3");
        Assert.assertEquals(sampleMap.getSampleNameToVcfPath().size(), 2);
        Assert.assertEquals(new ArrayList<>(sampleMap.getSampleNameToVcfPath().keySet()), Arrays.asList("Sample1", "Sample3"));
        Assert.assertEquals(sampleMap.getVCFForSample("Sample1").toString(), "file1");
        Assert.assertEquals(sampleMap.getVCFForSampleAsPath("Sample1").toString(), "file1");
        Assert.assertEquals(sampleMap.getVCFForSample("Sample3").toString(), "file3");
        Assert.assertEquals(sampleMap.getVCFForSampleAsPath("Sample3").toString(), "file3");
        Assert.assertTrue(sampleMap.indicesSpecified());
        Assert.assertEquals(sampleMap.getVCFIndexForSample("Sample1").toString(), "index1");
        Assert.assertEquals(sampleMap.getVCFIndexForSampleAsPath("Sample1").toString(), "index1");
        Assert.assertNull(sampleMap.getVCFIndexForSample("Sample3"));
        Assert.assertNull(sampleMap.getVCFIndexForSampleAsPath("Sample3"));

        sampleMap.addSample("Sample2", new URI("file2"));
        Assert.assertEquals(sampleMap.getNumSamples(), 3);
        Assert.assertEquals(sampleMap.getSampleNamesInSortedOrder().size(), 3);
        Assert.assertEquals(sampleMap.getSampleNamesInSortedOrder().get(0), "Sample1");
        Assert.assertEquals(sampleMap.getSampleNamesInSortedOrder().get(1), "Sample2");
        Assert.assertEquals(sampleMap.getSampleNamesInSortedOrder().get(2), "Sample3");
        Assert.assertEquals(sampleMap.getSampleNameToVcfPath().size(), 3);
        Assert.assertEquals(new ArrayList<>(sampleMap.getSampleNameToVcfPath().keySet()), Arrays.asList("Sample1", "Sample2", "Sample3"));
        Assert.assertEquals(sampleMap.getVCFForSample("Sample1").toString(), "file1");
        Assert.assertEquals(sampleMap.getVCFForSampleAsPath("Sample1").toString(), "file1");
        Assert.assertEquals(sampleMap.getVCFForSample("Sample2").toString(), "file2");
        Assert.assertEquals(sampleMap.getVCFForSampleAsPath("Sample2").toString(), "file2");
        Assert.assertEquals(sampleMap.getVCFForSample("Sample3").toString(), "file3");
        Assert.assertEquals(sampleMap.getVCFForSampleAsPath("Sample3").toString(), "file3");
        Assert.assertTrue(sampleMap.indicesSpecified());
        Assert.assertEquals(sampleMap.getVCFIndexForSample("Sample1").toString(), "index1");
        Assert.assertEquals(sampleMap.getVCFIndexForSampleAsPath("Sample1").toString(), "index1");
        Assert.assertNull(sampleMap.getVCFIndexForSample("Sample3"));
        Assert.assertNull(sampleMap.getVCFIndexForSampleAsPath("Sample3"));
        Assert.assertNull(sampleMap.getVCFIndexForSample("Sample2"));
        Assert.assertNull(sampleMap.getVCFIndexForSampleAsPath("Sample2"));
    }

    @DataProvider
    public Object[][] badInputsToAddSample() {
        return new Object[][] {
                { " Sample1", "vcf1" },
                { "Sample1 ", "vcf1" },
                { " Sample1 ", "vcf1" },
                { "", "vcf1" },
                { " ", "vcf1" },
                { null, "vcf1" },
                { "Sample1", null}
        };
    }

    @Test(dataProvider = "badInputsToAddSample", expectedExceptions = UserException.BadInput.class)
    public void testBadInputToAddSample(final String sampleName, final String vcf) throws URISyntaxException {
        final URI vcfURI = vcf != null ? new URI(vcf) : null;
        final SampleNameMap sampleMap = new SampleNameMap();
        sampleMap.addSample(sampleName, vcfURI);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testAddDuplicateSample() throws URISyntaxException {
        final SampleNameMap sampleMap = new SampleNameMap();
        sampleMap.addSample("Sample1", new URI("vcf1"));
        sampleMap.addSample("Sample1", new URI("vcf1alt"));
    }

    @Test(expectedExceptions = UserException.class)
    public void testCheckVcfIsCompressedAndIndexed() {
        final File sampleFile = IOUtils.writeTempFile(ORDERED_SAMPLE_MAP, "goodSampleMapping", ".txt");
        final SampleNameMap sampleMap = new SampleNameMap(sampleFile.toPath(), true);
    }
}
