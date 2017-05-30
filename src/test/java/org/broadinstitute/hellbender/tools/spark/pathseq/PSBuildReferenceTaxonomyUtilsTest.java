package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.*;
import java.util.*;

public final class PSBuildReferenceTaxonomyUtilsTest extends BaseTest {

    @Test
    public void testParseReferenceRecords() {

        final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();

        final List<SAMSequenceRecord> dictList = new ArrayList<>();

        dictList.add(new SAMSequenceRecord("record1|test", 500));
        dictList.add(new SAMSequenceRecord("record2|taxid|1", 5000));
        dictList.add(new SAMSequenceRecord("record3|taxid|2|", 1000));
        dictList.add(new SAMSequenceRecord("record4|ref|NC_1", 2000));
        dictList.add(new SAMSequenceRecord("record5|ref|NC_2|taxid|1", 1000));

        Map<String, Tuple2<String, Long>> result = PSBuildReferenceTaxonomyUtils.parseReferenceRecords(dictList, taxIdToProperties);

        Assert.assertNotNull(result);
        Assert.assertEquals(result.size(), 2);
        Assert.assertTrue(result.containsKey("record1"));
        Assert.assertEquals(result.get("record1")._1, "record1|test");
        Assert.assertEquals((long) result.get("record1")._2, 500);
        Assert.assertTrue(result.containsKey("NC_1"));
        Assert.assertEquals(result.get("NC_1")._1, "record4|ref|NC_1");
        Assert.assertEquals((long) result.get("NC_1")._2, 2000);

        Assert.assertEquals(taxIdToProperties.size(), 2);
        Assert.assertTrue(taxIdToProperties.containsKey("1"));
        Assert.assertEquals(taxIdToProperties.get("1").accessions.size(), 2);
        Assert.assertTrue(taxIdToProperties.get("1").accessions.contains("record2|taxid|1"));
        Assert.assertTrue(taxIdToProperties.get("1").accessions.contains("record5|ref|NC_2|taxid|1"));
        Assert.assertEquals(taxIdToProperties.get("1").length, 6000);
        Assert.assertTrue(taxIdToProperties.containsKey("2"));
        Assert.assertEquals(taxIdToProperties.get("2").accessions.size(), 1);
        Assert.assertTrue(taxIdToProperties.get("2").accessions.contains("record3|taxid|2|"));
        Assert.assertEquals(taxIdToProperties.get("2").length, 1000);
    }

    @Test(expectedExceptions = Exception.class)
    public void testParseCatalog() {

        final String input =    "2\tx\tacc_A\tx\tx\tx\tx\n"
                            +   "3\tx\tacc_B\tx\tx\tx\tx\n";
        final BufferedReader reader = new BufferedReader(new StringReader(input));

        final Map<String, Tuple2<String, Long>> accessionToNameAndLength = new HashMap<>();
        accessionToNameAndLength.put("acc_A", new Tuple2<>("ref|acc_A", 2000L));
        accessionToNameAndLength.put("acc_C", new Tuple2<>("ref|acc_C", 3000L));

        final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();
        final PSPathogenReferenceTaxonProperties taxonProperties = new PSPathogenReferenceTaxonProperties();
        taxonProperties.accessions.add("ref|acc_B|taxid|3");
        taxonProperties.length = 1000;
        taxIdToProperties.put("3", taxonProperties);

        final Collection<String> result = PSBuildReferenceTaxonomyUtils.parseCatalog(reader, accessionToNameAndLength, taxIdToProperties, false, null);

        //Accessions not found in the catalog
        Assert.assertNotNull(result);
        Assert.assertEquals(result.size(), 1);
        final Iterator<String> resultIter = result.iterator();
        Assert.assertEquals(resultIter.next(), "acc_C");

        //Information expected to be added to taxIdToProperties
        final PSPathogenReferenceTaxonProperties expectedProperties = new PSPathogenReferenceTaxonProperties();
        expectedProperties.length = 2000L;
        expectedProperties.accessions.add("ref|acc_A");

        Assert.assertEquals(taxIdToProperties.size(), 2);
        Assert.assertTrue(taxIdToProperties.containsKey("2"));
        Assert.assertEquals(taxIdToProperties.get("2").length, expectedProperties.length);
        Assert.assertEquals(taxIdToProperties.get("2").name, expectedProperties.name);
        Assert.assertEquals(taxIdToProperties.get("2").rank, expectedProperties.rank);
        Assert.assertEquals(taxIdToProperties.get("2").accessions.size(), expectedProperties.accessions.size());
        Assert.assertTrue(taxIdToProperties.get("2").accessions.containsAll(expectedProperties.accessions));

        //Test bad catalog input
        final String badInput = "2\tx\n";
        final BufferedReader badReader = new BufferedReader(new StringReader(badInput));
        PSBuildReferenceTaxonomyUtils.parseCatalog(badReader, accessionToNameAndLength, taxIdToProperties, false, null);
    }

    @Test(expectedExceptions = Exception.class)
    public void testParseNamesFile() {

        final String input = "1\t|\tname A\t|\t-\t|\taccording to xyz\t|\n"
                + "1\t|\tname A\t|\t-\t|\tscientific name\t|\n"
                + "2\t|\tname B\t|\t-\t|\tscientific name\t|\n"
                + "3\t|\tname C\t|\t-\t|\tscientific name\t|\n";
        final BufferedReader reader = new BufferedReader(new StringReader(input));

        final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();
        taxIdToProperties.put("1", new PSPathogenReferenceTaxonProperties());
        taxIdToProperties.put("2", new PSPathogenReferenceTaxonProperties());
        taxIdToProperties.put("4", new PSPathogenReferenceTaxonProperties());

        PSBuildReferenceTaxonomyUtils.parseNamesFile(reader, taxIdToProperties);

        Assert.assertEquals(taxIdToProperties.get("1").name, "name A");
        Assert.assertEquals(taxIdToProperties.get("2").name, "name B");
        Assert.assertEquals(taxIdToProperties.get("3").name, "name C");
        Assert.assertEquals(taxIdToProperties.get("4").name, null);

        //Throw exception when the input has wrong number of columns
        final String badInput = input + "1\t|\tname A\t|\t-\t|\n";
        final BufferedReader badReader = new BufferedReader(new StringReader(badInput));
        PSBuildReferenceTaxonomyUtils.parseNamesFile(badReader, taxIdToProperties);
    }

    @Test(expectedExceptions = Exception.class)
    public void testParseNodesFile() {

        final String input = "1\t|\t0\t|\troot\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\n"
                + "2\t|\t1\t|\tkingdom\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\n"
                + "3\t|\t1\t|\tkingdom\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\n"
                + "4\t|\t1\t|\tkingdom\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\n";
        final BufferedReader reader = new BufferedReader(new StringReader(input));

        final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();

        final PSPathogenReferenceTaxonProperties rootProperties = new PSPathogenReferenceTaxonProperties();
        rootProperties.name = "root";
        taxIdToProperties.put("1", rootProperties);

        final PSPathogenReferenceTaxonProperties taxPropertiesA = new PSPathogenReferenceTaxonProperties();
        taxPropertiesA.accessions.add("ref|acc_A|taxid|2");
        taxIdToProperties.put("2", taxPropertiesA);

        final PSPathogenReferenceTaxonProperties taxPropertiesB = new PSPathogenReferenceTaxonProperties();
        taxPropertiesB.accessions.add("ref|acc_B|taxid|3");
        taxIdToProperties.put("3", taxPropertiesB);

        final PSPathogenReferenceTaxonProperties taxPropertiesC = new PSPathogenReferenceTaxonProperties();
        taxPropertiesC.accessions.add("ref|acc_C|taxid|5");
        taxIdToProperties.put("5", taxPropertiesC);

        final Collection<String> result = PSBuildReferenceTaxonomyUtils.parseNodesFile(reader, taxIdToProperties);

        rootProperties.parentTaxId = null;
        rootProperties.rank = "root";
        Assert.assertEquals(taxIdToProperties.get("1"), rootProperties);

        taxPropertiesA.parentTaxId = "1";
        taxPropertiesA.rank = "kingdom";
        Assert.assertEquals(taxIdToProperties.get("2"), taxPropertiesA);

        taxPropertiesB.parentTaxId = "1";
        taxPropertiesB.rank = "kingdom";
        Assert.assertEquals(taxIdToProperties.get("3"), taxPropertiesB);

        Assert.assertNotNull(result);
        Assert.assertEquals(result.size(), 1);
        final Iterator<String> resultIter = result.iterator();
        Assert.assertEquals(resultIter.next(), "4");

        //Throw exception when the input has wrong number of columns
        final String badInput = input + "5\t|\t1\t|\tkingdom\t|\t-\t|\t-\t|\t-\t|\t-\t|\n";
        final BufferedReader badReader = new BufferedReader(new StringReader(badInput));
        PSBuildReferenceTaxonomyUtils.parseNodesFile(badReader, taxIdToProperties);
    }

    @Test
    public void testBuildReferenceNameToTaxMap() {
        final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();

        //Empty input
        Map<String, String> result = PSBuildReferenceTaxonomyUtils.buildAccessionToTaxIdMap(taxIdToProperties);
        Assert.assertTrue(result.isEmpty());
        Assert.assertTrue(taxIdToProperties.isEmpty());

        //Main test case
        final PSPathogenReferenceTaxonProperties taxPropertiesRoot = new PSPathogenReferenceTaxonProperties();
        taxPropertiesRoot.name = "root";
        taxIdToProperties.put("1", taxPropertiesRoot);

        final PSPathogenReferenceTaxonProperties taxPropertiesA = new PSPathogenReferenceTaxonProperties();
        taxPropertiesA.accessions.add("ref|acc_A|taxid|2");
        taxPropertiesA.accessions.add("ref|acc_B|taxid|2");
        taxIdToProperties.put("2", taxPropertiesA);

        final PSPathogenReferenceTaxonProperties taxPropertiesB = new PSPathogenReferenceTaxonProperties();
        taxPropertiesB.accessions.add("ref|acc_C|taxid|3");
        taxIdToProperties.put("3", taxPropertiesB);

        result = PSBuildReferenceTaxonomyUtils.buildAccessionToTaxIdMap(taxIdToProperties);

        Assert.assertEquals(result.size(), 3);
        Assert.assertTrue(result.containsKey("ref|acc_A|taxid|2"));
        Assert.assertEquals(result.get("ref|acc_A|taxid|2"), "2");
        Assert.assertTrue(result.containsKey("ref|acc_B|taxid|2"));
        Assert.assertEquals(result.get("ref|acc_B|taxid|2"), "2");
        Assert.assertTrue(result.containsKey("ref|acc_C|taxid|3"));
        Assert.assertEquals(result.get("ref|acc_C|taxid|3"), "3");
    }

    //Empty input should throw exception
    @Test(expectedExceptions = Exception.class)
    public void testEmptyBuildTaxonomicTree() {
        final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();

        PSBuildReferenceTaxonomyUtils.buildTaxonomicTree(taxIdToProperties);
    }

    @Test
    public void testBuildTaxonomicTree() {
        final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();

        //Tree that reduces to an empty tree because no nodes are associated with references
        final PSPathogenReferenceTaxonProperties taxPropertiesRoot = new PSPathogenReferenceTaxonProperties();
        taxIdToProperties.put("1", taxPropertiesRoot);

        //Note taxPropertiesC.accessions is empty
        final PSPathogenReferenceTaxonProperties taxPropertiesC = new PSPathogenReferenceTaxonProperties();
        taxPropertiesC.parentTaxId = "1";
        taxPropertiesC.name = "species Y";
        taxPropertiesC.rank = "species";
        taxIdToProperties.put("4", taxPropertiesC);

        final PSPathogenReferenceTaxonProperties taxPropertiesD = new PSPathogenReferenceTaxonProperties();
        taxIdToProperties.put("5", taxPropertiesD);

        //Main test: add nodes associated with references
        final PSPathogenReferenceTaxonProperties taxPropertiesA = new PSPathogenReferenceTaxonProperties();
        taxPropertiesA.parentTaxId = "1";
        taxPropertiesA.accessions.add("ref|acc_A|taxid|2");
        taxPropertiesA.accessions.add("ref|acc_B|taxid|2");
        taxPropertiesA.name = "species X";
        taxPropertiesA.rank = "species";
        taxIdToProperties.put("2", taxPropertiesA);

        final PSPathogenReferenceTaxonProperties taxPropertiesB = new PSPathogenReferenceTaxonProperties();
        taxPropertiesB.parentTaxId = "1";
        taxPropertiesB.accessions.add("ref|acc_C|taxid|3");
        taxPropertiesB.name = "species Z";
        taxPropertiesB.rank = "species";
        taxIdToProperties.put("3", taxPropertiesB);

        final PSTree tree = PSBuildReferenceTaxonomyUtils.buildTaxonomicTree(taxIdToProperties);

        Assert.assertEquals(tree.getNodeIDs().size(), 3);
        Assert.assertTrue(tree.hasNode("1"));
        Assert.assertTrue(tree.hasNode("2"));
        Assert.assertTrue(tree.hasNode("3"));
    }

    @Test(expectedExceptions = Exception.class)
    public void testBuildEmptyTaxonomicTree() {
        final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();

        //Tree that reduces to an empty tree because no nodes are associated with references
        final PSPathogenReferenceTaxonProperties taxPropertiesRoot = new PSPathogenReferenceTaxonProperties();
        taxIdToProperties.put("1", taxPropertiesRoot);

        final PSPathogenReferenceTaxonProperties taxPropertiesC = new PSPathogenReferenceTaxonProperties();
        taxPropertiesC.parentTaxId = "1";
        taxPropertiesC.name = "species Y";
        taxPropertiesC.rank = "species";
        taxIdToProperties.put("4", taxPropertiesC);

        final PSPathogenReferenceTaxonProperties taxPropertiesD = new PSPathogenReferenceTaxonProperties();
        taxIdToProperties.put("5", taxPropertiesD);

        //Note taxPropertiesC.accessions and taxPropertiesD.accessions are empty
        PSBuildReferenceTaxonomyUtils.buildTaxonomicTree(taxIdToProperties);
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testGetFileReaderGz() throws Exception {
        final File testFile = getTestFile("test_gz.txt.gz");
        final BufferedReader bufferedReader = PSBuildReferenceTaxonomyUtils.getBufferedReaderGz(testFile.getPath());
        final String expectedString = "you read my test file.";
        final String firstLine = bufferedReader.readLine();
        Assert.assertEquals(firstLine, expectedString);
        Assert.assertFalse(bufferedReader.ready());
    }

    @Test(expectedExceptions = Exception.class)
    public void testGetTarGzReaderBadArchivedFile() throws Exception {
        final File testFile = getTestFile("test_gz.txt.gz");
        PSBuildReferenceTaxonomyUtils.getBufferedReaderTarGz(testFile.getPath(), "invalid_file.txt");
    }

    @Test(expectedExceptions = Exception.class)
    public void testGetTarGzReaderBadTar() throws Exception {
        PSBuildReferenceTaxonomyUtils.getBufferedReaderTarGz(getSafeNonExistentFile("bad.tar.gz").getAbsolutePath(), "invalid_file.txt");
    }

    @Test
    public void testCloseReader() {
        final BufferedReader r;
        try {
            r = new BufferedReader(new FileReader(hg19MiniReference));
        } catch (IOException e) {
            throw new TestException(e);
        }
    }
}