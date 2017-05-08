package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.*;
import java.util.*;

public final class PSTaxonomyBuilderUtilsTest extends BaseTest {

    @Test
    public void testParseReferenceRecords() {

        final Map<String, PSTaxInfo> taxToInfo = new HashMap<>();

        final List<SAMSequenceRecord> dictList = new ArrayList<>();

        dictList.add(new SAMSequenceRecord("record1|test", 500));
        dictList.add(new SAMSequenceRecord("record2|taxid|1", 5000));
        dictList.add(new SAMSequenceRecord("record3|taxid|2|", 1000));
        dictList.add(new SAMSequenceRecord("record4|ref|NC_1", 2000));
        dictList.add(new SAMSequenceRecord("record5|ref|NC_2|taxid|1", 1000));

        Map<String, Tuple2<String, Long>> result = PSTaxonomyBuilderUtils.parseReferenceRecords(dictList, taxToInfo);

        Assert.assertNotNull(result);
        Assert.assertEquals(result.size(), 2);
        Assert.assertTrue(result.containsKey("record1"));
        Assert.assertEquals(result.get("record1")._1, "record1|test");
        Assert.assertEquals((long) result.get("record1")._2, 500);
        Assert.assertTrue(result.containsKey("NC_1"));
        Assert.assertEquals(result.get("NC_1")._1, "record4|ref|NC_1");
        Assert.assertEquals((long) result.get("NC_1")._2, 2000);

        Assert.assertEquals(taxToInfo.size(), 2);
        Assert.assertTrue(taxToInfo.containsKey("1"));
        Assert.assertEquals(taxToInfo.get("1").ref_names.size(), 2);
        Assert.assertTrue(taxToInfo.get("1").ref_names.contains("record2|taxid|1"));
        Assert.assertTrue(taxToInfo.get("1").ref_names.contains("record5|ref|NC_2|taxid|1"));
        Assert.assertEquals(taxToInfo.get("1").length, 6000);
        Assert.assertTrue(taxToInfo.containsKey("2"));
        Assert.assertEquals(taxToInfo.get("2").ref_names.size(), 1);
        Assert.assertTrue(taxToInfo.get("2").ref_names.contains("record3|taxid|2|"));
        Assert.assertEquals(taxToInfo.get("2").length, 1000);
    }

    @Test
    public void testParseCatalog() {


        final String input = "2\tx\tacc_A\tx\tx\tx\tx\n"
                + "3\tx\tacc_B\tx\tx\tx\tx\n";
        final BufferedReader reader = new BufferedReader(new StringReader(input));

        final Map<String, Tuple2<String, Long>> accToRefInfo = new HashMap<>();
        accToRefInfo.put("acc_A", new Tuple2<>("ref|acc_A", 2000L));
        accToRefInfo.put("acc_C", new Tuple2<>("ref|acc_C", 3000L));

        final Map<String, PSTaxInfo> taxToInfo = new HashMap<>();
        final PSTaxInfo taxInfo = new PSTaxInfo();
        taxInfo.ref_names.add("ref|acc_B|taxid|3");
        taxInfo.length = 1000;
        taxToInfo.put("3", taxInfo);

        final Collection<String> result = PSTaxonomyBuilderUtils.parseCatalog(reader, accToRefInfo, taxToInfo, false, null);

        //Accessions not found in the catalog
        Assert.assertNotNull(result);
        Assert.assertEquals(result.size(), 1);
        final Iterator<String> resultItr = result.iterator();
        Assert.assertEquals(resultItr.next(), "acc_C");

        //Information expected to be added to taxToInfo
        final PSTaxInfo expectedInfo2 = new PSTaxInfo();
        expectedInfo2.length = 2000L;
        expectedInfo2.ref_names.add("ref|acc_A");

        Assert.assertEquals(taxToInfo.size(), 2);
        Assert.assertTrue(taxToInfo.containsKey("2"));
        Assert.assertEquals(taxToInfo.get("2"), expectedInfo2);

        //Throw exception when the input has wrong number of columns
        try {
            final String badInput = "2\tx\n";
            final BufferedReader badReader = new BufferedReader(new StringReader(badInput));
            PSTaxonomyBuilderUtils.parseCatalog(badReader, accToRefInfo, taxToInfo, false, null);
            Assert.fail("Did not throw an exception for bad input (too few columns)");
        } catch (Exception e) {
        }
    }

    @Test
    public void testParseNamesFile() {

        final String input = "1\t|\tname A\t|\t-\t|\taccording to xyz\t|\n"
                + "1\t|\tname A\t|\t-\t|\tscientific name\t|\n"
                + "2\t|\tname B\t|\t-\t|\tscientific name\t|\n"
                + "3\t|\tname C\t|\t-\t|\tscientific name\t|\n";
        final BufferedReader reader = new BufferedReader(new StringReader(input));

        final Map<String, PSTaxInfo> taxToInfo = new HashMap<>();
        taxToInfo.put("1", new PSTaxInfo());
        taxToInfo.put("2", new PSTaxInfo());
        taxToInfo.put("4", new PSTaxInfo());

        PSTaxonomyBuilderUtils.parseNamesFile(reader, taxToInfo);

        Assert.assertEquals(taxToInfo.get("1").name, "name A");
        Assert.assertEquals(taxToInfo.get("2").name, "name B");
        Assert.assertEquals(taxToInfo.get("3").name, "name C");
        Assert.assertEquals(taxToInfo.get("4").name, null);

        //Throw exception when the input has wrong number of columns
        try {
            final String badInput = input + "1\t|\tname A\t|\t-\t|\n";
            final BufferedReader badReader = new BufferedReader(new StringReader(badInput));
            PSTaxonomyBuilderUtils.parseNamesFile(badReader, taxToInfo);
            Assert.fail("Did not throw exception with bad input format (missing column)");
        } catch (Exception e) {
        }
    }

    @Test
    public void testParseNodesFile() {

        final String input = "1\t|\t0\t|\troot\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\n"
                + "2\t|\t1\t|\tkingdom\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\n"
                + "3\t|\t1\t|\tkingdom\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\n"
                + "4\t|\t1\t|\tkingdom\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\t-\t|\n";
        final BufferedReader reader = new BufferedReader(new StringReader(input));

        final Map<String, PSTaxInfo> taxToInfo = new HashMap<>();

        final PSTaxInfo taxInfoRoot = new PSTaxInfo();
        taxInfoRoot.name = "root";
        taxToInfo.put("1", taxInfoRoot);

        final PSTaxInfo taxInfoA = new PSTaxInfo();
        taxInfoA.ref_names.add("ref|acc_A|taxid|2");
        taxToInfo.put("2", taxInfoA);

        final PSTaxInfo taxInfoB = new PSTaxInfo();
        taxInfoB.ref_names.add("ref|acc_B|taxid|3");
        taxToInfo.put("3", taxInfoB);

        final PSTaxInfo taxInfoC = new PSTaxInfo();
        taxInfoC.ref_names.add("ref|acc_C|taxid|5");
        taxToInfo.put("5", taxInfoC);

        final Collection<String> result = PSTaxonomyBuilderUtils.parseNodesFile(reader, taxToInfo);

        taxInfoRoot.parent_tax = null;
        taxInfoRoot.rank = "root";
        Assert.assertEquals(taxToInfo.get("1"), taxInfoRoot);

        taxInfoA.parent_tax = "1";
        taxInfoA.rank = "kingdom";
        Assert.assertEquals(taxToInfo.get("2"), taxInfoA);

        taxInfoB.parent_tax = "1";
        taxInfoB.rank = "kingdom";
        Assert.assertEquals(taxToInfo.get("3"), taxInfoB);

        Assert.assertNotNull(result);
        Assert.assertEquals(result.size(), 1);
        final Iterator<String> resultItr = result.iterator();
        Assert.assertEquals(resultItr.next(), "4");

        //Throw exception when the input has wrong number of columns
        try {
            final String badInput = input + "5\t|\t1\t|\tkingdom\t|\t-\t|\t-\t|\t-\t|\t-\t|\n";
            final BufferedReader badReader = new BufferedReader(new StringReader(badInput));
            PSTaxonomyBuilderUtils.parseNodesFile(badReader, taxToInfo);
            Assert.fail("Did not throw exception with bad input format (missing columns)");
        } catch (Exception e) {
        }
    }

    @Test
    public void testBuildReferenceNameToTaxMap() {
        final Map<String, PSTaxInfo> taxToInfo = new HashMap<>();

        //Empty input
        Map<String, String> result = PSTaxonomyBuilderUtils.buildReferenceNameToTaxMap(taxToInfo);
        Assert.assertTrue(result.isEmpty());
        Assert.assertTrue(taxToInfo.isEmpty());

        //Main test case
        final PSTaxInfo taxInfoRoot = new PSTaxInfo();
        taxInfoRoot.name = "root";
        taxToInfo.put("1", taxInfoRoot);

        final PSTaxInfo taxInfoA = new PSTaxInfo();
        taxInfoA.ref_names.add("ref|acc_A|taxid|2");
        taxInfoA.ref_names.add("ref|acc_B|taxid|2");
        taxToInfo.put("2", taxInfoA);

        final PSTaxInfo taxInfoB = new PSTaxInfo();
        taxInfoB.ref_names.add("ref|acc_C|taxid|3");
        taxToInfo.put("3", taxInfoB);

        result = PSTaxonomyBuilderUtils.buildReferenceNameToTaxMap(taxToInfo);

        Assert.assertEquals(result.size(), 3);
        Assert.assertTrue(result.containsKey("ref|acc_A|taxid|2"));
        Assert.assertEquals(result.get("ref|acc_A|taxid|2"), "2");
        Assert.assertTrue(result.containsKey("ref|acc_B|taxid|2"));
        Assert.assertEquals(result.get("ref|acc_B|taxid|2"), "2");
        Assert.assertTrue(result.containsKey("ref|acc_C|taxid|3"));
        Assert.assertEquals(result.get("ref|acc_C|taxid|3"), "3");
    }

    @Test
    public void testBuildTaxonomicTree() {
        final Map<String, PSTaxInfo> taxToInfo = new HashMap<>();

        //Empty input
        try {
            PSTaxonomyBuilderUtils.buildTaxonomicTree(taxToInfo);
            Assert.fail("Did not throw exception when the resulting tree is empty");
        } catch (Exception e) {
        }

        //Tree that reduces to an empty tree because no nodes are associated with references
        final PSTaxInfo taxInfoRoot = new PSTaxInfo();
        taxToInfo.put("1", taxInfoRoot);

        //Note taxInfoC.ref_names is empty
        final PSTaxInfo taxInfoC = new PSTaxInfo();
        taxInfoC.parent_tax = "1";
        taxInfoC.name = "species Y";
        taxInfoC.rank = "species";
        taxToInfo.put("4", taxInfoC);

        final PSTaxInfo taxInfoD = new PSTaxInfo();
        taxToInfo.put("5", taxInfoD);

        try {
            PSTaxonomyBuilderUtils.buildTaxonomicTree(taxToInfo);
            Assert.fail("Did not throw exception when the resulting tree is empty");
        } catch (Exception e) {
        }

        //Main test: add nodes associated with references
        final PSTaxInfo taxInfoA = new PSTaxInfo();
        taxInfoA.parent_tax = "1";
        taxInfoA.ref_names.add("ref|acc_A|taxid|2");
        taxInfoA.ref_names.add("ref|acc_B|taxid|2");
        taxInfoA.name = "species X";
        taxInfoA.rank = "species";
        taxToInfo.put("2", taxInfoA);

        final PSTaxInfo taxInfoB = new PSTaxInfo();
        taxInfoB.parent_tax = "1";
        taxInfoB.ref_names.add("ref|acc_C|taxid|3");
        taxInfoB.name = "species Z";
        taxInfoB.rank = "species";
        taxToInfo.put("3", taxInfoB);

        final PSTree tree = PSTaxonomyBuilderUtils.buildTaxonomicTree(taxToInfo);

        Assert.assertEquals(tree.getNodeIDs().size(), 3);
        Assert.assertTrue(tree.hasNode("1"));
        Assert.assertTrue(tree.hasNode("2"));
        Assert.assertTrue(tree.hasNode("3"));
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testGetFileReaderGz() throws Exception {
        final File testFile = getTestFile("test_gz.txt.gz");
        final BufferedReader bufferedReader = PSTaxonomyBuilderUtils.getBufferedReaderGz(testFile.getPath());
        final String expectedString = "you read my test file.";
        final String firstLine = bufferedReader.readLine();
        Assert.assertEquals(firstLine, expectedString);
        Assert.assertFalse(bufferedReader.ready());

        try {
            PSTaxonomyBuilderUtils.getBufferedReaderGz("non_existent.gz");
            Assert.fail("Did not throw UserException for bad file name");
        } catch (UserException e) {
        }
    }

    @Test
    public void testGetTarGzReader() throws Exception {
        final File testFile = getTestFile("test.tar.gz");
        final BufferedReader bufferedReader = PSTaxonomyBuilderUtils.getBufferedReaderTarGz(testFile.getPath(), "test_gz.txt");
        final String expectedString = "you read my test file.";
        final String firstLine = bufferedReader.readLine();
        Assert.assertEquals(firstLine, expectedString);
        Assert.assertFalse(bufferedReader.ready());

        try {
            PSTaxonomyBuilderUtils.getBufferedReaderTarGz("non_existent.tar.gz", "invalid_file.txt");
            Assert.fail("Did not throw UserException for bad archive name");
        } catch (UserException e) {
        }

        try {
            PSTaxonomyBuilderUtils.getBufferedReaderTarGz(testFile.getPath(), "invalid_file.txt");
            Assert.fail("Did not throw UserException for bad file name");
        } catch (UserException e) {
        }
    }

    @Test
    public void testCloseReader() {
        final BufferedReader r;
        try {
            r = new BufferedReader(new FileReader(hg19MiniReference));
        } catch (IOException e) {
            throw new TestException(e);
        }
        PSTaxonomyBuilderUtils.closeReader(r);
        try {
            r.ready();
            Assert.fail("Reader stream still open");
        } catch (IOException e) {
        }
    }
}