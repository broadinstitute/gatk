package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.*;
import java.util.*;

public final class PSBuildReferenceTaxonomyUtilsTest extends GATKBaseTest {

    @Test
    public void testParseReferenceRecords() {

        final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();

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
        Assert.assertTrue(taxIdToProperties.containsKey(1));
        Assert.assertEquals(taxIdToProperties.get(1).getAccessions().size(), 2);
        Assert.assertTrue(taxIdToProperties.get(1).getAccessions().contains("record2|taxid|1"));
        Assert.assertTrue(taxIdToProperties.get(1).getAccessions().contains("record5|ref|NC_2|taxid|1"));
        Assert.assertEquals(taxIdToProperties.get(1).getTotalLength(), 6000);
        Assert.assertTrue(taxIdToProperties.containsKey(2));
        Assert.assertEquals(taxIdToProperties.get(2).getAccessions().size(), 1);
        Assert.assertTrue(taxIdToProperties.get(2).getAccessions().contains("record3|taxid|2|"));
        Assert.assertEquals(taxIdToProperties.get(2).getTotalLength(), 1000);
    }

    @Test(expectedExceptions = Exception.class)
    public void testParseCatalog() {

        final String input =    "2\tx\tacc_A\tx\tx\tx\tx\n"
                            +   "3\tx\tacc_B\tx\tx\tx\tx\n";
        final BufferedReader reader = new BufferedReader(new StringReader(input));

        final Map<String, Tuple2<String, Long>> accessionToNameAndLength = new HashMap<>();
        accessionToNameAndLength.put("acc_A", new Tuple2<>("ref|acc_A", 2000L));
        accessionToNameAndLength.put("acc_C", new Tuple2<>("ref|acc_C", 3000L));

        final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();
        final PSPathogenReferenceTaxonProperties taxonProperties = new PSPathogenReferenceTaxonProperties();
        taxonProperties.addAccession("ref|acc_B|taxid|3", 1000);
        taxIdToProperties.put(3, taxonProperties);

        final Collection<String> result = PSBuildReferenceTaxonomyUtils.parseCatalog(reader, accessionToNameAndLength, taxIdToProperties, false, null);

        //Accessions not found in the catalog
        Assert.assertNotNull(result);
        Assert.assertEquals(result.size(), 1);
        final Iterator<String> resultIter = result.iterator();
        Assert.assertEquals(resultIter.next(), "acc_C");

        //Information expected to be added to taxIdToProperties
        final PSPathogenReferenceTaxonProperties expectedProperties = new PSPathogenReferenceTaxonProperties();
        expectedProperties.addAccession("ref|acc_A", 2000);

        Assert.assertEquals(taxIdToProperties.size(), 2);
        Assert.assertTrue(taxIdToProperties.containsKey(2));
        Assert.assertEquals(taxIdToProperties.get(2).getTotalLength(), expectedProperties.getTotalLength());
        Assert.assertEquals(taxIdToProperties.get(2).getName(), expectedProperties.getName());
        Assert.assertEquals(taxIdToProperties.get(2).getRank(), expectedProperties.getRank());
        Assert.assertEquals(taxIdToProperties.get(2).getAccessions().size(), expectedProperties.getAccessions().size());
        Assert.assertTrue(taxIdToProperties.get(2).getAccessions().containsAll(expectedProperties.getAccessions()));

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

        final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();
        taxIdToProperties.put(1, new PSPathogenReferenceTaxonProperties());
        taxIdToProperties.put(2, new PSPathogenReferenceTaxonProperties());
        taxIdToProperties.put(4, new PSPathogenReferenceTaxonProperties());

        PSBuildReferenceTaxonomyUtils.parseNamesFile(reader, taxIdToProperties);

        Assert.assertEquals(taxIdToProperties.get(1).getName(), "name A");
        Assert.assertEquals(taxIdToProperties.get(2).getName(), "name B");
        Assert.assertEquals(taxIdToProperties.get(3).getName(), "name C");
        Assert.assertEquals(taxIdToProperties.get(4).getName(), null);

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

        final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();

        final PSPathogenReferenceTaxonProperties rootProperties = new PSPathogenReferenceTaxonProperties("root");
        taxIdToProperties.put(1, rootProperties);

        final PSPathogenReferenceTaxonProperties taxPropertiesA = new PSPathogenReferenceTaxonProperties();
        taxPropertiesA.addAccession("ref|acc_A|taxid|2", 100);
        taxIdToProperties.put(2, taxPropertiesA);

        final PSPathogenReferenceTaxonProperties taxPropertiesB = new PSPathogenReferenceTaxonProperties();
        taxPropertiesB.addAccession("ref|acc_B|taxid|3", 100);
        taxIdToProperties.put(3, taxPropertiesB);

        final PSPathogenReferenceTaxonProperties taxPropertiesC = new PSPathogenReferenceTaxonProperties();
        taxPropertiesC.addAccession("ref|acc_C|taxid|5", 100);
        taxIdToProperties.put(5, taxPropertiesC);

        final Collection<Integer> result = PSBuildReferenceTaxonomyUtils.parseNodesFile(reader, taxIdToProperties);

        rootProperties.setRank("root");
        Assert.assertEquals(taxIdToProperties.get(1), rootProperties);

        taxPropertiesA.setParent(1);
        taxPropertiesA.setRank("kingdom");
        Assert.assertEquals(taxIdToProperties.get(2), taxPropertiesA);

        taxPropertiesB.setParent(1);
        taxPropertiesB.setRank("kingdom");
        Assert.assertEquals(taxIdToProperties.get(3), taxPropertiesB);

        Assert.assertNotNull(result);
        Assert.assertEquals(result.size(), 1);
        final Iterator<Integer> resultIter = result.iterator();
        Assert.assertEquals(resultIter.next().intValue(), 4);

        //Throw exception when the input has insufficient number of columns
        final String badInput = input + "5\t|\t1\t|\n";
        final BufferedReader badReader = new BufferedReader(new StringReader(badInput));
        PSBuildReferenceTaxonomyUtils.parseNodesFile(badReader, taxIdToProperties);
    }

    @Test
    public void testBuildReferenceNameToTaxMap() {
        final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();
        final PSTree tree = new PSTree(1);
        tree.addNode(PSTaxonomyConstants.VIRUS_ID,"Virus",1,0,"kingdom");
        tree.addNode(3,"Bacteria",1,0,"kingdom");
        tree.addNode(4,"A", PSTaxonomyConstants.VIRUS_ID,600,"species");
        tree.addNode(5,"B",3,1000,"species");

        //Empty input
        Map<String, Integer> result = PSBuildReferenceTaxonomyUtils.buildAccessionToTaxIdMap(taxIdToProperties, tree, 0);
        Assert.assertTrue(result.isEmpty());
        Assert.assertTrue(taxIdToProperties.isEmpty());

        //Main test case
        final PSPathogenReferenceTaxonProperties taxPropertiesRoot = new PSPathogenReferenceTaxonProperties("root");
        taxIdToProperties.put(1, taxPropertiesRoot);

        final String accA = "ref|acc_A|taxid|4";
        final String accB = "ref|acc_B|taxid|4";
        final String accC = "ref|acc_C|taxid|5";
        final String accD = "ref|acc_D|taxid|5";

        final PSPathogenReferenceTaxonProperties taxPropertiesA = new PSPathogenReferenceTaxonProperties();
        taxPropertiesA.addAccession(accA, 100);
        taxPropertiesA.addAccession(accB, 500);
        taxIdToProperties.put(4, taxPropertiesA);

        final PSPathogenReferenceTaxonProperties taxPropertiesB = new PSPathogenReferenceTaxonProperties();
        taxPropertiesB.addAccession(accC, 1000);
        taxPropertiesB.addAccession(accD, 100);
        taxIdToProperties.put(5, taxPropertiesB);

        result = PSBuildReferenceTaxonomyUtils.buildAccessionToTaxIdMap(taxIdToProperties, tree, 100);

        Assert.assertEquals(result.size(), 4);
        Assert.assertTrue(result.containsKey(accA));
        Assert.assertEquals(result.get(accA).intValue(), 4);
        Assert.assertTrue(result.containsKey(accB));
        Assert.assertEquals(result.get(accB).intValue(), 4);
        Assert.assertTrue(result.containsKey(accC));
        Assert.assertEquals(result.get(accC).intValue(), 5);
        Assert.assertTrue(result.containsKey(accD));
        Assert.assertEquals(result.get(accD).intValue(), 5);

        result = PSBuildReferenceTaxonomyUtils.buildAccessionToTaxIdMap(taxIdToProperties, tree, 101);

        Assert.assertEquals(result.size(), 3);
        Assert.assertTrue(result.containsKey(accA));
        Assert.assertEquals(result.get(accA).intValue(), 4);
        Assert.assertTrue(result.containsKey(accB));
        Assert.assertEquals(result.get(accB).intValue(), 4);
        Assert.assertTrue(result.containsKey(accC));
        Assert.assertEquals(result.get(accC).intValue(), 5);
    }

    //Empty input should throw exception
    @Test(expectedExceptions = Exception.class)
    public void testEmptyBuildTaxonomicTree() {
        final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();

        PSBuildReferenceTaxonomyUtils.buildTaxonomicTree(taxIdToProperties);
    }

    @Test
    public void testBuildTaxonomicTree() {
        final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();

        //Tree that reduces to an empty tree because no nodes are associated with references
        final PSPathogenReferenceTaxonProperties taxPropertiesRoot = new PSPathogenReferenceTaxonProperties();
        taxIdToProperties.put(1, taxPropertiesRoot);

        //Note taxPropertiesC.accessions is empty
        final PSPathogenReferenceTaxonProperties taxPropertiesC = new PSPathogenReferenceTaxonProperties();
        taxPropertiesC.setParent(1);
        taxPropertiesC.setName("species Y");
        taxPropertiesC.setRank("species");
        taxIdToProperties.put(4, taxPropertiesC);

        final PSPathogenReferenceTaxonProperties taxPropertiesD = new PSPathogenReferenceTaxonProperties();
        taxIdToProperties.put(5, taxPropertiesD);

        //Main test: add nodes associated with references
        final PSPathogenReferenceTaxonProperties taxPropertiesA = new PSPathogenReferenceTaxonProperties();
        taxPropertiesA.setParent(1);
        taxPropertiesA.addAccession("ref|acc_A|taxid|2", 100);
        taxPropertiesA.addAccession("ref|acc_B|taxid|2", 100);
        taxPropertiesA.setName("species X");
        taxPropertiesA.setRank("species");
        taxIdToProperties.put(2, taxPropertiesA);

        final PSPathogenReferenceTaxonProperties taxPropertiesB = new PSPathogenReferenceTaxonProperties();
        taxPropertiesB.setParent(1);
        taxPropertiesB.addAccession("ref|acc_C|taxid|3", 100);
        taxPropertiesB.setName("species Z");
        taxPropertiesB.setRank("species");
        taxIdToProperties.put(3, taxPropertiesB);

        final PSTree tree = PSBuildReferenceTaxonomyUtils.buildTaxonomicTree(taxIdToProperties);

        Assert.assertEquals(tree.getNodeIDs().size(), 3);
        Assert.assertTrue(tree.hasNode(1));
        Assert.assertTrue(tree.hasNode(2));
        Assert.assertTrue(tree.hasNode(3));
    }

    @Test(expectedExceptions = Exception.class)
    public void testBuildEmptyTaxonomicTree() {
        final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();

        //Tree that reduces to an empty tree because no nodes are associated with references
        final PSPathogenReferenceTaxonProperties taxPropertiesRoot = new PSPathogenReferenceTaxonProperties();
        taxIdToProperties.put(1, taxPropertiesRoot);

        final PSPathogenReferenceTaxonProperties taxPropertiesC = new PSPathogenReferenceTaxonProperties();
        taxPropertiesC.setParent(1);
        taxPropertiesC.setName("species Y");
        taxPropertiesC.setRank("species");
        taxIdToProperties.put(4, taxPropertiesC);

        final PSPathogenReferenceTaxonProperties taxPropertiesD = new PSPathogenReferenceTaxonProperties();
        taxIdToProperties.put(5, taxPropertiesD);

        //Note taxPropertiesC.accessions and taxPropertiesD.accessions are empty
        PSBuildReferenceTaxonomyUtils.buildTaxonomicTree(taxIdToProperties);
    }

    @Test
    public void testRemoveUnusedTaxIds() {
        final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();
        taxIdToProperties.put(1, new PSPathogenReferenceTaxonProperties());
        taxIdToProperties.put(2, new PSPathogenReferenceTaxonProperties());
        taxIdToProperties.put(3, new PSPathogenReferenceTaxonProperties());
        taxIdToProperties.put(4, new PSPathogenReferenceTaxonProperties());
        final PSTree tree = new PSTree(1);
        tree.addNode(2, "node2", 1, 0, "child");
        tree.addNode(4, "node2", 1, 0, "child");
        PSBuildReferenceTaxonomyUtils.removeUnusedTaxIds(taxIdToProperties, tree);
        Assert.assertEquals(taxIdToProperties.size(), 3);
        Assert.assertTrue(taxIdToProperties.containsKey(1));
        Assert.assertTrue(taxIdToProperties.containsKey(2));
        Assert.assertTrue(taxIdToProperties.containsKey(4));
        Assert.assertTrue(!taxIdToProperties.containsKey(3));
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