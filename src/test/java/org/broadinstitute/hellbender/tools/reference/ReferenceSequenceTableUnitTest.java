package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class ReferenceSequenceTableUnitTest extends GATKBaseTest {

    private ReferenceSequenceTable tableGenerator(List<GATKPath> references, CompareReferences.MD5CalculationMode md5CalculationMode){
        Map<GATKPath, ReferenceDataSource> referenceSources = new LinkedHashMap<>();
        for(GATKPath path : references){
            referenceSources.put(path, ReferenceDataSource.of(path.toPath()));
        }

        ReferenceSequenceTable table = new ReferenceSequenceTable(referenceSources, md5CalculationMode);
        table.build();

        return table;
    }

    @DataProvider(name = "testQueryByMD5Data")
    public Object[][] testQueryByMD5Data() {
        return new Object[][]{
                // references, md5, length, table
                new Object[]{ "8c0c38e352d8f3309eabe4845456f274",
                        16000,
                        tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta")),
                                CompareReferences.MD5CalculationMode.USE_DICT),
                },
                new Object[]{ "8c0c38e352d8f3309eabe4845456f274",
                        16000,
                        tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_1renamed.fasta")),
                                CompareReferences.MD5CalculationMode.USE_DICT),
                },
        };
    }

    public boolean compareRows(ReferenceSequenceTable.TableRow row1, ReferenceSequenceTable.TableRow row2){
        return row1.getMd5().equals(row2.getMd5()) && row1.getLength() == row2.getLength();
    }

    @Test(dataProvider = "testQueryByMD5Data")
    public void testQueryByMD5(String md5, int length, ReferenceSequenceTable table){
        ReferenceSequenceTable.TableRow row = table.queryByMD5(md5);

        Assert.assertEquals(row.getMd5(), md5);
        Assert.assertEquals(row.getLength(), length);

    }

    @DataProvider(name = "testQueryBySequenceNameData")
    public Object[][] testQueryBySequenceNameData() {
        return new Object[][]{
                // sequence name, expected number of rows, table
                new Object[]{ "2", 2,
                        tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                        new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_chr2snp.fasta")),
                        CompareReferences.MD5CalculationMode.USE_DICT),
                },
                new Object[]{ "1", 1,
                        tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                                        new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_chr2snp.fasta")),
                                CompareReferences.MD5CalculationMode.USE_DICT),
                },
                new Object[]{ "5", 0,
                        tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                                        new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_chr2snp.fasta")),
                                CompareReferences.MD5CalculationMode.USE_DICT),
                },
        };
    }


    @Test(dataProvider = "testQueryBySequenceNameData")
    public void testQueryBySequenceName(String expectedSequenceName, int expectedRows, ReferenceSequenceTable table){
        Set<ReferenceSequenceTable.TableRow> rows = table.queryBySequenceName(expectedSequenceName);

        if(expectedRows == 0){
            Assert.assertTrue(rows.isEmpty());
        } else {
            Assert.assertEquals(rows.size(), expectedRows);

            boolean contigPresent = false;
            for (ReferenceSequenceTable.TableRow row : rows) {
                contigPresent = false;
                for (ReferenceSequenceTable.TableEntry entry : row.getEntries()) {
                    if (entry.getColumnValue().equals(expectedSequenceName)) {
                        contigPresent = true;
                        break;
                    }
                }
                Assert.assertTrue(contigPresent, "Contig " + expectedSequenceName + " not found in all rows returned by queryBySequenceName(expectedSequenceName)");
            }
        }
    }

    @DataProvider(name = "testGenerateReferencePairsData")
    public Object[][] testGenerateReferencePairsData() {
        return new Object[][]{
                // references, expected number of pairs
                new Object[]{ Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                        new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_chr2snp.fasta")),
                        1
                },
                new Object[]{ Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                        new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_chr2snp.fasta"),
                        new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_1renamed.fasta")),
                        3
                },
        };
    }

    @Test(dataProvider = "testGenerateReferencePairsData")
    public void testGenerateReferencePairs(List<GATKPath> references, int expectedPairs){
        ReferenceSequenceTable table = tableGenerator(references, CompareReferences.MD5CalculationMode.USE_DICT);
        List<ReferencePair> pairs = table.generateReferencePairs();
        Assert.assertEquals(pairs.size(), expectedPairs);
    }

    @DataProvider(name = "testAnalyzeTableTwoRefsData")
    public Object[][] testAnalyzeTableTwoReferencesData(){
        return new Object[][]{
                // table, set of expected statuses
                new Object[]{ tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_chr2snp.fasta")),
                        CompareReferences.MD5CalculationMode.USE_DICT),
                        new HashSet<>(Arrays.asList(ReferencePair.Status.DIFFER_IN_SEQUENCE))
                },
                new Object[]{ tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_1renamed.fasta")),
                        CompareReferences.MD5CalculationMode.USE_DICT),
                        new HashSet<>(Arrays.asList(ReferencePair.Status.DIFFER_IN_SEQUENCE_NAMES))
                },
                new Object[]{ tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_missingchr3.fasta")),
                        CompareReferences.MD5CalculationMode.USE_DICT),
                        new HashSet<>(Arrays.asList(ReferencePair.Status.SUPERSET))
                },
                new Object[]{ tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_missingchr3.fasta"),
                                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta")),
                        CompareReferences.MD5CalculationMode.USE_DICT),
                        new HashSet<>(Arrays.asList(ReferencePair.Status.SUBSET))
                },
                new Object[]{ tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_missingchr1.fasta"),
                                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_missingchr3.fasta")),
                        CompareReferences.MD5CalculationMode.USE_DICT),
                        new HashSet<>(Arrays.asList(ReferencePair.Status.DIFFER_IN_SEQUENCES_PRESENT))
                },
                new Object[]{ tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta")),
                        CompareReferences.MD5CalculationMode.USE_DICT),
                        new HashSet<>(Arrays.asList(ReferencePair.Status.EXACT_MATCH))
                },
                new Object[]{ tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_chr2snp.fasta"),
                                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_missingchr3.fasta")),
                        CompareReferences.MD5CalculationMode.USE_DICT),
                        new HashSet<>(Arrays.asList(ReferencePair.Status.DIFFER_IN_SEQUENCES_PRESENT, ReferencePair.Status.DIFFER_IN_SEQUENCE))
                },
                new Object[]{ tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_missingchr1.fasta"),
                                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_missingchr1_renamedchr2.fasta")),
                        CompareReferences.MD5CalculationMode.USE_DICT),
                        new HashSet<>(Arrays.asList(ReferencePair.Status.DIFFER_IN_SEQUENCE_NAMES))
                },
        };
    }

    @Test(dataProvider = "testAnalyzeTableTwoRefsData")
    public void testAnalyzeTableTwoReferences(ReferenceSequenceTable table, Set<ReferencePair.Status> expectedStatus){
        List<ReferencePair> refPairs = table.analyzeTable();
        for(ReferencePair pair : refPairs){
            Assert.assertEquals(pair.getAnalysis(), expectedStatus);
        }
    }

    private Map<ReferencePair, Set<ReferencePair.Status>> manuallySetReferencePairStatus(){
        ReferenceSequenceTable table = tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                        new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_1renamed.fasta"),
                        new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_chr2snp.fasta")),
                CompareReferences.MD5CalculationMode.USE_DICT);

        ReferencePair refPair1 = new ReferencePair(table, new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_1renamed.fasta"));
        refPair1.removeStatus(ReferencePair.Status.EXACT_MATCH);
        refPair1.addStatus(ReferencePair.Status.DIFFER_IN_SEQUENCE_NAMES);

        ReferencePair refPair2 = new ReferencePair(table, new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_chr2snp.fasta"));
        refPair2.removeStatus(ReferencePair.Status.EXACT_MATCH);
        refPair2.addStatus(ReferencePair.Status.DIFFER_IN_SEQUENCE);

        ReferencePair refPair3 = new ReferencePair(table, new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_1renamed.fasta"),
                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_chr2snp.fasta"));
        refPair3.removeStatus(ReferencePair.Status.EXACT_MATCH);
        refPair3.addStatus(ReferencePair.Status.DIFFER_IN_SEQUENCE_NAMES);
        refPair3.addStatus(ReferencePair.Status.DIFFER_IN_SEQUENCE);

        Map<ReferencePair, Set<ReferencePair.Status>> referencePairStatuses = new HashMap<>();
            referencePairStatuses.put(refPair1, refPair1.getAnalysis());
            referencePairStatuses.put(refPair2, refPair2.getAnalysis());
            referencePairStatuses.put(refPair3, refPair3.getAnalysis());

        return referencePairStatuses;
    }

    @DataProvider(name = "testAnalyzeTableMultipleRefsData")
    public Object[][] testAnalyzeTableMultipleReferencesData(){
        return new Object[][]{
                new Object[]{tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_1renamed.fasta"),
                                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_chr2snp.fasta")),
                        CompareReferences.MD5CalculationMode.USE_DICT), manuallySetReferencePairStatus()},
        };
    }

    @Test(dataProvider = "testAnalyzeTableMultipleRefsData")
    public void testAnalyzeTableMultipleReferences(ReferenceSequenceTable table, Map<ReferencePair, Set<ReferencePair.Status>> expectedStatus){
        List<ReferencePair> refPairs = table.analyzeTable();
        for(ReferencePair pair : refPairs){ ;
            Assert.assertEquals(pair.getAnalysis(), expectedStatus.get(pair));
        }
    }

}
