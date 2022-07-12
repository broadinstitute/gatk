package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class ReferenceSequenceTableUnitTest {

    public ReferenceSequenceTable tableGenerator(List<GATKPath> references, CompareReferences.MD5CalculationMode md5CalculationMode){
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
            Assert.assertNull(rows);
        } else {
            Assert.assertEquals(rows.size(), expectedRows);

            boolean contigPresent = false;
            for (ReferenceSequenceTable.TableRow row : rows) {
                contigPresent = false;
                for (ReferenceSequenceTable.TableEntry entry : row.getEntries()) {
                    if (entry.getColumnValue().equals(expectedSequenceName)) {
                        contigPresent = true;
                        continue;
                    }
                }
                Assert.assertTrue(contigPresent);
            }
        }
    }

}
