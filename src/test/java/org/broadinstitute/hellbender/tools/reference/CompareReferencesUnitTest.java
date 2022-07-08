package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import java.io.File;
import java.util.*;

public class CompareReferencesUnitTest {

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
                new Object[]{ Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta")),
                        "8c0c38e352d8f3309eabe4845456f274",
                        16000,
                        tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta")),
                                CompareReferences.MD5CalculationMode.USE_DICT),
                },
                new Object[]{ Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                        new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_1renamed.fasta")),
                        "8c0c38e352d8f3309eabe4845456f274",
                        16000,
                        tableGenerator(Arrays.asList(new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini.fasta"),
                                new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_1renamed.fasta")),
                                CompareReferences.MD5CalculationMode.USE_DICT),
                },
        };
    }

    @Test(dataProvider = "testQueryByMD5Data")
    public void testQueryByMD5(List<GATKPath> references, String md5, int length, ReferenceSequenceTable table){
        Map<GATKPath, ReferenceDataSource> referenceSources = new LinkedHashMap<>();
        for(GATKPath path : references){
            referenceSources.put(path, ReferenceDataSource.of(path.toPath()));
        }

        ReferenceSequenceTable.TableRow expected = table.new TableRow(md5, length, references);
        ReferenceSequenceTable.TableRow actual = table.queryByMD5(md5);
        Assert.assertEquals(actual, expected);
    }

    @Test
    public void testQueryBySequenceName(){
        GATKPath refPath = new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/hg19mini_chr2snp.fasta");
        Map<GATKPath, ReferenceDataSource> referenceSources = new LinkedHashMap<>();
        referenceSources.put(refPath, ReferenceDataSource.of(refPath.toPath()));
        CompareReferences.MD5CalculationMode md5CalculationMode = CompareReferences.MD5CalculationMode.USE_DICT;
        ReferenceSequenceTable table = new ReferenceSequenceTable(referenceSources, md5CalculationMode);
        table.build();

        String sequenceName = "2";
        int expectedLength = 16000;
        List<GATKPath> references = new ArrayList<>();
        references.add(refPath);

        Set<ReferenceSequenceTable.TableRow> expected = new HashSet<>();
        expected.add(table.new TableRow("600ccbfe836f04b26c0eea36dd0485ee", expectedLength, references));

        Set<ReferenceSequenceTable.TableRow> actual = table.queryBySequenceName(sequenceName);
        Assert.assertEquals(actual, expected);
    }

}
