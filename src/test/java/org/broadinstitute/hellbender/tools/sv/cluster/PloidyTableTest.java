package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

public class PloidyTableTest extends CommandLineProgramTest {

    @Test
    public void testCreateTableFromFile() {
        final Path tablePath = Paths.get( getTestDataDir() + "/walkers/sv/svcluster/1kgp.batch1.ploidy.tsv");
        final PloidyTable table = new PloidyTable(tablePath);

        Assert.assertEquals(table.get("HG00096", "chr1").intValue(), 2);
        Assert.assertEquals(table.get("HG00096", "chrX").intValue(), 1);
        Assert.assertEquals(table.get("HG00096", "chrY").intValue(), 1);

        Assert.assertEquals(table.get("NA19143", "chr1").intValue(), 2);
        Assert.assertEquals(table.get("NA19143", "chrX").intValue(), 2);
        Assert.assertEquals(table.get("NA19143", "chrY").intValue(), 0);

        Assert.assertEquals(table.get("NA21133", "chr1").intValue(), 2);
        Assert.assertEquals(table.get("NA21133", "chrX").intValue(), 1);
        Assert.assertEquals(table.get("NA21133", "chrY").intValue(), 1);
    }

    @Test
    public void testCreateTableFromMap() {
        final Map<String, Map<String, Integer>> map = new HashMap<>();

        final Map<String, Integer> sampleMap1 = new HashMap<>();
        sampleMap1.put("chr1", 2);
        sampleMap1.put("chrX", 2);
        map.put("sample1", sampleMap1);

        final Map<String, Integer> sampleMap2 = new HashMap<>();
        sampleMap2.put("chr1", 2);
        sampleMap2.put("chrX", 1);
        map.put("sample2", sampleMap2);

        final PloidyTable table = new PloidyTable(map);

        Assert.assertEquals(table.get("sample1", "chr1").intValue(), 2);
        Assert.assertEquals(table.get("sample1", "chrX").intValue(), 2);
        Assert.assertEquals(table.get("sample2", "chr1").intValue(), 2);
        Assert.assertEquals(table.get("sample2", "chrX").intValue(), 1);
    }

}