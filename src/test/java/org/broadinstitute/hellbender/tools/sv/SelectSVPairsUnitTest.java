package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Map;

public class SelectSVPairsUnitTest extends GATKBaseTest {

    private static final String VALID_TABLE = String.join("\n",
            "VID_A\tVID_B\tSCORE",
            "A1\tB1\t10.0",
            "A1\tB2\t9.0",
            "A2\tB1\t8.0",
            "A3\tB3\t7.5",
            "A3\tB3\t7.5",
            "A4\tB4\t5.0",
            "");

    @Test
    public void testSelectHighestScoringNonConflictingPairs() throws IOException {
        final File table = createTempFile("select-sv-pairs", ".tsv");
        Files.writeString(table.toPath(), VALID_TABLE, StandardCharsets.UTF_8);
        final SelectSVPairs selector = new SelectSVPairs(new GATKPath(table.getAbsolutePath()));

        Assert.assertEquals(selector.getVidAToBMap(), Map.of("A1", "B1", "A3", "B3", "A4", "B4"));
        Assert.assertEquals(selector.getVidBToAMap(), Map.of("B1", "A1", "B3", "A3", "B4", "A4"));
    }

    @Test
    public void testRejectMissingColumn() throws IOException {
        final File table = createTempFile("select-sv-pairs", ".tsv");
        Files.writeString(table.toPath(), String.join("\n",
                "VID_A\tVID_B",
                "A1\tB1",
            ""), StandardCharsets.UTF_8);

        Assert.expectThrows(RuntimeException.class, () -> new SelectSVPairs(new GATKPath(table.getAbsolutePath())));
    }

    @Test
    public void testRejectExtraColumn() throws IOException {
        final File table = createTempFile("select-sv-pairs", ".tsv");
        Files.writeString(table.toPath(), String.join("\n",
                "VID_A\tVID_B\tSCORE\tEXTRA",
                "A1\tB1\t1.0\textra",
            ""), StandardCharsets.UTF_8);

        Assert.expectThrows(RuntimeException.class, () -> new SelectSVPairs(new GATKPath(table.getAbsolutePath())));
    }

    @Test
    public void testSVPairEquality() {
        final SelectSVPairs.SVPair first = new SelectSVPairs.SVPair("A1", "B1", 1.0);
        final SelectSVPairs.SVPair same = new SelectSVPairs.SVPair("A1", "B1", 1.0);
        final SelectSVPairs.SVPair differentScore = new SelectSVPairs.SVPair("A1", "B1", 2.0);

        Assert.assertEquals(first, same);
        Assert.assertEquals(first.hashCode(), same.hashCode());
        Assert.assertNotEquals(first, differentScore);
    }
}