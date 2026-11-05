package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class CollectSTRSlippageCountsIntegrationTest extends CommandLineProgramTest {
    private static final File DREAM_4_NORMAL = new File(largeFileTestDir + "mutect/dream_synthetic_bams/", "normal_4.bam");

    @Test
    public void testDreamNormal() throws IOException {
        final File output = createTempFile("slippage-counts", ".table");
        runCommandLine(new ArgumentsBuilder()
                .addInput(DREAM_4_NORMAL)
                .addReference(b37Reference)
                .addInterval("20:10000000-12000000")
                .addOutput(output));

        final List<String> lines = Files.readAllLines(output.toPath()).stream()
                .filter(line -> !line.startsWith("#")).collect(Collectors.toList());
        final List<String> columns = Arrays.asList(lines.get(0).split("\t"));
        Assert.assertEquals(columns.subList(5, columns.size()), CollectSTRSlippageCounts.COUNT_COLUMNS);
        Assert.assertTrue(lines.size() > 100);

        int previousEnd = 0;
        long unchanged = 0;
        long changed = 0;
        for (final String line : lines.subList(1, lines.size())) {
            final String[] fields = line.split("\t");
            final int start = Integer.parseInt(fields[1]);
            final int end = Integer.parseInt(fields[2]);
            final String unit = fields[3];
            final int units = Integer.parseInt(fields[4]);

            // STRs are in order, don't overlap, and satisfy the default definition
            Assert.assertTrue(start > previousEnd);
            Assert.assertEquals(end - start + 1, unit.length() * units);
            Assert.assertTrue(unit.matches("[ACGT]+"));
            Assert.assertTrue(end - start + 1 >= CollectSTRSlippageCounts.DEFAULT_MIN_REPEAT_LENGTH);
            Assert.assertTrue(units >= CollectSTRSlippageCounts.DEFAULT_MIN_REPEAT_UNITS);
            Assert.assertTrue(unit.length() <= CollectSTRSlippageCounts.DEFAULT_MAX_REPEAT_UNIT_LENGTH);
            previousEnd = end;

            final int[] counts = Arrays.stream(fields, 5, fields.length).mapToInt(Integer::parseInt).toArray();
            Assert.assertTrue(Arrays.stream(counts).sum() > 0);
            unchanged += counts[CollectSTRSlippageCounts.MAX_BINNED_CHANGE];
            changed += Arrays.stream(counts).sum() - counts[CollectSTRSlippageCounts.MAX_BINNED_CHANGE];
        }

        // slippage is the exception
        Assert.assertTrue(changed > 0);
        Assert.assertTrue(changed < 0.05 * unchanged);
    }
}
