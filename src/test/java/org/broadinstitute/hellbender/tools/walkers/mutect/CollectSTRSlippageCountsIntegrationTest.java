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
import java.util.Map;
import java.util.stream.Collectors;

public class CollectSTRSlippageCountsIntegrationTest extends CommandLineProgramTest {
    private static final File DREAM_4_NORMAL = new File(largeFileTestDir + "mutect/dream_synthetic_bams/", "normal_4.bam");
    private static final File GNOMAD = new File(largeFileTestDir, "very-small-gnomad.vcf");

    @Test
    public void testDreamNormal() throws IOException {
        final File output = createTempFile("slippage-counts", ".table");
        runCommandLine(new ArgumentsBuilder()
                .addInput(DREAM_4_NORMAL)
                .addReference(b37Reference)
                .addInterval("20:10000000-12000000")
                .addOutput(output));

        final List<String> lines = readLines(output);
        final List<String> columns = Arrays.asList(lines.get(0).split("\t"));
        final int altColumnsStart = 5 + CollectSTRSlippageCounts.COUNT_COLUMNS.size();
        Assert.assertEquals(columns.subList(5, altColumnsStart), CollectSTRSlippageCounts.COUNT_COLUMNS);
        Assert.assertEquals(columns.subList(altColumnsStart, columns.size()), CollectSTRSlippageCounts.ALT_COLUMNS);
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

            final int[] counts = Arrays.stream(fields, 5, altColumnsStart).mapToInt(Integer::parseInt).toArray();
            Assert.assertTrue(Arrays.stream(counts).sum() > 0);
            unchanged += counts[CollectSTRSlippageCounts.MAX_BINNED_CHANGE];
            changed += Arrays.stream(counts).sum() - counts[CollectSTRSlippageCounts.MAX_BINNED_CHANGE];

            // without a germline resource there are no allele frequencies
            Assert.assertEquals(fields[columns.indexOf("top_alt_af")], "NaN");
            Assert.assertEquals(fields[columns.indexOf("second_alt_af")], "NaN");
        }

        // slippage is the exception
        Assert.assertTrue(changed > 0);
        Assert.assertTrue(changed < 0.05 * unchanged);
    }

    @Test
    public void testGermlineResource() throws IOException {
        final File output = createTempFile("slippage-counts", ".table");
        runCommandLine(new ArgumentsBuilder()
                .addInput(new File(NA12878_20_21_WGS_bam))
                .addReference(b37Reference)
                .addInterval("20:10000000-11000000")
                .add(M2ArgumentCollection.GERMLINE_RESOURCE_LONG_NAME, GNOMAD)
                .addOutput(output));

        final List<String> lines = readLines(output);
        final List<String> columns = Arrays.asList(lines.get(0).split("\t"));
        final Map<Integer, String[]> rowsByStart = lines.subList(1, lines.size()).stream()
                .map(line -> line.split("\t"))
                .collect(Collectors.toMap(fields -> Integer.parseInt(fields[1]), fields -> fields));

        int alleleFrequencies = 0;
        for (final String[] fields : rowsByStart.values()) {
            for (int n = 0; n < CollectSTRSlippageCounts.ANNOTATED_ALTS; n++) {
                final int change = Integer.parseInt(fields[columns.indexOf(CollectSTRSlippageCounts.ALT_COLUMNS.get(2 * n))]);
                final double alleleFrequency = Double.parseDouble(fields[columns.indexOf(CollectSTRSlippageCounts.ALT_COLUMNS.get(2 * n + 1))]);
                if (change == 0) {
                    Assert.assertTrue(Double.isNaN(alleleFrequency));
                } else if (!Double.isNaN(alleleFrequency)) {
                    Assert.assertTrue(alleleFrequency > 0 && alleleFrequency <= 1);
                    alleleFrequencies++;
                }
            }
        }
        Assert.assertTrue(alleleFrequencies > 100);

        // every fragment of this GATA x 3 that changes it by whole units has one more GATA
        Assert.assertEquals(alts(columns, rowsByStart.get(10037111)), List.of("1", "0.796", "0", "NaN"));

        // fragments of this AC x 17 have 5 or 3 fewer ACs, which gnomAD has in a multiallelic record extending past the STR
        Assert.assertEquals(alts(columns, rowsByStart.get(10241337)), List.of("-5", "0.452", "-3", "0.331"));
    }

    private static List<String> readLines(final File table) throws IOException {
        return Files.readAllLines(table.toPath()).stream().filter(line -> !line.startsWith("#")).collect(Collectors.toList());
    }

    // a row's most common changes and their allele frequencies
    private static List<String> alts(final List<String> columns, final String[] fields) {
        return CollectSTRSlippageCounts.ALT_COLUMNS.stream().map(column -> fields[columns.indexOf(column)]).collect(Collectors.toList());
    }
}
