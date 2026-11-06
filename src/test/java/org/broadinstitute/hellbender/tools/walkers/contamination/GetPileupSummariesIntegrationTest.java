package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Created by David Benjamin on 2/16/17.
 */
public class GetPileupSummariesIntegrationTest extends CommandLineProgramTest {
    private static final File NA12878 = new File(largeFileTestDir, "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam");
    private static final File GNOMAD = new File(largeFileTestDir, "very-small-gnomad.vcf");

    @Test
    public void test() {
        final File output = createTempFile("output", ".table");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInput(NA12878)
                .addVCF(new File(thousandGenomes))
                .addIntervals(new File(thousandGenomes))
                .addOutput(output)
                .add(GetPileupSummaries.MAX_SITE_AF_SHORT_NAME, 0.9);

        runCommandLine(args);

        final ImmutablePair<String, List<PileupSummary>> sampleAndResult = PileupSummary.readFromFile(output);

        final List<PileupSummary> result = sampleAndResult.getRight();
        final String sample = sampleAndResult.getLeft();
        Assert.assertEquals(sample, "NA12878");

        // compare to IGV manual inspection
        final PileupSummary ps1 = result.get(0);
        Assert.assertEquals(ps1.getContig(), "20");
        Assert.assertEquals(ps1.getStart(), 10000117);
        Assert.assertEquals(ps1.getRefCount(), 35);
        Assert.assertEquals(ps1.getAltCount(), 28);
        Assert.assertEquals(ps1.getOtherAltCount(), 0);
        Assert.assertEquals(ps1.getAlleleFrequency(), 0.605);
        Assert.assertEquals(ps1.getRefBase(), Nucleotide.C);
        Assert.assertEquals(ps1.getAltBase(), Nucleotide.T);

        final PileupSummary ps2 = result.get(1);
        Assert.assertEquals(ps2.getStart(), 10000211);
        Assert.assertEquals(ps2.getRefCount(), 27);
        Assert.assertEquals(ps2.getAltCount(), 28);
        Assert.assertEquals(ps2.getAlleleFrequency(), 0.603);

        final PileupSummary ps3 = result.get(2);
        Assert.assertEquals(ps3.getStart(), 10000439);
        Assert.assertEquals(ps3.getRefCount(), 0);
        Assert.assertEquals(ps3.getAltCount(), 76);
        Assert.assertEquals(ps3.getAlleleFrequency(), 0.81);
        Assert.assertEquals(ps3.getRefBase(), Nucleotide.T);
        Assert.assertEquals(ps3.getAltBase(), Nucleotide.G);

        final PileupSummary ps4 = result.get(8);
        Assert.assertEquals(ps4.getStart(), 10001298);
        Assert.assertEquals(ps4.getRefCount(), 0);
        Assert.assertEquals(ps4.getAltCount(), 53);
        Assert.assertEquals(ps4.getOtherAltCount(), 0);
        Assert.assertEquals(ps4.getAlleleFrequency(), 0.809);

    }

    // the bam only covers some of the gnomAD sites in 20:10-11Mb
    @Test
    public void testIncludeUncoveredSites() {
        final List<PileupSummary> coveredSites = runOnGnomadSites(false);
        final List<PileupSummary> allSites = runOnGnomadSites(true);

        // the same sites with reads, with the same counts, plus sites without reads
        final Map<Integer, PileupSummary> allSitesByPosition = allSites.stream()
                .collect(Collectors.toMap(PileupSummary::getStart, Function.identity()));
        Assert.assertTrue(allSites.size() > coveredSites.size());
        for (final PileupSummary site : coveredSites) {
            final PileupSummary sameSite = allSitesByPosition.get(site.getStart());
            Assert.assertEquals(sameSite.getRefCount(), site.getRefCount());
            Assert.assertEquals(sameSite.getAltCount(), site.getAltCount());
            Assert.assertEquals(sameSite.getOtherAltCount(), site.getOtherAltCount());
        }

        final Set<Integer> coveredPositions = coveredSites.stream().map(PileupSummary::getStart).collect(Collectors.toSet());
        allSites.stream().filter(site -> !coveredPositions.contains(site.getStart()))
                .forEach(site -> Assert.assertEquals(site.getTotalCount(), 0));
    }

    private List<PileupSummary> runOnGnomadSites(final boolean includeUncoveredSites) {
        final File output = createTempFile("output", ".table");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInput(NA12878)
                .addVCF(GNOMAD)
                .addIntervals(GNOMAD)
                .addOutput(output)
                .add(GetPileupSummaries.MAX_SITE_AF_SHORT_NAME, 0.9);
        if (includeUncoveredSites) {
            args.addFlag(GetPileupSummaries.INCLUDE_UNCOVERED_SITES_LONG_NAME);
        }
        runCommandLine(args);
        return PileupSummary.readFromFile(output).getRight();
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNoAFFieldInHeader() {
        final File vcfWithoutAF = new File(publicTestDir, "empty.vcf");
        final File output = createTempFile("output", ".table");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInput(NA12878)
                .addVCF(vcfWithoutAF)
                .addIntervals(vcfWithoutAF)
                .addOutput(output);

        runCommandLine(args);
    }

}