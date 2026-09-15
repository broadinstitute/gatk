package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.downsampling.ChainedReadsDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.MaxDepthDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.PositionalDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class AssemblyRegionArgumentCollectionUnitTest extends GATKBaseTest {

    private final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();

    @DataProvider(name = "downsamplerCombinations")
    public Object[][] downsamplerCombinations() {
        return new Object[][] {
                // max reads per alignment start, max effective depth, expected downsampler class (null for none)
                { 0, 0, null },
                { 50, 0, PositionalDownsampler.class },
                { 0, 200, MaxDepthDownsampler.class },
                { 50, 200, ChainedReadsDownsampler.class },
        };
    }

    @Test(dataProvider = "downsamplerCombinations")
    public void testCreateReadsDownsamplerCombinations(final int maxReadsPerAlignmentStart, final int maxEffectiveDepth, final Class<?> expected) {
        final AssemblyRegionArgumentCollection args = new AssemblyRegionArgumentCollection();
        args.maxReadsPerAlignmentStart = maxReadsPerAlignmentStart;
        args.maxEffectiveDepth = maxEffectiveDepth;

        final ReadsDownsampler downsampler = args.createReadsDownsampler(header, false);

        if (expected == null) {
            Assert.assertNull(downsampler);
        } else {
            Assert.assertEquals(downsampler.getClass(), expected);
        }
    }

    @Test
    public void testDefaultsProduceOnlyThePerStartDownsampler() {
        final ReadsDownsampler downsampler = new AssemblyRegionArgumentCollection().createReadsDownsampler(header, false);
        Assert.assertEquals(downsampler.getClass(), PositionalDownsampler.class);
    }

    @DataProvider(name = "invalidDepthArguments")
    public Object[][] invalidDepthArguments() {
        return new Object[][] {
                // max effective depth, window
                { -1, 1000 },
                { 100, 0 },
                { 100, -10 },
        };
    }

    @Test(dataProvider = "invalidDepthArguments", expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void testValidateRejectsInvalidDepthArguments(final int maxEffectiveDepth, final int window) {
        final AssemblyRegionArgumentCollection args = new AssemblyRegionArgumentCollection();
        args.maxEffectiveDepth = maxEffectiveDepth;
        args.maxEffectiveDepthWindow = window;
        args.validate();
    }

    @Test
    public void testValidateAcceptsDisabledAndEnabledDepthCap() {
        final AssemblyRegionArgumentCollection args = new AssemblyRegionArgumentCollection();
        args.validate();
        args.maxEffectiveDepth = 200;
        args.validate();
    }
}
