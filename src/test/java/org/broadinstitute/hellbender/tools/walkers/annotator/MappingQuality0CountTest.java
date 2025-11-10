package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.testng.Assert.*;

public class MappingQuality0CountTest extends GATKBaseTest {
    private static final Allele REF = Allele.create("A", true);
    private static final Allele ALT = Allele.create("C");
    private static final List<Allele> ALLELES = Arrays.asList(REF, ALT);
    private static final String SAMPLE = "sample1";

    private GATKRead makeRead(final byte readBase, final byte mapQuality) {
        final int readLength = 10;
        return ArtificialAnnotationUtils.makeRead(readBase, (byte)30, mapQuality);
    }

    @Test
    public void testDescription(){
        String[] constants = {GATKVCFConstants.MQ0_COUNT_KEY};
        VCFFormatHeaderLine[] hlines    = {GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.MQ0_COUNT_KEY)};
        Assert.assertEquals(new MappingQuality0Count().getKeyNames(), new ArrayList<>(Arrays.asList(constants)));
        Assert.assertEquals(new MappingQuality0Count().getDescriptions(), new ArrayList<>(Arrays.asList(hlines)));
    }

    @DataProvider(name = "readData")
    public Object[][] readDepthData() {
        return new Object[][]{{0,17}};
    }

    @Test(dataProvider = "readData")
    public void testUsingReads(final int expectedRefDepth, final int expectedAltDepth){


        final int dpDepth = 50; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(SAMPLE, ALLELES).DP(dpDepth).make();

        final double log10PError = -5;

        final List<GATKRead> refReads = IntStream.range(0, expectedRefDepth).mapToObj(i -> makeRead(REF.getBases()[0], (byte)30)).collect(Collectors.toList());
        final List<GATKRead> altReads = IntStream.range(0, expectedAltDepth).mapToObj(i -> makeRead(ALT.getBases()[0], (byte)0)).collect(Collectors.toList());
        final AlleleLikelihoods<GATKRead, Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(SAMPLE, refReads, altReads, -100.0, -100.0, REF, ALT);

        final VariantContext vc = new VariantContextBuilder("test", "1", 10005, 10005, ALLELES).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        new MappingQuality0Count().annotate(new ReferenceContext(new ReferenceFileSource(Path.of(b37_reference_20_21)),
                new SimpleInterval("20", 10005, 10005)), vc, gAC, gb, likelihoods);
        final int[] mq0 = (int [])gb.make().getExtendedAttribute("MQ0C");
        Assert.assertEquals(mq0, new int[] {expectedRefDepth,expectedAltDepth});
    }

}