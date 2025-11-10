package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.testng.Assert.*;

public class MismatchCountTest extends GATKBaseTest {
    private static final int READ_LENGTH = 50;
    private static final Allele REF = Allele.create("A", true);
    private static final Allele ALT = Allele.create("C");
    private static final List<Allele> ALLELES = Arrays.asList(REF, ALT);
    private static final String SAMPLE = "sample1";

    private GATKRead makeRead(final byte[] readBases, final int nMismatches) {
        for (int i = 0 ; i < nMismatches; i++){
            if (readBases[i] == 'A'){
                readBases[i] = 'T';
            } else {
                readBases[i] = 'A';
            }
        }

        byte[] qual = Utils.dupBytes((byte)30, READ_LENGTH);
        String cigar = String.format("%dM", READ_LENGTH);
        GATKRead read =  ArtificialReadUtils.createArtificialRead(readBases, qual, cigar);
        read.setMappingQuality(20);
        read.setAttribute(ReadUtils.NUM_MISMATCH_TAG, nMismatches);
        return read;
    }


    @Test
    public void testDescription(){
        String[] constants = {GATKVCFConstants.NM_COUNT_KEY};
        VCFFormatHeaderLine[] hlines    = {GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.NM_COUNT_KEY)};
        Assert.assertEquals(new MismatchCount().getKeyNames(), new ArrayList<>(Arrays.asList(constants)));
        Assert.assertEquals(new MismatchCount().getDescriptions(), new ArrayList<>(Arrays.asList(hlines)));
    }

    @BeforeClass

    @DataProvider(name = "MismatchCountDataProvider")
    public Object[][] makeMismatchCountDataProvider() {
        ReferenceFileSource refSource = new ReferenceFileSource(Path.of(b37_reference_20_21));
        ReferenceContext    rc        = new ReferenceContext(refSource, new SimpleInterval("20", 1, 10205));

        final List<GATKRead> refReads = IntStream.range(0, 10).mapToObj(i ->
                makeRead(Arrays.copyOfRange(rc.getBases(),10000,10000+READ_LENGTH),
                0)).collect(Collectors.toList());
        final List<GATKRead> altReads = IntStream.range(0, 10).mapToObj(i ->
                makeRead(Arrays.copyOfRange(rc.getBases(),10000,10000+READ_LENGTH),
                        5)).collect(Collectors.toList());
        final AlleleLikelihoods<GATKRead, Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(SAMPLE, refReads, altReads, -100.0, -100.0, REF, ALT);
        final int dpDepth = 50; //Note: using a different value on purpose so that we can check that reads are preferred over DP

        final Genotype gAC = new GenotypeBuilder(SAMPLE, ALLELES).DP(dpDepth).make();
        final double log10PError = -5;
        final VariantContext vc = new VariantContextBuilder("test", "1", 10005, 10005, ALLELES)
                .log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();
        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        int[] result = {0,5};
        return new Object[][]{{likelihoods, vc, gb, result}};
    }

    @Test(dataProvider = "MismatchCountDataProvider")
    public void testUsingReads(final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                               final VariantContext vc, final GenotypeBuilder gb, int[] result){

        new MismatchCount().annotate(new ReferenceContext(new ReferenceFileSource(Path.of(b37_reference_20_21)),
                new SimpleInterval("20", 1, 10205)), vc, vc.getGenotype(SAMPLE), gb, likelihoods);
        final int[] mismatchCounts = (int [])gb.make().getExtendedAttribute("NMC");

        Assert.assertEquals(mismatchCounts, result);
    }

}