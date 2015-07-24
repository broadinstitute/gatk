package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public final class ActivityProfileStateUnitTest {
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() {
        // sequence
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 100);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    @DataProvider(name = "ActiveProfileResultProvider")
    public Object[][] makeActiveProfileResultProvider() {
        final List<Object[]> tests = new LinkedList<>();

        final String chr = genomeLocParser.getSequenceDictionary().getSequence(0).getSequenceName();
        for ( final GenomeLoc loc : Arrays.asList(
                genomeLocParser.createGenomeLoc(chr, 10, 10),
                genomeLocParser.createGenomeLoc(chr, 100, 100))) {
            for ( final double prob : Arrays.asList(0.0, 0.5, 1.0) ) {
                for ( final ActivityProfileState.Type state : ActivityProfileState.Type.values() ) {
                    for ( final Number value : Arrays.asList(1, 2, 4) ) {
                        tests.add(new Object[]{ loc, prob, state, value});
                    }
                }
                tests.add(new Object[]{ loc, prob, null, null});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ActiveProfileResultProvider")
    public void testActiveProfileResultProvider(GenomeLoc loc, final double prob, ActivityProfileState.Type maybeState, final Number maybeNumber) {
        final ActivityProfileState apr = maybeState == null
                ? new ActivityProfileState(loc, prob)
                : new ActivityProfileState(loc, prob, maybeState, maybeNumber);

        Assert.assertEquals(apr.getLoc(), loc);
        Assert.assertEquals(apr.getOffset(loc), 0);
        Assert.assertNotNull(apr.toString());
        Assert.assertEquals(apr.isActiveProb(), prob);
        Assert.assertEquals(apr.getResultState(), maybeState == null ? ActivityProfileState.Type.NONE : maybeState);
        Assert.assertEquals(apr.getResultValue(), maybeState == null ? null : maybeNumber);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testError1(){
        final String chr = genomeLocParser.getSequenceDictionary().getSequence(0).getSequenceName();
        new ActivityProfileState(genomeLocParser.createGenomeLoc(chr, 10, 10), 0.1, null, -1.0);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testError2(){
        final String chr = genomeLocParser.getSequenceDictionary().getSequence(0).getSequenceName();
        new ActivityProfileState(genomeLocParser.createGenomeLoc(chr, 10, 11), 0.1, null, 1.0);
    }
}
