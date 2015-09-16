package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AlignmentContextUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import java.util.List;

public final class ConsensusAlleleCounterUnitTest extends BaseTest {
    
    final int nSamples = 1;
    final int[] numReadsPerAllele = {10,10};
    final String SAMPLE_PREFIX = "sample";

    @Test
    public void testBasicConsensusCounts() {
        final ArtificialReadPileupTestProvider pp = new ArtificialReadPileupTestProvider(nSamples, SAMPLE_PREFIX);

        // 4 inserted bases, min cnt = 10
        String altBases = "CCTC";
        int eventLength = 4;
        List<Allele> alleles = getConsensusAlleles(pp,  eventLength,true,10,0.1, altBases);
        
        Assert.assertEquals(alleles.size(), 2);
        Assert.assertEquals(alleles.get(1).getBaseString().substring(1), altBases.substring(0, eventLength));


        // test deletions
        eventLength = 3;
        alleles = getConsensusAlleles(pp,  eventLength,false,10,0.1, altBases);
        Assert.assertEquals(alleles.size(), 2);
        Assert.assertEquals(alleles.get(0).getBaseString().substring(1, eventLength), new String(pp.getReferenceContext().getForwardBases()).substring(1, eventLength));

        // same with min Reads = 11
        alleles = getConsensusAlleles(pp,  eventLength,false,11,0.1, altBases);
        Assert.assertEquals(alleles.size(), 0);

        // increase required fraction per sample to just below threshold
        alleles = getConsensusAlleles(pp,  eventLength,false,10,0.49999, altBases);
        Assert.assertEquals(alleles.size(), 2);
        alleles = getConsensusAlleles(pp,  eventLength,false,10,0.5001, altBases);
        Assert.assertEquals(alleles.size(), 0);

        // test N's in insertions
        altBases = "CCTC";
        eventLength = 4;
        alleles = getConsensusAlleles(pp,  eventLength,true,10,0.1, altBases);

        Assert.assertEquals(alleles.size(), 2);
        Assert.assertEquals(alleles.get(1).getBaseString().substring(1, eventLength + 1), altBases);

        altBases = "CCTCN";
        eventLength = 5;
        alleles = getConsensusAlleles(pp,  eventLength,true,10,0.1, altBases);

        Assert.assertEquals(alleles.size(), 2);

    }
    
    private List<Allele> getConsensusAlleles(final ArtificialReadPileupTestProvider pp, int eventLength, boolean isInsertion, int minCnt, double minFraction, String altBases) {

        final ConsensusAlleleCounter counter = new ConsensusAlleleCounter(minCnt, minFraction);
        return counter.computeConsensusAlleles(pp.referenceContext,
                pp.getAlignmentContextFromAlleles(isInsertion?eventLength:-eventLength,altBases,numReadsPerAllele),
                AlignmentContextUtils.ReadOrientation.COMPLETE);

    }
}
