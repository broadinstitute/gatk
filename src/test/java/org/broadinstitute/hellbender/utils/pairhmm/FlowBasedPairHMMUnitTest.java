package org.broadinstitute.hellbender.utils.pairhmm;

import com.google.common.base.Strings;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.LikelihoodEngineArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PairHMMLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParamUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.FlowBasedHMMEngine;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;

import java.util.*;

public class FlowBasedPairHMMUnitTest extends GATKBaseTest {

    // at this point, this constant MUST be 4, unless otherwise tested
    @Test
    public void testFlowSize4() {
        Assert.assertEquals(FlowBasedPairHMM.FLOW_SIZE, 4);
    }


    // TODO NOTE this test is disabled. That is because currently it is difficult to generate ArtificialFlowBasedReads from scratch...
    @Test (enabled = false)
    public void testComputeLikelihoods(){
        final LikelihoodEngineArgumentCollection LEAC = new LikelihoodEngineArgumentCollection();

        final FlowBasedArgumentCollection defaultArgs = new FlowBasedArgumentCollection();
        final FlowBasedHMMEngine lce = new FlowBasedHMMEngine(defaultArgs, (byte) LEAC.gcpHMM, QualityUtils.qualToErrorProbLog10(LEAC.phredScaledGlobalReadMismappingRate), LEAC.expectedErrorRatePerBase, PairHMMLikelihoodCalculationEngine.PCRErrorModel.CONSERVATIVE,
                LEAC.dontUseDragstrPairHMMScores ? null : DragstrParamUtils.parse(LEAC.dragstrParams), LEAC.enableDynamicReadDisqualification, LEAC.readDisqualificationThresholdConstant,
                LEAC.minUsableIndelScoreToUse, (byte) LEAC.flatDeletionPenalty, (byte) LEAC.flatInsertionPenatly);


        final Map<String, List<GATKRead>> perSampleReadList = new HashMap<>();
        final int n = 10;
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode(n + "M"));
        read1.setMappingQuality(60);
        final String sample1 = "sample1";
        perSampleReadList.put(sample1, Arrays.asList(read1));

        final SampleList samples = new IndexedSampleList(sample1);

        final List<Haplotype> haplotypes = new ArrayList<>();
        final byte[] bases = Strings.repeat("A", n + 1).getBytes();
        final Haplotype hap1 = new Haplotype(bases, true);
        hap1.setGenomeLocation(read1);
        haplotypes.add(hap1);

        final byte[] basesModified = bases;
        basesModified[5] = 'C';//different bases
        final Haplotype hap2 = new Haplotype(basesModified, false);
        hap2.setGenomeLocation(read1);//use same loc
        haplotypes.add(hap2);

        final SAMReadGroupRecord rg = new SAMReadGroupRecord("x");
        rg.setAttribute(SAMReadGroupRecord.FLOW_ORDER_TAG, "ACTG");
        final SAMFileHeader hdr = new SAMFileHeader();
        hdr.addReadGroup(rg);

        final AlleleLikelihoods<GATKRead, Haplotype> likes = lce.computeReadLikelihoods(haplotypes, hdr, samples, perSampleReadList, true);
        final LikelihoodMatrix<GATKRead, Haplotype> mtx = likes.sampleMatrix(0);

        Assert.assertEquals(mtx.numberOfAlleles(), 2);
        Assert.assertEquals(mtx.evidenceCount(), 1);
        final double v1 = mtx.get(0, 0);
        final double v2 = mtx.get(1, 0);

        Assert.assertTrue(v1 > v2, "matching haplotype should have a higher likelihood");
        lce.close();
    }
}