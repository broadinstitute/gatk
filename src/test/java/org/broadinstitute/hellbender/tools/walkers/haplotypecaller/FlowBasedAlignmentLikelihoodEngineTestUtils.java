package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class FlowBasedAlignmentLikelihoodEngineTestUtils {

    public static AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(
            final List<Haplotype> haplotypeList, final List<GATKRead> reads,
            final boolean filterPoorly, final SAMFileHeader hdr,
            final FlowBasedAlignmentLikelihoodEngine engine) {


        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        final ArrayList<String> _sampList = new ArrayList<>();
        _sampList.add("HG001");
        final SampleList samples = new IndexedSampleList(_sampList);

        // Add likelihoods for each sample's reads to our result
        final HashMap<String, List<GATKRead>> perSampleReadList = new HashMap<>();
        perSampleReadList.put("HG001", reads);

        final AlleleLikelihoods<GATKRead, Haplotype> result = new AlleleLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            engine.computeReadLikelihoods(result.sampleMatrix(i), hdr);
        }

        result.normalizeLikelihoods(engine.getLog10globalReadMismappingRate(), engine.isSymmetricallyNormalizeAllelesToReference());
        if ( filterPoorly ) {
            result.filterPoorlyModeledEvidence(engine.log10MinTrueLikelihood(engine.getExpectedErrorRatePerBase(), false));
        }

        return result;
    }
}
