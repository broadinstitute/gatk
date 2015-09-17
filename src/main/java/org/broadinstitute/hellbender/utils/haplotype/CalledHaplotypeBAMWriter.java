package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.util.Locatable;

import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.Set;


/**
 * Writes a BAM containing just the reads in stratifiedReadMap aligned to their
 * most likely haplotype among all of the called haplotypes.
 *
 * Primarily useful for users of the HaplotypeCaller who want to better understand the
 * support of their calls w.r.t. the reads.
 *
 */
final class CalledHaplotypeBAMWriter extends HaplotypeBAMWriter {
    public CalledHaplotypeBAMWriter(final HaplotypeBAMDestination destination) {
        super(destination);
    }

    /**
     * Write out a BAM representing for the haplotype caller at this site
     *
     * @param haplotypes a list of all possible haplotypes at this loc, must not be null
     * @param paddedReferenceLoc the span of the based reference here, must not be null
     * @param bestHaplotypes a list of the best (a subset of all) haplotypes that actually went forward into genotyping,
     *                       must not be null
     * @param calledHaplotypes a list of the haplotypes at where actually called as non-reference, must not be null
     * @param readLikelihoods a map from sample -> likelihoods for each read for each of the best haplotypes
     */
    @Override
    public void writeReadsAlignedToHaplotypes(final Collection<Haplotype> haplotypes,
                                              final Locatable paddedReferenceLoc,
                                              final Collection<Haplotype> bestHaplotypes,
                                              final Set<Haplotype> calledHaplotypes,
                                              final ReadLikelihoods<Haplotype> readLikelihoods) {
        Utils.nonNull(haplotypes, "haplotypes cannot be null");
        Utils.nonNull(paddedReferenceLoc, "paddedReferenceLoc cannot be null");
        Utils.nonNull(calledHaplotypes, "calledHaplotypes cannot be null");
        Utils.nonNull(readLikelihoods, "readLikelihoods cannot be null");

        if (calledHaplotypes.isEmpty()) { // only write out called haplotypes
            return;
        }

        writeHaplotypesAsReads(calledHaplotypes, calledHaplotypes, paddedReferenceLoc);

        final int sampleCount = readLikelihoods.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            for (final GATKRead read : readLikelihoods.sampleReads(i)) {
                writeReadAgainstHaplotype(read);
            }
        }
    }
}
