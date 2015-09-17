package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * A haplotype bam writer that writes out all haplotypes as reads and then
 * the alignment of reach read to its best match among the best haplotypes.
 *
 * Primarily useful for people working on the HaplotypeCaller method itself
 */
final class AllHaplotypeBAMWriter extends HaplotypeBAMWriter {

    public AllHaplotypeBAMWriter(final HaplotypeBAMDestination destination) {
        super(destination);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void writeReadsAlignedToHaplotypes(final Collection<Haplotype> haplotypes,
                                              final Locatable paddedReferenceLoc,
                                              final Collection<Haplotype> bestHaplotypes,
                                              final Set<Haplotype> calledHaplotypes,
                                              final ReadLikelihoods<Haplotype> readLikelihoods) {

        Utils.nonNull(haplotypes, "haplotypes cannot be null");
        Utils.nonNull(paddedReferenceLoc, "paddedReferenceLoc cannot be null");
        Utils.nonNull(readLikelihoods, "readLikelihoods cannot be null");

        writeHaplotypesAsReads(haplotypes, new HashSet<>(bestHaplotypes), paddedReferenceLoc);

        final int sampleCount = readLikelihoods.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            for (final GATKRead read : readLikelihoods.sampleReads(i)) {
                writeReadAgainstHaplotype(read);
            }
        }
    }
}
