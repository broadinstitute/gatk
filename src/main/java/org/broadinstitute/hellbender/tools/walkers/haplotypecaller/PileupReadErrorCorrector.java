package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class PileupReadErrorCorrector implements ReadErrorCorrector {
    private final double logOddsThreshold;
    private final SAMFileHeader header;
    private final List<Byte> altQualBuffer = new ArrayList<>();

    private static final byte GOOD_QUAL = 30;

    // if there are {@code INDEL_SPAN} mismatches with {@code INDEL_SPAN} bases at the end of the read, there may be an indel
    // and we don't correct bases
    private static final int INDEL_SPAN = 15;
    private static final int INDEL_MISMATCHES = 3;

    public PileupReadErrorCorrector(final double logOddsThreshold, final SAMFileHeader header) {
        this.logOddsThreshold = logOddsThreshold;
        this.header = header;
    }

    @Override
    public final List<GATKRead> correctReads(final Collection<GATKRead> originalReads) {
        final List<GATKRead> reads = originalReads.stream().map(GATKRead::deepCopy).collect(Collectors.toList());

        final Iterator<AlignmentContext> locusIterator = new LocusIteratorByState(reads.iterator(), DownsamplingMethod.NONE, ReadUtils.getSamplesFromHeader(header), header, false);

        final Map<GATKRead, List<Pair<Integer, Byte>>> potentialCorrections = reads.stream().collect(Collectors.toMap(read -> read, read -> new ArrayList<>()));

        Utils.stream(locusIterator).map(AlignmentContext::getBasePileup).forEach(pileup -> {
            final Nucleotide.Counter counter = new Nucleotide.Counter();
            pileup.forEach(pe -> counter.add(pe.getBase()));

            final Optional<Nucleotide> pluralityBase = Nucleotide.STANDARD_BASES.stream().max(Comparator.comparingLong(counter::get));
            if (!pluralityBase.isPresent()) {
                return;
            }
            // TODO: not really the ref
            final byte ref = pluralityBase.get().encodeAsByte();
            altQualBuffer.clear();

            int refCount = 0;
            for (final PileupElement pe : pileup) {
                // TODO: pe.getBase() == ref
                if (pe.getBase() == ref) {
                    refCount++;
                } else {
                    altQualBuffer.add(pe.getQual());
                }
            }

            // TODO: skip if altQual buffer is empty?
            final double logOdds = Mutect2Engine.logLikelihoodRatio(refCount, altQualBuffer);
            // TODO: what if there is a good variant and an error?
            if (logOdds < logOddsThreshold) {
                for (final PileupElement pe : pileup) {
                    if (pe.getBase() != ref && !(pe.isDeletion() || pe.isBeforeInsertion() || pe.isAfterDeletionEnd() || pe.isBeforeDeletionStart() || pe.isAfterInsertion() || pe.isAfterSoftClip())) {
                        potentialCorrections.get(pe.getRead()).add(ImmutablePair.of(pe.getOffset(), ref));
                    }
                }
            }

        });

        potentialCorrections.entrySet().forEach(entry -> {
            final GATKRead read = entry.getKey();
            final byte[] bases = read.getBasesNoCopy();
            final int length = bases.length;
            final byte[] quals = read.getBaseQualitiesNoCopy();
            final List<Pair<Integer, Byte>> edits = entry.getValue();
            final int size = edits.size();

            int firstEdit = 0;
            for (int n = 0; n + INDEL_MISMATCHES < size && edits.get(n + INDEL_MISMATCHES - 1).getLeft() - edits.get(n).getLeft() < INDEL_SPAN; n++) {
                firstEdit = n + INDEL_MISMATCHES;
            }

            int lastEdit = size - 1;
            for (int n = size - 1; n >= INDEL_MISMATCHES - 1 && edits.get(n).getLeft() - edits.get(n - INDEL_MISMATCHES + 1).getLeft()  < INDEL_SPAN; n--) {
                lastEdit = n - INDEL_MISMATCHES;
            }

            for (int n = firstEdit; n <= lastEdit; n++) {
                bases[edits.get(n).getLeft()] = edits.get(n).getRight();
                quals[edits.get(n).getLeft()] = GOOD_QUAL;
            }
        });

        return reads;
    }
}
