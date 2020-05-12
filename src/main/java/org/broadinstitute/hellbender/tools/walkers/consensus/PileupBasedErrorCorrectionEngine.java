package org.broadinstitute.hellbender.tools.walkers.consensus;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.LocationAndAlleles;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.util.*;
import java.util.stream.Collectors;

public class PileupBasedErrorCorrectionEngine extends ErrorCorrectionEngineBase {
    public PileupBasedErrorCorrectionEngine(SAMFileHeader header, ReferenceInputArgumentCollection referenceArguments, SAMFileGATKReadWriter outputWriter) {
        super(header, referenceArguments, outputWriter);

    }

    public GATKRead getIndelConsensusRead(final List<GATKRead> reads, final ReferenceContext referenceContext) {
        // Can I just go from reads to haplotypes to events? Yes.
        final int regionStart = reads.stream().mapToInt(GATKRead::getStart).min().getAsInt();
        final int regionEnd = reads.stream().mapToInt(GATKRead::getEnd).max().getAsInt();

        // We need padding because in the event of insertion we might need bases in the reference that falls outside of the tight window
        final int REFERENCE_PADDING = 50;

        final SimpleInterval refInterval = new SimpleInterval(referenceContext.getContig(),
                Math.max(regionStart - REFERENCE_PADDING, 0),
                Math.min(regionEnd + REFERENCE_PADDING, referenceContext.getContig().length()));
        final SimpleInterval justCurious = referenceContext.getInterval();

        final List<ReadAndHaplotype> readAndHaplotypes = getReadAndHaplotypes(reads, refInterval);

        /*** BOOM ***/
        final List<Pair<LocationAndAlleles, Double>> indelEvents = findIndelEventsAndVariances(readAndHaplotypes, refInterval);

        final int readLength = reads.get(0).getLength();
        final byte[] newInsertionQualities = new byte[readLength];
        final byte[] newDeletionQualities = new byte[readLength];
        Arrays.fill(newInsertionQualities, ReadUtils.DEFAULT_INSERTION_DELETION_QUAL);
        Arrays.fill(newDeletionQualities, ReadUtils.DEFAULT_INSERTION_DELETION_QUAL);
        final Optional<GATKRead> consensusReadOption = findRepresentativeReadWithDesiredIndel(readAndHaplotypes, indelEvents);

        if (!consensusReadOption.isPresent()){
            // TODO: handle this case
            throw new UserException("WHAAAA");
        }

        boolean updatedInsertionQuality = false;
        boolean updatedDeletionQuality = false;

        /*** Update the insertion/deletion qualities at the position of insertion/deletion ***/
        final GATKRead consensusRead = consensusReadOption.get();
        final int readStart = consensusRead.getStart();
        for (Pair<LocationAndAlleles, Double> pair : indelEvents){
            final LocationAndAlleles event = pair.getLeft();
            final double variance = pair.getRight();
            final int offset = event.getLoc() - readStart; // off by one?
            if (event.isInsertion(1)){
                newInsertionQualities[offset] = varianceToIndelQuality(variance);
                updatedInsertionQuality = true;
            } else {
                newDeletionQualities[offset] = varianceToIndelQuality(variance);
                updatedDeletionQuality = true;
            }
        }

        if (updatedDeletionQuality) consensusRead.setAttribute(ReadUtils.BQSR_BASE_DELETION_QUALITIES, newDeletionQualities);
        if (updatedInsertionQuality) consensusRead.setAttribute(ReadUtils.BQSR_BASE_INSERTION_QUALITIES, newInsertionQualities);

        return consensusRead;
    }

    private List<Pair<LocationAndAlleles, Double>> findIndelEventsAndVariances(final List<ReadAndHaplotype> readAndHaplotypes, SimpleInterval refInterval){
        final List<Haplotype> haplotypes = ReadAndHaplotype.getHaplotypeList(readAndHaplotypes);

        // TODO: use referenceContext Instead.
        final ReferenceSequenceFile referenceReader = AssemblyBasedCallerUtils.createReferenceReader(Utils.nonNull(referenceArguments.getReferenceFileName()));
        final byte[] referenceBases = AssemblyRegion.getReference(referenceReader, 0, refInterval);

        // I'd have to write down ref, ref position, and what they are set to in the case of with assembly.
        // And then identify what kind of padding, refloc, etc., needs to be set to by hand.

        // Identify all the events using buildEventMapsForHaplotypes: why is this necessary?
        final TreeSet<Integer> startPositions = EventMap.buildEventMapsForHaplotypes(haplotypes, referenceBases, refInterval, true, 1);

        final List<Pair<LocationAndAlleles, Double>> indelEvents = new ArrayList<>();
        for (int start : startPositions){
            final List<VariantContext> eventsAtThisLoc = AssemblyBasedCallerUtils.getVariantContextsFromActiveHaplotypes(start, haplotypes, false);
            final List<VariantContext> indelEventsAtThisLoc = eventsAtThisLoc.stream().filter(vc -> vc.isIndel()).collect(Collectors.toList());
            if (indelEventsAtThisLoc.isEmpty()){
                continue;
            }

            final Map<LocationAndAlleles, List<Haplotype>> haplotypesByAllele = CorrectErrorsUtils.groupHaplotypesByAllelesAtThisLoc(start, haplotypes, false);

            Pair<LocationAndAlleles, Integer> mostCommonAlleleAndCount = new ImmutablePair<>(null, -1);
            for (Map.Entry<LocationAndAlleles, List<Haplotype>> pairs : haplotypesByAllele.entrySet()){
                final LocationAndAlleles locationAndAlleles = pairs.getKey();
                final int count = pairs.getValue().size();
                if (count > mostCommonAlleleAndCount.getRight()){
                    mostCommonAlleleAndCount = new ImmutablePair<>(locationAndAlleles, count);
                }
            }

            final double variance = CorrectErrorsUtils.computeVarianceAroundMostCommon(haplotypesByAllele, mostCommonAlleleAndCount.getLeft());
            indelEvents.add(new ImmutablePair<>(mostCommonAlleleAndCount.getLeft(), variance));
            int f = 3; // GOOD!
        }
        return indelEvents;
    }

    public static byte varianceToIndelQuality(final double variance) {
        // This is completely adhoc for now, and that's OK.
        return (byte) Math.max(ReadUtils.DEFAULT_INSERTION_DELETION_QUAL - variance, 10);
    }

    /**
     * Creating the new read by hand seems error-prone, and perhaps not necessary for prototyping.
     */
    private Optional<GATKRead> findRepresentativeReadWithDesiredIndel(final List<ReadAndHaplotype> readAndHaplotypes,
                                                                      final List<Pair<LocationAndAlleles, Double>> indelEvents){
        Utils.validateArg(!indelEvents.isEmpty(), "indelEvents may not be empty");
        for (final ReadAndHaplotype rah : readAndHaplotypes) {
            final Collection<VariantContext> vcs = rah.getHaplotype().getEventMap().getVariantContexts();
            final List<LocationAndAlleles> laas = vcs.stream().filter(VariantContext::isIndel)
                    .map(vc -> new LocationAndAlleles(vc.getStart(), vc.getAlleles()))
                    .collect(Collectors.toList());
            if (indelEvents.stream().allMatch(pair -> laas.contains(pair.getLeft()))){
                final GATKRead consensusRead = rah.getRead();
                return Optional.of(consensusRead);
            }
        }
        return Optional.empty();
    }
}
