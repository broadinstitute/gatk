package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.positional;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

/**
 * EXPERIMENTAL: a local assembler whose kmer graph is anchored to reference positions instead of a de Bruijn graph.
 *
 * <p>Read kmers join graph nodes only in increasing order of rank, roughly the reference position, so the graph is
 * acyclic for any kmer size and needs no larger kmers in repeats.  Haplotypes are not paths through the graph:
 * the graph only turns reads into events, and the haplotypes are a small set of combinations of supported events
 * that is consistent with what fragments carry together (see {@link PositionalHaplotypeFinder}).</p>
 */
public final class PositionalAssembler {
    private static final Logger logger = LogManager.getLogger(PositionalAssembler.class);

    public static final int DEFAULT_KMER_SIZE = 10;
    /** as in {@link org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler} */
    private static final int MIN_HAPLOTYPE_REFERENCE_LENGTH = 30;

    private final PositionalAssemblyParameters params;
    private final int maxHaplotypes;
    private final boolean debug;

    /**
     * @param kmerSize kmer size of the graph
     * @param minPruning reads an event needs, and fragments a combination of alleles needs to make a new haplotype
     * @param adaptivePruning whether an event also needs the log odds against sequencing error that adaptive pruning
     *                        requires of a graph chain
     * @param initialErrorRateForPruning initial error rate of adaptive pruning
     * @param pruningLogOddsThreshold natural log odds threshold of adaptive pruning
     * @param minBaseQuality bases of lower quality are not used in kmers
     * @param minDanglingBranchLength a read's dangling end with fewer kmers than this only supports events that other reads make
     * @param maxHaplotypes maximum number of haplotypes, including the reference
     * @param eventScores alignment scores for turning reads into events: those of haplotype-to-reference alignment
     */
    public PositionalAssembler(final int kmerSize, final int minPruning, final boolean adaptivePruning, final double initialErrorRateForPruning,
                               final double pruningLogOddsThreshold, final int minBaseQuality, final int minDanglingBranchLength,
                               final int maxHaplotypes, final SWParameters eventScores, final boolean debug) {
        Utils.validateArg(maxHaplotypes > 0, "maxHaplotypes must be positive");
        this.params = PositionalAssemblyParameters.defaults(kmerSize, minPruning, adaptivePruning, initialErrorRateForPruning,
                pruningLogOddsThreshold, minBaseQuality, minDanglingBranchLength, eventScores);
        this.maxHaplotypes = maxHaplotypes;
        this.debug = debug;
    }

    /**
     * Assemble the reads of an assembly region whose reads have been finalized for assembly.
     *
     * @param refHaplotype the reference haplotype of the region's padded span
     * @param fullReferenceWithPadding reference bases of {@code refLoc}, which contains the padded span
     * @return the reference haplotype and the assembled haplotypes of the padded span, aligned to the reference
     */
    public AssemblyResultSet runLocalAssembly(final AssemblyRegion region, final Haplotype refHaplotype, final byte[] fullReferenceWithPadding,
                                              final SimpleInterval refLoc, final SAMFileHeader header, final SmithWatermanAligner aligner,
                                              final SWParameters haplotypeToReferenceSWParameters) {
        Utils.validateArg(fullReferenceWithPadding.length == refLoc.size(), "Reference bases and reference loc must be the same size.");
        final AssemblyResultSet resultSet = new AssemblyResultSet();
        resultSet.setRegionForGenotyping(region);
        resultSet.setFullReferenceWithPadding(fullReferenceWithPadding);
        resultSet.setPaddedReferenceLoc(refLoc);
        refHaplotype.setGenomeLocation(region.getPaddedSpan());

        final List<PositionalHaplotypeFinder.AssembledHaplotype> assembled = assemble(region.getReads(), fullReferenceWithPadding, refLoc.getStart(), header);

        final List<Haplotype> haplotypes = new ArrayList<>();
        haplotypes.add(refHaplotype);
        final byte[] refBases = refHaplotype.getBases();
        final int refLength = refHaplotype.getCigar().getReferenceLength();
        // window offset of the reference haplotype's first base
        final int offset = refHaplotype.getAlignmentStartHapwrtRef();
        int failed = 0;
        for (final PositionalHaplotypeFinder.AssembledHaplotype assembledHaplotype : assembled) {
            if (haplotypes.size() == maxHaplotypes) {
                break;
            }
            final byte[] bases = applyWithinSpan(refBases, offset, assembledHaplotype.events());
            if (Arrays.equals(bases, refBases)) {
                continue;
            }
            final Cigar cigar = CigarUtils.calculateCigar(refBases, bases, aligner, haplotypeToReferenceSWParameters, SWOverhangStrategy.SOFTCLIP);
            if (cigar == null || cigar.isEmpty() || cigar.getReferenceLength() != refLength || cigar.getReferenceLength() < MIN_HAPLOTYPE_REFERENCE_LENGTH
                    || cigar.getCigarElements().stream().anyMatch(ce -> ce.getOperator() == CigarOperator.N)) {
                failed++;
                continue;
            }
            final Haplotype haplotype = new Haplotype(bases, false);
            haplotype.setCigar(cigar);
            haplotype.setAlignmentStartHapwrtRef(refHaplotype.getAlignmentStartHapwrtRef());
            haplotype.setGenomeLocation(region.getPaddedSpan());
            haplotype.setKmerSize(params.k());
            if (!haplotypes.contains(haplotype)) {
                haplotypes.add(haplotype);
            }
        }
        if (debug) {
            logger.info(String.format(Locale.ROOT, "Positional assembly of %s: %d reads, %d haplotypes (%d failed alignment to the reference)",
                    region.getSpan(), region.getReads().size(), haplotypes.size(), failed));
        }
        resultSet.replaceAllHaplotypes(haplotypes);
        return resultSet;
    }

    /**
     * The haplotypes of reads over a reference window, as events in window coordinates, most supported first.
     *
     * @param window reference bases
     * @param windowStart 1-based reference position of the window's first base
     */
    List<PositionalHaplotypeFinder.AssembledHaplotype> assemble(final List<GATKRead> reads, final byte[] window, final int windowStart,
                                                                final SAMFileHeader header) {
        final List<PositionalRead> prepared = new ArrayList<>(reads.size());
        for (final GATKRead read : reads) {
            final PositionalRead positionalRead = PositionalRead.create(read, ReadUtils.getSampleName(read, header), windowStart,
                    params.k(), params.minBaseQuality());
            if (positionalRead != null) {
                prepared.add(positionalRead);
            }
        }
        final byte[] upperCaseWindow = new String(window, StandardCharsets.US_ASCII).toUpperCase(Locale.ROOT).getBytes(StandardCharsets.US_ASCII);
        return new PositionalHaplotypeFinder(params, upperCaseWindow, prepared).findHaplotypes();
    }

    /**
     * The haplotype's events applied to the reference haplotype, which starts at window offset {@code offset}.  Events
     * that don't leave a reference base on each side within the reference haplotype are dropped, as assembled
     * haplotypes always start and end with reference bases.
     */
    private static byte[] applyWithinSpan(final byte[] refBases, final int offset, final List<PositionalEvent> events) {
        final StringBuilder bases = new StringBuilder(refBases.length + 16);
        int next = 0;
        for (final PositionalEvent event : events) {
            final int start = event.position() - offset, end = start + event.ref().length();
            if (start < 1 || end > refBases.length - 1 || start < next) {
                continue;
            }
            if (!new String(refBases, start, end - start, StandardCharsets.US_ASCII).equalsIgnoreCase(event.ref())) {
                throw new IllegalStateException("Positional assembly event " + event + " does not match the reference");
            }
            bases.append(new String(refBases, next, start - next, StandardCharsets.US_ASCII)).append(event.alt());
            next = end;
        }
        bases.append(new String(refBases, next, refBases.length - next, StandardCharsets.US_ASCII));
        return bases.toString().getBytes(StandardCharsets.US_ASCII);
    }
}
