package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.walkers.contamination.PileupSummaryCollector;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine.PileupQualBuffer;
import org.broadinstitute.hellbender.tools.walkers.readorientation.F1R2CountsCollector;
import org.broadinstitute.hellbender.tools.walkers.readorientation.F1R2CountsCollector.PileupStatistics;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Deque;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

/**
 * What Mutect2 needs to know about the pileups at the loci of a contig, accumulated one read at a time.  Rather than
 * building a pileup at every locus, Mutect2 walks each read's alignment once, when the read arrives, and adds what its
 * pileup elements tell it to what we know of its sample's pileup at each locus: the number of elements (bases and
 * deletions), the qualities of the elements that are evidence of an alt allele, as isActive used to accumulate them, and,
 * when Mutect2 collects F1R2 counts, the {@link PileupStatistics} of each tumor sample.  When Mutect2 collects pileup
 * summaries it also counts the bases of each table's pileup at the common variant sites, which it looks up before the
 * reads that cover them arrive.
 *
 * The elements are those of the pileups of {@link org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState}
 * with deletions: one for each aligned base and each deleted reference base, except inside the adaptor
 * ({@link ReadUtils#isBaseInsideAdaptor}).  Reads must be added in order of their starts, as they are to pileups, so the
 * alt qualities of each sample's pileup come in the pileup's order.  Once every read starting at or before a locus has
 * been added, the locus's evidence is complete.
 */
final class LocusEvidence {
    private static final int INITIAL_CAPACITY = 1 << 12;

    // the reference bases are fetched in blocks of at least this many bases, and the pileup summary sites are looked up
    // in blocks of this many loci
    private static final int REFERENCE_BLOCK_SIZE = 100_000;
    private static final int SITE_BLOCK_SIZE = 10_000;

    private final SAMFileHeader header;
    private final String contig;
    private final int contigLength;
    private final ReferenceDataSource reference;
    private final FeatureManager features;

    // the samples, in the order of LocusIteratorByState's pileups, the index of each, and which are normal
    private final List<String> samples;
    private final Map<String, Integer> sampleIndices;
    private final boolean[] isNormal;

    private final int halfPcrSnvQual;
    private final double multipleSubstitutionBaseQualCorrection;

    // the reference bases of [referenceStart, referenceStart + referenceBases.length)
    private int referenceStart = 1;
    private byte[] referenceBases = new byte[0];

    // What each sample's pileup at each locus tells us, at [sample][locus & (capacity - 1)], for the loci from start, the
    // first that hasn't been discarded, to start + capacity, a power of 2: its size, its alt qualities, if it has any, and
    // its F1R2 statistics, if we collect the sample's F1R2 counts.  Most loci have no alt evidence, and their arrays
    // of alt qualities come from a pool of spares, so that the dense depths are all most loci touch.
    private int start;
    private int capacity = INITIAL_CAPACITY;
    private final int[][] depths;
    private final PileupQualBuffer[][] altQuals;
    private final PileupStatistics[][] f1r2Statistics;
    private final Deque<PileupQualBuffer> spareAltQuals = new ArrayDeque<>();

    // no alt qualities, and the alt qualities of several samples' pileups at a locus, combined
    private final PileupQualBuffer noAltQuals;
    private final PileupQualBuffer combinedAltQuals;

    private final F1R2CountsCollector f1r2CountsCollector;
    private final int f1r2MinBaseQuality;

    // pileup summaries: the collector or null, and the candidate sites from the first at or after start, looked up through
    // sitesKnownThrough, with the base counts of each table's pileup at each
    private final PileupSummaryCollector pileupSummaryCollector;
    private final NavigableMap<Integer, Site> sites = new TreeMap<>();
    private int sitesKnownThrough;

    /**
     * @param start the first locus whose evidence is recorded
     * @param normalSamples the normal samples; the others are tumor samples
     * @param f1r2CountsCollector the collector of the tumor samples' F1R2 counts, or null if Mutect2 doesn't collect them
     * @param pileupSummaryCollector the collector of pileup summaries, or null if Mutect2 doesn't collect them
     */
    LocusEvidence(final SAMFileHeader header, final String contig, final int start, final ReferenceDataSource reference,
                  final FeatureManager features, final Set<String> normalSamples, final int pcrSnvQual,
                  final double multipleSubstitutionBaseQualCorrection, final F1R2CountsCollector f1r2CountsCollector,
                  final PileupSummaryCollector pileupSummaryCollector) {
        this.header = Utils.nonNull(header);
        this.contig = Utils.nonNull(contig);
        this.contigLength = header.getSequenceDictionary().getSequence(contig).getSequenceLength();
        this.start = start;
        this.reference = Utils.nonNull(reference);
        this.features = features;

        this.halfPcrSnvQual = pcrSnvQual / 2;
        this.multipleSubstitutionBaseQualCorrection = multipleSubstitutionBaseQualCorrection;
        noAltQuals = new PileupQualBuffer(multipleSubstitutionBaseQualCorrection);
        combinedAltQuals = new PileupQualBuffer(multipleSubstitutionBaseQualCorrection);
        this.f1r2CountsCollector = f1r2CountsCollector;
        f1r2MinBaseQuality = f1r2CountsCollector == null ? 0 : f1r2CountsCollector.getMinBaseQuality();

        samples = new ArrayList<>(ReadUtils.getSamplesFromHeader(header));
        sampleIndices = samples.stream().collect(Collectors.toMap(sample -> sample, samples::indexOf));
        isNormal = new boolean[samples.size()];
        depths = new int[samples.size()][capacity];
        altQuals = new PileupQualBuffer[samples.size()][capacity];
        f1r2Statistics = new PileupStatistics[samples.size()][];
        for (int sample = 0; sample < samples.size(); sample++) {
            isNormal[sample] = normalSamples.contains(samples.get(sample));
            if (f1r2CountsCollector != null && !isNormal[sample]) {
                f1r2Statistics[sample] = new PileupStatistics[capacity];
                for (int slot = 0; slot < capacity; slot++) {
                    f1r2Statistics[sample][slot] = f1r2CountsCollector.makePileupStatistics();
                }
            }
        }
        this.pileupSummaryCollector = pileupSummaryCollector;
        sitesKnownThrough = start - 1;
    }

    /**
     * A pileup summary site and the counts of A, C, G and T and of all elements in each table's pileup there
     */
    private record Site(VariantContext variant, int[][] baseCounts, int[] numElements) { }

    // ------------------------------------------------------------------------------------------------------------------
    // accumulating
    // ------------------------------------------------------------------------------------------------------------------

    /**
     * Add the pileup elements of a read, which must not start before any read added before it, at loci from {@link #start} on
     */
    void add(final GATKRead read) {
        final int readEnd = read.getEnd();
        if (readEnd < start) {
            return;
        }
        Utils.validateArg(read.getContig().equals(contig), () -> "read " + read + " is not on " + contig);
        ensureCapacity(readEnd);

        final Integer sample = sampleIndices.get(ReadUtils.getSampleName(read, header));
        Utils.validate(sample != null, () -> "read " + read + " has a sample that is not in the header");
        final int mask = capacity - 1;
        final int[] sampleDepths = depths[sample];
        final PileupQualBuffer[] sampleAltQuals = altQuals[sample];
        final PileupStatistics[] sampleF1R2Statistics = f1r2Statistics[sample];
        final boolean isF1R2 = ReadUtils.isF1R2(read);
        final int mappingQuality = read.getMappingQuality();

        // the loci of the read in pileups: from start on, and outside the adaptor, which, if the fragment size is well
        // defined and at most the adaptor size, covers the loci at most its boundary for reverse-strand reads and at least
        // it for forward-strand reads
        final int adaptorBoundary = read.getAdaptorBoundary();
        final boolean hasAdaptor = adaptorBoundary != ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY && ReadUtils.isBaseInsideAdaptor(read, adaptorBoundary);
        final int firstLocus = Math.max(Math.max(start, read.getStart()), hasAdaptor && read.isReverseStrand() ? adaptorBoundary + 1 : 0);
        final int lastLocus = hasAdaptor && !read.isReverseStrand() ? Math.min(readEnd, adaptorBoundary - 1) : readEnd;
        if (firstLocus > lastLocus) {
            return;
        }
        if (firstLocus < referenceStart || lastLocus >= referenceStart + referenceBases.length) {
            referenceStart = firstLocus;
            referenceBases = reference.queryAndPrefetch(contig, firstLocus, Math.min(contigLength, Math.max(lastLocus, firstLocus + REFERENCE_BLOCK_SIZE - 1))).getBases();
        }

        // isActive caps the qualities of substitutions where the read overlaps its mate, which it reckons by the mate's
        // start and this read's length, to half the PCR qual, so that the pair's evidence adds up to the PCR qual
        final int mateStart = (!read.isProperlyPaired() || read.mateIsUnmapped()) ? Integer.MAX_VALUE : read.getMateStart();
        final int readLength = read.getLength();

        final byte[] bases = read.getBasesNoCopy();
        final byte[] quals = read.getBaseQualitiesNoCopy();
        final List<CigarElement> cigar = read.getCigarElements();
        int elementStart = read.getStart();     // the locus of the current element's first on-genome base
        int elementOffset = 0;                  // the offset into the read of the current element's first base

        for (int k = 0; k < cigar.size(); k++) {
            final CigarOperator operator = cigar.get(k).getOperator();
            final int length = cigar.get(k).getLength();
            final int from = Math.max(elementStart, firstLocus);
            final int to = Math.min(elementStart + length - 1, lastLocus);
            switch (operator) {
                case M, EQ, X -> {
                    // what this element's first and last bases are next to: PileupElement's isAfterSoftClip and
                    // isAfterInsertion, isBeforeSoftClip and isBeforeDeletionStart, and the indel whose length
                    // getLengthOfImmediatelyFollowingIndel gives
                    final CigarOperator previous = k > 0 ? cigar.get(k - 1).getOperator() : null;
                    final CigarOperator next = k + 1 < cigar.size() ? cigar.get(k + 1).getOperator() : null;
                    final int followingDeletion = followingDeletionLength(cigar, k);
                    final int followingIndel = followingDeletion > 0 ? followingDeletion : (next == CigarOperator.I ? cigar.get(k + 1).getLength() : 0);
                    final int elementEnd = elementStart + length - 1;
                    for (int locus = from; locus <= to; locus++) {
                        sampleDepths[locus & mask]++;
                    }
                    for (int locus = from; locus <= to; locus++) {
                        // without F1R2 statistics to accumulate, skip the bases inside the element that match the
                        // reference, which are no evidence, comparing many bases at once
                        if (sampleF1R2Statistics == null && locus > elementStart && locus < elementEnd) {
                            final int insideEnd = Math.min(to, elementEnd - 1);
                            final int mismatch = Arrays.mismatch(bases, elementOffset + locus - elementStart, elementOffset + insideEnd - elementStart + 1,
                                    referenceBases, locus - referenceStart, insideEnd - referenceStart + 1);
                            if (mismatch < 0) {
                                locus = insideEnd;
                                continue;
                            }
                            locus += mismatch;
                        }
                        final int offset = elementOffset + locus - elementStart;
                        final boolean atStart = locus == elementStart;
                        final boolean atEnd = locus == elementEnd;
                        final byte qual = quals[offset];
                        final byte base = bases[offset];
                        final int slot = locus & mask;

                        // the alt evidence of isActive: an indel right after this base, a soft clip next to it whose base
                        // isn't junk, or a substitution
                        if (atEnd && followingIndel > 0) {
                            altQualsToAddTo(sampleAltQuals, slot).accumulateIndel(Mutect2Engine.indelQual(followingIndel));
                        } else if (qual > Mutect2Engine.MINIMUM_BASE_QUALITY &&
                                ((atEnd && next == CigarOperator.S && quals[offset + 1] > Mutect2Engine.MINIMUM_BASE_QUALITY) ||
                                        (atStart && previous == CigarOperator.S && quals[offset - 1] > Mutect2Engine.MINIMUM_BASE_QUALITY))) {
                            altQualsToAddTo(sampleAltQuals, slot).accumulateIndel(Mutect2Engine.indelQual(1));
                        } else if (base != referenceBases[locus - referenceStart] && qual > Mutect2Engine.MINIMUM_BASE_QUALITY) {
                            final boolean overlapsMate = mateStart <= locus && locus < mateStart + readLength;
                            altQualsToAddTo(sampleAltQuals, slot).accumulateSubstitution(base, overlapsMate ? (byte) FastMath.min(qual, halfPcrSnvQual) : qual);
                        }

                        if (sampleF1R2Statistics != null && qual > f1r2MinBaseQuality) {
                            sampleF1R2Statistics[slot].add(false, base, isF1R2,
                                    (atStart && previous == CigarOperator.I) || (atEnd && followingDeletion > 0), mappingQuality);
                        }
                    }
                    elementStart += length;
                    elementOffset += length;
                }
                case D -> {
                    for (int locus = from; locus <= to; locus++) {
                        final int slot = locus & mask;
                        sampleDepths[slot]++;
                        altQualsToAddTo(sampleAltQuals, slot).accumulateIndel(Mutect2Engine.indelQual(length));
                        if (sampleF1R2Statistics != null && PileupElement.DELETION_QUAL > f1r2MinBaseQuality) {
                            sampleF1R2Statistics[slot].add(true, PileupElement.DELETION_BASE, isF1R2, true, mappingQuality);
                        }
                    }
                    elementStart += length;
                }
                case N -> elementStart += length;  // pileups leave out reference skips
                case I, S -> elementOffset += length;
                default -> { }    // hard clips and padding
            }
        }

        if (pileupSummaryCollector != null) {
            countSiteElements(read, firstLocus, lastLocus);
        }
    }

    /**
     * The length of the deletion that is the next on-genome element after an element of a cigar, or 0 if the next
     * on-genome element isn't a deletion, as {@link PileupElement#isBeforeDeletionStart} finds it
     */
    private static int followingDeletionLength(final List<CigarElement> cigar, final int k) {
        for (int n = k + 1; n < cigar.size(); n++) {
            final CigarOperator operator = cigar.get(n).getOperator();
            if (operator == CigarOperator.D) {
                return cigar.get(n).getLength();
            } else if (operator.isAlignment()) {
                return 0;
            }
        }
        return 0;
    }

    /**
     * The alt qualities of a sample's pileup at a slot, from the spares if it has none yet
     */
    private PileupQualBuffer altQualsToAddTo(final PileupQualBuffer[] sampleAltQuals, final int slot) {
        if (sampleAltQuals[slot] == null) {
            sampleAltQuals[slot] = spareAltQuals.isEmpty() ? new PileupQualBuffer(multipleSubstitutionBaseQualCorrection) : spareAltQuals.pop();
        }
        return sampleAltQuals[slot];
    }

    /**
     * Make room for the loci through a locus, moving what we know of the loci we have to their slots in larger arrays
     */
    private void ensureCapacity(final int locus) {
        if (locus < start + capacity) {
            return;
        }
        final int newCapacity = Math.max(2 * capacity, Integer.highestOneBit(locus - start) << 1);
        for (int sample = 0; sample < samples.size(); sample++) {
            final int[] newDepths = new int[newCapacity];
            final PileupQualBuffer[] newAltQuals = new PileupQualBuffer[newCapacity];
            final PileupStatistics[] newF1R2Statistics = f1r2Statistics[sample] == null ? null : new PileupStatistics[newCapacity];
            for (int l = start; l < start + newCapacity; l++) {
                final int from = l & (capacity - 1);
                final int to = l & (newCapacity - 1);
                final boolean have = l < start + capacity;
                newDepths[to] = have ? depths[sample][from] : 0;
                newAltQuals[to] = have ? altQuals[sample][from] : null;
                if (newF1R2Statistics != null) {
                    newF1R2Statistics[to] = have ? f1r2Statistics[sample][from] : f1r2CountsCollector.makePileupStatistics();
                }
            }
            depths[sample] = newDepths;
            altQuals[sample] = newAltQuals;
            f1r2Statistics[sample] = newF1R2Statistics;
        }
        capacity = newCapacity;
    }

    /**
     * Free the loci before a locus, which must not be after the loci of reads still to come
     */
    void discardLociBefore(final int locus) {
        if (locus <= start) {
            return;
        }
        for (int l = start; l < Math.min(locus, start + capacity); l++) {
            final int slot = l & (capacity - 1);
            for (int sample = 0; sample < samples.size(); sample++) {
                depths[sample][slot] = 0;
                if (altQuals[sample][slot] != null) {
                    altQuals[sample][slot].clear();
                    spareAltQuals.push(altQuals[sample][slot]);
                    altQuals[sample][slot] = null;
                }
                if (f1r2Statistics[sample] != null) {
                    f1r2Statistics[sample][slot].reset();
                }
            }
        }
        start = locus;
        while (!sites.isEmpty() && sites.firstKey() < locus) {
            sites.pollFirstEntry();
        }
        // don't look up sites in a gap between intervals
        sitesKnownThrough = Math.max(sitesKnownThrough, locus - 1);
    }

    // ------------------------------------------------------------------------------------------------------------------
    // pileup summary sites
    // ------------------------------------------------------------------------------------------------------------------

    /**
     * Count the elements of a read's pileups, from firstLocus to lastLocus, at pileup summary sites in the tables that include it
     */
    private void countSiteElements(final GATKRead read, final int firstLocus, final int lastLocus) {
        ensureSitesThrough(lastLocus);
        final Collection<Site> readSites = sites.subMap(firstLocus, true, lastLocus, true).values();
        if (readSites.isEmpty()) {
            return;
        }
        final boolean[] tablesIncludingRead = new boolean[pileupSummaryCollector.numberOfTables()];
        boolean anyTable = false;
        for (int table = 0; table < tablesIncludingRead.length; table++) {
            tablesIncludingRead[table] = pileupSummaryCollector.tableIncludes(table, read);
            anyTable |= tablesIncludingRead[table];
        }
        if (!anyTable) {
            return;
        }

        final int softStart = read.getSoftStart();
        final Cigar cigar = read.getCigar();
        for (final Site site : readSites) {
            final Pair<Integer, CigarOperator> offsetAndOperator = ReadUtils.getReadIndexForReferenceCoordinate(softStart, cigar, site.variant().getStart());
            if (offsetAndOperator.getRight() == CigarOperator.N) {
                continue;   // pileups leave out reference skips
            }
            final int baseIndex = offsetAndOperator.getRight() == CigarOperator.D ? -1 : BaseUtils.simpleBaseToBaseIndex(read.getBase(offsetAndOperator.getLeft()));
            for (int table = 0; table < tablesIncludingRead.length; table++) {
                if (tablesIncludingRead[table]) {
                    site.numElements()[table]++;
                    if (baseIndex != -1) {
                        site.baseCounts()[table][baseIndex]++;
                    }
                }
            }
        }
    }

    /**
     * Look up the pileup summary sites of the loci through a locus, if we haven't
     */
    private void ensureSitesThrough(final int locus) {
        while (sitesKnownThrough < Math.min(locus, contigLength)) {
            final SimpleInterval block = new SimpleInterval(contig, sitesKnownThrough + 1, Math.min(contigLength, sitesKnownThrough + SITE_BLOCK_SIZE));
            final int numTables = pileupSummaryCollector.numberOfTables();
            for (final VariantContext variant : pileupSummaryCollector.candidateSites(new FeatureContext(features, block))) {
                sites.put(variant.getStart(), new Site(variant, new int[numTables][4], new int[numTables]));
            }
            sitesKnownThrough = block.getEnd();
        }
    }

    /**
     * Write the pileup summaries of a complete locus if it is a site with pileup summaries
     */
    void writePileupSummaries(final int locus) {
        if (pileupSummaryCollector == null) {
            return;
        }
        ensureSitesThrough(locus);
        final Site site = sites.get(locus);
        if (site == null || !pileupSummaryCollector.summarizes(site.variant())) {
            return;
        }
        for (int table = 0; table < site.numElements().length; table++) {
            if (site.numElements()[table] > 0) {
                pileupSummaryCollector.writeRecord(table, site.variant(), site.baseCounts()[table]);
            }
        }
    }

    // ------------------------------------------------------------------------------------------------------------------
    // reading the evidence of complete loci
    // ------------------------------------------------------------------------------------------------------------------

    /**
     * The total size of the tumor or normal samples' pileups at a locus
     */
    int depth(final int locus, final boolean normal) {
        int result = 0;
        for (int sample = 0; sample < samples.size(); sample++) {
            result += isNormal[sample] == normal ? depths[sample][locus & (capacity - 1)] : 0;
        }
        return result;
    }

    /**
     * The alt qualities of the tumor or normal samples' pileups at a locus, in the order of the samples and then of their
     * pileups' elements.  The result must not be modified, and is valid until the next call.
     */
    PileupQualBuffer altQuals(final int locus, final boolean normal) {
        PileupQualBuffer only = noAltQuals;
        int count = 0;
        for (int sample = 0; sample < samples.size(); sample++) {
            final PileupQualBuffer sampleAltQuals = altQuals[sample][locus & (capacity - 1)];
            if (isNormal[sample] == normal && sampleAltQuals != null) {
                only = sampleAltQuals;
                count++;
            }
        }
        if (count <= 1) {
            return only;
        }
        combinedAltQuals.clear();
        for (int sample = 0; sample < samples.size(); sample++) {
            final PileupQualBuffer sampleAltQuals = altQuals[sample][locus & (capacity - 1)];
            if (isNormal[sample] == normal && sampleAltQuals != null) {
                combinedAltQuals.addAll(sampleAltQuals);
            }
        }
        return combinedAltQuals;
    }

    /**
     * The F1R2 statistics of a tumor sample's pileup at a locus
     */
    @VisibleForTesting
    PileupStatistics f1r2Statistics(final String sample, final int locus) {
        return f1r2Statistics[sampleIndices.get(sample)][locus & (capacity - 1)];
    }

    /**
     * Collect the F1R2 counts of a complete locus from the tumor samples' statistics there
     *
     * @param referenceContext a reference context on the contig, whose data source supplies reference contexts
     */
    void collectF1R2Counts(final int locus, final ReferenceContext referenceContext) {
        for (int sample = 0; sample < samples.size(); sample++) {
            if (f1r2Statistics[sample] != null) {
                f1r2CountsCollector.process(samples.get(sample), f1r2Statistics[sample][locus & (capacity - 1)], contig, locus, referenceContext);
            }
        }
    }
}
