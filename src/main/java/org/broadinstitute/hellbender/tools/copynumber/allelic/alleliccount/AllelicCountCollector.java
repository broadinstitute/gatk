package org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount;

import htsjdk.samtools.*;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.NotPrimaryAlignmentFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Collects reference/alternate allele counts at sites.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicCountCollector {

    private static final Logger logger = LogManager.getLogger(AllelicCountCollector.class);

    public static final List<Nucleotide> BASES = Collections.unmodifiableList(Arrays.asList(Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.T));

    private static final int MAX_INTERVALS_FOR_INDEX = 25000;  //set interval threshold for indexing for SamLocusIterator
    private static final int NUMBER_OF_SITES_PER_LOGGED_STATUS_UPDATE = 10000;

    private final SamReaderFactory readerFactory;
    private final ReferenceSequenceFileWalker referenceWalker;

    /**
     * Constructs a {@link AllelicCountCollector} object for calculating {@link AllelicCountCollection} objects from
     * BAMs based on a reference genome.
     * @param referenceFile         file containing the reference
     * @param validationStringency  validation stringency to use for reading BAM files
     */
    public AllelicCountCollector(final File referenceFile,
                                 final ValidationStringency validationStringency) {
        Utils.regularReadableUserFile(referenceFile);
        readerFactory = SamReaderFactory.makeDefault().validationStringency(validationStringency).referenceSequence(referenceFile);
        referenceWalker = new ReferenceSequenceFileWalker(referenceFile);
    }

    /**
     * Returns an {@link AllelicCountCollection} based on the pileup at sites (specified by an interval list)
     * in a sorted BAM file.  Reads and bases below the specified mapping quality and base quality, respectively,
     * are filtered out of the pileup.  The alt count is defined as the total count minus the ref count, and the
     * alt nucleotide is defined as the non-ref base with the highest count, with ties broken by the order of the
     * bases in {@link AllelicCountCollector#BASES}.
     * @param bamFile           sorted BAM file
     * @param siteIntervals     interval list of sites
     * @param minMappingQuality minimum mapping quality required for reads to be included in pileup
     * @param minBaseQuality    minimum base quality required for bases to be included in pileup
     * @return                  AllelicCountCollection of ref/alt counts at sites in BAM file
     */
    public AllelicCountCollection collect(final File bamFile,
                                          final IntervalList siteIntervals,
                                          final int minMappingQuality,
                                          final int minBaseQuality) {
        try (final SamReader reader = readerFactory.open(bamFile)) {
            ParamUtils.isPositiveOrZero(minMappingQuality, "Minimum mapping quality must be nonnegative.");
            ParamUtils.isPositiveOrZero(minBaseQuality, "Minimum base quality must be nonnegative.");
            if (reader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new UserException.BadInput("BAM file " + bamFile.toString() + " must be coordinate sorted.");
            }

            final int numberOfSites = siteIntervals.size();
            final boolean useIndex = numberOfSites < MAX_INTERVALS_FOR_INDEX;
            final SamLocusIterator locusIterator = new SamLocusIterator(reader, siteIntervals, useIndex);

            //set read and locus filters [note: read counts match IGV, but off by a few from pysam.mpileup]
            final List<SamRecordFilter> samFilters = Arrays.asList(new NotPrimaryAlignmentFilter(), new DuplicateReadFilter());
            locusIterator.setSamFilters(samFilters);
            locusIterator.setEmitUncoveredLoci(true);
            locusIterator.setIncludeNonPfReads(false);
            locusIterator.setMappingQualityScoreCutoff(minMappingQuality);
            locusIterator.setQualityScoreCutoff(minBaseQuality);

            logger.info("Examining " + numberOfSites + " sites in total...");
            int locusCount = 0;
            final AllelicCountCollection counts = new AllelicCountCollection();
            for (final SamLocusIterator.LocusInfo locus : locusIterator) {
                if (locusCount % NUMBER_OF_SITES_PER_LOGGED_STATUS_UPDATE == 0) {
                    logger.info("Examined " + locusCount + " sites.");
                }
                locusCount++;

                final Nucleotide refBase = Nucleotide.valueOf(referenceWalker.get(locus.getSequenceIndex()).getBases()[locus.getPosition() - 1]);
                if (!BASES.contains(refBase)) {
                    logger.warn(String.format("The reference position at %d has an unknown base call (value: %s). Skipping...",
                            locus.getPosition(), refBase.toString()));
                    continue;
                }

                final Nucleotide.Counter baseCounts = getPileupBaseCounts(locus);
                final int totalBaseCount = BASES.stream().mapToInt(b -> (int) baseCounts.get(b)).sum(); //only include total ACGT counts in binomial test (exclude N, etc.)
                final int refReadCount = (int) baseCounts.get(refBase);
                final int altReadCount = totalBaseCount - refReadCount;                                 //we take alt = total - ref instead of the actual alt count
                final Nucleotide altBase = inferAltFromPileupBaseCounts(baseCounts, refBase);

                counts.add(new AllelicCount(
                        new SimpleInterval(locus.getSequenceName(), locus.getPosition(), locus.getPosition()),
                        refReadCount, altReadCount, refBase, altBase));
            }
            logger.info(locusCount + " sites out of " + numberOfSites + " total sites were examined.");
            return counts;
        } catch (final IOException | SAMFormatException e) {
            throw new UserException("Unable to collect allelic counts from " + bamFile);
        }
    }

    /**
     * Returns base-pair counts at a given locus.
     */
    private static Nucleotide.Counter getPileupBaseCounts(final SamLocusIterator.LocusInfo locus) {
        final Nucleotide.Counter result = new Nucleotide.Counter();
        locus.getRecordAndOffsets().stream().forEach(r -> result.add(r.getReadBase()));
        return result;
    }

    /**
     * Returns the non-ref base with highest count (if there is a tie, the first base in the order given in
     * {@link AllelicCountCollector#BASES} will be returned).
     */
    private static Nucleotide inferAltFromPileupBaseCounts(final Nucleotide.Counter baseCounts,
                                                           final Nucleotide refNucleotide) {
        return BASES.stream()
                .filter(b -> b != refNucleotide)
                .sorted((b1, b2) -> Long.compare(baseCounts.get(b1), baseCounts.get(b2)))
                .findFirst().get();
    }
}
