package org.broadinstitute.hellbender.tools.exome;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.NotPrimaryAlignmentFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Gets heterozygous SNP pulldown for normal and tumor samples.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class HetPulldownCalculator {
    private final Logger logger = LogManager.getLogger(HetPulldownCalculator.class);

    private final File refFile;
    private final IntervalList snpIntervals;

    private final int minMappingQuality;
    private final int minBaseQuality;
    private final ValidationStringency validationStringency;

    private static final int NUMBER_OF_LOG_UPDATES = 20;    //sets (approximate) number of status updates printed to log
    private static final double HET_ALLELE_FRACTION = 0.5;

    //set interval threshold for indexing for SamLocusIterator
    private static final int MAX_INTERVALS_FOR_INDEX = 25000;

    /**
     * Constructs a {@link HetPulldownCalculator} object for calculating {@link Pulldown} objects from files
     * containing a reference genome and an interval list of common SNP sites.  Reads and bases below the specified
     * mapping quality and base quality, respectively, are filtered out of the pileup.
     * @param refFile           file containing the reference
     * @param snpFile           file containing the interval list of common SNP sites
     * @param minMappingQuality minimum mapping quality required for reads to be included in pileup
     * @param minBaseQuality    minimum base quality required for bases to be included in pileup
     * @param validationStringency  validation stringency to use for reading BAM files
     */
    public HetPulldownCalculator(final File refFile, final File snpFile,
                                 final int minMappingQuality, final int minBaseQuality,
                                 final ValidationStringency validationStringency) {
        ParamUtils.isPositiveOrZero(minMappingQuality, "Minimum mapping quality must be nonnegative.");
        ParamUtils.isPositiveOrZero(minBaseQuality, "Minimum base quality must be nonnegative.");
        this.refFile = refFile;
        this.snpIntervals = IntervalList.fromFile(snpFile);
        this.minMappingQuality = minMappingQuality;
        this.minBaseQuality = minBaseQuality;
        this.validationStringency = validationStringency;
    }

    /**
     * Provides flags for running getHetPulldown based on sample type (normal or tumor).
     */
    private enum SampleType {
        NORMAL, TUMOR
    }

    /**
     * Returns map of base-pair counts at a given locus.  Only includes ACGT counts.
     * @param locus locus
     * @return      map of base-pair counts
     */
    public static Map<Character,Integer> getPileupBaseCounts(final SamLocusIterator.LocusInfo locus) {
        int aCount = 0, cCount = 0, gCount = 0, tCount = 0;

        for (final SamLocusIterator.RecordAndOffset rec : locus.getRecordAndPositions()) {
            switch (Character.toUpperCase(rec.getReadBase())) {
                case 'A':
                    ++aCount;
                    break;
                case 'C':
                    ++cCount;
                    break;
                case 'G':
                    ++gCount;
                    break;
                case 'T':
                    ++tCount;
                    break;
                default:
                    break;
            }
        }

        Map<Character,Integer> baseCounts = new HashMap<>();
        baseCounts.put('A', aCount);
        baseCounts.put('C', cCount);
        baseCounts.put('G', gCount);
        baseCounts.put('T', tCount);
        return baseCounts;
    }

    /**
     * Returns true if the distribution of major and other base-pair counts from a pileup at a locus is compatible with
     * allele fraction of 0.5.
     *
     * <p>
     *     Compatibility is defined by a p-value threshold.  That is, compute the two-sided p-value of observing
     *     a number of major read counts out of a total number of reads, assuming the given heterozygous
     *     allele fraction.  If the p-value is less than the given threshold, then reject the null hypothesis
     *     that the heterozygous allele fraction is 0.5 (i.e., SNP is likely to be homozygous) and return false,
     *     otherwise return true.
     * </p>
     * @param baseCounts        map of base-pair counts
     * @param totalBaseCount    total base-pair counts (excluding N, etc.)
     * @param pvalThreshold     p-value threshold for two-sided binomial test (should be in [0, 1], but no check is performed)
     * @return                  boolean compatibility with heterozygous allele fraction
     */
    @VisibleForTesting
    protected static boolean isPileupHetCompatible(final Map<Character, Integer> baseCounts, final int totalBaseCount,
                                                   final double pvalThreshold) {
        final int majorReadCount = Collections.max(baseCounts.values());

        if (majorReadCount == 0 || totalBaseCount - majorReadCount == 0) {
            return false;
        }

        final double pval = new BinomialTest().binomialTest(totalBaseCount, majorReadCount, HET_ALLELE_FRACTION,
                AlternativeHypothesis.TWO_SIDED);

        return pval >= pvalThreshold;
    }

    /**
     * Calls {@link HetPulldownCalculator#getHetPulldown} with flags set for a normal sample.
     */
    public Pulldown getNormal(final File normalBAMFile, final double pvalThreshold, final int minReadCount) {
        ParamUtils.inRange(pvalThreshold, 0., 1., "p-value threshold must be in [0, 1].");
        return getHetPulldown(normalBAMFile, this.snpIntervals, SampleType.NORMAL, pvalThreshold, minReadCount);
    }

    /**
     * Calls {@link HetPulldownCalculator#getHetPulldown} with flags set for a tumor sample.
     */
    public Pulldown getTumor(final File tumorBAMFile, final IntervalList normalHetIntervals, final int minReadCount) {
        return getHetPulldown(tumorBAMFile, normalHetIntervals, SampleType.TUMOR, -1, minReadCount);
    }

    /**
     * For a normal or tumor sample, returns a data structure giving (intervals, reference counts, alternate counts),
     * where intervals give positions of likely heterozygous SNP sites.
     *
     * <p>
     *     For a normal sample:
     *     <ul>
     *         The IntervalList snpIntervals gives common SNP sites in 1-based format.
     *     </ul>
     *     <ul>
     *         The p-value threshold must be specified for a two-sided binomial test,
     *         which is used to determine SNP sites from snpIntervals that are
     *         compatible with a heterozygous SNP, given the sample.  Only these sites are output.
     *     </ul>
     * </p>
     * <p>
     *     For a tumor sample:
     *     <ul>
     *         The IntervalList snpIntervals gives heterozygous SNP sites likely to be present in the normal sample.
     *         This should be from {@link HetPulldownCalculator#getNormal} in 1-based format.
     *         Only these sites are output.
     *     </ul>
     * </p>
     * @param bamFile           sorted BAM file for sample
     * @param snpIntervals      IntervalList of SNP sites
     * @param sampleType        flag indicating type of sample (SampleType.NORMAL or SampleType.TUMOR)
     *                          (determines whether to perform binomial test)
     * @param pvalThreshold     p-value threshold for two-sided binomial test, used for normal sample
     * @param minimumRawReads   minimum number of total reads that must be present at a het site
     * @return                  Pulldown of heterozygous SNP sites in 1-based format
     */
    private Pulldown getHetPulldown(final File bamFile, final IntervalList snpIntervals, final SampleType sampleType,
                                    final double pvalThreshold, final int minimumRawReads) {
        try (final SamReader bamReader = SamReaderFactory.makeDefault().validationStringency(validationStringency)
                .referenceSequence(refFile).open(bamFile);
             final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(this.refFile)) {
            if (bamReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new UserException.BadInput("BAM file " + bamFile.toString() + " must be coordinate sorted.");
            }

            final Pulldown hetPulldown = new Pulldown(bamReader.getFileHeader());

            final int totalNumberOfSNPs = snpIntervals.size();
            final SamLocusIterator locusIterator = new SamLocusIterator(bamReader, snpIntervals,
                    totalNumberOfSNPs < MAX_INTERVALS_FOR_INDEX);

            //set read and locus filters [note: read counts match IGV, but off by a few from pysam.mpileup]
            final List<SamRecordFilter> samFilters = Arrays.asList(new NotPrimaryAlignmentFilter(),
                    new DuplicateReadFilter());
            locusIterator.setSamFilters(samFilters);
            locusIterator.setEmitUncoveredLoci(false);
            locusIterator.setIncludeNonPfReads(false);
            locusIterator.setMappingQualityScoreCutoff(minMappingQuality);
            locusIterator.setQualityScoreCutoff(minBaseQuality);

            logger.info("Examining " + totalNumberOfSNPs + " sites...");
            final int iterationsPerStatus =
                    Math.max((int) Math.floor((float) totalNumberOfSNPs / NUMBER_OF_LOG_UPDATES), 1);
            int locusCount = 1;
            for (final SamLocusIterator.LocusInfo locus : locusIterator) {
                if (locusCount % iterationsPerStatus == 0) {
                    logger.info("Examined " + locusCount + " out of " + totalNumberOfSNPs + " sites.");
                }
                locusCount++;

                //include N, etc. reads here
                final int totalReadCount = locus.getRecordAndPositions().size();
                if (totalReadCount < minimumRawReads) {
                    continue;
                }

                final Map<Character, Integer> baseCounts = getPileupBaseCounts(locus);
                //only include total ACGT counts in binomial test (exclude N, etc.)
                final int totalBaseCount = baseCounts.values().stream().mapToInt(Number::intValue).sum();

                if (sampleType == SampleType.NORMAL &&
                        !isPileupHetCompatible(baseCounts, totalBaseCount, pvalThreshold)) {
                    continue;
                }

                final char refBase = (char) refWalker.get(locus.getSequenceIndex()).getBases()[locus.getPosition() - 1];
                final int refReadCount = baseCounts.get(refBase);
                final int altReadCount = totalBaseCount - refReadCount;

                hetPulldown.add(new SimpleInterval(locus.getSequenceName(), locus.getPosition(), locus.getPosition()),
                        refReadCount, altReadCount);
            }
            logger.info("Examined " + totalNumberOfSNPs + " sites.");
            return hetPulldown;
        } catch (final IOException | SAMFormatException e) {
            throw new UserException(e.getMessage());
        }
    }
}
