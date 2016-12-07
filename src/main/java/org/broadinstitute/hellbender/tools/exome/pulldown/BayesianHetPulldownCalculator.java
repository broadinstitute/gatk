package org.broadinstitute.hellbender.tools.exome.pulldown;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import htsjdk.samtools.*;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.NotPrimaryAlignmentFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountTableColumn;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * A Bayesian heterozygous SNP pulldown calculator. Base qualities are taken into account
 * to increase precision (see CNV-methods.pdf for details).
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */

public final class BayesianHetPulldownCalculator {

    private final Logger logger = LogManager.getLogger(BayesianHetPulldownCalculator.class);

    /**
     * A simple class to handle base read and mapping error probabilitites
     */
    @VisibleForTesting
    static final class BaseQuality {
        private final double readErrorProbability, mappingErrorProbability;

        BaseQuality(final double readErrorProbability, final double mappingErrorProbability) {
            this.readErrorProbability = readErrorProbability;
            this.mappingErrorProbability = mappingErrorProbability;
        }

        double getReadErrorProbability() { return readErrorProbability; }
        double getMappingErrorProbability() { return mappingErrorProbability; }
    }

    private final HeterozygousPileupPriorModel hetPrior;

    private static final Nucleotide[] PROPER_BASES = {Nucleotide.A, Nucleotide.C, Nucleotide.T, Nucleotide.G};

    private final File refFile;
    private final IntervalList snpIntervals;

    private final int readDepthThreshold;
    private final int minMappingQuality;
    private final int minBaseQuality;
    private final ValidationStringency validationStringency;

    /* experimental */
    private final double errorProbabilityAdjustmentFactor;

    /* default priors */
    private static final double DEFAULT_PRIOR_REF_HOM = 0.5; /* a homozygous site being the ref allele */
    private static final double DEFAULT_PRIOR_HET = 0.5; /* a site being heterozygous */

    /* approximate number of status updates printed to log */
    private static final int NUMBER_OF_SITES_PER_LOGGED_STATUS_UPDATE = 10000;

    /**
     * Constructor of {@link BayesianHetPulldownCalculator} object
     *
     * @param refFile the reference genome file
     * @param snpIntervals {@link IntervalList} of common SNPs
     * @param minMappingQuality minimum phred mapping quality
     * @param minBaseQuality minimum phred base quality
     * @param readDepthThreshold minimum read depth
     * @param validationStringency validation stringency
     * @param errorProbabilityAdjustmentFactor (experimental) multiplicative factor for read and mapping error
     *                                         probabilities
     * @param hetPrior the prior model for heterzygous pileups
     */
    public BayesianHetPulldownCalculator(final File refFile, final IntervalList snpIntervals,
                                         final int minMappingQuality, final int minBaseQuality,
                                         final int readDepthThreshold, final ValidationStringency validationStringency,
                                         final double errorProbabilityAdjustmentFactor,
                                         final HeterozygousPileupPriorModel hetPrior) {
        ParamUtils.isPositiveOrZero(minMappingQuality, "Minimum mapping quality must be nonnegative.");
        ParamUtils.isPositiveOrZero(minBaseQuality, "Minimum base quality must be nonnegative.");

        this.refFile = Utils.nonNull(refFile);
        this.snpIntervals = Utils.nonNull(snpIntervals);
        this.minMappingQuality = ParamUtils.isPositive(minMappingQuality, "Minimum mapping quality must be a positive integer");
        this.minBaseQuality = ParamUtils.isPositive(minBaseQuality, "Minimum base quality must be a positive integer");
        this.readDepthThreshold = ParamUtils.isPositive(readDepthThreshold, "Read depth threshold must be a positive integer");
        this.validationStringency = Utils.nonNull(validationStringency);
        this.errorProbabilityAdjustmentFactor = ParamUtils.isPositive(errorProbabilityAdjustmentFactor,
                "Error adjustment factor must be positive.");
        this.hetPrior = Utils.nonNull(hetPrior);
    }

    /**
     * Calculate the log likelihood of a read pileup is composed of a single nucelotide {@code expectedNucleotide}
     * @param baseQualities map of bases to list of their calling error probabilities
     * @param expectedNucleotide the expected nucleotide
     * @return the log likelihood
     */
    private double getSingletonLogLikelihood(final Map<Nucleotide, List<BaseQuality>> baseQualities,
                                             final Nucleotide expectedNucleotide) {
        return baseQualities.get(expectedNucleotide).stream()
                .mapToDouble(bq -> FastMath.log(1 - bq.getReadErrorProbability() -
                        3 * bq.getMappingErrorProbability() / 4)).sum() +
                baseQualities.keySet().stream()
                        .filter(nucl -> nucl != expectedNucleotide)
                        .map(baseQualities::get).mapToDouble(bqlist -> bqlist.stream()
                        .mapToDouble(bq -> FastMath.log(bq.getReadErrorProbability() / 3 +
                                bq.getMappingErrorProbability() / 4)).sum()).sum();
    }

    /**
     * Calculate the log likelihood of a SNP site being homozygous for a given read pileup
     * (see CNV-method.pdf for details)
     * @param baseQualities map of bases to list of their calling error probabilities at the SNP site
     * @param refNucleotide the ref allele nucleotide
     * @param altNucleotide the alt allele nucleotide
     * @param homRefPrior the prior probability of the ref allele given that the site is homozygous
     * @return the log likelihood
     */
    @VisibleForTesting
    double getHomLogLikelihood(final Map<Nucleotide, List<BaseQuality>> baseQualities,
                               final Nucleotide refNucleotide, final Nucleotide altNucleotide,
                               final double homRefPrior) {
        /* return the sum of |hom,ref) and |hom,alt) likelihoods */
        return GATKProtectedMathUtils.logSumExp(
                FastMath.log(homRefPrior) + getSingletonLogLikelihood(baseQualities, refNucleotide),
                FastMath.log(1 - homRefPrior) + getSingletonLogLikelihood(baseQualities, altNucleotide));
    }

    /**
     * Calculate the log likelihood of a SNP site being heterozygous for a given read pileup
     * (see CNV-method.pdf for details)
     *
     * @param baseQualities map of bases to list of their calling error probabilities at the SNP site
     * @param refNucleotide the ref allele nucleotide
     * @param altNucleotide the alt allele nucleotide
     * @return the log likelihood
     */
    @VisibleForTesting
    double getHetLogLikelihood(final Map<Nucleotide, List<BaseQuality>> baseQualities,
                               final Nucleotide refNucleotide, final Nucleotide altNucleotide) {
        /* non-Ref-Alt entries */
        final double errorLogLikelihood = Sets.difference(baseQualities.keySet(),
                new HashSet<>(Arrays.asList(refNucleotide, altNucleotide))).stream()
                .map(baseQualities::get).mapToDouble(bqlist -> bqlist.stream()
                        .mapToDouble(bq -> FastMath.log(bq.getReadErrorProbability() / 3 +
                                bq.getMappingErrorProbability() / 4))
                        .sum())
                .sum();

        final List<ImmutablePair<Double, Double>> coeffs = new ArrayList<>();

        /* Ref entries in the pile up */
        baseQualities.get(refNucleotide).stream().forEach(bq ->coeffs.add(new ImmutablePair<>(
                        bq.getReadErrorProbability() / 3 + bq.getMappingErrorProbability() / 4,
                        1 - 4 * bq.getReadErrorProbability() / 3 - bq.getMappingErrorProbability() / 4)));

        /* Alt entries in the pile up */
        baseQualities.get(altNucleotide).stream().forEach(bq -> coeffs.add(new ImmutablePair<>(
                        1 - bq.getReadErrorProbability() - 3 * bq.getMappingErrorProbability() / 4,
                        -1 + 4 * bq.getReadErrorProbability() / 3 + bq.getMappingErrorProbability())));

        return hetPrior.getHetLogLikelihood(coeffs) + errorLogLikelihood;
    }

    /**
     * Checks of a given base is in PROPER_BASES
     * @param base a nucleotide
     * @return boolean
     */
    private boolean isProperBase(final Nucleotide base) {
        for (final Nucleotide properBase : PROPER_BASES) {
            if (base.equals(properBase)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns map of base-pair to error probabilities at a given locus. All reads are considered (not just ACTG)
     * @param locus locus
     * @return map of base-pair to error probabilities
     */
    private Map<Nucleotide, List<BaseQuality>> getPileupBaseQualities(final SamLocusIterator.LocusInfo locus) {
        final Map<Nucleotide, List<BaseQuality>> baseQualities = locus.getRecordAndPositions().stream()
                .map(rp -> new ImmutablePair<>(
                        Nucleotide.valueOf(rp.getReadBase()),
                        new BaseQuality(
                                errorProbabilityAdjustmentFactor * QualityUtils.qualToErrorProb(rp.getBaseQuality()),
                                errorProbabilityAdjustmentFactor * QualityUtils.qualToErrorProb(rp.getRecord().getMappingQuality())
                        )
                ))
                .filter(rp -> isProperBase(rp.getLeft()))
                .collect(Collectors.groupingBy(
                        ImmutablePair::getLeft,
                        Collectors.mapping(ImmutablePair::getRight, Collectors.toList()))
                );

        /* make sure that the main bases {A, C, T, G} are included in the map */
        for (final Nucleotide base : PROPER_BASES) {
            if (!baseQualities.containsKey(base)) {
                baseQualities.put(base, new ArrayList<>());
            }
        }

        return baseQualities;
    }

    /**
     * Returns base-pair counts at a given locus.
     * @param locus locus
     * @return      base-pair counts
     */
    private static Nucleotide.Counter getPileupBaseCounts(final SamLocusIterator.LocusInfo locus) {
        final Nucleotide.Counter result = new Nucleotide.Counter();
        for (final SamLocusIterator.RecordAndOffset rec : locus.getRecordAndPositions()) {
            result.add(rec.getReadBase());
        }
        return result;
    }

    /**
     * Infer the alt allele base from the pileup: we pick the base with the highest frequency that is not Ref as Alt.
     * This is the maximum likelihood estimation. <br>
     *
     * <b>Remark</b>: While this is definitely not the perfect way to do it and it is desirable to have both Ref/Alt
     * allele information, it has no serious pitfalls either:<BR><BR>
     * <ul>
     *     <li>
     *         If the reads are balanced between between Ref and Alt, it correctly infers the Alt base.
     *     </li>
     *     <li>
     *         If the input is all Ref, then Alt is chosen randomly, but it's fine since there are no Alt counts in
     *         the pileup anyway and the likelihoods are not affected. If the input is all Alt, it still works properly.
     *     </li>
     *     <li>
     *         If the input in all Ref + a few read/mapping errors, then Alt will be chosen as the most frequent
     *         erroneous base; again, this is fine because Hom likelihood is far greater than Het likelihood in this
     *         case anyway. Likewise, if the input is all Alt + a few errors, it still works properly.
     *     </li>
     * </ul>
     *
     * @param baseQualities map from bases to error probabilities
     * @param refBase the ref allele base
     * @return the likely alt allele
     */
    @VisibleForTesting
    static Nucleotide inferAltFromPileup(final Map<Nucleotide, List<BaseQuality>> baseQualities,
                                                final Nucleotide refBase) {
        /* sort the bases in the descending order by their frequency */
        final Nucleotide[] bases = PROPER_BASES.clone();
        Arrays.sort(bases, (L, R) -> Integer.compare(baseQualities.get(R).size(), baseQualities.get(L).size()));
        /* pick the base with highest frequency, skip over ref */
        for (Nucleotide base : bases) {
            if (base != refBase) {
                return base;
            }
        }

        /* we shouldn't be here unless the baseErrorProbabilities is malformed */
        throw new GATKException.ShouldNeverReachHereException("The Alt base can not be inferred from the " +
                "pileup. The size of the pileup is: " + baseQualities.size());
    }

    /**
     * Returns a {@link SamLocusIterator} object for a given {@link SamReader} and {@link IntervalList} with filters
     * on minimum base quality and minimum mapping quality
     *
     * @param samReader a SamReader object
     * @return a SamLocusIterator object
     */
    private SamLocusIterator getSamLocusIteratorWithDefaultFilters(final SamReader samReader) {
        final SamLocusIterator locusIterator = new SamLocusIterator(samReader, snpIntervals, false);

        /* set read and locus filters */
        final List<SamRecordFilter> samFilters = Arrays.asList(new NotPrimaryAlignmentFilter(),
                new DuplicateReadFilter());
        locusIterator.setSamFilters(samFilters);
        locusIterator.setEmitUncoveredLoci(false);
        locusIterator.setIncludeNonPfReads(false);
        locusIterator.setMappingQualityScoreCutoff(minMappingQuality);
        locusIterator.setQualityScoreCutoff(minBaseQuality);

        return locusIterator;
    }

    /**
     * For a given normal or tumor BAM file, walks through the list of common SNPs,
     * {@link BayesianHetPulldownCalculator#snpIntervals}), detects heterozygous sites, and returns
     * a {@link Pulldown} containing detailed information on the called heterozygous SNP sites.
     *
     * The {@code hetCallingStrigency} parameters sets the threshold posterior for calling a Het SNP site:
     *
     *      hetPosteriorThreshold = 1 - 10^{-hetCallingStringency}
     *      hetThresholdLogOdds = log(hetPosteriorThreshold/(1-hetPosteriorThreshold))
     *                          = log(10^{hetCallingStringency} - 1)
     *
     * (see CNV-methods.pdf for details)
     *
     * @param bamFile sorted BAM file for sample
     * @param hetCallingStringency strigency for calling a Het site
     * @return Pulldown of heterozygous SNP sites in 1-based format
     */
    public Pulldown getHetPulldown(final File bamFile, final double hetCallingStringency) {
        /* log odds from stringency */
        final double hetThresholdLogOdds = FastMath.log(FastMath.pow(10, hetCallingStringency) - 1);

        try (final SamReader bamReader = SamReaderFactory.makeDefault().validationStringency(validationStringency)
                .referenceSequence(refFile).open(bamFile);
             final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(refFile)) {
            if (bamReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new UserException.BadInput("BAM file " + bamFile.toString() + " must be coordinate sorted.");
            }

            final Pulldown hetPulldown = new Pulldown(bamReader.getFileHeader());
            final SamLocusIterator locusIterator = getSamLocusIteratorWithDefaultFilters(bamReader);

            final int totalNumberOfSNPs = snpIntervals.size();
            logger.info("Examining " + totalNumberOfSNPs + " sites in total...");
            int locusCount = 0;
            for (final SamLocusIterator.LocusInfo locus : locusIterator) {
                if (locusCount % NUMBER_OF_SITES_PER_LOGGED_STATUS_UPDATE == 0) {
                    logger.info("Examined " + locusCount + " covered sites.");
                }
                locusCount++;

                final int totalReadCount = locus.getRecordAndPositions().size();
                if (totalReadCount <= readDepthThreshold) {
                    continue;
                }
                final Nucleotide refBase = Nucleotide.valueOf(refWalker.get(locus.getSequenceIndex())
                        .getBases()[locus.getPosition() - 1]);
                if (!isProperBase(refBase)) {
                    logger.warn(String.format("The reference position at %d has an unknown base call (value: %s). Even though" +
                            " this position is indicated to be a possible heterozygous SNP in the provided SNP interval list," +
                            " no inference can be made. Continuing ...", locus.getPosition(), refBase.toString()));
                    continue;
                }

                final Map<Nucleotide, List<BaseQuality>> baseQualities = getPileupBaseQualities(locus);
                final Nucleotide altBase = inferAltFromPileup(baseQualities, refBase);

                /* calculate Het log odds */
                final double hetLogLikelihood = getHetLogLikelihood(baseQualities, refBase, altBase);
                final double homLogLikelihood = getHomLogLikelihood(baseQualities, refBase, altBase,
                        DEFAULT_PRIOR_REF_HOM);
                final double hetLogOdds = (hetLogLikelihood + FastMath.log(DEFAULT_PRIOR_HET)) -
                        (homLogLikelihood + FastMath.log(1 - DEFAULT_PRIOR_HET));

                if (hetLogOdds > hetThresholdLogOdds) {
                    hetPulldown.add(new AllelicCount(
                            new SimpleInterval(locus.getSequenceName(), locus.getPosition(), locus.getPosition()),
                            baseQualities.get(refBase).size(), baseQualities.get(altBase).size(),
                            refBase, altBase, totalReadCount, hetLogOdds));
                }
            }

            logger.info(locusCount + " covered sites out of " + totalNumberOfSNPs + " total sites were examined.");

            return hetPulldown;

        } catch (final IOException | SAMFormatException e) {
            throw new UserException(e.getMessage());
        }
    }

    /**
     * Calculates the Het pulldown from a tumor file, given the tumor BAM file and the pulldown from a matched
     * normal BAM file.
     *
     * Note: this method does not perform any statistical inference. The Het SNP sites are directly carried over
     * from the matched normal pulldown. Here, we only collect statistics (ref count, alt count, read depth) and
     * save to a pulldown. The verbosity level of the pulldown is INTERMEDIATE (see {@link AllelicCountTableColumn}).
     *
     * @param tumorBamFile the tumor BAM file
     * @param normalHetPulldown the matched normal Het pulldown
     * @return tumor Het pulldown
     */
    public Pulldown getTumorHetPulldownFromNormalPulldown(final File tumorBamFile, final Pulldown normalHetPulldown) {
        try (final SamReader bamReader = SamReaderFactory.makeDefault().validationStringency(validationStringency)
                .referenceSequence(refFile).open(tumorBamFile)) {
            if (bamReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new UserException.BadInput("BAM file " + tumorBamFile.toString() + " must be coordinate sorted.");
            }

            final Pulldown tumorHetPulldown = new Pulldown(bamReader.getFileHeader());
            final SamLocusIterator locusIterator = getSamLocusIteratorWithDefaultFilters(bamReader);

            /* get a map of SimpleIntervals in the pulldown to their index */
            final Map<SimpleInterval, Integer> normalPulldownIndexMap = normalHetPulldown.getSimpleIntervalToIndexMap();

            final int totalNumberOfSNPs = snpIntervals.size();
            logger.info("Examining " + totalNumberOfSNPs + " sites in total...");
            int locusCount = 0;
            for (final SamLocusIterator.LocusInfo locus : locusIterator) {
                if (locusCount % NUMBER_OF_SITES_PER_LOGGED_STATUS_UPDATE == 0) {
                    logger.info("Examined " + locusCount + " covered sites.");
                }
                locusCount++;

                final int totalReadCount = locus.getRecordAndPositions().size();
                if (totalReadCount <= readDepthThreshold) {
                    continue;
                }

                /* find the AllelicCount from the normal pulldown */
                int indexInNormalPulldown;
                try {
                    indexInNormalPulldown = normalPulldownIndexMap.get(
                            new SimpleInterval(locus.getSequenceName(), locus.getPosition(), locus.getPosition()));
                } catch (NullPointerException e) {
                    throw new GATKException.ShouldNeverReachHereException("Can not find the required AllelicCount " +
                            "object in the normal pulldown. Stopping.");
                }

                /* just count the alt and ref nucleotide and add to the tumor pulldown */
                final Nucleotide.Counter baseCounts = getPileupBaseCounts(locus);
                tumorHetPulldown.add(new AllelicCount(
                        new SimpleInterval(locus.getSequenceName(), locus.getPosition(), locus.getPosition()),
                        (int) baseCounts.get(normalHetPulldown.getCounts().get(indexInNormalPulldown).getRefNucleotide()),
                        (int) baseCounts.get(normalHetPulldown.getCounts().get(indexInNormalPulldown).getAltNucleotide()),
                        normalHetPulldown.getCounts().get(indexInNormalPulldown).getRefNucleotide(),
                        normalHetPulldown.getCounts().get(indexInNormalPulldown).getAltNucleotide(),
                        totalReadCount)
                );
            }

            logger.info(locusCount + " covered sites out of " + totalNumberOfSNPs + " total sites were examined.");

            return tumorHetPulldown;

        } catch (final IOException | SAMFormatException e) {
            throw new UserException(e.getMessage());
        }
    }
}
