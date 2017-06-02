package org.broadinstitute.hellbender.tools.pon.allelic;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.CreateAllelicPanelOfNormals;
import org.broadinstitute.hellbender.tools.exome.CreatePanelOfNormals;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionInitializer;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionLikelihoods;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * <p>
 *     Represents the panel of normals used for reference-bias correction.
 * </p>
 * <p>
 *     The global distribution of reference bias across all sites is assumed to be Gamma(alpha, beta).  Hence,
 *     MLE values for the global hyperparameters alpha and beta are calculated using {@link AlleleFractionInitializer}
 *     and stored, along with the corresponding mean and variance of the distribution.
 * </p>
 * <p>
 *     Further, after these global hyperparameters are fit, posteriors for the local reference bias at each site
 *     in the panel are calculated.  These posteriors can also be represented as Gamma(alpha_i, beta_i), so we
 *     store them as a map with sites as keys and the local hyperparameters alpha_i and beta_i for each site i as values.
 *     These posteriors are then used as priors for the local reference bias downstream.
 * </p>
 * <p>
 *     Reading/writing to/from both TSV and HDF5 are implemented.
 * </p>
 *
 * See docs/CNVs/CNV-methods.pdf for details.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicPanelOfNormals {
    private static final Logger logger = LogManager.getLogger(AllelicPanelOfNormals.class);

    public static final AllelicPanelOfNormals EMPTY_PON = new AllelicPanelOfNormals();

    //comment strings for global hyperparameters in TSV output
    private static final String GLOBAL_ALPHA_COMMENT_STRING = "GLOBAL_ALPHA=";
    private static final String GLOBAL_BETA_COMMENT_STRING = "GLOBAL_BETA=";

    //local (per site) hyperparameter values
    private final Map<SimpleInterval, HyperparameterValues> siteToHyperparameterValuesMap = new HashMap<>();

    //global hyperparameter values
    private final HyperparameterValues globalHyperparameterValues;
    private final double globalMeanBias;
    private final double globalBiasVariance;

    //constructor for the empty PoN
    private AllelicPanelOfNormals() {
        globalHyperparameterValues = new HyperparameterValues(Double.NaN, Double.NaN);
        globalMeanBias = Double.NaN;
        globalBiasVariance = Double.NaN;
    }

    /**
     * Constructs an allelic panel of normals with the specified values.
     * @param globalHyperparameterValues    MLE global hyperparameter values
     * @param siteToHyperparameterValuesMap site-to-hyperparameter values map
     */
    public AllelicPanelOfNormals(final HyperparameterValues globalHyperparameterValues,
                                 final Map<SimpleInterval, HyperparameterValues> siteToHyperparameterValuesMap) {
        this.globalHyperparameterValues = globalHyperparameterValues;
        globalMeanBias = AlleleFractionLikelihoods.meanBias(globalHyperparameterValues.alpha, globalHyperparameterValues.beta);
        globalBiasVariance = AlleleFractionLikelihoods.biasVariance(globalHyperparameterValues.alpha, globalHyperparameterValues.beta);
        this.siteToHyperparameterValuesMap.putAll(siteToHyperparameterValuesMap);
    }

    /**
     * Constructs an allelic panel of normals from an {@link AllelicCountCollection} that contains
     * total alt and ref counts observed across all normals at each site.
     * @param counts    total alt and ref counts observed across all normals at each site
     */
    public AllelicPanelOfNormals(final AllelicCountCollection counts) {
        globalHyperparameterValues = calculateMLEGlobalHyperparameterValues(counts);
        globalMeanBias = AlleleFractionLikelihoods.meanBias(globalHyperparameterValues.alpha, globalHyperparameterValues.beta);
        globalBiasVariance = AlleleFractionLikelihoods.biasVariance(globalHyperparameterValues.alpha, globalHyperparameterValues.beta);
        initializeSiteToHyperparameterPairMap(counts);
    }

    /**
     * Gets the MLE global hyperparameter values alpha and beta.
     * @return  MLE global hyperparameter values alpha and beta
     */
    public HyperparameterValues getGlobalHyperparameterValues() {
        throwExceptionIfPoNIsEmpty();
        return globalHyperparameterValues;
    }

    /**
     * Gets the MLE global mean-bias hyperparameter.
     * @return  MLE global mean-bias hyperparameter
     */
    public double getGlobalMeanBias() {
        throwExceptionIfPoNIsEmpty();
        return globalMeanBias;
    }

    /**
     * Gets the MLE global bias-variance hyperparameter.
     * @return  MLE global bias-variance hyperparameter
     */
    public double getGlobalBiasVariance() {
        throwExceptionIfPoNIsEmpty();
        return globalBiasVariance;
    }

    /**
     * Returns a list of the site-hyperparameter map entries sorted by lexicographical contig order.
     * @return  list of the site-hyperparameter map entries sorted by lexicographical contig order
     */
    public List<Map.Entry<SimpleInterval, HyperparameterValues>> getSortedMapEntries() {
        return collectSortedMapEntries(siteToHyperparameterValuesMap);
    }

    /**
     * Gets the reference-bias alpha hyperparameter at a given SNP site if it is in the panel of normals
     * and the MLE global alpha hyperparameter if it is not.
     * @param site  SNP site
     * @return      reference-bias alpha hyperparameter if site is in panel of normals,
     *              MLE global alpha hyperparameter if it is not
     */
    public double getAlpha(final SimpleInterval site) {
        throwExceptionIfPoNIsEmpty();
        Utils.nonNull(site);
        return siteToHyperparameterValuesMap.getOrDefault(site, globalHyperparameterValues).alpha;
    }

    /**
     * Gets the reference-bias beta hyperparameter at a given SNP site if it is in the panel of normals
     * and the MLE global beta hyperparameter if it is not.
     * @param site  SNP site
     * @return      reference-bias beta hyperparameter if site is in panel of normals,
     *              MLE global beta hyperparameter if it is not
     */
    public double getBeta(final SimpleInterval site) {
        throwExceptionIfPoNIsEmpty();
        Utils.nonNull(site);
        return siteToHyperparameterValuesMap.getOrDefault(site, globalHyperparameterValues).beta;
    }

    /**
     * Reads an allelic panel of normals from an HDF5 or tab-separated file (file type is automatically detected).
     * Tab-separated files should have global hyperparameter values alpha and beta specified by comment lines
     * denoted by {@code GLOBAL_ALPHA_COMMENT_STRING} and {@code GLOBAL_BETA_COMMENT_STRING}:
     * <p>
     *     #GLOBAL_ALPHA=...<br>
     *     #GLOBAL_BETA=...
     * </p>
     * followed by lines specifying hyperparameter values at each site,
     * with column headers as in {@link AllelicPanelOfNormalsTableColumn}:
     * <p>
     *     CONTIG \t POSITION \t ALPHA \t BETA
     * </p>
     *
     * Note that we opt for a static read method as opposed to a constructor that takes a file parameter because
     * the former allows us to return the static {@code EMPTY_PON} if the allelic panel of normals is not present in
     * an HDF5 file.
     * @param inputFile    HDF5 file containing a coverage panel of normals created by {@link CreatePanelOfNormals}
     *                     and an allelic panel of normals created and set by {@link CreateAllelicPanelOfNormals}
     *                     ({@code EMPTY_PON} is returned if the latter was never set), or a
     *                     tab-separated file that contains global hyperparameters in comment lines and lines specifying hyperparameter values at each site
     */
    public static AllelicPanelOfNormals read(final File inputFile) {
        IOUtils.canReadFile(inputFile);

        if (isHDF5File(inputFile)) {
            //construct from HDF5 file
            try (final HDF5File hdf5File = new HDF5File(inputFile)) {
                final AllelicPanelOfNormals allelicPoN = HDF5AllelicPoNUtils.read(hdf5File);
                logger.info(String.format("Loaded allelic panel of normals from HDF5 file: %s.", inputFile));
                return allelicPoN;
            }
        } else {
            //construct from TSV file
            try (final AllelicPanelOfNormalsReader reader = new AllelicPanelOfNormalsReader(inputFile)) {
                //parse comment lines for global hyperparameter values
                final AllelicPanelOfNormals.HyperparameterValues globalHyperparameterValues = parseGlobalHyperparameterValues(inputFile);
                //parse data lines for local hyperparameter values
                final Map<SimpleInterval, HyperparameterValues> siteToHyperparameterValuesMap = new HashMap<>();
                reader.stream().forEach(s -> siteToHyperparameterValuesMap.put(s.getKey(), s.getValue()));
                final AllelicPanelOfNormals allelicPoN = new AllelicPanelOfNormals(globalHyperparameterValues, siteToHyperparameterValuesMap);
                logger.info(String.format("Loaded allelic panel of normals from TSV file: %s.", inputFile));
                return allelicPoN;
            } catch (final IOException | UncheckedIOException ex) {
                throw new UserException.CouldNotReadInputFile(inputFile, ex);
            }
        }
    }

    /**
     * Writes out the {@link AllelicPanelOfNormals} as tab-separated values to the specified file.
     * Output file will have global hyperparameter values alpha and beta specified by comment lines
     * denoted by {@code GLOBAL_ALPHA_COMMENT_STRING} and {@code GLOBAL_BETA_COMMENT_STRING}:
     * <p>
     *     #GLOBAL_ALPHA=...<br>
     *     #GLOBAL_BETA=...
     * </p>
     * followed by lines specifying hyperparameter values at each site,
     * with column headers as in {@link AllelicPanelOfNormalsTableColumn}:
     * <p>
     *     CONTIG \t POSITION \t ALPHA \t BETA
     * </p>
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     */
    public void write(final File outputFile) {
        Utils.nonNull(outputFile);
        final List<Map.Entry<SimpleInterval, HyperparameterValues>> sortedMapEntries = collectSortedMapEntries(siteToHyperparameterValuesMap);
        try (final AllelicPanelOfNormalsWriter writer = new AllelicPanelOfNormalsWriter(outputFile)) {
            writer.writeComment(GLOBAL_ALPHA_COMMENT_STRING + AllelicPanelOfNormalsWriter.formatDouble(globalHyperparameterValues.alpha));
            writer.writeComment(GLOBAL_BETA_COMMENT_STRING + AllelicPanelOfNormalsWriter.formatDouble(globalHyperparameterValues.beta));
            writer.writeAllRecords(sortedMapEntries);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    /**
     * Writes out the {@link AllelicPanelOfNormals} to the specified HDF5 file with the specified {@link HDF5File.OpenMode}.
     * See {@link HDF5AllelicPoNUtils#write}.
     * @param outputHDF5File    HDF5 file to write to (if it exists, it will be overwritten)
     * @param openMode          desired {@link HDF5File.OpenMode} (if {@code HDF5File.OpenMode.READ_ONLY}, an exception will be thrown)
     */
    public void write(final File outputHDF5File, final HDF5File.OpenMode openMode) {
        Utils.nonNull(outputHDF5File);
        try (final HDF5File hdf5File = new HDF5File(outputHDF5File, openMode)) {
            HDF5AllelicPoNUtils.write(hdf5File, this);
        }
    }

    public static class HyperparameterValues {
        private final double alpha;
        private final double beta;

        public HyperparameterValues(final double alpha, final double beta) {
            this.alpha = alpha;
            this.beta = beta;
        }

        /**
         * Initializes the hyperparameter values at a site given the total counts observed across all normals.
         * @param alpha MLE global hyperparameter value for alpha
         * @param beta  MLE global hyperparameter value for beta
         * @param a     total alt counts observed across all normals at site
         * @param r     total ref counts observed across all normals at site
         */
        private HyperparameterValues(final double alpha, final double beta, final int a, final int r) {
            final double f = 0.5;
            final int n = a + r;
            final double lambda0 = AlleleFractionLikelihoods.biasPosteriorMode(alpha, beta, f, a, r);
            final double kappa = AlleleFractionLikelihoods.biasPosteriorCurvature(alpha, f, r, n, lambda0);
            this.alpha = AlleleFractionLikelihoods.biasPosteriorEffectiveAlpha(lambda0, kappa);
            this.beta = AlleleFractionLikelihoods.biasPosteriorEffectiveBeta(lambda0, kappa);
        }

        public double getAlpha() {
            return alpha;
        }

        public double getBeta() {
            return beta;
        }
    }

    private static boolean isHDF5File(final File inputFile) {
        try (final HDF5File hdf5File = new HDF5File(inputFile)) {
            logger.info(String.format("Opening %s as an HDF5 file...", hdf5File.getFile()));
            return true;
        } catch (final HDF5LibException e) {
            return false;
        }
    }

    //parses global hyperparameters from comment lines in a TSV
    private static HyperparameterValues parseGlobalHyperparameterValues(final File inputFile) {
        try {
            final List<String> commentLines = FileUtils.readLines(inputFile).stream()
                    .filter(l -> l.contains(TableUtils.COMMENT_PREFIX)).collect(Collectors.toList());
            final List<String> globalAlphaCommentLines = commentLines.stream()
                    .filter(l -> l.contains(GLOBAL_ALPHA_COMMENT_STRING)).collect(Collectors.toList());
            final List<String> globalBetaCommentLines = commentLines.stream()
                    .filter(l -> l.contains(GLOBAL_BETA_COMMENT_STRING)).collect(Collectors.toList());
            if (globalAlphaCommentLines.size() != 1 || globalBetaCommentLines.size() != 1) {
                throw new UserException.BadInput("MLE global hyperparameter values for allelic panel of normals " +
                        "were not specified correctly in comment lines of file " + inputFile.getAbsolutePath() + ".");
            }
            final double globalAlpha = Double.valueOf(globalAlphaCommentLines.get(0)
                    .replace(TableUtils.COMMENT_PREFIX + GLOBAL_ALPHA_COMMENT_STRING, ""));
            final double globalBeta = Double.valueOf(globalBetaCommentLines.get(0)
                    .replace(TableUtils.COMMENT_PREFIX + GLOBAL_BETA_COMMENT_STRING, ""));
            ParamUtils.isPositive(globalAlpha, "MLE global hyperparameter alpha must be positive.");
            ParamUtils.isPositive(globalBeta, "MLE global hyperparameter alpha must be positive.");
            return new HyperparameterValues(globalAlpha, globalBeta);
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        } catch (final NumberFormatException e) {
            throw new UserException.BadInput("MLE global hyperparameter values for allelic panel of normals " +
                    "were not specified correctly in comment lines of file " + inputFile.getAbsolutePath() + ".");
        }
    }

    /**
     * Find MLE hyperparameter values from counts in panel of normals.  See analogous code in {@link AlleleFractionInitializer}.
     */
    private HyperparameterValues calculateMLEGlobalHyperparameterValues(final AllelicCountCollection counts) {
        double meanBias = AlleleFractionInitializer.INITIAL_MEAN_BIAS;
        double biasVariance = AlleleFractionInitializer.INITIAL_BIAS_VARIANCE;
        double previousIterationLogLikelihood;
        double nextIterationLogLikelihood = Double.NEGATIVE_INFINITY;
        logger.info(String.format("Initializing MLE hyperparameter values for allelic panel of normals.  Iterating until log likelihood converges to within %.3f.",
                AlleleFractionInitializer.LOG_LIKELIHOOD_CONVERGENCE_THRESHOLD));
        int iteration = 1;
        do {
            previousIterationLogLikelihood = nextIterationLogLikelihood;
            meanBias = estimateMeanBias(meanBias, biasVariance, counts);
            biasVariance = estimateBiasVariance(meanBias, biasVariance, counts);
            nextIterationLogLikelihood = AlleleFractionLikelihoods.logLikelihoodForAllelicPanelOfNormals(meanBias, biasVariance, counts);
            logger.info(String.format("Iteration %d, model log likelihood = %.3f.", iteration, nextIterationLogLikelihood));
            iteration++;
        } while (iteration < AlleleFractionInitializer.MAX_ITERATIONS &&
                nextIterationLogLikelihood - previousIterationLogLikelihood > AlleleFractionInitializer.LOG_LIKELIHOOD_CONVERGENCE_THRESHOLD);

        final double alpha = AlleleFractionLikelihoods.alpha(meanBias, biasVariance);
        final double beta = AlleleFractionLikelihoods.beta(meanBias, biasVariance);
        logger.info("MLE hyperparameter values for allelic panel of normals found:");
        logger.info("alpha = " + alpha);
        logger.info("beta = " + beta);
        return new HyperparameterValues(alpha, beta);
    }

    private static double estimateMeanBias(final double meanBias, final double biasVariance, final AllelicCountCollection counts) {
        final Function<Double, Double> objective = proposedMeanBias ->
                AlleleFractionLikelihoods.logLikelihoodForAllelicPanelOfNormals(proposedMeanBias, biasVariance, counts);
        return OptimizationUtils.argmax(objective, 0.0, AlleleFractionInitializer.MAX_REASONABLE_MEAN_BIAS, meanBias);
    }

    private static double estimateBiasVariance(final double meanBias, final double biasVariance, final AllelicCountCollection counts) {
        final Function<Double, Double> objective = proposedBiasVariance ->
                AlleleFractionLikelihoods.logLikelihoodForAllelicPanelOfNormals(meanBias, proposedBiasVariance, counts);
        return OptimizationUtils.argmax(objective, 0.0, AlleleFractionInitializer.MAX_REASONABLE_BIAS_VARIANCE, biasVariance);
    }

    //transforms ref/alt counts at each site to hyperparameters, see docs/CNVs/CNV-methods.pdf for details
    private void initializeSiteToHyperparameterPairMap(final AllelicCountCollection counts) {
        logger.info("Initializing allelic panel of normals...");
        for (final AllelicCount count : counts.getCounts()) {
            final SimpleInterval site = count.getInterval();
            final HyperparameterValues hyperparameterValues = new HyperparameterValues(
                    globalHyperparameterValues.alpha, globalHyperparameterValues.beta,
                    count.getAltReadCount(), count.getRefReadCount());
            if (siteToHyperparameterValuesMap.containsKey(site)) {
                throw new UserException.BadInput("Input AllelicCountCollection for allelic panel of normals contains duplicate sites.");
            } else {
                siteToHyperparameterValuesMap.put(site, hyperparameterValues);
            }
        }
        logger.info("Allelic panel of normals initialized.");
    }

    private void throwExceptionIfPoNIsEmpty() {
        if (equals(EMPTY_PON)) {
            throw new UnsupportedOperationException("Cannot get MLE hyperparameters for empty panel of normals.");
        }
    }

    //returns a list of the sites in siteToHyperarameterPairMap that are sorted by SimpleInterval
    //(contigs are sorted by lexicographical order)
    private static List<Map.Entry<SimpleInterval, HyperparameterValues>> collectSortedMapEntries(
            final Map<SimpleInterval, HyperparameterValues> siteToHyperparameterPairMap) {
        final List<SimpleInterval> sortedMapKeys = new ArrayList<>(siteToHyperparameterPairMap.keySet());
        Collections.sort(sortedMapKeys, IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
        return sortedMapKeys.stream().map(si -> new AbstractMap.SimpleEntry<>(si, siteToHyperparameterPairMap.get(si)))
                .collect(Collectors.toList());
    }
}
