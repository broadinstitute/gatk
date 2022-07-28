package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.primitives.Doubles;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.broadinstitute.hellbender.tools.walkers.vqsr.VariantRecalibrator;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.BGMMVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.PythonSklearnVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsModelBackend;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsScorer;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Scores variant calls in a VCF file based on site-level annotations using a previously trained model.
 *
 * <p>
 *     This tool is intended to be used as the last step in a variant-filtering workflow that supersedes the
 *     {@link VariantRecalibrator} workflow. Using a previously trained model produced by {@link TrainVariantAnnotationsModel},
 *     this tool assigns a score to each call (with a lower score indicating that a call is more likely to be an artifact).
 *     Each score can also be converted to a corresponding sensitivity to a calibration set, if the latter is available.
 *     Each VCF record can also be annotated with additional resource labels and/or hard filtered based on its
 *     calibration-set sensitivity, if desired.
 * </p>
 *
 * <p>
 *     Note that annotations and metadata are collected in memory during traversal until they are written to HDF5 files
 *     upon completion of the traversal. Memory requirements thus roughly scale linearly with both the number of sites
 *     scored and the number of annotations. For large callsets, this tool may be run in parallel over separate
 *     genomic shards using the {@value StandardArgumentDefinitions#INTERVALS_LONG_NAME} argument as usual.
 * </p>
 *
 * <p>
 *     Scores and annotations are also output to HDF5 files, which may be viewed using
 *     <a href="https://support.hdfgroup.org/products/java/hdfview/">hdfview</a> or loaded in Python using
 *     <a href="http://www.pytables.org/">PyTables</a> or <a href="http://www.h5py.org/">h5py</a>.
 * </p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Input VCF file. Site-level annotations will be extracted from the contained variants (or alleles,
 *         if the {@value USE_ALLELE_SPECIFIC_ANNOTATIONS_LONG_NAME} argument is specified).
 *     </li>
 *     <li>
 *         Annotations to use for scoring. These should be identical to those used in the {@link ExtractVariantAnnotations}
 *         step to create the training set.
 *     </li>
 *     <li>
 *         Variant types (i.e., SNP and/or INDEL) to score. Logic for determining variant type was retained from
 *         {@link VariantRecalibrator}; see {@link VariantType}. To use different models for SNPs and INDELs
 *         (e.g., if it is desired to use different sets of annotations for each variant type), one can first run
 *         this tool to score SNPs and then again on the resulting output to score INDELs.
 *     </li>
 *     <li>
 *         Model prefix. This should denote the path of model files produced by {@link TrainVariantAnnotationsModel}.
 *     </li>
 *     <li>
 *         (Optional) Model backend. This should be identical to that specified in {@link TrainVariantAnnotationsModel}.
 *         The default Python IsolationForest implementation requires either the GATK Python environment
 *         or that certain Python packages (argparse, h5py, numpy, sklearn, and dill) are otherwise available.
 *         A custom backend can also be specified in conjunction with the {@value PYTHON_SCRIPT_LONG_NAME} argument.
 *     </li>
 *     <li>
 *         (Optional) Resource VCF file(s). See the corresponding documentation in {@link ExtractVariantAnnotations}.
 *         In typical usage, the same resource VCFs and tags provided to that tool should also be provided here.
 *         In addition, the sites-only VCF that is produced by that tool can also be provided here and used to
 *         mark those labeled sites that were extracted, which can be useful if these are a subset of the resource sites.
 *     </li>
 *     <li>
 *         (Optional) Calibration-set sensitivity thresholds for SNPs and INDELs. If the corresponding SNP or INDEL
 *         calibration-set scores are available in the provided model files, sites that have a calibration-set
 *         sensitivity falling above the corresponding threshold (i.e., a score falling below the corresponding
 *         score threshold) will have a filter applied.
 *     </li>
 *     <li>
 *         Output prefix.
 *         This is used as the basename for output files.
 *     </li>
 * </ul>
 *
 * <h3>Outputs</h3>
 *
 * <ul>
 *     <li>
 *         Scored VCF file and index. The VCF will not be gzipped if the {@value DO_NOT_GZIP_VCF_OUTPUT_LONG_NAME}
 *         argument is set to true. The INFO field in each VCF record will be annotated with:
 *
 *         <p>
 *             1) a score (with a key as given by the {@value SCORE_KEY_LONG_NAME} argument,
 *             which has a default value of {@value DEFAULT_SCORE_KEY}),
 *         </p>
 *         <p>
 *             2) if resources are provided, flags corresponding to the labels (e.g.,
 *             {@value LabeledVariantAnnotationsData#TRAINING_LABEL}, {@value LabeledVariantAnnotationsData#CALIBRATION_LABEL}, etc.)
 *             of resources containing the record,
 *         </p>
 *         <p>
 *             3) if the {@value SNP_KEY_LONG_NAME} argument (which has a default value of {@value DEFAULT_SNP_KEY})
 *             is non-null, a flag corresponding to whether a site is treated as a SNP,
 *         </p>
 *         <p>
 *             4) if {@value SNP_CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME} and/or
 *             {@value INDEL_CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME} are provided, a filter (with name given by
 *             the {@value LOW_SCORE_FILTER_NAME_LONG_NAME} argument, which has a default value of
 *             {@value DEFAULT_LOW_SCORE_FILTER_NAME}) will be applied if a record has a calibration-set sensitivity
 *             falling above the appropriate threshold (i.e., if it has a score falling below the corresponding
 *             score threshold).
 *         </p>
 *         <p>
 *             If {@value USE_ALLELE_SPECIFIC_ANNOTATIONS_LONG_NAME} is true, the score, SNP flag, calibration sensitivity,
 *             and filter appropriate for the highest scoring allele are used; however, the resource labels for all alleles
 *             are applied.
 *         </p>
 *
 *     </li>
 *     <li>
 *         (Optional) Annotations HDF5 file (.annot.hdf5). Annotation data and metadata for all scored sites
 *         (labeled and unlabeled) are stored in the HDF5 directory structure given in the documentation for the
 *         {@link ExtractVariantAnnotations} tool. This file will only be produced if the number of scored sites
 *         is nonzero.
 *         </p>
 *
 *     </li>
 *     <li>
 *         (Optional) Scores HDF5 file (.scores.hdf5). Scores for all scored sites are stored in the
 *         HDF5 path {@value VariantAnnotationsScorer#SCORES_PATH}. Scores are given in the same order as records
 *         in both the VCF and the annotations HDF5 file. This file will only be produced if the number of scored sites
 *         is nonzero.
 *         </p>
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <p>
 *     Score sites using a model (produced by {@link TrainVariantAnnotationsModel} using the default
 *     {@link VariantAnnotationsModelBackend#PYTHON_IFOREST} model backend and contained in the directory
 *     {@code model_dir}), producing the outputs 1) {@code output.vcf.gz}, 2) {@code output.vcf.gz.tbi},
 *     3) {@code output.annot.hdf5}, and 4) {@code output.scores.hdf5}. Note that {@code extract.vcf.gz} is
 *     produced by {@link ExtractVariantAnnotations}. Records will be filtered according to the values provided to the
 *     {@value SNP_CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME} and {@value INDEL_CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME}
 *     arguments; the values below are only meant to be illustrative and should be set as appropriate for a given analysis.
 *
 * <pre>
 *     gatk ScoreVariantAnnotations \
 *          -V input.vcf \
 *          -A annotation_1 \
 *          ...
 *          -A annotation_N \
 *          --model-prefix model_dir \
 *          --mode SNP \
 *          --resource snp-training,training=true snp-training.vcf \
 *          --resource snp-calibration,calibration=true snp-calibration.vcf \
 *          --mode INDEL \
 *          --resource indel-training,training=true indel-training.vcf \
 *          --resource indel-calibration,calibration=true indel-calibration.vcf \
 *          --resource extracted,extracted=true extract.vcf.gz \
 *          --snp-calibration-sensitivity-threshold 0.99 \
 *          --indel-calibration-sensitivity-threshold 0.99 \
 *          -O output
 * </pre>
 *
 * <p>
 *     One may chain together two runs of this tool to score SNPs and INDELs using different models
 *     (note that SNP and INDEL models have "snp" and "indel" tags in their respective filenames, so these
 *     models can still be contained in the same {@code model_dir} directory).
 *     This may have implications for mixed SNP/INDEL sites, especially if filters are applied; see also the
 *     {@value IGNORE_ALL_FILTERS_LONG_NAME} and {@value IGNORE_FILTER_LONG_NAME} arguments.
 *
 * <pre>
 *     gatk ScoreVariantAnnotations \
 *          -V input.vcf \
 *          -A snp_annotation_1 \
 *          ...
 *          -A snp_annotation_N \
 *          --model-prefix model_dir \
 *          --mode SNP \
 *          --resource snp-training,training=true snp-training.vcf \
 *          --resource snp-calibration,calibration=true snp-calibration.vcf \
 *          --resource extracted,extracted=true snp-extract.vcf.gz \
 *          --snp-calibration-sensitivity-threshold 0.99 \
 *          -O intermediate-output
 *
 *     gatk ScoreVariantAnnotations \
 *          -V intermediate-output.vcf \
 *          -A indel_annotation_1 \
 *          ...
 *          -A indel_annotation_M \
 *          --model-prefix model_dir \
 *          --mode INDEL \
 *          --resource indel-training,training=true indel-training.vcf \
 *          --resource indel-calibration,calibration=true indel-calibration.vcf \
 *          --resource extracted,extracted=true indel-extract.vcf.gz \
 *          --indel-calibration-sensitivity-threshold 0.99 \
 *          -O output
 * </pre>
 *
 * <h3>Custom modeling/scoring backends (ADVANCED)</h3>
 *
 * <p>
 *     The primary scoring functionality performed by this tool is accomplished by a "scoring backend"
 *     whose fundamental contract is to take an input annotation matrix and to output corresponding scores,
 *     with both input and output given as HDF5 files. Rather than using one of the available, implemented backends,
 *     advanced users may provide their own backend via the {@value PYTHON_SCRIPT_LONG_NAME} argument.
 *     See documentation in the modeling and scoring interfaces ({@link VariantAnnotationsModel} and
 *     {@link VariantAnnotationsScorer}, respectively), as well as the default Python IsolationForest implementation at
 *     {@link PythonSklearnVariantAnnotationsScorer} and
 *     org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/isolation-forest.py.
 * </p>
 *
 * DEVELOPER NOTE: See documentation in {@link LabeledVariantAnnotationsWalker}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Scores variant calls in a VCF file based on site-level annotations using a previously trained model.",
        oneLineSummary = "Scores variant calls in a VCF file based on site-level annotations using a previously trained model",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class ScoreVariantAnnotations extends LabeledVariantAnnotationsWalker {

    public static final String MODEL_PREFIX_LONG_NAME = "model-prefix";
    public static final String MODEL_BACKEND_LONG_NAME = TrainVariantAnnotationsModel.MODEL_BACKEND_LONG_NAME;
    public static final String PYTHON_SCRIPT_LONG_NAME = "python-script";
    public static final String SNP_CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME = "snp-calibration-sensitivity-threshold";
    public static final String INDEL_CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME = "indel-calibration-sensitivity-threshold";

    public static final String SNP_KEY_LONG_NAME = "snp-key";
    public static final String SCORE_KEY_LONG_NAME = "score-key";
    public static final String CALIBRATION_SENSITIVITY_KEY_LONG_NAME = "calibration-sensitivity-key";
    public static final String LOW_SCORE_FILTER_NAME_LONG_NAME = "low-score-filter-name";
    public static final String DOUBLE_FORMAT_LONG_NAME = "double-format";

    public static final String DEFAULT_SNP_KEY = LabeledVariantAnnotationsData.SNP_LABEL;
    public static final String DEFAULT_SCORE_KEY = "SCORE";
    public static final String DEFAULT_CALIBRATION_SENSITIVITY_KEY = "CALIBRATION_SENSITIVITY";
    public static final String DEFAULT_LOW_SCORE_FILTER_NAME = "LOW_SCORE";
    public static final String DEFAULT_DOUBLE_FORMAT = "%.4f";

    public static final String SCORES_HDF5_SUFFIX = ".scores.hdf5";

    @Argument(
            fullName = MODEL_PREFIX_LONG_NAME)
    private String modelPrefix;

    @Argument(
            fullName = MODEL_BACKEND_LONG_NAME,
            doc = "Backend to use for scoring. " +
                    "JAVA_BGMM will use a pure Java implementation (ported from Python scikit-learn) of the Bayesian Gaussian Mixture Model. " +
                    "PYTHON_IFOREST will use the Python scikit-learn implementation of the IsolationForest method and " +
                    "will require that the corresponding Python dependencies are present in the environment. " +
                    "PYTHON_SCRIPT will use the script specified by the " + PYTHON_SCRIPT_LONG_NAME + " argument. " +
                    "See the tool documentation for more details." )
    private VariantAnnotationsModelBackend modelBackend = VariantAnnotationsModelBackend.PYTHON_IFOREST;

    @Argument(
            fullName = PYTHON_SCRIPT_LONG_NAME,
            doc = "Python script used for specifying a custom scoring backend. If provided, " + MODEL_BACKEND_LONG_NAME + " must also be set to PYTHON_SCRIPT.",
            optional = true)
    private File pythonScriptFile;

    @Argument(
            fullName = SNP_CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME,
            doc = "If specified, SNPs with scores corresponding to a calibration sensitivity that is greater than or equal to this threshold will be hard filtered.",
            optional = true,
            minValue = 0.,
            maxValue = 1.)
    private Double snpCalibrationSensitivityThreshold;

    @Argument(
            fullName = INDEL_CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME,
            doc = "If specified, indels with scores corresponding to a calibration sensitivity that is greater than or equal to this threshold will be hard filtered.",
            optional = true,
            minValue = 0.,
            maxValue = 1.)
    private Double indelCalibrationSensitivityThreshold;

    @Argument(
            fullName = SNP_KEY_LONG_NAME,
            doc = "Annotation flag to use for labeling sites as SNPs in output. " +
                    "Set this to \"null\" to omit these labels.")
    private String snpKey = DEFAULT_SNP_KEY;

    @Argument(
            fullName = SCORE_KEY_LONG_NAME,
            doc = "Annotation key to use for score values in output.")
    private String scoreKey = DEFAULT_SCORE_KEY;

    @Argument(
            fullName = CALIBRATION_SENSITIVITY_KEY_LONG_NAME,
            doc = "Annotation key to use for calibration-sensitivity values in output.")
    private String calibrationSensitivityKey = DEFAULT_CALIBRATION_SENSITIVITY_KEY;

    @Argument(
            fullName = LOW_SCORE_FILTER_NAME_LONG_NAME,
            doc = "Name to use for low-score filter in output.")
    private String lowScoreFilterName = DEFAULT_LOW_SCORE_FILTER_NAME;

    @Argument(
            fullName = DOUBLE_FORMAT_LONG_NAME,
            doc = "Format string to use for formatting score and calibration-sensitivity values in output.")
    private String doubleFormat = DEFAULT_DOUBLE_FORMAT;

    private File outputScoresFile;
    private Iterator<Double> scoresIterator;
    private Iterator<Boolean> isSNPIterator;

    private VariantAnnotationsScorer snpScorer;
    private VariantAnnotationsScorer indelScorer;

    private Function<Double, Double> snpCalibrationSensitivityConverter;
    private Function<Double, Double> indelCalibrationSensitivityConverter;

    @Override
    protected int numberOfPasses() {
        return 2;
    }

    @Override
    public void afterOnTraversalStart() {

        Utils.nonNull(scoreKey);
        Utils.nonNull(calibrationSensitivityKey);
        Utils.nonNull(lowScoreFilterName);
        Utils.nonNull(doubleFormat);

        switch (modelBackend) {
            case JAVA_BGMM:
                Utils.validateArg(pythonScriptFile == null,
                        "Python script should not be provided when using JAVA_BGMM backend.");
                logger.info("Running in JAVA_BGMM mode...");
                snpScorer = deserializeScorerFromSerFiles(VariantType.SNP);
                indelScorer = deserializeScorerFromSerFiles(VariantType.INDEL);
                break;
            case PYTHON_IFOREST:
                Utils.validateArg(pythonScriptFile == null,
                        "Python script should not be provided when using PYTHON_IFOREST backend.");

                pythonScriptFile = IOUtils.writeTempResource(new Resource(TrainVariantAnnotationsModel.ISOLATION_FOREST_PYTHON_SCRIPT, TrainVariantAnnotationsModel.class));
                PythonScriptExecutor.checkPythonEnvironmentForPackage("argparse");
                PythonScriptExecutor.checkPythonEnvironmentForPackage("h5py");
                PythonScriptExecutor.checkPythonEnvironmentForPackage("numpy");
                PythonScriptExecutor.checkPythonEnvironmentForPackage("sklearn");
                PythonScriptExecutor.checkPythonEnvironmentForPackage("dill");
                logger.info("Running in PYTHON_IFOREST mode...");
                snpScorer = deserializeScorerFromPklFiles(VariantType.SNP);
                indelScorer = deserializeScorerFromPklFiles(VariantType.INDEL);
                break;
            case PYTHON_SCRIPT:
                IOUtils.canReadFile(pythonScriptFile);
                logger.info("Running in PYTHON_SCRIPT mode...");
                snpScorer = deserializeScorerFromPklFiles(VariantType.SNP);
                indelScorer = deserializeScorerFromPklFiles(VariantType.INDEL);
                break;
            default:
                throw new GATKException.ShouldNeverReachHereException("Unknown model-backend mode.");
        }

        if (snpScorer == null && indelScorer == null) {
            throw new UserException.BadInput(String.format("At least one serialized scorer must be present " +
                    "in the model files with the prefix %s.", modelPrefix));
        }
        if (variantTypesToExtract.contains(VariantType.SNP) && snpScorer == null) {
            throw new UserException.BadInput(String.format("SNPs were indicated for extraction via the %s argument, " +
                    "but no serialized SNP scorer was available in the model files with the prefix.", MODE_LONG_NAME, modelPrefix));
        }
        if (variantTypesToExtract.contains(VariantType.INDEL) && indelScorer == null) {
            throw new UserException.BadInput(String.format("INDELs were indicated for extraction via the %s argument, " +
                    "but no serialized INDEL scorer was available in the model files with the prefix.", MODE_LONG_NAME, modelPrefix));
        }

        snpCalibrationSensitivityConverter = readCalibrationScoresAndCreateConverter(VariantType.SNP);
        indelCalibrationSensitivityConverter = readCalibrationScoresAndCreateConverter(VariantType.INDEL);

        if (snpCalibrationSensitivityConverter == null && snpCalibrationSensitivityThreshold != null) {
            throw new UserException.BadInput(String.format("The %s argument was specified, " +
                            "but no SNP calibration scores were provided in the model files with the prefix %s.",
                    SNP_CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME, modelPrefix));
        }
        if (indelCalibrationSensitivityConverter == null && indelCalibrationSensitivityThreshold != null) {
            throw new UserException.BadInput(String.format("The %s argument was specified, " +
                    "but no INDEL calibration scores were provided in the model files with the prefix %s.",
                    INDEL_CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME, modelPrefix));
        }

        outputScoresFile = new File(outputPrefix + SCORES_HDF5_SUFFIX);

        // TODO this validation method should perhaps be moved outside of the CNV code
        CopyNumberArgumentValidationUtils.validateOutputFiles(outputScoresFile);
    }

    @Override
    protected void nthPassApply(final VariantContext variant,
                                final ReadsContext readsContext,
                                final ReferenceContext referenceContext,
                                final FeatureContext featureContext,
                                final int n) {
        final List<Triple<List<Allele>, VariantType, TreeSet<String>>> metadata = extractVariantMetadata(variant, featureContext, true);
        final boolean isVariantExtracted = !metadata.isEmpty();
        if (n == 0 && isVariantExtracted) {
            addExtractedVariantToData(data, variant, metadata);
        }
        if (n == 1) {
            if (isVariantExtracted) {
                writeExtractedVariantToVCF(variant, metadata);
            } else {
                vcfWriter.add(variant);
            }
        }
    }

    @Override
    protected void afterNthPass(final int n) {
        if (n == 0) {
            // TODO if BGMM, preprocess annotations and write to HDF5 with BGMMVariantAnnotationsScorer.preprocessAnnotationsWithBGMMAndWriteHDF5
            writeAnnotationsToHDF5();
            if (data.size() > 0) {
                data.clear();
                readAnnotationsAndWriteScoresToHDF5();
                scoresIterator = Arrays.stream(VariantAnnotationsScorer.readScores(outputScoresFile)).iterator();
                isSNPIterator = LabeledVariantAnnotationsData.readLabel(outputAnnotationsFile, LabeledVariantAnnotationsData.SNP_LABEL).iterator();
            } else {
                scoresIterator = Collections.emptyIterator();
                isSNPIterator = Collections.emptyIterator();
            }
        }
        if (n == 1) {
            if (scoresIterator.hasNext()) {
                throw new IllegalStateException("Traversals of scores and variants " +
                        "(or alleles, in allele-specific mode) were not correctly synchronized.");
            }
            if (vcfWriter != null) {
                vcfWriter.close();
            }
        }
    }

    private VariantAnnotationsScorer deserializeScorerFromPklFiles(final VariantType variantType) {
        final String variantTypeTag = '.' + variantType.toString().toLowerCase();
        final File scorerPklFile = new File(
                modelPrefix + variantTypeTag + PythonSklearnVariantAnnotationsScorer.PYTHON_SCORER_PKL_SUFFIX);
        final File negativeScorerPklFile = new File(
                modelPrefix + variantTypeTag + TrainVariantAnnotationsModel.NEGATIVE_TAG + PythonSklearnVariantAnnotationsScorer.PYTHON_SCORER_PKL_SUFFIX);
        return scorerPklFile.canRead()
                ? negativeScorerPklFile.canRead()
                ? VariantAnnotationsScorer.combinePositiveAndNegativeScorer(
                new PythonSklearnVariantAnnotationsScorer(pythonScriptFile, scorerPklFile),
                new PythonSklearnVariantAnnotationsScorer(pythonScriptFile, negativeScorerPklFile))
                : new PythonSklearnVariantAnnotationsScorer(pythonScriptFile, scorerPklFile)
                : null;
    }

    private VariantAnnotationsScorer deserializeScorerFromSerFiles(final VariantType variantType) {
        final String variantTypeTag = '.' + variantType.toString().toLowerCase();
        final File scorerSerFile = new File(
                modelPrefix + variantTypeTag + BGMMVariantAnnotationsScorer.BGMM_SCORER_SER_SUFFIX);
        final File negativeScorerSerFile = new File(
                modelPrefix + variantTypeTag + TrainVariantAnnotationsModel.NEGATIVE_TAG + BGMMVariantAnnotationsScorer.BGMM_SCORER_SER_SUFFIX);
        return scorerSerFile.canRead()
                ? negativeScorerSerFile.canRead()
                ? VariantAnnotationsScorer.combinePositiveAndNegativeScorer(
                BGMMVariantAnnotationsScorer.deserialize(scorerSerFile),
                BGMMVariantAnnotationsScorer.deserialize(negativeScorerSerFile))
                : BGMMVariantAnnotationsScorer.deserialize(scorerSerFile)
                : null;
    }

    private Function<Double, Double> readCalibrationScoresAndCreateConverter(final VariantType variantType) {
        final String variantTypeTag = '.' + variantType.toString().toLowerCase();
        final File calibrationScores = new File(
                modelPrefix + variantTypeTag + TrainVariantAnnotationsModel.CALIBRATION_SCORES_HDF5_SUFFIX);
        return calibrationScores.canRead()
                ? VariantAnnotationsScorer.createScoreToCalibrationSensitivityConverter(VariantAnnotationsScorer.readScores(calibrationScores))
                : null;
    }

    private void readAnnotationsAndWriteScoresToHDF5() {
        final List<String> annotationNames = LabeledVariantAnnotationsData.readAnnotationNames(outputAnnotationsFile);
        final List<Boolean> isSNP = LabeledVariantAnnotationsData.readLabel(outputAnnotationsFile, LabeledVariantAnnotationsData.SNP_LABEL);
        final double[][] allAnnotations = LabeledVariantAnnotationsData.readAnnotations(outputAnnotationsFile);
        final int numAll = allAnnotations.length;
        final List<Double> allScores = new ArrayList<>(Collections.nCopies(numAll, Double.NaN));
        if (variantTypesToExtract.contains(VariantType.SNP)) {
            logger.info("Scoring SNP variants...");
            scoreVariantTypeAndSetElementsOfAllScores(annotationNames, allAnnotations, isSNP, snpScorer, allScores);
        }
        if (variantTypesToExtract.contains(VariantType.INDEL)) {
            logger.info("Scoring INDEL variants...");
            final List<Boolean> isIndel = isSNP.stream().map(x -> !x).collect(Collectors.toList());
            scoreVariantTypeAndSetElementsOfAllScores(annotationNames, allAnnotations, isIndel, indelScorer, allScores);
        }
        VariantAnnotationsScorer.writeScores(outputScoresFile, Doubles.toArray(allScores));
        logger.info(String.format("Scores written to %s.", outputScoresFile.getAbsolutePath()));
    }

    private static void scoreVariantTypeAndSetElementsOfAllScores(final List<String> annotationNames,
                                                                  final double[][] allAnnotations,
                                                                  final List<Boolean> isVariantType,
                                                                  final VariantAnnotationsScorer variantTypeScorer,
                                                                  final List<Double> allScores) {
        final File variantTypeAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, allAnnotations, isVariantType);
        final File variantTypeScoresFile = IOUtils.createTempFile("temp", ".scores.hdf5");
        variantTypeScorer.score(variantTypeAnnotationsFile, variantTypeScoresFile); // TODO we do not fail until here in the case of mismatched annotation names; we could fail earlier
        final double[] variantTypeScores = VariantAnnotationsScorer.readScores(variantTypeScoresFile);
        final Iterator<Double> variantTypeScoresIterator = Arrays.stream(variantTypeScores).iterator();
        IntStream.range(0, allScores.size()).filter(isVariantType::get).forEach(i -> allScores.set(i, variantTypeScoresIterator.next()));
    }

    @Override
    void writeExtractedVariantToVCF(final VariantContext vc,
                                    final List<Allele> altAlleles,
                                    final Set<String> labels) {
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        labels.forEach(l -> builder.attribute(l, true)); // labels should already be sorted as a TreeSet

        final List<Double> scores = useASAnnotations
                ? altAlleles.stream().map(a -> scoresIterator.next()).collect(Collectors.toList())
                : Collections.singletonList(scoresIterator.next());
        final double score = Collections.max(scores);
        final int scoreIndex = scores.indexOf(score);
        builder.attribute(scoreKey, formatDouble(score));

        final List<Boolean> isSNP = useASAnnotations
                ? altAlleles.stream().map(a -> isSNPIterator.next()).collect(Collectors.toList())
                : Collections.singletonList(isSNPIterator.next());
        final boolean isSNPMax = isSNP.get(scoreIndex);

        if (snpKey != null) {
            builder.attribute(snpKey, isSNPMax);
        }

        final Function<Double, Double> calibrationSensitivityConverter = isSNPMax ? snpCalibrationSensitivityConverter : indelCalibrationSensitivityConverter;
        if (calibrationSensitivityConverter != null) {
            final double calibrationSensitivity = calibrationSensitivityConverter.apply(score);
            builder.attribute(calibrationSensitivityKey, formatDouble(calibrationSensitivity));
            final Double calibrationSensitivityThreshold = isSNPMax ? snpCalibrationSensitivityThreshold : indelCalibrationSensitivityThreshold;
            if (calibrationSensitivityThreshold != null && calibrationSensitivity >= calibrationSensitivityThreshold) {
                builder.filter(lowScoreFilterName); // TODO does this sufficiently cover the desired behavior when dealing with previously filtered sites, etc.?
            }
        }

        vcfWriter.add(builder.make());
    }

    private String formatDouble(final double x) {
        return String.format(doubleFormat, x);
    }

    /**
     * Copies the header from the input VCF and adds info lines for the score, calibration-sensitivity, and label keys,
     * as well as the filter line.
     */
    @Override
    VCFHeader constructVCFHeader(final List<String> sortedLabels) {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();

        final Set<VCFHeaderLine> hInfo = new HashSet<>(inputHeaders);
        hInfo.add(new VCFInfoHeaderLine(scoreKey, 1, VCFHeaderLineType.Float,
                "Score according to the model applied by ScoreVariantAnnotations"));
        hInfo.add(new VCFInfoHeaderLine(calibrationSensitivityKey, 1, VCFHeaderLineType.Float,
                String.format("Calibration sensitivity corresponding to the value of %s", scoreKey)));
        hInfo.add(new VCFFilterHeaderLine(lowScoreFilterName, "Low score (corresponding to high calibration sensitivity)"));

        hInfo.addAll(getDefaultToolVCFHeaderLines());
        if (snpKey != null) {
            hInfo.add(new VCFInfoHeaderLine(snpKey, 1, VCFHeaderLineType.Flag, "This site was considered a SNP during filtering"));
        }
        hInfo.addAll(sortedLabels.stream()
                .map(l -> new VCFInfoHeaderLine(l, 1, VCFHeaderLineType.Flag, String.format(RESOURCE_LABEL_INFO_HEADER_LINE_FORMAT_STRING, l)))
                .collect(Collectors.toList()));

        return new VCFHeader(hInfo, inputHeader.getGenotypeSamples());
    }

    @Override
    public Object onTraversalSuccess() {

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }
}