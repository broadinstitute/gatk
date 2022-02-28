package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.primitives.Doubles;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.BGMMVariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.BGMMVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.PythonSklearnVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsScorer;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
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
 * TODO
 *
 * DEVELOPER NOTE: See documentation in {@link LabeledVariantAnnotationsWalker}.
 */
@CommandLineProgramProperties(
        // TODO
        summary = "",
        oneLineSummary = "",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public class ScoreVariantAnnotations extends LabeledVariantAnnotationsWalker {

    private static final String SCORE_KEY = GATKVCFConstants.VQS_LOD_KEY;
    private static final String TRUTH_SENSITIVITY_KEY = "TRUTH_SENSITIVITY";
    private static final String SCORE_AND_TRUTH_SENSITIVITY_FORMAT = "%.4f";

    private static final String SCORES_HDF5_SUFFIX = ".scores.hdf5";

    @Argument(
            fullName = "python-script",
            optional = true)
    private File pythonScriptFile;

    @Argument(
            fullName = "model-prefix")
    private String modelPrefix;

    private File outputScoresFile;
    private Iterator<Double> scoresIterator;
    private Iterator<Boolean> isSNPIterator;

    private VariantAnnotationsScorer snpScorer;
    private VariantAnnotationsScorer indelScorer;

    private Function<Double, Double> snpTruthSensitivityConverter;
    private Function<Double, Double> indelTruthSensitivityConverter;

    // TODO document, make enum (extract labeled vs. extract all)
    @Override
    public boolean isExtractOnlyLabeledVariants() {
        return false;
    }

    @Override
    protected int numberOfPasses() {
        return 2;
    }

    @Override
    public void beforeOnTraversalStart() {

        if (pythonScriptFile != null) {
            logger.info("Python script was provided, running in PYTHON mode...");

            IOUtils.canReadFile(pythonScriptFile);

            PythonScriptExecutor.checkPythonEnvironmentForPackage("argparse");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("h5py");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("numpy");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("sklearn");
            PythonScriptExecutor.checkPythonEnvironmentForPackage("dill");

            // TODO extract method and constants
            final File snpScorerPklFile = new File(modelPrefix + ".snp.scorer.pkl");
            snpScorer = snpScorerPklFile.canRead()
                    ? new PythonSklearnVariantAnnotationsScorer(pythonScriptFile, snpScorerPklFile)
                    : null;
            final File indelScorerPklFile = new File(modelPrefix + ".indel.scorer.pkl");
            indelScorer = indelScorerPklFile.canRead()
                    ? new PythonSklearnVariantAnnotationsScorer(pythonScriptFile, indelScorerPklFile)
                    : null;
        } else {
            logger.info("Python script was not provided, running in BGMM mode...");
            // TODO validate BGMM model inputs

            // TODO extract method and constants
            final File snpScorerSerFile = new File(modelPrefix + ".snp" + BGMMVariantAnnotationsModel.BGMM_SCORER_SER_SUFFIX);
            snpScorer = snpScorerSerFile.canRead()
                    ? BGMMVariantAnnotationsScorer.deserialize(snpScorerSerFile)
                    : null;
            final File indelScorerSerFile = new File(modelPrefix + ".indel" + BGMMVariantAnnotationsModel.BGMM_SCORER_SER_SUFFIX);
            indelScorer = indelScorerSerFile.canRead()
                    ? BGMMVariantAnnotationsScorer.deserialize(indelScorerSerFile)
                    : null;
        }

        // TODO validate modes
        if (snpScorer == null && indelScorer == null) {
            throw new UserException.BadInput("At least one serialized scorer must be present in model files.");
        }

        // TODO extract method and constants
        final File snpTruthScores = new File(modelPrefix + ".snp.truthScores.hdf5");
        snpTruthSensitivityConverter = snpTruthScores.canRead()
                ? VariantAnnotationsScorer.createScoreToTruthSensitivityConverter(VariantAnnotationsScorer.readScores(snpTruthScores))
                : null;
        final File indelTruthScores = new File(modelPrefix + ".indel.truthScores.hdf5");
        indelTruthSensitivityConverter = indelTruthScores.canRead()
                ? VariantAnnotationsScorer.createScoreToTruthSensitivityConverter(VariantAnnotationsScorer.readScores(indelTruthScores))
                : null;

        outputScoresFile = new File(outputPrefix + SCORES_HDF5_SUFFIX);

        for (final File outputFile : Collections.singletonList(outputScoresFile)) {
            if ((outputFile.exists() && !outputFile.canWrite()) ||
                    (!outputFile.exists() && !outputFile.getAbsoluteFile().getParentFile().canWrite())) {
                throw new UserException(String.format("Cannot create output file at %s.", outputFile));
            }
        }
    }

    @Override
    protected void nthPassApply(final VariantContext variant,
                                final ReadsContext readsContext,
                                final ReferenceContext referenceContext,
                                final FeatureContext featureContext,
                                final int n) {
        final List<Triple<List<Allele>, VariantType, TreeSet<String>>> metadata = extractVariantMetadata(variant, featureContext);
        final boolean isVariantExtracted = !metadata.isEmpty();
        if (n == 0 && isVariantExtracted) {
            addExtractedVariantToData(variant, metadata);
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
            // TODO if BGMM, preprocess annotations and write to HDF5
//            // write preprocessed annotations
//            // TODO clean this up
//            final List<String> annotationNames = this.dataBatch.getSortedAnnotationKeys();
//
//            final File outputPreprocessedAnnotationsFile = new File(outputPrefix + ".annot.pre.hdf5");
//            try (final HDF5File hdf5File = new HDF5File(outputPreprocessedAnnotationsFile, HDF5File.OpenMode.CREATE)) { // TODO allow appending
//                IOUtils.canReadFile(hdf5File.getFile());
//
//                hdf5File.makeStringArray("/data/annotation_names", annotationNames.toArray(new String[0]));
//                HDF5Utils.writeChunkedDoubleMatrix(hdf5File, "/data/annotations", preprocessedData, maximumChunkSize);
//                hdf5File.makeDoubleArray("/data/is_training", this.dataBatch.getFlattenedData().stream().mapToDouble(x -> x.labels.contains("training") ? 1 : 0).toArray());
//            } catch (final HDF5LibException exception) {
//                throw new GATKException(String.format("Exception encountered during writing of preprocessed annotations (%s). Output file at %s may be in a bad state.",
//                        exception, outputPreprocessedAnnotationsFile.getAbsolutePath()));
//            }
//            logger.info(String.format("Preprocessed annotations written to %s.", outputPreprocessedAnnotationsFile.getAbsolutePath()));
            writeAnnotationsToHDF5AndClearData();
            readAnnotationsAndWriteScoresToHDF5();
            scoresIterator = Arrays.stream(VariantAnnotationsScorer.readScores(outputScoresFile)).iterator();
            isSNPIterator = LabeledVariantAnnotationsData.readLabel(outputAnnotationsFile, "snp").iterator();
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

    private void readAnnotationsAndWriteScoresToHDF5() {
        final List<String> annotationNames = LabeledVariantAnnotationsData.readAnnotationNames(outputAnnotationsFile);
        final List<Boolean> isSNP = LabeledVariantAnnotationsData.readLabel(outputAnnotationsFile, "snp");
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
        variantTypeScorer.score(variantTypeAnnotationsFile, variantTypeScoresFile);
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
        builder.attribute(SCORE_KEY, formatDouble(score));

        final List<Boolean> isSNP = useASAnnotations
                ? altAlleles.stream().map(a -> isSNPIterator.next()).collect(Collectors.toList())
                : Collections.singletonList(isSNPIterator.next());
        final boolean isSNPMax = isSNP.get(scoreIndex);
        builder.attribute("snp", isSNPMax);

        final Function<Double, Double> truthSensitivityConverter = isSNPMax ? snpTruthSensitivityConverter : indelTruthSensitivityConverter;
        if (truthSensitivityConverter != null) {
            final double truthSensitivity = truthSensitivityConverter.apply(score);
            builder.attribute(TRUTH_SENSITIVITY_KEY, formatDouble(truthSensitivity));
        }

        vcfWriter.add(builder.make());
    }

    private static String formatDouble(final double x) {
        return String.format(SCORE_AND_TRUTH_SENSITIVITY_FORMAT, x);
    }

    /**
     * Copies the header from the input VCF and adds info lines for the score and label keys.
     */
    @Override
    VCFHeader constructVCFHeader(final List<String> sortedLabels) {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();

        final Set<VCFHeaderLine> hInfo = new HashSet<>(inputHeaders);
        hInfo.add(GATKVCFHeaderLines.getInfoLine(SCORE_KEY));
        hInfo.add(new VCFInfoHeaderLine(TRUTH_SENSITIVITY_KEY, 1, VCFHeaderLineType.Float,
                String.format("Truth sensitivity corresponding to the value of %s", SCORE_KEY)));

        hInfo.addAll(getDefaultToolVCFHeaderLines());
        // TODO extract
        hInfo.add(new VCFInfoHeaderLine("snp", 1, VCFHeaderLineType.Flag, "This site was considered a SNP during filtering"));
        hInfo.addAll(sortedLabels.stream()
                .map(l -> new VCFInfoHeaderLine(l, 1, VCFHeaderLineType.Flag, String.format("This site was labeled as %s according to resources", l)))
                .collect(Collectors.toList()));

        return new VCFHeader(hInfo, inputHeader.getGenotypeSamples());
    }

    @Override
    public Object onTraversalSuccess() {

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }
}