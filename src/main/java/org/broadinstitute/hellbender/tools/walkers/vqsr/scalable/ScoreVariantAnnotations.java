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
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * TODO
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

    private VariantAnnotationsScorer snpScorer;
    private VariantAnnotationsScorer indelScorer;

    // TODO document, make enum (extract labeled vs. extract all)
    @Override
    public boolean isExtractOnlyLabeledVariants() {
        return false;
    }

    @Override
    public boolean isOutputSitesOnlyVCF() {
        return outputSitesOnlyVCFs;
    }

    @Override
    protected int numberOfPasses() {
        return pythonScriptFile == null ? 1 : 2;
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

            final File snpScorerSerFile = new File(modelPrefix + ".snp" + BGMMVariantAnnotationsModel.BGMM_SCORER_SER_SUFFIX);
            snpScorer = snpScorerSerFile.canRead()
                    ? BGMMVariantAnnotationsScorer.deserialize(snpScorerSerFile)
                    : null;
            final File indelScorerSerFile = new File(modelPrefix + ".indel" + BGMMVariantAnnotationsModel.BGMM_SCORER_SER_SUFFIX);
            indelScorer = indelScorerSerFile.canRead()
                    ? BGMMVariantAnnotationsScorer.deserialize(indelScorerSerFile)
                    : null;
        }

        if (snpScorer == null && indelScorer == null) {
            throw new UserException.BadInput("At least one serialized scorer must be present in model files.");
        }

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
        final List<Triple<List<Allele>, VariantType, Set<String>>> metadata = extractMetadata(variant, featureContext);
        if (!metadata.isEmpty()) {
            if (pythonScriptFile == null) {
                if (n == 0) {
                    data.add(variant,
                            metadata.stream().map(Triple::getLeft).collect(Collectors.toList()),
                            metadata.stream().map(Triple::getMiddle).collect(Collectors.toList()),
                            metadata.stream().map(Triple::getRight).collect(Collectors.toList()));
                    final List<Double> scores = data.getData().get(data.size() - 1).stream()
                            .map(d -> d.getVariantType() == VariantType.SNP
                                    ? snpScorer.scoreSamples(Stream.of(d.getAnnotations()).toArray(double[][]::new))[0]
                                    : indelScorer.scoreSamples(Stream.of(d.getAnnotations()).toArray(double[][]::new))[0])
                            .collect(Collectors.toList());
                    scoresIterator = scores.listIterator();
                    writeVariantToVCF(variant,
                            metadata.stream().map(Triple::getLeft).flatMap(List::stream).collect(Collectors.toList()),
                            metadata.stream().map(Triple::getRight).flatMap(Set::stream).collect(Collectors.toSet()));
                }
            } else {
                if (n == 0) {
                    data.add(variant,
                            metadata.stream().map(Triple::getLeft).collect(Collectors.toList()),
                            metadata.stream().map(Triple::getMiddle).collect(Collectors.toList()),
                            metadata.stream().map(Triple::getRight).collect(Collectors.toList()));
                } else if (n == 1) {
                    writeVariantToVCF(variant,
                            metadata.stream().map(Triple::getLeft).flatMap(List::stream).collect(Collectors.toList()),
                            metadata.stream().map(Triple::getRight).flatMap(Set::stream).collect(Collectors.toSet()));
                }
            }
        } else if (n == numberOfPasses()) {
            vcfWriter.add(variant);
        }
    }

    @Override
    protected void afterNthPass(final int n) {
        if (pythonScriptFile == null) {
            writeAnnotationsToHDF5();
            writeScoresToHDF5();
        } else {
            if (n == 0) {
                writeAnnotationsToHDF5();
                data.clear();
                writeScoresToHDF5();
                scoresIterator = Arrays.stream(VariantAnnotationsScorer.readScores(outputScoresFile)).iterator();
            }
        }
        if (n == numberOfPasses()) {
            if (vcfWriter != null) {
                vcfWriter.close();
            }
        }
    }

    private void writeScoresToHDF5() {
        final List<String> annotationNames = LabeledVariantAnnotationsData.readAnnotationNames(outputAnnotationsFile);
        final List<Boolean> isSNP = LabeledVariantAnnotationsData.readLabel(outputAnnotationsFile, "snp");
        final double[][] allAnnotations = LabeledVariantAnnotationsData.readAnnotations(outputAnnotationsFile);
        final int numAll = allAnnotations.length;
        final List<Double> allScores = new ArrayList<>(Collections.nCopies(numAll, Double.NaN));
        if (variantTypesToExtract.contains(VariantType.SNP)) {
            final File snpAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, allAnnotations, isSNP);
            final File snpScoresFile = IOUtils.createTempFile("snp", ".scores.hdf5");
            snpScorer.scoreSamples(snpAnnotationsFile, snpScoresFile);
            final double[] snpScores = VariantAnnotationsScorer.readScores(snpScoresFile);
            final Iterator<Double> snpScoresIterator = Arrays.stream(snpScores).iterator();
            IntStream.range(0, numAll).filter(isSNP::get).forEach(i -> allScores.set(i, snpScoresIterator.next()));
        }
        if (variantTypesToExtract.contains(VariantType.INDEL)) {
            final List<Boolean> isIndel = isSNP.stream().map(x -> !x).collect(Collectors.toList());
            final File indelAnnotationsFile = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(annotationNames, allAnnotations, isIndel);
            final File indelScoresFile = IOUtils.createTempFile("indel", ".scores.hdf5");
            indelScorer.scoreSamples(indelAnnotationsFile, indelScoresFile);
            final double[] indelScores = VariantAnnotationsScorer.readScores(indelScoresFile);
            final Iterator<Double> indelScoresIterator = Arrays.stream(indelScores).iterator();
            IntStream.range(0, numAll).filter(isIndel::get).forEach(i -> allScores.set(i, indelScoresIterator.next()));
        }
        VariantAnnotationsScorer.writeScores(outputScoresFile, Doubles.toArray(allScores));
        logger.info(String.format("Scores written to %s.", outputScoresFile.getAbsolutePath()));
    }

    @Override
    void writeVariantToVCF(final VariantContext vc,
                           final List<Allele> altAlleles,
                           final Set<String> labels) {
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        final List<String> sortedLabels = labels.stream().sorted().collect(Collectors.toList());
        sortedLabels.forEach(l -> builder.attribute(l, true));
        final double bestScore = useASAnnotations
                ? altAlleles.stream().mapToDouble(a -> scoresIterator.next()).max().getAsDouble()
                : scoresIterator.next();
        builder.attribute(SCORE_KEY, bestScore);
        vcfWriter.add(builder.make());
    }

    @Override
    VCFHeader constructVCFHeader(final List<String> sortedLabels) {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();

        final Set<VCFHeaderLine> hInfo = new HashSet<>(inputHeaders);
        hInfo.add(GATKVCFHeaderLines.getInfoLine(SCORE_KEY));

        hInfo.addAll(getDefaultToolVCFHeaderLines());
        hInfo.addAll(sortedLabels.stream()
                .map(l -> new VCFInfoHeaderLine(l, 1, VCFHeaderLineType.Flag, String.format("This site was labeled as %s according to resources", l)))
                .collect(Collectors.toList()));

        return new VCFHeader(hInfo, inputHeader.getGenotypeSamples());
    }

    @Override
    public Object onTraversalSuccess() {
//            final VariantAnnotationUtils.Scorer scorer = VariantAnnotationUtils.deserialize(
//                    new File(modelPrefix + SCORER_SER_SUFFIX), // TODO clean up
//                    VariantAnnotationUtils.Scorer.class);
//            final double[][] data = this.dataBatch.getFlattenedData().stream().map(vd -> vd.annotations).toArray(double[][]::new);
//            final Pair<double[][], double[]> preprocessedDataAndScores = scorer.preprocessAndScoreSamples(data);
//            final double[][] preprocessedData = preprocessedDataAndScores.getLeft();
//            scores = preprocessedDataAndScores.getRight();
//            VariantAnnotationUtils.writeScores(outputScoresFile, scores);
//
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

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

//    @Override
//    void writeVariantToVCF(final VariantContext vc,
//                           final Triple<List<Allele>, VariantType, Set<String>> metadata) {
//        // TODO validate
//
//        final VariantContextBuilder builder = new VariantContextBuilder(SCORE_KEY, datum.loc.getContig(), datum.loc.getStart(), datum.loc.getEnd(), alleles);
//        builder.attribute(VCFConstants.END_KEY, datum.loc.getEnd());
//
//        if (writeScores) {
//            builder.attribute(SCORE_KEY, String.format("%.4f", datum.score));
//        }
//
//        if (datum.labels.contains(VariantLabeledAnnotationsData.TRAINING_LABEL)) {
//            builder.attribute(GATKVCFConstants.POSITIVE_LABEL_KEY, true);
//        }
//        vcfWriter.add(builder.make());
//    }
}