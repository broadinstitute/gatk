package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.BGMMVariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.BGMMVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.PythonSklearnVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsScorer;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.Collections;
import java.util.Set;

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
public class ScoreVariantAnnotations extends LabeledVariantAnnotationsBatchWalker {

    private static final String SCORE_KEY = GATKVCFConstants.VQS_LOD_KEY;

    private static final String SCORER_PKL_SUFFIX = ".scorer.pkl";
    private static final String SCORER_SER_SUFFIX = ".scorer.ser";
    private static final String SCORES_HDF5_SUFFIX = ".scores.hdf5";

    @Argument(
            fullName = "python-script",
            optional = true)
    private File pythonScriptFile;

    @Argument(
            fullName = "model-prefix")
    private String modelPrefix;

    private File outputScoresFile;
    private File outputVCFFile;
    private VariantContextWriter vcfWriter;
    private VariantAnnotationsScorer snpScorer;
    private VariantAnnotationsScorer indelScorer;

    @Override
    public boolean isExtractVariantsNotOverlappingResources() {
        return true;
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

        outputScoresFile = new File(outputPrefix + SCORES_HDF5_SUFFIX);
        outputVCFFile = new File(outputPrefix + ".vcf"); // TODO add -O VCF output

        for (final File outputFile : Collections.singletonList(outputScoresFile)) {
            if ((outputFile.exists() && !outputFile.canWrite()) ||
                    (!outputFile.exists() && !outputFile.getAbsoluteFile().getParentFile().canWrite())) {
                throw new UserException(String.format("Cannot create output file at %s.", outputFile));
            }
        }

        vcfWriter = createVCFWriter(outputVCFFile);
        vcfWriter.writeHeader(constructVCFHeader());
    }

    @Override
    void doBatchWork() {
//        final List<Boolean> isSNP = dataBatch.getFlattenedData().flatstream().map(datum -> datum.is);
//        if (variantTypesToExtract.contains(VariantType.SNP)) {
//        }
//        if (variantTypesToExtract.contains(VariantType.INDEL)) {
//        }
//        final File outputScoresFile = new File(outputPrefix + outputPrefixTag + outputSuffix);
//        scorer.scoreSamples(annotationsFile, outputScoresFile);
    }

    @Override
    public void afterOnTraversalSuccess() {

//        logger.info(String.format("Extracted annotations for %s total variants.", dataBatch.getFlattenedData().size()));

//        logger.info("Scoring...");
//        final double[] scores;
//        if (pythonScriptFile != null) {
//
//            scores = VariantAnnotationUtils.readScores(outputScoresFile);
//        } else {
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
//        }
//
//        // TODO decouple setting of scores and VCF writing
//        LabeledVariantAnnotationsData.setScores(dataBatch.getFlattenedData(), scores);
//        logger.info("Scoring complete.");
//
//        logger.info("Writing VCF...");
//        writeVCF(false, false,true);

        vcfWriter.close();

        logger.info(String.format("%s complete.", getClass().getSimpleName()));
    }

    private VCFHeader constructVCFHeader() {
        //TODO: this should be refactored/consolidated as part of
        // https://github.com/broadinstitute/gatk/issues/2112
        // https://github.com/broadinstitute/gatk/issues/121,
        // https://github.com/broadinstitute/gatk/issues/1116 and
        // Initialize VCF header lines
        Set<VCFHeaderLine> hInfo = getDefaultToolVCFHeaderLines();
        hInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(SCORE_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.POSITIVE_LABEL_KEY));
        hInfo.add(GATKVCFHeaderLines.getFilterLine(VCFConstants.PASSES_FILTERS_v4));
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        if (hasReference()) {
            hInfo = VcfUtils.updateHeaderContigLines(
                    hInfo, referenceArguments.getReferencePath(), sequenceDictionary, true);
        } else if (null != sequenceDictionary) {
            hInfo = VcfUtils.updateHeaderContigLines(hInfo, null, sequenceDictionary, true);
        }
        return new VCFHeader(hInfo);
    }
}