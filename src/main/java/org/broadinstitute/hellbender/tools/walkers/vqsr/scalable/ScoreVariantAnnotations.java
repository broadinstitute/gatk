//package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;
//
//import htsjdk.samtools.SAMSequenceDictionary;
//import htsjdk.variant.variantcontext.writer.VariantContextWriter;
//import htsjdk.variant.vcf.VCFConstants;
//import htsjdk.variant.vcf.VCFHeader;
//import htsjdk.variant.vcf.VCFHeaderLine;
//import htsjdk.variant.vcf.VCFStandardHeaderLines;
//import org.apache.commons.lang3.tuple.Pair;
//import org.broadinstitute.barclay.argparser.Argument;
//import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
//import org.broadinstitute.barclay.help.DocumentedFeature;
//import org.broadinstitute.hdf5.HDF5File;
//import org.broadinstitute.hdf5.HDF5LibException;
//import org.broadinstitute.hellbender.exceptions.GATKException;
//import org.broadinstitute.hellbender.exceptions.UserException;
//import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
//import org.broadinstitute.hellbender.utils.io.IOUtils;
//import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
//import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
//import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
//import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
//import org.broadinstitute.hellbender.utils.variant.VcfUtils;
//import picard.cmdline.programgroups.VariantFilteringProgramGroup;
//
//import java.io.File;
//import java.io.IOException;
//import java.util.ArrayList;
//import java.util.Arrays;
//import java.util.Collections;
//import java.util.List;
//import java.util.Set;
//
///**
// * TODO
// */
//@CommandLineProgramProperties(
//        // TODO
//        summary = "",
//        oneLineSummary = "",
//        programGroup = VariantFilteringProgramGroup.class
//)
//@DocumentedFeature
//public class ScoreVariantAnnotations extends VariantLabeledAnnotationsWalker {
//
//    private static final String SCORE_KEY = GATKVCFConstants.VQS_LOD_KEY;
//
//    private static final String SCORER_PKL_SUFFIX = ".scorer.pkl";
//    private static final String SCORER_SER_SUFFIX = ".scorer.ser";
//    private static final String SCORES_HDF5_SUFFIX = ".scores.hdf5";
//    private static final String RECALIBRATION_VCF_SUFFIX = ".recal.vcf";
//
//    @Argument(
//            fullName = "python-script",
//            optional = true)
//    private File pythonScriptFile;
//
//    @Argument(
//            fullName = "model-prefix")
//    private String modelPrefix;
//
//    private File inputScorerPklFile;
//    private File outputScoresFile;
//    private File outputVCFFile;
//    private VariantContextWriter vcfWriter;
//
//    @Override
//    public boolean isExtractVariantsNotOverlappingResources() {
//        return true;
//    }
//
//    @Override
//    public void beforeOnTraversalStart() {
//
//        if (pythonScriptFile != null) {
//            logger.info("Python script was provided, running in PYTHON mode...");
//
//            inputScorerPklFile = new File(modelPrefix + SCORER_PKL_SUFFIX);
//
//            IOUtils.canReadFile(pythonScriptFile);
//            IOUtils.canReadFile(inputScorerPklFile);
//
//            PythonScriptExecutor.checkPythonEnvironmentForPackage("sklearn");
//            PythonScriptExecutor.checkPythonEnvironmentForPackage("dill");
//        } else {
//            // TODO validate BGMM model inputs
//        }
//
//        outputScoresFile = new File(outputPrefix + SCORES_HDF5_SUFFIX);
//        outputVCFFile = new File(outputPrefix + ".vcf"); // TODO add -O VCF output
//
//        for (final File outputFile : Collections.singletonList(outputScoresFile)) {
//            if ((outputFile.exists() && !outputFile.canWrite()) ||
//                    (!outputFile.exists() && !outputFile.getAbsoluteFile().getParentFile().canWrite())) {
//                throw new UserException(String.format("Cannot create output file at %s.", outputFile));
//            }
//        }
//
//        vcfWriter = createVCFWriter(outputVCFFile);
//        vcfWriter.writeHeader(constructVCFHeader());
//    }
//
//    @Override
//    public void afterOnTraversalSuccess() {
//
//        logger.info(String.format("Extracted annotations for %s total variants.", dataBatch.getFlattenedData().size()));
//
//        logger.info("Scoring...");
//        final double[] scores;
//        if (pythonScriptFile != null) {
//            final PythonScriptExecutor executor = new PythonScriptExecutor(true);
//            final ProcessOutput pythonProcessOutput = executor.executeScriptAndGetOutput(
//                    pythonScriptFile.getAbsolutePath(),
//                    null,
//                    composePythonArguments(outputAnnotationsFile, inputScorerPklFile, outputScoresFile));
//
//            if (pythonProcessOutput.getExitValue() != 0) {
//                throw executor.getScriptException(executor.getExceptionMessageFromScriptError(pythonProcessOutput));
//            }
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
//        VariantLabeledAnnotationsData.setScores(dataBatch.getFlattenedData(), scores);
//        logger.info("Scoring complete.");
//
//        logger.info("Writing VCF...");
//        writeVCF(false, false,true);
//
//        vcfWriter.close();
//
//        logger.info(String.format("%s complete.", getClass().getSimpleName()));
//    }
//
//    private static List<String> composePythonArguments(final File rawAnnotationsFile,
//                                                       final File scorerPklFile,
//                                                       final File outputScoresFile) {
//        try {
//            return new ArrayList<>(Arrays.asList(
//                    "--raw_annotations_file=" + rawAnnotationsFile.getCanonicalPath(),
//                    "--scorer_pkl_file=" + scorerPklFile.getCanonicalPath(),
//                    "--output_scores_file=" + outputScoresFile.getCanonicalPath()));
//        } catch (final IOException e) {
//            throw new UserException.BadInput(String.format("Encountered exception resolving canonical file paths: %s", e));
//        }
//    }
//
//    private VCFHeader constructVCFHeader() {
//        //TODO: this should be refactored/consolidated as part of
//        // https://github.com/broadinstitute/gatk/issues/2112
//        // https://github.com/broadinstitute/gatk/issues/121,
//        // https://github.com/broadinstitute/gatk/issues/1116 and
//        // Initialize VCF header lines
//        Set<VCFHeaderLine> hInfo = getDefaultToolVCFHeaderLines();
//        hInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
//        hInfo.add(GATKVCFHeaderLines.getInfoLine(SCORE_KEY));
//        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.POSITIVE_LABEL_KEY));
//        hInfo.add(GATKVCFHeaderLines.getFilterLine(VCFConstants.PASSES_FILTERS_v4));
//        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
//        if (hasReference()) {
//            hInfo = VcfUtils.updateHeaderContigLines(
//                    hInfo, referenceArguments.getReferencePath(), sequenceDictionary, true);
//        } else if (null != sequenceDictionary) {
//            hInfo = VcfUtils.updateHeaderContigLines(hInfo, null, sequenceDictionary, true);
//        }
//        return new VCFHeader(hInfo);
//    }
//}