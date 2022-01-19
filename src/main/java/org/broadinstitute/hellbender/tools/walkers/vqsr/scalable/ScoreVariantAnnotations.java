package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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
public class ScoreVariantAnnotations extends VariantAnnotationWalker {

    private static final String SCORE_KEY = GATKVCFConstants.VQS_LOD_KEY;
    private static final String DUMMY_ALLELE = "<VQSR>";

    private static final String SCORER_PKL_SUFFIX = ".scorer.pkl";
    private static final String ANNOTATIONS_HDF5_SUFFIX = ".annot.hdf5";
    private static final String SCORES_HDF5_SUFFIX = ".scores.hdf5";
    private static final String RECALIBRATION_VCF_SUFFIX = ".recal.vcf";

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output prefix.")
    private String outputPrefix;

    @Argument(
            fullName = "python-script")
    private File pythonScriptFile;

    @Argument(
            fullName = "model-prefix")
    private String modelPrefix;

    private VariantContextWriter recalWriter;
    private File scorerPklFile;
    private File outputAnnotationsHDF5File;
    private File outputScoresHDF5File;
    private File outputRecalVCFFile;

    @Override
    public void beforeOnTraversalStart() {
        isExtractTrainingAndTruthOnly = false;

        scorerPklFile = new File(modelPrefix + SCORER_PKL_SUFFIX);
        outputAnnotationsHDF5File = new File(outputPrefix + ANNOTATIONS_HDF5_SUFFIX);
        outputScoresHDF5File = new File(outputPrefix + SCORES_HDF5_SUFFIX);
        outputRecalVCFFile = new File(outputPrefix + RECALIBRATION_VCF_SUFFIX);

        IOUtils.canReadFile(pythonScriptFile);
        IOUtils.canReadFile(scorerPklFile);

        PythonScriptExecutor.checkPythonEnvironmentForPackage("sklearn");
        PythonScriptExecutor.checkPythonEnvironmentForPackage("dill");
    }

    @Override
    public void afterTraversalSuccess() {

        logger.info(String.format("Extracted annotations for %s total variants.", dataManager.getData().size()));

        logger.info("Writing annotations...");
        writeAnnotationsHDF5(outputAnnotationsHDF5File);
        logger.info(String.format("Annotations written to %s.", outputAnnotationsHDF5File.getAbsolutePath()));

        logger.info("Scoring...");
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final ProcessOutput pythonProcessOutput = executor.executeScriptAndGetOutput(
                pythonScriptFile.getAbsolutePath(),
                null,
                composePythonArguments(outputAnnotationsHDF5File, scorerPklFile, outputScoresHDF5File));

        if (pythonProcessOutput.getExitValue() != 0) {
            throw executor.getScriptException(executor.getExceptionMessageFromScriptError(pythonProcessOutput));
        }

        try (final HDF5File outputScoresFileHDF5File = new HDF5File(outputScoresHDF5File, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(outputScoresFileHDF5File.getFile());
            final double[] scores = outputScoresFileHDF5File.readDoubleArray("/scores");
            dataManager.setScores(dataManager.getData(), scores);
        } catch (final RuntimeException exception) {
            throw new GATKException(String.format("Exception encountered during reading of scores from %s: %s",
                    outputScoresHDF5File.getAbsolutePath(), exception));
        }
        logger.info("Scoring complete.");

        logger.info("Writing out recalibration VCF...");
        recalWriter = createVCFWriter(outputRecalVCFFile);
        recalWriter.writeHeader(constructVCFHeader());
        writeOutRecalibrationTable(recalWriter, getBestAvailableSequenceDictionary());
        logger.info(String.format("Recalibration VCF written to %s.", outputRecalVCFFile.getAbsolutePath()));

        logger.info(String.format("%s complete.", getClass().getSimpleName()));
    }

    @Override
    public void closeTool(){
        if (recalWriter != null) {
            recalWriter.close();
        }
    }

    private static List<String> composePythonArguments(final File rawAnnotationsFile,
                                                       final File scorerPklFile,
                                                       final File outputScoresFile) {
        try {
            return new ArrayList<>(Arrays.asList(
                    "--raw_annotations_file=" + rawAnnotationsFile.getCanonicalPath(),
                    "--scorer_pkl_file=" + scorerPklFile.getCanonicalPath(),
                    "--output_scores_file=" + outputScoresFile.getCanonicalPath()));
        } catch (final IOException e) {
            throw new UserException.BadInput(String.format("Encountered exception resolving canonical file paths: %s", e));
        }
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

    private void writeOutRecalibrationTable(final VariantContextWriter recalWriter,
                                            final SAMSequenceDictionary seqDictionary) {
        final List<VariantDatum> data = dataManager.getData();

        // we need to sort in coordinate order in order to produce a valid VCF
        data.sort((vd1, vd2) -> IntervalUtils.compareLocatables(vd1.loc, vd2.loc, seqDictionary));

        // create dummy alleles to be used
        List<Allele> alleles = Arrays.asList(Allele.create("N", true), Allele.create(DUMMY_ALLELE, false));

        for (final VariantDatum datum : data) {
            if (useASannotations) {
                alleles = Arrays.asList(datum.referenceAllele, datum.alternateAllele); //use the alleles to distinguish between multiallelics in AS mode
            }
            final VariantContextBuilder builder = new VariantContextBuilder(SCORE_KEY, datum.loc.getContig(), datum.loc.getStart(), datum.loc.getEnd(), alleles);
            builder.attribute(VCFConstants.END_KEY, datum.loc.getEnd());
            builder.attribute(SCORE_KEY, String.format("%.4f", datum.score));

            if (datum.atTrainingSite) {
                builder.attribute(GATKVCFConstants.POSITIVE_LABEL_KEY, true);
            }

            recalWriter.add(builder.make());
        }
    }
}