package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.acnvconversion.ACNVModeledSegmentConversionUtils;
import org.broadinstitute.hellbender.tools.exome.acsconversion.ACSModeledSegmentUtils;
import org.broadinstitute.hellbender.tools.exome.allelicbalancecaller.AllelicSplitCall;
import org.broadinstitute.hellbender.tools.exome.allelicbalancecaller.CNLOHCaller;
import org.broadinstitute.hellbender.tools.exome.titanconversion.TitanFileConverter;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * CLI for calling ACNV segments as CNLoH and Balanced.  Additionally, provides some useful conversions of the ACNV
 * output file:
 *  - GATK CNV (total copy ratio)
 *  - Broad CGA Allelic CapSeg
 *  - TITAN
 * REQUIRES Apache Spark
 *
 */
@CommandLineProgramProperties(

        summary = "Call whether the segments are balanced (MAF = 0.5).  For detailed explanation, please see the docs.  " +
                "NOTE: This tool only works with diploid organisms.  This tool also converts files into TITAN and " +
                "Broad CGA Allelic CapSeg formats.  This tool uses spark, though running locally is fine.\n" +
                "\nIMPORTANT:  The additional CNLoH calls from this tool should be treated with a lot of skepticism.  Preliminary results " +
                "indicated very poor performance.",
        oneLineSummary = "Call whether segments are balanced (i.e. MAF=0.5) and write ACS and TITAN files converted from ACNV.",
        programGroup = CopyNumberProgramGroup.class
)
public class CallAllelicSplits extends SparkCommandLineProgram {

    static final long serialVersionUID = 42123132L;

    protected static final String RHO_SHORT_NAME="r";
    protected static final String RHO_LONG_NAME="rhoThreshold";
    protected static final int DEFAULT_NUM_ITERATIONS = 10;
    protected static final String NUM_ITERATIONS_SHORT_NAME="n";
    protected static final String NUM_ITERATIONS_LONG_NAME="numIterations";
    protected static final String OUTPUT_DIR_SHORT_NAME="od";
    protected static final String OUTPUT_DIR_LONG_NAME="outputDir";
    protected static final String GATK_SEG_FILE_TAG = "cnv";
    protected static final String CGA_ACS_SEG_FILE_TAG = "acs";
    protected static final String BALANCED_SEG_FILE_TAG = "cnb_called";
    protected static final String TITAN_TN_FILE_TAG = "titan.tn";
    protected static final String TITAN_HET_FILE_TAG = "titan.het";


    @Argument(
            doc = "Input file for tumor-sample ref/alt read counts at normal-sample heterozygous-SNP sites (output of GetHetCoverage tool).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpCountsFile;

    @Argument(
            doc = "File of segmented regions of the genome, produced by AllelicCNV.",
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            fullName =  ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            optional = false
    )
    protected String segmentsFile;

    @Argument(
            doc = "(Advanced) Minimum rho (CCF*purity) to allow in calls.  Decreasing this value can yield false positives for CN LoH and false negatives for balanced segments.",
            shortName = CallAllelicSplits.RHO_SHORT_NAME,
            fullName =  CallAllelicSplits.RHO_LONG_NAME,
            optional = true
    )
    protected double rhoThreshold = CNLOHCaller.RHO_THRESHOLD_DEFAULT;

    @Argument(
            doc = "(Advanced) Number of iterations to perform for EM step in calling..",
            shortName = CallAllelicSplits.NUM_ITERATIONS_SHORT_NAME,
            fullName =  CallAllelicSplits.NUM_ITERATIONS_LONG_NAME,
            optional = true
    )
    protected int numIterations = DEFAULT_NUM_ITERATIONS;

    @Argument(
            doc = "Output directory to save all files.  This will be created if it does not already exist.",
            shortName = CallAllelicSplits.OUTPUT_DIR_SHORT_NAME,
            fullName =  CallAllelicSplits.OUTPUT_DIR_LONG_NAME,
            optional = false
    )
    protected File outputDir;

    @Argument(
            doc = "Input file for tumor-sample tangent-normalized target coverages (.tn.tsv output of GATK CNV tool).",
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File targetCoveragesFile;


    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        final String originalLogLevel =
                (ctx.getLocalProperty("logLevel") != null) ? ctx.getLocalProperty("logLevel") : "INFO";
        ctx.setLogLevel("WARN");

        final CNLOHCaller cnlohCaller = new CNLOHCaller();
        cnlohCaller.setRhoThreshold(rhoThreshold);

        final List<ACNVModeledSegment> segments = SegmentUtils.readACNVModeledSegmentFile(new File(segmentsFile));
        String sampleName = determineSampleName(new File(segmentsFile));

        // Create the outputdir
        try {
            FileUtils.forceMkdir(outputDir);
        } catch (final IOException ioe) {
            throw new UserException("Cannot create " + outputDir +".  Does a file (not directory) exist with the same name?  Do you have access to create?", ioe);
        }

        final Genome genome = new Genome(targetCoveragesFile, snpCountsFile);

        // Make the calls
        logger.info("Making the balanced-segment (and CNLoH) calls...");
        final List<AllelicSplitCall> calls = cnlohCaller.makeCalls(segments, numIterations, ctx);

        // Write updated ACNV file with calls
        logger.info("Writing updated ACNV file with calls ...");
        final File finalACNVModeledSegmentsFile = new File(outputDir, getSegmentsBaseFilename() + "." + BALANCED_SEG_FILE_TAG + ".seg");
        SegmentUtils.writeCnLoHACNVModeledSegmentFile(finalACNVModeledSegmentsFile, calls, genome);

        // write file for GATK CNV formatted seg file
        logger.info("Writing file with same output format as GATK CNV...");
        final File finalModeledSegmentsFileAsGatkCNV = new File(outputDir, getSegmentsBaseFilename() + "." + GATK_SEG_FILE_TAG + ".seg");
        SegmentUtils.writeModeledSegmentFile(finalModeledSegmentsFileAsGatkCNV,
                ACNVModeledSegmentConversionUtils.convertACNVModeledSegmentsToModeledSegments(segments, genome), sampleName);

        //write file for ACS-compatible output to help Broad CGA
        logger.info("Writing file with same output format as Broad CGA Allelic CapSeg ...");
        final File finalACSModeledSegmentsFile = new File(outputDir, getSegmentsBaseFilename() + "." + CGA_ACS_SEG_FILE_TAG + ".seg");
        ACSModeledSegmentUtils.writeACNVModeledSegmentFileAsAllelicCapSegFile(finalACSModeledSegmentsFile, calls, genome);

        //write files for TITAN-compatible output to help Broad CGA
        logger.info("Writing het file with input format for TITAN ...");
        final File finalTitanHetFile = new File(outputDir, getSegmentsBaseFilename() + "." + TITAN_HET_FILE_TAG + ".tsv");
        TitanFileConverter.convertHetPulldownToTitanHetFile(snpCountsFile, finalTitanHetFile);

        logger.info("Writing tangent-normalized target CR estimates with input format for TITAN ...");
        final File finalTitanTNFile = new File(outputDir, getSegmentsBaseFilename() + "." + TITAN_TN_FILE_TAG + ".tsv");
        TitanFileConverter.convertCRToTitanCovFile(targetCoveragesFile, finalTitanTNFile);

        ctx.setLogLevel(originalLogLevel);
        logger.info("SUCCESS: CNLoH and splits called for sample " + sampleName + ".");
    }

    private String determineSampleName(final File segmentsFile) {
        final List<String> sampleNames = SegmentUtils.readSampleNamesFromSegmentFile(segmentsFile);
        if (sampleNames.size() != 1) {
            throw new UserException.BadInput("This tool only supports exactly one sample in the input.  More than one found in " + segmentsFile.getAbsolutePath());
        }
        return sampleNames.get(0);
    }

    private String getSegmentsBaseFilename() {
        return FilenameUtils.removeExtension(new File(segmentsFile).getAbsoluteFile().getName());
    }


}
