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
import org.broadinstitute.hellbender.tools.exome.cnlohcaller.CNLOHCall;
import org.broadinstitute.hellbender.tools.exome.cnlohcaller.CNLOHCaller;
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

        summary = "NOTE: This tool only works with diploid organisms.",
        oneLineSummary = "",
        programGroup = CopyNumberProgramGroup.class
)
public class CallCNLoHAndSplits extends SparkCommandLineProgram {

    static final long serialVersionUID = 42123132L;

    protected final static double DEFAULT_RHO_THRESHOLD=0.2;
    protected final static String RHO_SHORT_NAME="r";
    protected final static String RHO_LONG_NAME="rhoThreshold";
    protected final static int DEFAULT_NUM_ITERATIONS = 10;
    protected final static String NUM_ITERATIONS_SHORT_NAME="n";
    protected final static String NUM_ITERATIONS_LONG_NAME="numIterations";
    protected final static String OUTPUT_DIR_SHORT_NAME="od";
    protected final static String OUTPUT_DIR_LONG_NAME="outputDir";
    protected final static String GATK_SEG_FILE_TAG = "cnv";
    protected final static String CGA_ACS_SEG_FILE_TAG = "acs";
    protected final static String CNLOH_BALANCED_SEG_FILE_TAG = "cnb_called";
    protected final static String TITAN_TN_FILE_TAG = "titan.tn";
    protected final static String TITAN_HET_FILE_TAG = "titan.het";


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
            doc = "(Advanced) Minimum rho (CCF*purity) to allow in calls.  Decreasing this value can yield false alarms for CN LoH and false negatives for balanced segments.",
            shortName = CallCNLoHAndSplits.RHO_SHORT_NAME,
            fullName =  CallCNLoHAndSplits.RHO_LONG_NAME,
            optional = true
    )
    protected double rhoThreshold = DEFAULT_RHO_THRESHOLD;

    @Argument(
            doc = "(Advanced) Number of iterations to perform for EM step in calling..",
            shortName = CallCNLoHAndSplits.NUM_ITERATIONS_SHORT_NAME,
            fullName =  CallCNLoHAndSplits.NUM_ITERATIONS_LONG_NAME,
            optional = true
    )
    protected int numIterations = DEFAULT_NUM_ITERATIONS;

    @Argument(
            doc = "Output directory to save all files.  This will be created if it does not already exist.",
            shortName = CallCNLoHAndSplits.OUTPUT_DIR_SHORT_NAME,
            fullName =  CallCNLoHAndSplits.OUTPUT_DIR_LONG_NAME,
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
    protected void runPipeline(JavaSparkContext ctx) {
        final CNLOHCaller cnlohCaller = new CNLOHCaller();
        cnlohCaller.setRhoThreshold(rhoThreshold);

        final List<ACNVModeledSegment> segs = SegmentUtils.readACNVModeledSegmentFile(new File(segmentsFile));
        final List<String> sampleNames = SegmentUtils.readSampleNamesFromSegmentFile(new File(segmentsFile));
        if (sampleNames.size() != 1) {
            throw new UserException("This tool only supports exactly one sample in the input.");
        }
        String sampleName = sampleNames.get(0);

        // Create the outputdir
        try {
            FileUtils.forceMkdir(outputDir);
        } catch (final IOException ioe) {
            throw new UserException("Cannot create " + outputDir +".  Does a file (not directory) exist with the same name?  Do you have access to create?", ioe);
        }

        // Create a Genome
        final Genome genome = new Genome(targetCoveragesFile, snpCountsFile, sampleName);

        // Make the calls
        logger.info("Making the CNLoH and balanced-segment calls...");
        final List<CNLOHCall> calls = cnlohCaller.makeCalls(segs, numIterations, ctx);

        // write file for GATK CNV formatted seg file
        logger.info("Writing file with same output format as GATK CNV...");
        final File finalModeledSegmentsFileAsGatkCNV = new File(outputDir + "/" + getSegmentsBasefilename() + "." + GATK_SEG_FILE_TAG + ".seg");
        SegmentUtils.writeModeledSegmentFile(finalModeledSegmentsFileAsGatkCNV,
                ACNVModeledSegmentConversionUtils.convertACNVModeledSegmentsToModeledSegments(segs, genome), sampleName);

        //write file for ACS-compatible output to help Broad CGA
        logger.info("Writing file with same output format as Broad CGA Allelic CapSeg ...");
        final File finalACSModeledSegmentsFile = new File(outputDir + "/" + getSegmentsBasefilename() + "." + CGA_ACS_SEG_FILE_TAG + ".seg");
        ACSModeledSegmentUtils.writeACNVModeledSegmentFileAsAllelicCapSegFile(finalACSModeledSegmentsFile, calls, genome);

        //write files for TITAN-compatible output to help Broad CGA
        logger.info("Writing het file with input format for TITAN ...");
        final File finalTitanHetFile = new File(outputDir + "/" + getSegmentsBasefilename() + "." + TITAN_HET_FILE_TAG + ".tsv");
        TitanFileConverter.convertHetPulldownToTitanHetFile(snpCountsFile, finalTitanHetFile);

        logger.info("Writing tangent-normalized target CR estimates with input format for TITAN ...");
        final File finalTitanTNFile = new File(outputDir + "/" + getSegmentsBasefilename() + "." + TITAN_TN_FILE_TAG + ".tsv");
        TitanFileConverter.convertCRToTitanCovFile(targetCoveragesFile, finalTitanTNFile);

        // Write updated ACNV file with calls
        logger.info("Writing updated ACNV file with calls ...");
        final File finalACNVModeledSegmentsFile = new File(outputDir + "/" + getSegmentsBasefilename() + "." + CNLOH_BALANCED_SEG_FILE_TAG + ".seg");
        SegmentUtils.writeCnLoHACNVModeledSegmentFile(finalACNVModeledSegmentsFile, calls, genome);
    }

    private String getSegmentsBasefilename() {
        return FilenameUtils.removeExtension(new File(segmentsFile).getAbsoluteFile().getName());
    }


}
