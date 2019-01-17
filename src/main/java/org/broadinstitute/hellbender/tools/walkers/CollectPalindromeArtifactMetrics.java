package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.transformers.PalindromeArtifactClipReadTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import picard.analysis.MergeableMetricBase;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

@CommandLineProgramProperties(
        summary = "Prints reads from the provided file(s) with corresponding reference bases (if a reference is provided) to the specified output file (or STDOUT if none specified)",
        oneLineSummary = "Print reads with reference context",
        programGroup = DiagnosticsAndQCProgramGroup.class,
        omitFromCommandLine = true
)
public class CollectPalindromeArtifactMetrics extends ReadWalker {

    private PalindromicArtifactMetric metric = new PalindromicArtifactMetric();
    private static final String CLIPPED_KEY = "XC";
    private static final String UMI_KEY = "RX";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File OUTPUT_FILE = null;

    @Argument(fullName = "metrics-file", shortName = "M", doc = "Metrics file, if provided", common = false, optional = true)
    private File METRICS_FILE = null;

    private PrintStream outputStream = null;
    private String sample;

    @Override
    public void onTraversalStart() {
        try {
            outputStream = OUTPUT_FILE != null ? new PrintStream(OUTPUT_FILE) : System.out;
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotReadInputFile(OUTPUT_FILE, e);
        }
    }

    @Override
    public ReadTransformer makePostReadFilterTransformer() {
        return super.makePostReadFilterTransformer().andThen(new PalindromeArtifactClipReadTransformer(new ReferenceFileSource(referenceArguments.getReferencePath()), Mutect2Engine.MIN_PALINDROME_SIZE));
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        if(read.failsVendorQualityCheck()){
            return;
        }
        if (metric.SAMPLE == null) {
            metric.SAMPLE = ReadUtils.getSampleName(read, getHeaderForReads());
            logger.info("got sample name from first read: ", metric.SAMPLE);
        }

        final boolean isArtifact = read.hasAttribute(CLIPPED_KEY);
        if (isArtifact) {
            metric.ARTIFACT_READS++;
            metric.CLIPPED_BASES += read.getAttributeAsInteger(CLIPPED_KEY);
        }
        metric.TOTAL_READS++;
        if (read.hasAttribute(UMI_KEY)){
            String[] umi = read.getAttributeAsString(UMI_KEY).split("-");
            if (isArtifact){
                metric.COUNT_ARTIFACT_UMIS++;
            }
            metric.TOTAL_UMIS++;
            if(umi[0].equals(umi[1])) {
                if (isArtifact) {
                    metric.ARTIFACT_UMIS_EQUAL++;
                }
                metric.TOTAL_UMIS_EQUAL++;
            } else if (StringUtil.hammingDistance(umi[0], umi[1]) == 1) {
                if (isArtifact) {
                    metric.ARTIFACT_UMIS_OFF_BY_ONE++;
                }
                metric.TOTAL_UMIS_OFF_BY_ONE++;
            }
        }
    }

    @Override
    public void closeTool() {

        final PalindromicArtifactMetric metric = new PalindromicArtifactMetric();
        metric.calculateDerivedFields();

        final MetricsFile<PalindromicArtifactMetric, ?> metricsFile = getMetricsFile();
        metricsFile.addMetric(metric);
        metricsFile.write(METRICS_FILE);

        outputStream.println("Transformed Reads: " + metric.ARTIFACT_READS);
        outputStream.println("Transformed Bases: " + metric.CLIPPED_BASES);
        outputStream.println("Total Reads: " + metric.TOTAL_READS);

        outputStream.println("Total matching UMIs: " + metric.TOTAL_UMIS_EQUAL);
        outputStream.println("Total off by one UMIs: " + metric.TOTAL_UMIS_OFF_BY_ONE);
        outputStream.println("Artifact matching UMIs: " + metric.ARTIFACT_UMIS_EQUAL);
        outputStream.println("Artifact off by one UMIs: " + metric.ARTIFACT_UMIS_OFF_BY_ONE);
        if (outputStream != null) {
            outputStream.close();
        }
    }

    static public class PalindromicArtifactMetric extends MergeableMetricBase {

        @MergeByAssertEquals
        String SAMPLE;

        @MergeByAdding
        Integer TOTAL_READS;

        @MergeByAdding
        Integer ARTIFACT_READS;

        @NoMergingIsDerived
        Double ARTIFACT_RATE;

        @MergeByAdding
        Integer CLIPPED_BASES;

        @MergeByAdding
        Integer TOTAL_BASES;

        @NoMergingIsDerived
        Double CLIPPING_RATE;

        @MergeByAdding
        Integer TOTAL_UMIS_EQUAL;

        @MergeByAdding
        Integer TOTAL_UMIS;

        @MergeByAdding
        Integer COUNT_ARTIFACT_UMIS;

        @MergeByAdding
        Integer ARTIFACT_UMIS_EQUAL;

        @MergeByAdding
        Double ARTIFACT_UMIS_EQUAL_RATE;

        @MergeByAdding
        Double UMIS_EQUAL_RATE;

        @MergeByAdding
        Integer TOTAL_UMIS_OFF_BY_ONE;

        @MergeByAdding
        Integer ARTIFACT_UMIS_OFF_BY_ONE;

        @MergeByAdding
        Double UMIS_OFF_BY_ONE_RATE;

        @MergeByAdding
        Double ARTIFACT_UMIS_OFF_BY_ONE_RATE;

        @Override
        public void calculateDerivedFields() {
            super.calculateDerivedFields();

            ARTIFACT_RATE = ARTIFACT_READS / (double) TOTAL_READS;

            ARTIFACT_UMIS_EQUAL_RATE = ARTIFACT_UMIS_EQUAL / (double) COUNT_ARTIFACT_UMIS;
            UMIS_EQUAL_RATE = TOTAL_UMIS_EQUAL / (double) TOTAL_UMIS;
            UMIS_OFF_BY_ONE_RATE = TOTAL_UMIS_OFF_BY_ONE / (double)TOTAL_UMIS ;

            ARTIFACT_UMIS_EQUAL_RATE = ARTIFACT_UMIS_EQUAL / (double) COUNT_ARTIFACT_UMIS;
            ARTIFACT_UMIS_OFF_BY_ONE_RATE = ARTIFACT_UMIS_OFF_BY_ONE / (double) COUNT_ARTIFACT_UMIS;

            CLIPPING_RATE =  CLIPPED_BASES / (double) TOTAL_BASES;
        }
    }
}
