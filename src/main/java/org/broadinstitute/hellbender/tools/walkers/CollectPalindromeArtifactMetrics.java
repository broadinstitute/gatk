package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
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
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@ExperimentalFeature
public class CollectPalindromeArtifactMetrics extends ReadWalker {

    private final PalindromicArtifactMetric metric = new PalindromicArtifactMetric();
    public static final String CLIPPED_KEY = "XC";
    private static final String UMI_KEY = "RX";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File OUTPUT_FILE = null;

    @Argument(fullName = "metrics-file", shortName = "M", doc = "Metrics file, if provided", common = false, optional = true)
    private File METRICS_FILE = null;

    private PrintStream outputStream = null;

    @Override
    public boolean requiresReference() {
        return true;
    }

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
        if (read.failsVendorQualityCheck()) {
            return;
        }
        if (metric.SAMPLE == null) {
            metric.SAMPLE = ReadUtils.getSampleName(read, getHeaderForReads());
            logger.info("got sample name from first read: " + metric.SAMPLE);
        }

        final boolean isArtifact = read.hasAttribute(CLIPPED_KEY);
        if (isArtifact) {
            metric.ARTIFACT_READS++;
            metric.CLIPPED_BASES += read.getAttributeAsInteger(CLIPPED_KEY);
        }
        metric.TOTAL_READS++;
        if (read.hasAttribute(UMI_KEY)) {
            String[] umi = read.getAttributeAsString(UMI_KEY).split("-");
            if (isArtifact) {
                metric.COUNT_ARTIFACT_UMIS++;
            }
            metric.TOTAL_UMIS++;
            if (umi[0].equals(umi[1])) {
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
        public String SAMPLE;

        @MergeByAdding
        public Integer TOTAL_READS = 0;

        @MergeByAdding
        public Integer ARTIFACT_READS = 0;

        @NoMergingIsDerived
        public Double ARTIFACT_RATE = 0D;

        @MergeByAdding
        public Integer CLIPPED_BASES = 0;

        @MergeByAdding
        public Integer TOTAL_BASES = 0;

        @NoMergingIsDerived
        public Double CLIPPING_RATE = 0D;

        @MergeByAdding
        public Integer TOTAL_UMIS_EQUAL = 0;

        @MergeByAdding
        public Integer TOTAL_UMIS = 0;

        @MergeByAdding
        public Integer COUNT_ARTIFACT_UMIS = 0;

        @MergeByAdding
        public Integer ARTIFACT_UMIS_EQUAL = 0;

        @MergeByAdding
        public Double ARTIFACT_UMIS_EQUAL_RATE = 0D;

        @MergeByAdding
        public Double UMIS_EQUAL_RATE = 0D;

        @MergeByAdding
        public Integer TOTAL_UMIS_OFF_BY_ONE = 0;

        @MergeByAdding
        public Integer ARTIFACT_UMIS_OFF_BY_ONE = 0;

        @MergeByAdding
        public Double UMIS_OFF_BY_ONE_RATE = 0D;

        @MergeByAdding
        public Double ARTIFACT_UMIS_OFF_BY_ONE_RATE = 0D;

        @Override
        public void calculateDerivedFields() {
            super.calculateDerivedFields();

            ARTIFACT_RATE = ARTIFACT_READS / (double) TOTAL_READS;

            ARTIFACT_UMIS_EQUAL_RATE = ARTIFACT_UMIS_EQUAL / (double) COUNT_ARTIFACT_UMIS;
            UMIS_EQUAL_RATE = TOTAL_UMIS_EQUAL / (double) TOTAL_UMIS;
            UMIS_OFF_BY_ONE_RATE = TOTAL_UMIS_OFF_BY_ONE / (double) TOTAL_UMIS;

            ARTIFACT_UMIS_EQUAL_RATE = ARTIFACT_UMIS_EQUAL / (double) COUNT_ARTIFACT_UMIS;
            ARTIFACT_UMIS_OFF_BY_ONE_RATE = ARTIFACT_UMIS_OFF_BY_ONE / (double) COUNT_ARTIFACT_UMIS;

            CLIPPING_RATE = CLIPPED_BASES / (double) TOTAL_BASES;
        }
    }
}
