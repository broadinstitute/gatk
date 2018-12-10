package org.broadinstitute.hellbender.tools.walkers;

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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

@CommandLineProgramProperties(
        summary = "Prints reads from the provided file(s) with corresponding reference bases (if a reference is provided) to the specified output file (or STDOUT if none specified)",
        oneLineSummary = "Print reads with reference context",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public class PalindromeArtifactReadWalker extends ReadWalker {
    public int COUNT = 0;
    public int TOTAL_COUNT = 0;
    public long COUNT_BASES = 0;
    public int COUNT_ARTIFACT_UMIS_EQUAL = 0;
    public int COUNT_ARTIFACT_UMIS_OFF_BY_ONE = 0;
    public int TOTAL_ARTIFACT_UMIS_EQUAL = 0;
    public int TOTAL_ARTIFACT_UMIS_OFF_BY_ONE = 0;
    public static final String CLIPPED_KEY = "XC";
    public static final String UMI_KEY = "RX";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File OUTPUT_FILE = null;

    private PrintStream outputStream = null;

    @Override
    public void onTraversalStart() {
        try {
            outputStream = OUTPUT_FILE != null ? new PrintStream(OUTPUT_FILE) : System.out;
        }
        catch ( FileNotFoundException e ) {
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
        if(read.hasAttribute(CLIPPED_KEY)) {
            COUNT++;
            COUNT_BASES = COUNT_BASES + read.getAttributeAsInteger(CLIPPED_KEY);
        }
        TOTAL_COUNT++;
        if (read.hasAttribute(UMI_KEY)){
            String[] umi = read.getAttributeAsString(UMI_KEY).split("-");
            if(umi[0].equals(umi[1])) {
                if (read.hasAttribute(CLIPPED_KEY)) {
                    COUNT_ARTIFACT_UMIS_EQUAL++;
                }
                TOTAL_ARTIFACT_UMIS_EQUAL++;
            } else if (StringUtil.hammingDistance(umi[0], umi[1]) == 1) {
                if (read.hasAttribute(CLIPPED_KEY)) {
                    COUNT_ARTIFACT_UMIS_OFF_BY_ONE++;
                }
                TOTAL_ARTIFACT_UMIS_OFF_BY_ONE++;
            }
        }
    }

    @Override
    public void closeTool(){
        outputStream.println("Transformed Reads: " + COUNT);
        outputStream.println("Transformed Bases: " + COUNT_BASES);
        outputStream.println("Total Reads: " + TOTAL_COUNT);
        outputStream.println("Total matching UMIs: " + TOTAL_ARTIFACT_UMIS_EQUAL);
        outputStream.println("Total off by one UMIs: " + TOTAL_ARTIFACT_UMIS_OFF_BY_ONE);
        outputStream.println("Artifact matching UMIs: " + COUNT_ARTIFACT_UMIS_EQUAL);
        outputStream.println("Artifact off by one UMIs: " + COUNT_ARTIFACT_UMIS_OFF_BY_ONE);
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
