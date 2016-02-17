package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.utils.segmenter.RCBSSegmenter;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Segment genomic data into regions of constant copy-ratio.  Only supports one sample input.",
        oneLineSummary = "Segment genomic data into regions of constant copy-ratio",
        programGroup = CopyNumberProgramGroup.class
)
public final class PerformSegmentation extends CommandLineProgram {

    public static final String TARGET_WEIGHT_FILE_LONG_NAME= "targetWeights";
    public static final String TARGET_WEIGHT_FILE_SHORT_NAME = "tw";

    public final static String ALPHA_LONG_NAME="alpha";
    public final static String ALPHA_SHORT_NAME="alpha";

    public final static String NPERM_LONG_NAME="nperm";
    public final static String NPERM_SHORT_NAME="nperm";

    public final static String PMETHOD_LONG_NAME="pmethod";
    public final static String PMETHOD_SHORT_NAME="pmethod";

    public final static String MINWIDTH_LONG_NAME="minWidth";
    public final static String MINWIDTH_SHORT_NAME="minWidth";

    public final static String KMAX_LONG_NAME="kmax";
    public final static String KMAX_SHORT_NAME="kmax";

    public final static String NMIN_LONG_NAME="nmin";
    public final static String NMIN_SHORT_NAME="nmin";

    public final static String ETA_LONG_NAME="eta";
    public final static String ETA_SHORT_NAME="eta";

    public final static String TRIM_LONG_NAME="trim";
    public final static String TRIM_SHORT_NAME="trim";

    public final static String UNDOSPLITS_LONG_NAME="undoSplits";
    public final static String UNDOSPLITS_SHORT_NAME="undoSplits";

    public final static String UNDOPRUNE_LONG_NAME="undoPrune";
    public final static String UNDOPRUNE_SHORT_NAME="undoPrune";

    public final static String UNDOSD_LONG_NAME="undoSD";
    public final static String UNDOSD_SHORT_NAME="undoSD";

    @Argument(
            doc = "Genomic targets file",
            shortName = ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME,
            fullName =  ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME,
            optional = false
    )
    protected String tangentFile;

    @Argument(
            doc = "Full path to the outputted segment file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false
    )
    protected String outFile;

    @Argument(
            doc = "If input data has had a log2 transform applied",
            shortName = ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.LOG2_LONG_NAME,
            optional = true
    )
    protected Boolean log = false;

    @Argument(
            doc = "File with target weights.  This is the 1/var(post-projected targets for each normal).  " +
                    "Listed one value per line in plain text.  Values of zero or less, Nan, Inf, and -Inf are not " +
                    "acceptable.  Must have the same number of values as there are in the tangentFile.",
            shortName = TARGET_WEIGHT_FILE_SHORT_NAME,
            fullName = TARGET_WEIGHT_FILE_LONG_NAME,
            optional = true
    )
    protected File weightFile = null;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            shortName = ALPHA_SHORT_NAME,
            fullName = ALPHA_LONG_NAME,
            optional = true
    )
    protected Double alpha = 0.01;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            shortName = NPERM_SHORT_NAME,
            fullName = NPERM_LONG_NAME,
            optional = true
    )
    protected Integer nperm = 10000;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            shortName = PMETHOD_SHORT_NAME,
            fullName = PMETHOD_LONG_NAME,
            optional = true
    )
    protected RCBSSegmenter.PMethod pmethod = RCBSSegmenter.PMethod.HYBRID;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            shortName = MINWIDTH_SHORT_NAME,
            fullName = MINWIDTH_LONG_NAME,
            optional = true
    )
    protected Integer minWidth = 2;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            shortName = KMAX_SHORT_NAME,
            fullName = KMAX_LONG_NAME,
            optional = true
    )
    protected Integer kmax = 25;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            shortName = NMIN_SHORT_NAME,
            fullName = NMIN_LONG_NAME,
            optional = true
    )
    protected Integer nmin = 200;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            shortName = ETA_SHORT_NAME,
            fullName = ETA_LONG_NAME,
            optional = true
    )
    protected Double eta = 0.05;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            shortName = TRIM_SHORT_NAME,
            fullName = TRIM_LONG_NAME,
            optional = true
    )
    protected Double trim = 0.025;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            shortName = UNDOSPLITS_SHORT_NAME,
            fullName = UNDOSPLITS_LONG_NAME,
            optional = true
    )
    protected RCBSSegmenter.UndoSplits undoSplits = RCBSSegmenter.UndoSplits.NONE;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            shortName = UNDOPRUNE_SHORT_NAME,
            fullName = UNDOPRUNE_LONG_NAME,
            optional = true
    )
    protected Double undoPrune = 0.05;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            shortName = UNDOSD_SHORT_NAME,
            fullName = UNDOSD_LONG_NAME,
            optional = true
    )
    protected Integer undoSD = 3;

    @Override
    protected Object doWork() {

        final String sampleName = TargetCoverageUtils.getSampleNameForCLIsFromTargetCoverageFile(new File(tangentFile));

        applySegmentation(sampleName, tangentFile, outFile);
        return "Success";
    }

    private void applySegmentation(final String sampleName, final String tangentFile, final String outFile) {
        RCBSSegmenter.writeSegmentFile(sampleName, tangentFile, outFile, log, weightFile, alpha, nperm, pmethod,
                minWidth, kmax, nmin, eta, trim, undoSplits, undoPrune, undoSD);
    }
}
