package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.utils.segmenter.RCBSSegmenter;

import java.io.File;
/**
 * Tool groups contiguous targets with the same copy-ratio per sample.
 *
 * <p>
 *     The tool requires R and R components. Install R from https://www.r-project.org/.
 *     For the components, download the script at
 *     https://github.com/broadinstitute/gatk-protected/blob/master/scripts/install_R_packages.R and
 *     use it with the following comm
 *     <pre>
 *         Rscript install_R_packages.R
 *     </pre>
 * </p>
 *
 * <h3>Examples</h3>
 * <p>
 *     The commands encompass empirically determined parameters for TCGA project data.
 *     You may obtain better results with different parameters.
 * </p>
 *
 * <p>For whole exome sequencing (WES) data: </p>
 *
 * <pre>
 * java -Xmx4g -jar $gatk_jar PerformSegmentation \
 *   --tangentNormalized tn_coverage.tn.tsv \
 *   --log2Input true \
 *   --undoSD 2 \
 *   --output entity_id.seg
 * </pre>
 *
 * <p>For whole genome sequencing (WGS) data: </p>
 *
 * <pre>
 * java -Xmx4g -jar $gatk_jar PerformSegmentation \
 *   --tangentNormalized tn_coverage.tn.tsv \
 *   --log2Input true \
 *   --undoSplits SDUNDO \
 *   --output entity_id.seg
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Segment genomic data into regions of constant copy-ratio.  Only supports one sample input.",
        oneLineSummary = "Segment genomic data into regions of constant copy-ratio",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class PerformSegmentation extends CommandLineProgram {

    public static final String TARGET_WEIGHT_FILE_LONG_NAME = "targetWeights";
    public static final String TARGET_WEIGHT_FILE_SHORT_NAME = "tw";

    public static final String ALPHA_LONG_NAME = "alpha";
    public static final String ALPHA_SHORT_NAME = "alpha";

    public static final String NPERM_LONG_NAME = "nperm";
    public static final String NPERM_SHORT_NAME = "nperm";

    public static final String PMETHOD_LONG_NAME = "pmethod";
    public static final String PMETHOD_SHORT_NAME = "pmethod";

    public static final String MINWIDTH_LONG_NAME = "minWidth";
    public static final String MINWIDTH_SHORT_NAME = "minWidth";

    public static final String KMAX_LONG_NAME = "kmax";
    public static final String KMAX_SHORT_NAME = "kmax";

    public static final String NMIN_LONG_NAME = "nmin";
    public static final String NMIN_SHORT_NAME = "nmin";

    public static final String ETA_LONG_NAME = "eta";
    public static final String ETA_SHORT_NAME = "eta";

    public static final String TRIM_LONG_NAME = "trim";
    public static final String TRIM_SHORT_NAME = "trim";

    public static final String UNDOSPLITS_LONG_NAME = "undoSplits";
    public static final String UNDOSPLITS_SHORT_NAME = "undoSplits";

    public static final String UNDOPRUNE_LONG_NAME = "undoPrune";
    public static final String UNDOPRUNE_SHORT_NAME = "undoPrune";

    public static final String UNDOSD_LONG_NAME = "undoSD";
    public static final String UNDOSD_SHORT_NAME = "undoSD";

    @Argument(
            doc = "Tangent-normalized read counts file",
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected String tangentNormalizedCoverageFile;

    @Argument(
            doc = "Full path to the outputted segment file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false
    )
    protected String outFile;

    @Argument(
            doc = "If input data has had a log2 transform applied",
            fullName = ExomeStandardArgumentDefinitions.LOG2_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
            optional = true
    )
    protected Boolean log = false;

    @Argument(
            doc = "File with target weights.  This is the 1/var(post-projected targets for each normal).  " +
                    "Listed one value per line in plain text.  Values of zero or less, Nan, Inf, and -Inf are not " +
                    "acceptable.  Must have the same number of values as there are in the tangentNormalizedCoverageFile.",
            fullName = TARGET_WEIGHT_FILE_LONG_NAME,
            shortName = TARGET_WEIGHT_FILE_SHORT_NAME,
            optional = true
    )
    protected File weightFile = null;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            fullName = ALPHA_LONG_NAME,
            shortName = ALPHA_SHORT_NAME,
            optional = true
    )
    protected Double alpha = 0.01;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            fullName = NPERM_LONG_NAME,
            shortName = NPERM_SHORT_NAME,
            optional = true
    )
    protected Integer nperm = 10000;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            fullName = PMETHOD_LONG_NAME,
            shortName = PMETHOD_SHORT_NAME,
            optional = true
    )
    protected RCBSSegmenter.PMethod pmethod = RCBSSegmenter.PMethod.HYBRID;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            fullName = MINWIDTH_LONG_NAME,
            shortName = MINWIDTH_SHORT_NAME,
            optional = true
    )
    protected Integer minWidth = 2;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            fullName = KMAX_LONG_NAME,
            shortName = KMAX_SHORT_NAME,
            optional = true
    )
    protected Integer kmax = 25;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            fullName = NMIN_LONG_NAME,
            shortName = NMIN_SHORT_NAME,
            optional = true
    )
    protected Integer nmin = 200;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            fullName = ETA_LONG_NAME,
            shortName = ETA_SHORT_NAME,
            optional = true
    )
    protected Double eta = 0.05;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            fullName = TRIM_LONG_NAME,
            shortName = TRIM_SHORT_NAME,
            optional = true
    )
    protected Double trim = 0.025;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            fullName = UNDOSPLITS_LONG_NAME,
            shortName = UNDOSPLITS_SHORT_NAME,
            optional = true
    )
    protected RCBSSegmenter.UndoSplits undoSplits = RCBSSegmenter.UndoSplits.NONE;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            fullName = UNDOPRUNE_LONG_NAME,
            shortName = UNDOPRUNE_SHORT_NAME,
            optional = true
    )
    protected Double undoPrune = 0.05;

    @Argument(
            doc = "(Advanced) Please see https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf",
            fullName = UNDOSD_LONG_NAME,
            shortName = UNDOSD_SHORT_NAME,
            optional = true
    )
    protected Integer undoSD = 3;

    @Override
    protected Object doWork() {
        final String sampleName = ReadCountCollectionUtils.getSampleNameForCLIsFromReadCountsFile(new File(tangentNormalizedCoverageFile));
        applySegmentation(sampleName, tangentNormalizedCoverageFile, outFile);
        return "SUCCESS";
    }

    private void applySegmentation(final String sampleName, final String tangentFile, final String outFile) {
        RCBSSegmenter.writeSegmentFile(sampleName, tangentFile, outFile, log, weightFile, alpha, nperm, pmethod,
                minWidth, kmax, nmin, eta, trim, undoSplits, undoPrune, undoSD);
    }
}
