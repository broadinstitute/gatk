package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;

/**
 * Outputs reference/alternate read counts for a tumor sample at heterozygous SNP sites present in a normal sample.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Output ref/alt counts for tumor sample at heterozygous SNPs in normal sample.",
        oneLineSummary = "Output ref/alt counts for tumor sample at heterozygous SNPs in normal sample.",
        programGroup = CopyNumberProgramGroup.class
)
public final class GetHetCoverage extends CommandLineProgram {

    @ArgumentCollection
    protected static final ReferenceInputArgumentCollection REF_ARGUMENTS =
            new RequiredReferenceInputArgumentCollection();

    protected static final String PVALUE_THRESHOLD_FULL_NAME = "pvalueThreshold";
    protected static final String PVALUE_THRESHOLD_SHORT_NAME = "p";

    protected static final String MINIMUM_MAPPING_QUALITY_SHORT_NAME = "minMQ";
    protected static final String MINIMUM_MAPPING_QUALITY_FULL_NAME = "minimumMappingQuality";

    protected static final String MINIMUM_BASE_QUALITY_SHORT_NAME = "minBQ";
    protected static final String MINIMUM_BASE_QUALITY_FULL_NAME = "minimumBaseQuality";

    @Argument(
            doc = "BAM file for normal sample.",
            fullName = ExomeStandardArgumentDefinitions.NORMAL_BAM_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.NORMAL_BAM_FILE_SHORT_NAME,
            optional = false
    )
    protected File normalBAMFile;

    @Argument(
            doc = "BAM file for tumor sample.",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_SHORT_NAME,
            optional = false
    )
    protected File tumorBAMFile;

    @Argument(
            doc = "Interval-list file of common SNPs.",
            fullName = ExomeStandardArgumentDefinitions.SNP_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpFile;

    @Argument(
            doc = "Output file for normal-sample ref/alt read counts (at heterozygous SNPs).",
            fullName = ExomeStandardArgumentDefinitions.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.NORMAL_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File normalHetOutputFile;

    @Argument(
            doc = "Output file for tumor-sample ref/alt read counts (at heterozygous SNPs in normal sample).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File tumorHetOutputFile;

    @Argument(
            doc = "p-value threshold for binomial test for heterozygous SNPs in normal sample.",
            fullName = PVALUE_THRESHOLD_FULL_NAME,
            shortName = PVALUE_THRESHOLD_SHORT_NAME,
            optional = false
    )
    protected double pvalThreshold = 0.05;

    @Argument(
            doc = "Minimum mapping quality; reads with lower quality will be filtered out of pileup.",
            shortName = MINIMUM_MAPPING_QUALITY_SHORT_NAME,
            fullName  = MINIMUM_MAPPING_QUALITY_FULL_NAME,
            optional = true
    )
    protected int minimumMappingQuality = 30;

    @Argument(
            doc = "Minimum base quality; base calls with lower quality will be filtered out of pileup.",
            shortName = MINIMUM_BASE_QUALITY_SHORT_NAME,
            fullName = MINIMUM_BASE_QUALITY_FULL_NAME,
            optional = true
    )
    protected int minimumBaseQuality = 20;

    @Override
    protected Object doWork() {
        if (pvalThreshold < 0 || pvalThreshold > 1) {
            throw new UserException.BadArgumentValue(PVALUE_THRESHOLD_FULL_NAME,
                    Double.toString(pvalThreshold),
                    "p-value threshold should be in the [0, 1] range.");
        }

        final HetPulldownCalculator hetPulldown = new HetPulldownCalculator(REF_ARGUMENTS.getReferenceFile(), snpFile,
                minimumMappingQuality, minimumBaseQuality);

        logger.info("Getting normal het pulldown...");
        final Pulldown normalHetPulldown = hetPulldown.getNormal(normalBAMFile, pvalThreshold);
        normalHetPulldown.write(normalHetOutputFile);
        logger.info("Normal het pulldown written to " + normalHetOutputFile.toString());

        final IntervalList normalHetIntervals = normalHetPulldown.getIntervals();

        logger.info("Getting tumor het pulldown...");
        final Pulldown tumorHetPulldown = hetPulldown.getTumor(tumorBAMFile, normalHetIntervals);
        tumorHetPulldown.write(tumorHetOutputFile);
        logger.info("Tumor het pulldown written to " + tumorHetOutputFile.toString());

        return "SUCCESS";
    }
}