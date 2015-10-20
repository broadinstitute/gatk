package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;

/**
 * Outputs reference/alternate read counts for a tumor sample at heterozygous SNP sites present in a normal sample.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Output ref/alt counts for tumor sample at heterozygous SNPs in normal sample.",
        oneLineSummary = "Output ref/alt counts for tumor sample at heterozygous SNPs in normal sample."
)
public final class GetHetCoverage extends CommandLineProgram {
    protected static final String NORMAL_BAM_FILE_FULL_NAME = "normal";
    protected static final String NORMAL_BAM_FILE_SHORT_NAME = "N";

    protected static final String TUMOR_BAM_FILE_FULL_NAME = "tumor";
    protected static final String TUMOR_BAM_FILE_SHORT_NAME = "T";

    protected static final String SNP_FILE_FULL_NAME = "snpIntervals";
    protected static final String SNP_FILE_SHORT_NAME = "snp";

    @ArgumentCollection
    protected static final ReferenceInputArgumentCollection REF_ARGUMENTS =
            new RequiredReferenceInputArgumentCollection();

    protected static final String NORMAL_HET_REF_ALT_COUNTS_FILE_FULL_NAME = "normalOut";
    protected static final String NORMAL_HET_REF_ALT_COUNTS_FILE_SHORT_NAME = "NOut";

    protected static final String TUMOR_HET_REF_ALT_COUNTS_FILE_FULL_NAME = "tumorOut";
    protected static final String TUMOR_HET_REF_ALT_COUNTS_FILE_SHORT_NAME = "TOut";

    protected static final String HET_ALLELE_FRACTION_FULL_NAME = "hetAlleleFraction";
    protected static final String HET_ALLELE_FRACTION_SHORT_NAME = "AF";

    protected static final String PVALUE_THRESHOLD_FULL_NAME = "pvalueThreshold";
    protected static final String PVALUE_THRESHOLD_SHORT_NAME = "p";

    @Argument(
            doc = "BAM file for normal sample.",
            fullName = NORMAL_BAM_FILE_FULL_NAME,
            shortName = NORMAL_BAM_FILE_SHORT_NAME,
            optional = false
    )
    protected File normalBAMFile;

    @Argument(
            doc = "BAM file for tumor sample.",
            fullName = TUMOR_BAM_FILE_FULL_NAME,
            shortName = TUMOR_BAM_FILE_SHORT_NAME,
            optional = false
    )
    protected File tumorBAMFile;

    @Argument(
            doc = "Interval-list file of common SNPs.",
            fullName = SNP_FILE_FULL_NAME,
            shortName = SNP_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpFile;

    @Argument(
            doc = "Output file for normal-sample ref/alt read counts (at heterozygous SNPs).",
            fullName = NORMAL_HET_REF_ALT_COUNTS_FILE_FULL_NAME,
            shortName = NORMAL_HET_REF_ALT_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File normalHetOutputFile;

    @Argument(
            doc = "Output file for tumor-sample ref/alt read counts (at heterozygous SNPs in normal sample).",
            fullName = TUMOR_HET_REF_ALT_COUNTS_FILE_FULL_NAME,
            shortName = TUMOR_HET_REF_ALT_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File tumorHetOutputFile;

    @Argument(
            doc = "Heterozygous allele fraction.",
            fullName = HET_ALLELE_FRACTION_FULL_NAME,
            shortName = HET_ALLELE_FRACTION_SHORT_NAME,
            optional = false
    )
    protected double hetAlleleFraction = 0.5;

    @Argument(
            doc = "p-value threshold for binomial test for heterozygous SNPs in normal sample.",
            fullName = PVALUE_THRESHOLD_FULL_NAME,
            shortName = PVALUE_THRESHOLD_SHORT_NAME,
            optional = false
    )
    protected double pvalThreshold = 0.05;

    @Override
    protected Object doWork() {
        if (pvalThreshold < 0 || pvalThreshold > 1) {
            throw new UserException.BadArgumentValue(PVALUE_THRESHOLD_FULL_NAME,
                    Double.toString(pvalThreshold),
                    "p-value threshold should be in the [0, 1] range.");
        }

        final HetPulldownCalculator hetPulldown = new HetPulldownCalculator(REF_ARGUMENTS.getReferenceFile(), snpFile);

        logger.info("Getting normal het pulldown...");
        final Pulldown normalHetPulldown = hetPulldown.getNormal(normalBAMFile, hetAlleleFraction, pvalThreshold);
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