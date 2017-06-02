package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountTableColumn.AllelicCountTableVerbosity;
import org.broadinstitute.hellbender.tools.exome.pulldown.HetPulldownCalculator;
import org.broadinstitute.hellbender.tools.exome.pulldown.Pulldown;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import java.io.File;

/**
 * Outputs reference/alternate read counts at heterozygous SNP sites present in a normal sample,
 * and, if specified, at the same sites in a tumor sample.
 *
 * <p>The non-Bayesian version of {@link GetBayesianHetCoverage}. See GetBayesianHetCoverage for modes.</p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Output ref/alt counts at heterozygous SNPs in normal sample (and at same sites in tumor sample, if specified).",
        oneLineSummary = "Output ref/alt counts at heterozygous SNPs in normal sample (and at same sites in tumor sample, if specified)",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class GetHetCoverage extends CommandLineProgram {

    protected static final String PVALUE_THRESHOLD_FULL_NAME = "pvalueThreshold";
    protected static final String PVALUE_THRESHOLD_SHORT_NAME = "p";

    protected static final String MINIMUM_MAPPING_QUALITY_SHORT_NAME = "minMQ";
    protected static final String MINIMUM_MAPPING_QUALITY_FULL_NAME = "minimumMappingQuality";

    protected static final String MINIMUM_BASE_QUALITY_SHORT_NAME = "minBQ";
    protected static final String MINIMUM_BASE_QUALITY_FULL_NAME = "minimumBaseQuality";

    protected static final String MINIMUM_READ_COUNT_SHORT_NAME = "minRC";
    protected static final String MINIMUM_READ_COUNT_FULL_NAME = "minimumReadCount";

    @ArgumentCollection
    protected static final ReferenceInputArgumentCollection REFERENCE_ARGUMENTS =
            new RequiredReferenceInputArgumentCollection();

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
            optional = true
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
            optional = true
    )
    protected File tumorHetOutputFile;

    @Argument(
            doc = "p-value threshold for binomial test for heterozygous SNPs in normal sample (must be in [0, 1]).",
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

    @Argument(
            doc = "Validation stringency for all BAM files read by this program.  Setting stringency to SILENT " +
            "can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) " +
            "do not otherwise need to be decoded.",
            common=true)
    protected ValidationStringency VALIDATION_STRINGENCY = ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY;

    @Argument(
            doc = "Minimum raw number of reads that must be present at a site for it to be considered for het compatibility.",
            shortName = MINIMUM_READ_COUNT_SHORT_NAME,
            fullName = MINIMUM_READ_COUNT_FULL_NAME,
            optional = true
    )
    protected int minimumRawReads = 15;

    @Override
    protected Object doWork() {
        //if tumor arguments are missing, throw exception (and do not even get normal pulldown)
        final boolean doTumorPulldown;
        if (tumorHetOutputFile != null && tumorBAMFile != null) {
            doTumorPulldown = true;
        } else if ((tumorHetOutputFile == null) != (tumorBAMFile == null)) {
            throw new CommandLineException("Must specify both BAM and output files for tumor pulldown.");
        } else {
            doTumorPulldown = false;
        }

        final HetPulldownCalculator hetPulldown = new HetPulldownCalculator(REFERENCE_ARGUMENTS.getReferenceFile(),
                snpFile, minimumMappingQuality, minimumBaseQuality, VALIDATION_STRINGENCY);

        logger.info("Getting normal het pulldown...");
        final Pulldown normalHetPulldown = hetPulldown.getNormal(normalBAMFile, pvalThreshold, minimumRawReads);
        normalHetPulldown.write(normalHetOutputFile, AllelicCountTableVerbosity.BASIC);
        logger.info("Normal het pulldown written to " + normalHetOutputFile.toString());

        if (doTumorPulldown) {
            final IntervalList normalHetIntervals = normalHetPulldown.getIntervals();

            logger.info("Getting tumor het pulldown...");
            final Pulldown tumorHetPulldown = hetPulldown.getTumor(tumorBAMFile, normalHetIntervals, minimumRawReads);
            tumorHetPulldown.write(tumorHetOutputFile, AllelicCountTableVerbosity.BASIC);
            logger.info("Tumor het pulldown written to " + tumorHetOutputFile.toString());
        }

        return "SUCCESS";
    }
}