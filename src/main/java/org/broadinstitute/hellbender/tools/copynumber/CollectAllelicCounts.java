package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollector;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import java.io.File;

/**
 * Collects reference/alternate allele counts at sites.  The alt count is defined as the total count minus the ref count,
 * and the alt nucleotide is defined as the non-ref base with the highest count, with ties broken by the order of the
 * bases in {@link AllelicCountCollector#BASES}.
 *
 * <h3>Example</h3>
 *
 * <pre>
 * java -Xmx4g -jar $gatk_jar CollectAllelicCounts \
 *   --input sample.bam \
 *   --reference ref_fasta.fa \
 *   --siteIntervals sites.interval_list \
 *   --output allelic_counts.tsv
 * </pre>
 *
 * <p>
 *     The --siteIntervals is a Picard-style interval list, e.g.:
 * </p>
 *
 * <pre>
 *     @HD	VN:1.4	SO:unsorted
 *     @SQ	SN:1	LN:16000	M5:8c0c38e352d8f3309eabe4845456f274
 *     @SQ	SN:2	LN:16000	M5:5f8388fe3fb34aa38375ae6cf5e45b89
 *      1	10736	10736	+	normal
 *      1	11522	11522	+	normal
 *      1	12098	12098	+	normal
 *      1	12444	12444	+	normal
 *      1	13059	13059	+	normal
 *      1	14630	14630	+	normal
 *      1	15204	15204	+	normal
 *      2	14689	14689	+	normal
 *      2	14982	14982	+	normal
 *      2	15110	15110	+	normal
 *      2	15629	15629	+	normal
 * </pre>
 *
 * <p>
 *     The resulting table counts the reference versus alternate allele and indicates the alleles. For example:
 * </p>
 *
 * <pre>
 *     CONTIG	POSITION	REF_COUNT	ALT_COUNT	REF_NUCLEOTIDE	ALT_NUCLEOTIDE
 *     1	10736	0	0	G	A
 *     1	11522	7	4	G	C
 *     1	12098	8	6	G	A
 *     1	12444	0	18	T	A
 *     1	13059	0	8	C	G
 *     1	14630	9	8	T	A
 *     1	15204	4	4	C	G
 *     2	14689	6	9	T	A
 *     2	14982	6	5	G	A
 *     2	15110	6	0	G	A
 *     2	15629	5	3	T	C
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Collects ref/alt counts at sites.",
        oneLineSummary = "Collects ref/alt counts at sites",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class CollectAllelicCounts extends CommandLineProgram {

    @ArgumentCollection
    protected final ReferenceInputArgumentCollection referenceArguments = new RequiredReferenceInputArgumentCollection();

    @Argument(
            doc = "Input BAM file.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            optional = false
    )
    protected File inputBAMFile;

    @Argument(
            doc = "Interval-list file of sites.",
            fullName = ExomeStandardArgumentDefinitions.SITES_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SITES_FILE_SHORT_NAME,
            optional = false
    )
    protected File inputSiteIntervalsFile;

    @Argument(
            doc = "Output allelic-counts file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false
    )
    protected File outputAllelicCountsFile;

    @Argument(
            doc = "Minimum mapping quality; reads with lower quality will be filtered out of pileup.",
            fullName  = "minimumMappingQuality",
            shortName = "minMQ",
            optional = true
    )
    protected int minimumMappingQuality = 30;

    @Argument(
            doc = "Minimum base quality; base calls with lower quality will be filtered out of pileup.",
            fullName = "minimumBaseQuality",
            shortName = "minBQ",
            optional = true
    )
    protected int minimumBaseQuality = 20;

    @Argument(
            doc = "Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT " +
                    "can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) " +
                    "do not otherwise need to be decoded.",
            fullName = StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_LONG_NAME,
            shortName = StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_SHORT_NAME,
            common = true,
            optional = true
    )
    protected ValidationStringency readValidationStringency = ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY;

    @Override
    protected Object doWork() {
        validateArguments();

        final File referenceFile = referenceArguments.getReferenceFile();
        final IntervalList siteIntervals = IntervalList.fromFile(inputSiteIntervalsFile);

        final AllelicCountCollector allelicCountCollector = new AllelicCountCollector(referenceFile, readValidationStringency);

        logger.info("Collecting allelic counts...");
        final AllelicCountCollection allelicCounts = allelicCountCollector.collect(inputBAMFile, siteIntervals, minimumMappingQuality, minimumBaseQuality);
        allelicCounts.write(outputAllelicCountsFile);
        logger.info("Allelic counts written to " + outputAllelicCountsFile.toString());

        return "SUCCESS";
    }

    private void validateArguments() {
        IOUtils.canReadFile(inputBAMFile);
        IOUtils.canReadFile(inputSiteIntervalsFile);
        IOUtils.canReadFile(referenceArguments.getReferenceFile());
        ParamUtils.isPositiveOrZero(minimumMappingQuality, "Mapping-quality threshold must be greater than or equal to zero.");
        ParamUtils.isPositiveOrZero(minimumBaseQuality, "Base-quality threshold must be greater than or equal to zero.");
    }
}
