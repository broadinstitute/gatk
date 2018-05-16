package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVLocalContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.SVLocalContext.InvBreakEndContext;

/**
 * (Internal) Given GATK-SV discovery pipeline output VCF files and fastq files of short reads sent for local assembly,
 * interpret inversions.
 *
 * <p>
 *     This is an experimental tool and should not be of interest to most researchers.
 *     It scans the input VCF file containing BND records,
 *     applies a primitive filtering step, then
 *     tries to make inversion calls on BND records that pass the filter.
 *     Note that BND records that don't pass the filtering step
 *     are NOT necessarily bad variants, it is simply an indication
 *     that this tool's logic is not applicable to them.
 * </p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>GATK-SV discovery pipeline output simple variant vcf</li>
 *     <li>GATK-SV discovery pipeline output simple variant vcf</li>
 *     <li>Path to directory containing FASTQ files used for local assembly in the discovery pipeline</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>BED file annotating BND records that are not suited for this tool</li>
 *     <li>VCF containing inversion calls, and when available, associated deletion and duplication calls</li>
 *     <li>(optional) intermediate files used in making the records in the above VCF file</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk InversionBreakendInterpreterSpark \
 *     -R reference.2bit \
 *     -O pathToOutput \
 *     --simple-vcf \
 *     --cpx-vcf \
 *     --intermediate-files # optional
 * </pre>
 *
 * <h3>Notes</h3>
 * <p>The reference is broadcast by Spark, and must therefore be a .2bit file due to current restrictions.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) Examines aligned contigs from local assemblies and calls structural variants or their breakpoints",
        summary =
                "This tool is used in development and should not be of interest to most researchers. It is a prototype of" +
                        " structural variant calling, and has been under active developments. For more stable version," +
                        " please see DiscoverVariantsFromContigAlignmentsSAMSpark." +
                        " This tool takes a file containing the alignments of assembled contigs" +
                        " (typically the output file produced by FindBreakpointEvidenceSpark) and searches for reads with" +
                        " split alignments or large gaps indicating the presence of structural variation breakpoints." +
                        " Variations' types are determined by analyzing the signatures of the split alignments," +
                        " and are written to VCF files in the designated output directory.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public final class InversionBreakendInterpreterSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(InversionBreakendInterpreterSpark.class);

    @Argument(doc = "input VCF containing simple variants, particularly BND variants that are inversion breakpoint suspects",
            fullName = "simple-vcf")
    private String nonComplexVCF;

    @Argument(doc = "file containing complex variants as output by GATK-SV",
            fullName = "cpx-vcf")
    private String complexVCF;

    @Argument(doc = "path to directory containing fastqs",
            fullName = "fastqDir")
    private String pathToFastqsDir;

    @Argument(doc = "path to a (temporary) path on a local FS to score bwa mem reference index files",
            fullName = "local-fs-temp-dir")
    private String localFileSystemTempDir;

    @Argument(doc = "prefix for where to dump results, including intermediate files (if requested via --intermediate-files)",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputPrefix;

    @Argument(doc = "flag to indicate intermediate files (alignments of short reads to artificial references) are requested",
            fullName = "intermediate-files", optional = true)
    @Advanced
    private Boolean dumpIntermediateFiles = false;

    @Argument(doc = "", fullName = "mate-dist-threshold")
    @Advanced
    private Integer mateDistanceThreshold = 100_000;

    @Argument(doc = "", fullName = "contig-mq-filter")
    @Advanced
    private Integer contigMQFilter = 30;


    @Override
    public boolean requiresReference() {
        return true;
    }


    @Override
    protected void runTool(final JavaSparkContext ctx) {

        validateParamns();

        final ReferenceMultiSource reference = getReference();
        final List<SimpleInterval> intervals = getIntervals();
        final VCFHeader vcfHeader = VariantsSparkSource.getHeader(nonComplexVCF);
        final SAMSequenceDictionary refDict = vcfHeader.getSequenceDictionary();
        final VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);

        final JavaRDD<VariantContext> complexVariants = variantsSparkSource.getParallelVariantContexts(complexVCF, intervals);
        final JavaRDD<InvBreakEndContext> inversionBreakends = variantsSparkSource.getParallelVariantContexts(nonComplexVCF, intervals)
                .filter(SVLocalContext::indicatesInversion).map(SVLocalContext.InvBreakEndContext::new);

        // note: we are doing this because we know there's going to be only ~500 mate pairs on a particular sample, so paralleling with RDD may hurt performance
        List<InversionBreakendPreFilter.OverlappingPair> preprocessingResult = InversionBreakendPreFilter
                .preprocess(inversionBreakends.collect(), complexVariants.collect(),
                        mateDistanceThreshold, contigMQFilter, outputPrefix, refDict, localLogger);

        final List<VariantContext> variantContexts = LinkedInversionBreakpointsInference
                .makeInterpretation(preprocessingResult, pathToFastqsDir, reference, localLogger);

        // TODO: 5/16/18 artifact that we currently don't have sample columns for discovery VCF (until genotyping code is in)
        final String sampleId = "sample";
        final String vcfName = outputPrefix + (outputPrefix.endsWith("/") ? "" : "/") + sampleId + "_inversions_from_bnds.vcf" ;
        SVVCFWriter.writeVCF(variantContexts, vcfName, refDict, localLogger);
    }

    private void validateParamns() {
        IOUtils.assertFileIsReadable(IOUtils.getPath(nonComplexVCF));
        IOUtils.assertFileIsReadable(IOUtils.getPath(complexVCF));
        ParamUtils.isPositiveOrZero(mateDistanceThreshold,
                "Invalid value provided to mateDistanceThreshold: " + mateDistanceThreshold);
        ParamUtils.isPositiveOrZero(contigMQFilter,
                "Invalid value provided to contigMQFilter: " + contigMQFilter);
    }
}
