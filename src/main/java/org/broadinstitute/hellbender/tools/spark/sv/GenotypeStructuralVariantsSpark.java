package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.commons.lang.StringUtils;
import org.apache.logging.log4j.Logger;
import org.apache.spark.Partitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.sql.Dataset;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredVariantInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.ShardRDD;
import org.broadinstitute.hellbender.engine.spark.SparkSharder;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVContextBuilder;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalLocator;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVContext;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SerializableBiFunction;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignmentUtils;
import org.broadinstitute.hellbender.utils.gcs.BamBucketIoUtils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import picard.util.MathUtil;
import scala.Tuple2;

import java.io.File;
import java.io.Serializable;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.function.IntFunction;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Genotypes structural variants as discovered using {@link StructuralVariationDiscoveryPipelineSpark}.
 *
 * <h2>Inputs</h2>
 * <p>
 *     This tools takes on the following required inputs:
 *     <ul>
 *         <li>a discovered SV vcf file or files (-V or --variant argument),</li>
 *         <li>an assembled contigs alignment (-assemblies or --assembled-contigs-file argument) file
 *         <li>and the name of the directory that contains the assemblies input fastq files (-fastq-dir or --fastq-assembly-directory argument).</li>
 *     </ul>
 * </p>
 *
 * <h2>Outputs</h2>
 * <p>
 *     The main output of this tool is a new enriched vcf that contains the information in the
 *     input vcf file plus genotypes calls for the supported structural variant types.
 *     Non-supported variants are filtered out.
 * </p>
 *
 * <p>
 *     By default, the output genotype would include the total likelihoods and allele depth counts.
 *     However you can ask for stratified likelihoods and allele depths to be emitted for
 *     split-read, interval-size and discordant read pair information
 *     (--emit-stratified-likelihoods and --emit-stratified-allele-depths arguments).
 * </p>
 *
 * <h3>Diagnosis output alignment</h3>
 * <p>
 *     In addition, you can request an output alignment file that encloses for each genotyped variant
 *     the reconstructed haplotypes, overlapping contigs and read pairs involved in the genotyping
 *     process (-bamout --output-diagnosis-alignment argument).
 * </p>
 * <p>
 *     Nonetheless, this output is provided for diagnosis/debugging purposes and may add to much
 *     run-time so is only recommended to be use when analyzing small regions.
 * </p>
 * <p>
 *     Records are annotated with VC:Z:xxx tag with a unique id that identifies the corresponding variant; this can be used to
 *     disambiguate when variants are close to each other.
 * </p>
 * <p>
 *     Records with names "ref" and "alt" represent the <i>reference haplotype</i> and the <i>alternative haplotype</i> for that variant.
 *     These are synthetic records. There is always one and only one of each for each variant. These are all part of the {@code "HAP"} read-group.
 * </p>
 * <p>
 *     Each variant can have a number of <i>assembly contigs</i> with ids of the form <i>asmMMMMMM:tigNNNNNN</i>. These are the ones that
 *     overlap any of the break-points of that variant in the input assemblies file.
 *     These are all part of the {@code "CTG"} read-group.
 * </p>
 * <p>
 *     Finally the actual input read/template records are output under the {@code "TMPL"} read group.
 * </p>
 *  <p>
 *    Contig and template record have a number of annotations indicating amongst other things what allele/haplotype they support:
 *    <dl>
 *         <dt>HP</dt><dd>the haplotype the seem to support, {@code 'ref'} for reference, {@code 'alt'} for alternative and {@code '.'} for neither. (e.g. {@code "HP:Z:ref", "HP:Z:alt", "HP:Z:."}).</dd>
 *         <dt>HQ</dt><dd>the confidence in the support for the haplotype in {@code "HP"} (e.g. {@code "HQ:Z:10.0"}).
 *                        This value is equal to the difference between the score for the supported and the other haplotype</dd></dt>
 *         <dt>RS</dt><dd>the reference support score screen (e.g. {@code "RS:Z:-100,55,1,2,30,0"}).</dd>
 *         <dt>XS</dt><dd>the alternative allele score screen (e.g. {@code "XS:Z:-1241,134,2,1,5,1}).</dd>
 *         <dt>RA</dt><dd>alignment versus the "refHaplotype" haplotype for contigs and vs the reference for haplotypes</dd>
 *         <dt>XA</dt><dd>alignment versus the "altHaplotype" haplotype for contigs. It does not apply to haplotypes</dd>
 *         <dt>VC</dt><dd>coordinate of the targeted variant {@code chr:pos}</dd>
 *         <dt>lT</dt><dd>likelihoods under the reference and alternative haplotypes (templates only)</dd>
 *         <dt>lI</dt><dd>insert size likelihoods under the reference and alternative haplotypes (templates only)</dd>
 *         <dt>lD</dt><dd>discordant pair orientation under the reference and alternative haplotypes (templates only)</dd>
 *         <dt>lR</dt><dd>split read likelihoods under the reference and alternative haplotypes (templates only)</dd>
 *     </dl>
 * </p>
 * <p>
 *     {@code RS} and {@code XS} score annotations follows this format:
 * </p>
 *     <pre>Score,matches,mismatches,indels,indelLength,reversals</pre>
 * <p>
 *     Where:
 *     <dl>
 *         <dt>score</dt><dd>score the contig given the happlotype (refHaplotype. for {@code RS}, altHaplotype. for {@code XS})</dd>
 *         <dt>matches</dt><dd>number of base call matches</dd>
 *         <dt>mismatches</dt><dd>number of base call mismatches</dd>
 *         <dt>indels</dt><dd>number of insertion or deletions (i.e. gap-openings)</dd>
 *         <dt>indelLength</dt><dd>total length of insertion and deletions (i.e. gap-openings + gap-extensions)</dd>
 *     </dl>
 * </p>
 */
@CommandLineProgramProperties(summary = "genotype SV variant call files",
        oneLineSummary = "genotype SV variant call files",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class GenotypeStructuralVariantsSpark extends GATKSparkTool {

    private static final Logger logger = LogManager.getLogger(GenotypeStructuralVariantsSpark.class);

    public static final String FASTQ_FILE_DIR_SHORT_NAME = "fastq-dir";
    public static final String FASTQ_FILE_DIR_FULL_NAME = "fastq-assembly-directory";
    public static final String CTG_FILE_SHORT_NAME = "assemblies";
    public static final String CTG_FILE_FULL_NAME = "assembled-contigs-file";
    public static final String OUTPUT_ALIGNMENT_SHORT_NAME = "bamout";
    public static final String OUTPUT_ALIGNMENT_FULL_NAME = "output-diagnosis-alignment";
    public static final String SHARD_SIZE_SHORT_NAME = "shard";
    public static final String SHARD_SIZE_FULL_NAME = "shard-size";
    public static final String INSERT_SIZE_DISTR_SHORT_NAME = "ins-size-dist";
    public static final String INSERT_SIZE_DISTR_FULL_NAME = "insert-size-distribution";
    public static final String EMIT_GENOTYPING_PERFORMANCE_STATS_SHORT_NAME = "emit-genotyping-stats";
    public static final String EMIT_GENOTYPING_PERFORMANCE_STATS_FULL_NAME = EMIT_GENOTYPING_PERFORMANCE_STATS_SHORT_NAME;
    public static final String EMIT_STRATIFIED_LIKELIHOODS_SHORT_NAME = "emit-stratified-likelihoods";
    public static final String EMIT_STRATIFIED_LIKELIHOODS_FULL_NAME = EMIT_STRATIFIED_LIKELIHOODS_SHORT_NAME;
    public static final String EMIT_STRATIFIED_ALLELE_DEPTHS_SHORT_NAME = "emit-stratified-allele-depths";
    public static final String EMIT_STRATIFIED_ALLELE_DEPTHS_FULL_NAME = EMIT_STRATIFIED_ALLELE_DEPTHS_SHORT_NAME;
    public static final String PADDING_SIZE_SHORT_NAME = "padding";
    public static final String PADDING_SIZE_FULL_NAME = "padding-size";
    public static final int DEFAULT_PADDING_SIZE = 300;

    public static final String HAPLOTYPE_READ_GROUP = "HAP";
    public static final String CONTIG_READ_GROUP = "CTG";
    public static final String TEMPLATE_READ_GROUP = "TMPL";
    public static final String HAPLOTYPE_CALL_TAG = "HP";
    public static final String HAPLOTYPE_QUAL_TAG = "HQ";
    public static final String REFERENCE_SCORE_TAG = "RS";
    public static final String REFERENCE_ALIGNMENT_TAG = "RA";
    public static final String ALTERNATIVE_SCORE_TAG = "XS";
    public static final String ALTERNATIVE_ALIGNMENT_TAG = "XA";
    public static final String VARIANT_CONTEXT_TAG = "VC";
    public static final String TOTAL_LIKELIHOOD_TAG = "lT";
    public static final String STRATIFIED_LIKELIHOOD_TAG_PREFIX = "l";

    /**
     * If the traversal intervals is larger than this fraction of the reference, we simply process all contig alignments
     * in search of primary alignments. Otherwise, we do it in two steps:
     *
     * 1. We scan contig alignments in targeted intervals and add to the
     *    targeted intervals the location of relevant primary alignments
     *    that are not covered by the initial target intervals.
     *
     * 2. We scan the extended target intervals for primary alignments.
     *
     * @see #composeContigTraversalParameters
     */
    private static final double MAXIMUM_REFERENCE_FRACTION_FOR_TARGETED_CONTIG_SEARCH = 0.50;
    private static final String MAXIMUM_REFERENCE_FRACTION_FOR_TARGETED_CONTIG_SEARCH_STR = String.format("%.2f", MAXIMUM_REFERENCE_FRACTION_FOR_TARGETED_CONTIG_SEARCH);
    public static final String UNKNOWN_SAMPLE_NAME = "<UNKNOWN>";
    public static final int MAX_NUMBER_OF_TEMPLATES_IN_CONTEXT = 5000;
    private static final Pattern ASSEMBLY_NAME_ALPHAS = Pattern.compile("[a-zA-Z_]+");

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private RequiredVariantInputArgumentCollection variantArguments = new RequiredVariantInputArgumentCollection();

    @ArgumentCollection
    private RealignmentScoreParameters realignmentScoreArguments = new RealignmentScoreParameters();

    @Argument(doc = "fastq files location",
            shortName = FASTQ_FILE_DIR_SHORT_NAME,
            fullName = FASTQ_FILE_DIR_FULL_NAME)
    private String fastqDir = null;

    @Argument(doc = "assemblies SAM/BAM file location",
            shortName = CTG_FILE_SHORT_NAME,
            fullName = CTG_FILE_FULL_NAME)
    private String contigsFile = null;

    @Argument(doc = "assemblies SAM/BAM file aligned",
            shortName = OUTPUT_ALIGNMENT_SHORT_NAME,
            fullName = OUTPUT_ALIGNMENT_FULL_NAME,
            optional = true)
    private String outputAlignmentFile = null;

    @Argument(doc = "output VCF file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    private String outputFile = null;

    @Argument(doc = "sample name. If not provided, by default is taking from the input vcf it at all present otherwise we use " + UNKNOWN_SAMPLE_NAME + " as a placeholder",
            shortName = StandardArgumentDefinitions.SAMPLE_NAME_SHORT_NAME,
            fullName = StandardArgumentDefinitions.SAMPLE_NAME_LONG_NAME,
            optional = true)
    private String sampleName = null;

    @Argument(doc = "include in the output VCF genotyping performance (cpu-time) relevant annotations such as number of templates, haplotype/contigs involved and genotyping elapse time",
              shortName = EMIT_GENOTYPING_PERFORMANCE_STATS_SHORT_NAME,
              fullName = EMIT_GENOTYPING_PERFORMANCE_STATS_FULL_NAME,
              optional = true)
    private boolean emitGenotypingPerformanceStats = false;

    @Argument(doc = "include in the output VCF stratified likelihoods based on mapping, insert-size and discordant pair-orientation",
              shortName = EMIT_STRATIFIED_LIKELIHOODS_SHORT_NAME,
              fullName = EMIT_STRATIFIED_LIKELIHOODS_FULL_NAME,
              optional = true)
    private boolean emitStratifiedLikelihoods = false;

    @Argument(doc = "include in the output VCF stratified allele depths based on mapping, insert-size and discordant pair-orientation",
              shortName = EMIT_STRATIFIED_ALLELE_DEPTHS_SHORT_NAME,
              fullName = EMIT_STRATIFIED_ALLELE_DEPTHS_FULL_NAME,
              optional = true)
    private boolean emitStratifiedAlleleDepths = false;

    @SuppressWarnings("FieldCanBeLocal")
    @Argument(doc = "shard size",
            shortName = SHARD_SIZE_SHORT_NAME,
            fullName = SHARD_SIZE_FULL_NAME,
            optional = true)
    private int shardSize = 100000;

    @SuppressWarnings("FieldCanBeLocal")
    @Argument(doc = "minimum likelihood (Phred) difference to consider that a template or read support an allele over any other",
            shortName = "infoTLD", fullName = "informativeTemplateLikelihoodDifference", optional = true)
    private double informativeTemplateDifferencePhred = 2.0;

    @Argument(doc = "insert size distribution",
            shortName = INSERT_SIZE_DISTR_SHORT_NAME,
            fullName = INSERT_SIZE_DISTR_FULL_NAME,
            optional = true)
    private InsertSizeDistribution insertSizeDistribution = new InsertSizeDistribution("N(309,149)");

    @Argument(doc ="padding size",
            shortName = PADDING_SIZE_SHORT_NAME,
            fullName = PADDING_SIZE_FULL_NAME,
            minValue = 0.0,
            optional = true)
    private int paddingSize = DEFAULT_PADDING_SIZE;

    private VariantsSparkSource variantsSource;
    private ReadsSparkSource contigsSource;

    @Override
    public boolean requiresReads() {
        return false;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    private void setUp(final JavaSparkContext ctx) {
        variantsSource = new VariantsSparkSource(ctx);
        contigsSource = new ReadsSparkSource(ctx);
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        setUp(ctx);

        final SAMSequenceDictionary dictionary = getReferenceSequenceDictionary();
        final List<SimpleInterval> intervals = this.hasUserSuppliedIntervals() ? getIntervals()
                : IntervalUtils.getAllIntervalsForReference(dictionary);

        try (final SparkSharder sharder = new SparkSharder(ctx, dictionary, intervals, shardSize, 0)) {

            final SVIntervalLocator svIntervalLocator = SVIntervalLocator.of(dictionary);
            final Broadcast<SAMSequenceDictionary> dictionaryBroadcast = ctx.broadcast(dictionary);
            final Broadcast<SVIntervalLocator> svLocatorBroadcast = ctx.broadcast(svIntervalLocator);
            final Broadcast<InsertSizeDistribution> insertSizeDistributionBroadcast = ctx.broadcast(insertSizeDistribution);

            final VCFHeader outputVCFHeader = composeOutputVCFHeader(VariantsSparkSource.getHeader(
                    variantArguments.variantFiles.get(0).getFeaturePath()), sampleName, emitGenotypingPerformanceStats,
                    emitStratifiedLikelihoods, emitStratifiedAlleleDepths);
            final String sampleName = outputVCFHeader.getSampleNamesInOrder().get(0);
            final String referenceFilePath = referenceArguments.getReferenceFileName();

            final TraversalParameters contigAlignmentsTraversalParameters = composeContigTraversalParameters(ctx, contigsSource,
                    contigsFile, referenceFilePath, intervals, dictionaryBroadcast, svLocatorBroadcast, 1000); // 1000bp padding.

            final JavaRDD<SVContext> variants = variantsSource.getParallelVariantContexts(
                    variantArguments.variantFiles.get(0).getFeaturePath(), getIntervals())
                    .map(SVContext::of).filter(GenotypeStructuralVariantsSpark::structuralVariantAlleleIsSupported);

            final JavaRDD<GATKRead> contigAlignments = contigsSource
                    .getParallelReads(contigsFile, referenceArguments.getReferenceFileName(), contigAlignmentsTraversalParameters);

            final JavaRDD<Tuple2<SVContext, List<SVContig>>> variantsAndContigs = composeOverlappingContigRecordsPerVariant(variants, contigAlignments, sharder, dictionaryBroadcast, paddingSize);

            final JavaRDD<Tuple2<SVContext, List<SVHaplotype>>> variantsAndHaplotypes = composeHaplotypes(variantsAndContigs, referenceFilePath, paddingSize);

            final JavaRDD<SVGenotypingContext> genotypingContexts = composeGenotypingContexts(dictionaryBroadcast, svLocatorBroadcast, insertSizeDistributionBroadcast, variantsAndHaplotypes, sampleName)
                    .filter(context -> context.numberOfTemplates > 0 && context.numberOfHaplotypes >= 2);

            final SAMFileHeader outputAlignmentHeader = outputAlignmentFile == null ? null : composeOutputHeader(getReferenceSequenceDictionary(), sampleName);
            final JavaRDD<Call> calls = makeCalls(genotypingContexts, outputAlignmentHeader, ctx);
            SVVCFWriter.writeVCF(outputFile, referenceArguments.getReferenceFileName(), calls.map(c -> c.context), outputVCFHeader, logger);
            if (outputAlignmentFile != null) {
                try (final SAMFileWriter outputAlignmentWriter = BamBucketIoUtils.makeWriter(outputAlignmentFile, outputAlignmentHeader, true)) {
                    calls.flatMap(call -> call.outputAlignmentRecords.iterator())
                            .toLocalIterator()
                            .forEachRemaining(outputAlignmentWriter::addAlignment);
                } catch (final Exception ex) {
                    throw new UserException.CouldNotCreateOutputFile(outputAlignmentFile, ex);
                }
            }
        }
        tearDown(ctx);
    }

    private JavaRDD<SVGenotypingContext> composeGenotypingContexts(final Broadcast<SAMSequenceDictionary> dictionaryBroadcast,
                                                                   final Broadcast<SVIntervalLocator> svLocatorBroadcast,
                                                                   final Broadcast<InsertSizeDistribution> insertSizeDistributionBroadcast,
                                                                   final JavaRDD<Tuple2<SVContext, List<SVHaplotype>>> variantsAndHaplotypes, String sampleName) {
        final String fastqDir = this.fastqDir;
        final String fastqFileNameFormat = "asm%06d.fastq";
        return variantsAndHaplotypes
                .filter(tuple -> tuple._2().size() >= 2)
                .mapPartitions(it -> {
                    final SAMSequenceDictionary localDictionary = dictionaryBroadcast.getValue();
                    final AssemblyCollection assemblyCollection = new AssemblyCollection(fastqDir, fastqFileNameFormat);
                    final SVIntervalLocator locator = svLocatorBroadcast.getValue();
                    final InsertSizeDistribution insertSizeDistribution = insertSizeDistributionBroadcast.getValue();
                    return Utils.map(it, tuple -> composeGenotypingContext(tuple._1, sampleName, tuple._2, localDictionary, assemblyCollection, locator, insertSizeDistribution));
                });
    }

    private static JavaRDD<Tuple2<SVContext, List<SVHaplotype>>> composeHaplotypes(final JavaRDD<Tuple2<SVContext,List<SVContig>>> variantsAndContigs,
                                                                               final String referencePath, final int paddingSize) {
        return variantsAndContigs
                .mapPartitions(it -> {
                    final ReferenceMultiSparkSource referenceSource = new ReferenceMultiSparkSource(referencePath, SimpleInterval::valueOf);
                    return Utils.map(it, tuple -> {
                        final SVContext variant = tuple._1();
                        final List<SVContig> contigs = tuple._2();
                        final List<SimpleInterval> mandatoryCoverage = contigs.stream()
                                .flatMap(ctg -> ctg.getReferenceAlignment().stream().map(ai -> ai.referenceSpan))
                                .collect(Collectors.toList());
                        final List<SVHaplotype> haplotypes = variant.composeHaplotypesBasedOnReference(paddingSize, referenceSource, mandatoryCoverage);
                        final List<SVHaplotype> combined = new ArrayList<>(contigs.size() + haplotypes.size());
                        combined.addAll(haplotypes);
                        combined.addAll(contigs);
                        return new Tuple2<>(variant, combined);
                    });
                });
    }

    /**
     * Composes the appropriate traversal parameters for aligned contigs.
     * <p>
     *     If our traversal does int include the whole reference is possible that for some
     *     aligned contigs the primary alignment that contains its full sequence is not covered
     *     by the targeted interval.
     * </p>
     * <p>
     *     So the portion of the reference that is being targeted is large enough, is simply easier
     *     to scan the whole reference for those alignments.
     * </p>
     * <p>
     *     But it the targeted region is very small, it makes more sense to do a first scan to see whether new
     *     regions need to be added for those missing primary alignments and the proceed with the extended targeted
     *     region.
     * </p>
     * <p>
     *     This method will take care of that returning a traversal parameters that will cover all relevant contig primary
     *     alignments.
     * </p>
     */
    private static TraversalParameters composeContigTraversalParameters(final JavaSparkContext ctx,
                                                                        final ReadsSparkSource contigsSource,
                                                                        final String contigsSourcePath,
                                                                        final String referenceFilePath,
                                                                        final List<SimpleInterval> intervals,
                                                                        final Broadcast<SAMSequenceDictionary> dictionaryBroadcast,
                                                                        final Broadcast<SVIntervalLocator> svIntervalLocatorBroadcast,
                                                                        final int padding) {
        final SAMSequenceDictionary dictionary = dictionaryBroadcast.getValue();
        final SVIntervalLocator svIntervalLocator = svIntervalLocatorBroadcast.getValue();
        final long intervalsTotalSize = intervals.stream().mapToLong(SimpleInterval::size).sum();
        final long referenceTotalSize = dictionary.getSequences().stream().mapToLong(SAMSequenceRecord::getSequenceLength).sum();
        final double targetedReferenceFraction = intervalsTotalSize / (double) referenceTotalSize;
        final String percentage = String.format("%.2f", targetedReferenceFraction * 100);

        if (targetedReferenceFraction > MAXIMUM_REFERENCE_FRACTION_FOR_TARGETED_CONTIG_SEARCH) {
            logger.info("We proceed to scan the entire contig assembly file since the %% of targeted reference (%s) is larger than the maximum percentage for targeted scanning (%s)",
                    percentage, MAXIMUM_REFERENCE_FRACTION_FOR_TARGETED_CONTIG_SEARCH_STR);
            return new TraversalParameters(IntervalUtils.getAllIntervalsForReference(dictionary), false);
        } else {
            final TraversalParameters firstPassTraversalParameters = new TraversalParameters(intervals, false);
            final List<SimpleInterval> additionalIntervals = firstScanForMissingPrimaryContigAlignments(ctx, contigsSource,
                    contigsSourcePath, referenceFilePath, intervals, percentage,
                    firstPassTraversalParameters, svIntervalLocatorBroadcast);
            if (additionalIntervals.isEmpty()) {
                logger.info("No additional interval are need to look for primary contig alignments");
                return firstPassTraversalParameters;
            } else {
                logger.info("There is total of " + additionalIntervals.size() + " contigs whose primary alignment are not covered by the targeted intervals");
                final List<SimpleInterval> intervalsToAdd = composeAdditionIntervalsToCoverMissingPrimaryContigAlignments(padding, svIntervalLocator, additionalIntervals);
                final long basesAdded = intervalsToAdd.stream().mapToLong(SimpleInterval::size).sum();
                logger.info("We will add a total of " + intervalsToAdd + " intervals expanding " + basesAdded + " bases");
                final List<SimpleInterval> newIntervals = Stream.concat(intervals.stream(), intervalsToAdd.stream()).collect(Collectors.toList());
                return new TraversalParameters(newIntervals, false);
            }
        }
    }

    /**
     * Compose a concise list of additional intervals to add given the missing primary alignment locations.
     * These intervals are padded and merged if they overlap so that we avoid redundancy:
     */
    private static List<SimpleInterval> composeAdditionIntervalsToCoverMissingPrimaryContigAlignments(final int padding,
                                                                                                      final SVIntervalLocator locator,
                                                                                                      final List<SimpleInterval> additionalIntervals) {
        final SVIntervalTree<SimpleInterval> mergeTree = new SVIntervalTree<>();
        for (final SimpleInterval interval : additionalIntervals) {
            final SVInterval svInterval = locator.toSVInterval(interval, padding);
            if (!mergeTree.hasOverlapper(svInterval)) {
                mergeTree.put(svInterval, interval);
            } else { // we need to merge the overlapping intervals.
                final List<SVInterval> overlappers = Utils.stream(mergeTree.overlappers(svInterval))
                        .map(SVIntervalTree.Entry::getInterval)
                        .collect(Collectors.toList());
                int start = svInterval.getStart();
                int end = svInterval.getEnd();
                for (final SVInterval overlapper : overlappers) {
                    if (overlapper.getStart() < start) start = overlapper.getStart();
                    if (overlapper.getEnd() > end) end = overlapper.getEnd();
                    mergeTree.remove(overlapper);
                }
                final SVInterval superInterval = new SVInterval(svInterval.getContig(), start, end);
                mergeTree.put(superInterval, locator.toSimpleInterval(superInterval));
            }
        }
        return Utils.stream(mergeTree.iterator())
                  .map(SVIntervalTree.Entry::getValue)
                  .collect(Collectors.toList());
    }

    /**
     * Performs a first scan looking for missing primary alignments.
     *
     * @return list of the reference-span for each missing primary alignment.
     */
    private static List<SimpleInterval> firstScanForMissingPrimaryContigAlignments(
            final JavaSparkContext ctx,
            final ReadsSparkSource contigsSource,
            final String contigsSourcePath,
            final String referenceFilePath,
            final List<SimpleInterval> intervals,
            final String percentage,
            final TraversalParameters firstPassTraversalParameters,
            final Broadcast<SVIntervalLocator> locatorBroadcast) {

        logger.info("We proceed to scan targeted intervals for missing contig primary alignments as the %% of targeted reference (%s) is less or equal to the maximum for interval scan (%s)",
                percentage, MAXIMUM_REFERENCE_FRACTION_FOR_TARGETED_CONTIG_SEARCH_STR);
        final long startTime = System.currentTimeMillis();

        final SVIntervalLocator locator = locatorBroadcast.getValue();
        final SVIntervalTree<SimpleInterval> intervalSearchTree = new SVIntervalTree<>();
        for (final SimpleInterval interval : intervals) {
            intervalSearchTree.put(locator.toSVInterval(interval), interval);
        }
        final Broadcast<SVIntervalTree<SimpleInterval>> intervalTreeBroadcast = ctx.broadcast(intervalSearchTree);

        final List<AlignmentInterval> missingAlignments = contigsSource.getParallelReads(contigsSourcePath, referenceFilePath, firstPassTraversalParameters)
                .mapPartitions( it -> {
                    final SVIntervalTree<SimpleInterval> searchTree = intervalTreeBroadcast.getValue();
                    final SVIntervalLocator localLocator = locatorBroadcast.getValue();
                    return Utils.stream(it)
                            .filter(GATKRead::isSupplementaryAlignment)
                            .map(sa -> {
                                final List<AlignmentInterval> otherAlignments = AlignmentInterval
                                        .decodeList(sa.getAttributeAsString(SAMTag.SA.name()));
                                if (otherAlignments.isEmpty()) {
                                    throw new GATKException("unexpected lack of records in SA tag for supplementary alignment " + sa.getName() + " at " + sa.getAssignedContig() + ":" + sa.getAssignedStart());
                                } else {
                                    final AlignmentInterval first = otherAlignments.get(0);
                                    if (first.cigarAlong5to3DirectionOfContig.containsOperator(CigarOperator.H)) {
                                        throw new GATKException("unexpected hard-clips in what is mean to be the primary alignment in SA tag for supplementary alignment " + sa.getName() + " at " + sa.getAssignedContig() + ":" + sa.getAssignedStart());
                                    }
                                    return first;
                                }
                            })
                            .filter(pa -> !searchTree.hasOverlapper(localLocator.toSVInterval(pa.referenceSpan)))
                            .iterator();
                }).distinct().collect();
        intervalTreeBroadcast.destroy(false);
        final long stopTime = System.currentTimeMillis();
        final long elapse = stopTime - startTime;
        final long elapseMinutes =  TimeUnit.MILLISECONDS.toMinutes(elapse);
        final long elapseSeconds = TimeUnit.MILLISECONDS.toSeconds(elapse) - TimeUnit.MINUTES.toSeconds(elapseMinutes);
        logger.info("The scan took " + elapse + " milliseconds (" + elapseMinutes + " minutes, " + elapseSeconds + " seconds).");
        return missingAlignments.stream().map(ai -> ai.referenceSpan).collect(Collectors.toList());
    }

    /**
     * Given an RDD of assembled contigs alignments as {@link GATKRead} it composes and RDD of {@link AlignedContig}.
     * <p>
     *     Each assembled contigs may have more than one alignment GATKRecord although only one, the primary may contain
     *     the full contig sequence.
     * </p>
     * <p>
     *     This method get rid of these duplicates, returning a single record per contig that contains all
     *     the supplementary alignment information.
     * </p>
     * <p>
     *     This code assumes that there is exactly only one primary alignment per contig
     *     that contains the full assembled sequence. Supplementary alignments may or may not
     *     have all the bases.
     * </p>
     * @param alignments the input alignments.
     * @return never {@code null}.
     */
    private static JavaRDD<AlignedContig> composeAlignedContigs(final JavaRDD<GATKRead> alignments) {
        return alignments.filter(GATKRead::isPrimaryAlignment)
                  .map(r ->
                  {  if (r.isUnmapped()) {
                        throw new GATKException("the input assemblies file contain unmapped records");
                     }
                     return new AlignedContig(r);
                  });
    }

    private static JavaRDD<Tuple2<SVContext, List<SVContig>>> composeOverlappingContigRecordsPerVariant(
            final JavaRDD<SVContext> variants,
            final JavaRDD<GATKRead> contigs,
            final SparkSharder sharder,
            final Broadcast<SAMSequenceDictionary> dictionaryBroadcast,
            final int padding) {

        // We group SVContext and GATKReads (assembled contigs) into common shards:

        // Notice that SVContext should be added to all shards where it has a break-point:
        final ShardRDD<SVContext> variantsSharded = sharder
                .shard(variants, () -> {
                    final SAMSequenceDictionary dictionary = dictionaryBroadcast.getValue();
                    return v -> v.getBreakPointIntervals(padding, dictionary, true);
                });


        // Equally aligned-contigs should be added to any shard it has a primary or supplementary alignment.
        final JavaRDD<AlignedContig> alignedContigs = composeAlignedContigs(contigs);
        final ShardRDD<AlignedContig> alignedContigsSharded = sharder
                .shard(alignedContigs, ac -> ac.getAlignments().stream().map(ai -> ai.referenceSpan).collect(Collectors.toList()));
        // We join together variants and contigs that map to the same shard.
        final JavaPairRDD<SVContext, List<AlignedContig>> variantsAndContigsShared =
                variantsSharded.groupRight(alignedContigsSharded, (ctxs, alcs) -> {
                    final SAMSequenceDictionary dictionary = dictionaryBroadcast.getValue();
                    final SVIntervalLocator locator = SVIntervalLocator.of(dictionary);
                    final SVIntervalTree<List<AlignedContig>> contigsByInterval = composeLookupSVIntervalTree(locator, alcs,
                            contig -> contig.getAlignments().stream().map(ai -> ai.referenceSpan).collect(Collectors.toList()));
                    return ctxs.stream()
                            .map(v -> {
                                final List<AlignedContig> overlappingContigs = lookupSVIntervalTree(contigsByInterval, locator,
                                        v.getBreakPointIntervals(padding, dictionary, true).stream());
                                return new Tuple2<>(v, overlappingContigs);
                            }).collect(Collectors.toList());
                }).toPairRDD()
                  .values()
                  .flatMapToPair(List::iterator)
                  .reduceByKey((list1, list2) -> Stream.concat(list1.stream(), list2.stream()).distinct()
                          .collect(Collectors.toList()));
        // transform AlignedContig into SVContigs
        return variantsAndContigsShared
                .mapPartitions(it -> {
                    final SAMSequenceDictionary dictionary = dictionaryBroadcast.getValue();
                    return Utils.map(it, tuple -> new Tuple2<>(tuple._1,
                                 tuple._2.stream()
                                         .map(ac -> SVContig.of(ac, tuple._1, dictionary, padding))
                                         .collect(Collectors.toList())));
                });
    }

    private static <V> SVIntervalTree<List<V>> composeLookupSVIntervalTree(final SVIntervalLocator locator, Collection<? extends V> subjects, Function<? super V, Iterable<? extends Locatable>> locationsOf) {
        final SVIntervalTree<List<V>> result = new SVIntervalTree<>();
        List<V> nextEmptyList = new ArrayList<>(10);
        for (final V subject : subjects) {
            for (final Locatable loc : locationsOf.apply(subject)) {
                final List<V> extantList = result.putIfAbsent(locator.toSVInterval(loc), nextEmptyList);
                if (extantList == null) {
                    nextEmptyList.add(subject);
                    nextEmptyList = new ArrayList<>(10);
                } else {
                    extantList.add(subject);
                }
            }
        }
        return result;
    }

    private static <V> List<V> lookupSVIntervalTree(final SVIntervalTree<List<V>> tree, final SVIntervalLocator locator, final Stream<? extends Locatable> loc) {
        return loc.map(locator::toSVInterval)
                  .map(tree::overlappers)
                  .flatMap(Utils::stream)
                  .map(SVIntervalTree.Entry::getValue)
                  .flatMap(Collection::stream)
                  .distinct()
                  .collect(Collectors.toList());
    }

    /**
     * Given a variant and the list of relevant haplotypes composes the genotyping context object.
     * <p>
     *     For that it needs to gather all the relevant templates (fragment/read pairs) that were
     *     use on the assemblies associated to the haplotypes.
     * </p>
     */
    private static SVGenotypingContext composeGenotypingContext(final SVContext variant, final String sampleName,
                                                                final List<SVHaplotype> haplotypes,
                                                                final SAMSequenceDictionary dictionary,
                                                                final AssemblyCollection assemblyCollection,
                                                                final SVIntervalLocator locator,
                                                                final InsertSizeDistribution insertSizeDistribution) {
        final SVIntervalTree<SimpleInterval> coveredReference = haplotypes.stream()
                .filter(h -> !h.isNeitherReferenceNorAlternative())
                .flatMap(h -> h.getReferenceAlignment().stream())
                .flatMap(ai -> ai.referenceCoveredIntervals().stream())
                .collect(locator.toSVIntervalTree(si -> si));
        final IntStream assemblyNumbers = haplotypes.stream()
                .filter(SVHaplotype::isNeitherReferenceNorAlternative)
                .map(SVHaplotype::getName)
                .map(n -> n.substring(0, n.indexOf(":")))
                .map(a -> ASSEMBLY_NAME_ALPHAS.matcher(a).replaceAll("0"))
                .mapToInt(Integer::parseInt)
                .sorted()
                .distinct();
        final List<Template> allTemplates = assemblyNumbers
                .boxed()
                .flatMap(i -> assemblyCollection.templates(i).stream())
                .distinct()
                .peek(template -> template.calculateMaximumMappingQualities(coveredReference, locator, insertSizeDistribution))
                .collect(Collectors.toList());
        return new SVGenotypingContext(variant, haplotypes, allTemplates, sampleName, dictionary);
    }

    private SAMFileHeader composeOutputHeader(final SAMSequenceDictionary dictionary, final String sampleName) {
        final SAMFileHeader outputHeader = new SAMFileHeader();
        outputHeader.setSequenceDictionary(dictionary);
        final SAMReadGroupRecord haplotypeReadGroup = new SAMReadGroupRecord(HAPLOTYPE_READ_GROUP);
        haplotypeReadGroup.setSample(sampleName);
        final SAMReadGroupRecord contigReadGroup = new SAMReadGroupRecord(CONTIG_READ_GROUP);
        contigReadGroup.setSample(sampleName);
        final SAMReadGroupRecord templateReadGroup = new SAMReadGroupRecord(TEMPLATE_READ_GROUP);
        templateReadGroup.setSample(sampleName);
        outputHeader.addReadGroup(haplotypeReadGroup);
        outputHeader.addReadGroup(contigReadGroup);
        outputHeader.addReadGroup(templateReadGroup);
        final SAMProgramRecord programRecord = new SAMProgramRecord(getProgramName());
        programRecord.setCommandLine(getCommandLine());
        outputHeader.addProgramRecord(programRecord);
        outputHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return outputHeader;
    }

    private VCFHeader composeOutputVCFHeader(final VCFHeader inputHeader, final String sampleName,
                                          final boolean includePerformanceAnnotations,
                                          final boolean includeStratifiedLikelihoods,
                                          final boolean includeStratifiedAlleleDepths) {
        final List<String> sampleNames = inputHeader.getSampleNamesInOrder();
        final List<String> sampleNameList = Collections.singletonList(sampleName == null && sampleNames.isEmpty() ? "<UNKNOWN>" : (sampleName == null) ?sampleNames.get(0) : sampleName);
        final Set<VCFHeaderLine> inputHeaderLines = inputHeader.getMetaDataInInputOrder();
        final Collection<VCFHeaderLine> genotypingLines = GATKSVVCFHeaderLines.genotypingHeaderLines(includePerformanceAnnotations, includeStratifiedLikelihoods, includeStratifiedAlleleDepths);
        final Set<VCFHeaderLine> outputHeaderLines = new LinkedHashSet<>(inputHeaderLines.size() + 20);
        outputHeaderLines.addAll(inputHeaderLines);
        outputHeaderLines.addAll(genotypingLines);
        outputHeaderLines.add(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1, VCFHeaderLineType.Integer, "last base position of the variant"));
        return new VCFHeader(outputHeaderLines, sampleNameList);
    }

    private static class Call implements Serializable {
        private static final long serialVersionUID = -1L;
        public final SVContext context;
        public final List<SAMRecord> outputAlignmentRecords;

        private Call(final SVContext context, final List<SAMRecord> outputAlignmentRecords) {
            this.context = context;
            this.outputAlignmentRecords = outputAlignmentRecords;
        }

        public static Call of(final SVContext context) {
            return new Call(context, Collections.emptyList());
        }

        public static Call of(final SVContext context, final List<SAMRecord> outputAlignmentRecords) {
            return new Call(context, outputAlignmentRecords);
        }
    }

    @SuppressWarnings("unused")
    private static Object[] debug00(final LikelihoodMatrix<?> left, final LikelihoodMatrix<?> right) {
        final Object[] result = new Object[left.numberOfReads()];
        for (int i = 0; i < result.length; i++) {
            final double one = left.get(0, i) + right.get(0, i);
            final double two = left.get(1, i) + right.get(1, i);
            result[i] = StringUtils.join(new Object[]{left.get(0, i), right.get(0, i), left.get(1, i), right.get(1, i), left.get(0, i) + right.get(0, i),
                    left.get(1, i) + right.get(1, i), one - two, (one > two) ? -Math.min(0, two - 2 * one) : Math.min(0, one - 2 * two)}, ",");
        }
        return result;
    }

    private JavaRDD<Call> makeCalls(final JavaRDD<SVGenotypingContext> input, final SAMFileHeader outputAlignmentHeader, final JavaSparkContext ctx) {
        final Broadcast<SAMSequenceDictionary> broadCastDictionary = ctx.broadcast(getReferenceSequenceDictionary());
        final InsertSizeDistribution insertSizeDistribution = this.insertSizeDistribution;
        final RealignmentScoreParameters realignmentScoreArguments = this.realignmentScoreArguments;
        final boolean emitGenotypingPerformanceStats = this.emitGenotypingPerformanceStats;
        final boolean emitStratifiedAlleleDepths = this.emitStratifiedAlleleDepths;
        final boolean emitStratifiedLikelihoods = this.emitStratifiedLikelihoods;
        final double informativePhredLikelihoodDifferenceThreshold = this.informativeTemplateDifferencePhred;
        final int paddingSize = this.paddingSize;
        return input.mapPartitions(it -> {
            final SAMSequenceDictionary dictionary = broadCastDictionary.getValue();
            final SAMFileHeader header = new SAMFileHeader();
            header.setSequenceDictionary(dictionary);
            final GenotypeLikelihoodCalculator genotypeCalculator =
                    new GenotypeLikelihoodCalculators().getInstance(2, 2);
            return Utils.map(it, context -> {
                        context.reduceNumberOfTemplatesTo(MAX_NUMBER_OF_TEMPLATES_IN_CONTEXT);
                        final List<Template> templates = context.templates;

                        final TemplateMappingTable scoreTable = remapTemplatesOnHaplotypes(context, realignmentScoreArguments);
                        setMappingQualityToZeroForFragmentsThatMapDifferentlyOnFullReference(context, scoreTable);
                        setMissingAlignmentScores(realignmentScoreArguments, templates, scoreTable);

                        // we check what templates are relevant toward genotyping, by default only those that map across a break point.

                        final ReadLikelihoods<SVGenotypingContext.Allele> splitsReadlikelihoods = calculateSplitReadLikelihoods(context, realignmentScoreArguments, scoreTable);
                        final ReadLikelihoods<SVGenotypingContext.Allele> insertSizeLikelihoods = calculateInsertSizeLikelihoods(context, realignmentScoreArguments, scoreTable, insertSizeDistribution);
                        final ReadLikelihoods<SVGenotypingContext.Allele> discordantOrientationLikelihoods = calculateDiscordantOrientationLikelihoods(context, realignmentScoreArguments, scoreTable);
                        final ReadLikelihoods<SVGenotypingContext.Allele> totalLikelihoods = ReadLikelihoods.sum(splitsReadlikelihoods, insertSizeLikelihoods, discordantOrientationLikelihoods);

                        return composeCall(outputAlignmentHeader, realignmentScoreArguments, emitGenotypingPerformanceStats,
                                emitStratifiedAlleleDepths, emitStratifiedLikelihoods, informativePhredLikelihoodDifferenceThreshold,
                                insertSizeDistribution, paddingSize, dictionary, genotypeCalculator, context, splitsReadlikelihoods,
                                insertSizeLikelihoods, discordantOrientationLikelihoods, totalLikelihoods);
                    });
        });
    }

    private static Call composeCall(final SAMFileHeader outputAlignmentHeader,
                                    final RealignmentScoreParameters realignmentScoreArguments,
                                    final boolean emitGenotypingPerformanceStats,
                                    final boolean emitStratifiedAlleleDepths,
                                    final boolean emitStratifiedLikelihoods,
                                    final double informativePhredLikelihoodDifferenceThreshold,
                                    final InsertSizeDistribution insertSizeDistribution,
                                    final int paddingSize,
                                    final SAMSequenceDictionary dictionary,
                                    final GenotypeLikelihoodCalculator genotypeCalculator,
                                    final SVGenotypingContext context,
                                    final ReadLikelihoods<SVGenotypingContext.Allele> splitsReadlikelihoods,
                                    final ReadLikelihoods<SVGenotypingContext.Allele> insertSizeLikelihoods,
                                    final ReadLikelihoods<SVGenotypingContext.Allele> discordantOrientationLikelihoods,
                                    final ReadLikelihoods<SVGenotypingContext.Allele> totalLikelihoods) {
        adjustLikelihoodCalculatorAlleleFrequencies(context, insertSizeDistribution, genotypeCalculator);
        final SVContextBuilder outputBuilder = new SVContextBuilder(context.variant);
        if (emitGenotypingPerformanceStats) {
            outputBuilder.setGenotypingContextSizes(context.numberOfTemplates, context.numberOfHaplotypes);
            outputBuilder.startRecordingProcessingTime();
        }
        if (genotypeCalculator.getRelativeAlleleFrequency() != null) {
            outputBuilder.attribute(GATKSVVCFConstants.EXPECTED_RELATIVE_ALLELE_FREQUENCIES, genotypeCalculator.getRelativeAlleleFrequency());
        }
        final Genotype genotype = composeGenotype(context, totalLikelihoods, genotypeCalculator, splitsReadlikelihoods, discordantOrientationLikelihoods, insertSizeLikelihoods, realignmentScoreArguments, informativePhredLikelihoodDifferenceThreshold, emitStratifiedAlleleDepths, emitStratifiedLikelihoods);
        outputBuilder.genotypes(genotype);
        if (outputAlignmentHeader == null) {
            return Call.of(outputBuilder.make());
        } else {
            final List<SimpleInterval> relevantIntervals = context.variant.getBreakPointIntervals(paddingSize, true);
            final Map<String, ReadLikelihoods<SVGenotypingContext.Allele>> stratifiedLikelihoods = new HashMap<>(3);
            stratifiedLikelihoods.put(GATKSVVCFConstants.TEMPLATE_MAPPING_LIKELIHOODS, splitsReadlikelihoods);
            stratifiedLikelihoods.put(GATKSVVCFConstants.INSERT_SIZE_LIKELIHOODS, insertSizeLikelihoods);
            stratifiedLikelihoods.put(GATKSVVCFConstants.DISCORDANT_PAIR_ORIENTATION_LIKELIHOODS, discordantOrientationLikelihoods);
            return Call.of(outputBuilder.make(), composeDiagnosisReads(context, outputAlignmentHeader, dictionary, realignmentScoreArguments, relevantIntervals, totalLikelihoods, stratifiedLikelihoods, informativePhredLikelihoodDifferenceThreshold));
        }
    }

    private static List<SAMRecord> composeDiagnosisReads(final SVGenotypingContext context, final SAMFileHeader outputAlignmentHeader, final SAMSequenceDictionary dictionary,
                                                         final RealignmentScoreParameters parameters,
                                                         final List<SimpleInterval> relevantIntervals,
                                                         final ReadLikelihoods<SVGenotypingContext.Allele> totalLikelihoods,
                                                         final Map<String, ReadLikelihoods<SVGenotypingContext.Allele>> stratifiedLikelihoods, final double informativePherdDiff) {
        final List<SAMRecord> result = new ArrayList<>(context.numberOfHaplotypes + context.numberOfTemplates);
        for (final SVHaplotype haplotype : context.haplotypes) {
            final List<SAMRecord> haplotypesReads = composeAlignedOutputHaplotypeSAMRecords(outputAlignmentHeader, context, haplotype, parameters);
            result.addAll(haplotypesReads);
        }
        for (final Template template : context.templates) {
            final List<SAMRecord> templateReads = composeAlignedOutputTemplateSAMRecords(outputAlignmentHeader, context, template, relevantIntervals, totalLikelihoods, stratifiedLikelihoods, informativePherdDiff);
            result.addAll(templateReads);
        }
        //We do the sort outside.
        //result.sort(Comparator.comparingInt(SAMRecord::getReferenceIndex).thenComparingInt(SAMRecord::getAlignmentStart));
        return result;
    }

    private static List<SAMRecord> composeAlignedOutputTemplateSAMRecords(final SAMFileHeader outputAlignmentHeader, final SVGenotypingContext context, final Template template, final List<SimpleInterval> relevantIntervals,
                                                                          final ReadLikelihoods<SVGenotypingContext.Allele> totalLikelihoods,
                                                                          final Map<String,ReadLikelihoods<SVGenotypingContext.Allele>> stratifiedLikelihoods, final double informativePhredDiff) {

        final boolean hasRelevantMapping = template.fragments().stream()
                .flatMap(fragment -> fragment.alignmentIntervals().stream())
                .anyMatch(interval -> relevantIntervals.stream().anyMatch(target -> target.overlaps(interval.referenceSpan)));

        if (!hasRelevantMapping && (totalLikelihoods.containsRead(0, template.name()) ||
                    stratifiedLikelihoods.values().stream().anyMatch(lk -> lk.containsRead(0, template.name())))) {
            return Collections.singletonList(composeAlignedOutputTemplateSAMRecordForUnmapped(outputAlignmentHeader, template, context, totalLikelihoods, stratifiedLikelihoods, informativePhredDiff));
        } else if (hasRelevantMapping) {
            return composeAlignedOutputTemplateSAMRecordsForMapped(outputAlignmentHeader, template, context, relevantIntervals, totalLikelihoods, stratifiedLikelihoods, informativePhredDiff);
        } else {
            return Collections.emptyList();
        }
    }

    private static List<SAMRecord> composeAlignedOutputTemplateSAMRecordsForMapped(final SAMFileHeader outputAlignmentHeader, final Template template, final SVGenotypingContext context, final List<SimpleInterval> relevantIntervals, final ReadLikelihoods<SVGenotypingContext.Allele> totalLikelihoods, final Map<String, ReadLikelihoods<SVGenotypingContext.Allele>> stratifiedLikelihoods, final double informativePhredDiff) {
        final Collection<SAMRecord.SAMTagAndValue> likelihoodAnnotations = composeLikelihoodsAnnotations(template.name(), totalLikelihoods, stratifiedLikelihoods, informativePhredDiff);

        return template.fragments().stream()
                .flatMap(fragment ->
                    fragment.alignmentIntervals().stream()
                            .filter(ai -> relevantIntervals.stream().anyMatch(target -> target.overlaps(ai.referenceSpan)))
                            .map(ai -> ai.toSAMRecord(outputAlignmentHeader, template.name(), fragment.bases(), false,
                                    fragment == template.fragments().get(0) ? SAMFlag.FIRST_OF_PAIR.intValue() : 0, likelihoodAnnotations))
                )
                .collect(Collectors.toList());
    }

    private static Collection<SAMRecord.SAMTagAndValue> composeLikelihoodsAnnotations(final String name, final ReadLikelihoods<SVGenotypingContext.Allele> totalLikelihoods, final Map<String,ReadLikelihoods<SVGenotypingContext.Allele>> stratifiedLikelihoods, final double informativePhredDiff) {
        final Collection<SAMRecord.SAMTagAndValue> result = new ArrayList<>(10);
        if (totalLikelihoods.containsRead(0, name)) {
            result.add(new SAMRecord.SAMTagAndValue(TOTAL_LIKELIHOOD_TAG, composeLikelihoodTagValue(totalLikelihoods, name)));
            result.add(new SAMRecord.SAMTagAndValue(HAPLOTYPE_CALL_TAG, composeLikelihoodCall(totalLikelihoods, name, informativePhredDiff)));
        } else {
            result.add(new SAMRecord.SAMTagAndValue(HAPLOTYPE_CALL_TAG, Allele.NO_CALL_STRING));
        }
        for (final Map.Entry<String, ReadLikelihoods<SVGenotypingContext.Allele>> entry : stratifiedLikelihoods.entrySet()) {
            if (entry.getValue().containsRead(0, name)) {
                result.add(new SAMRecord.SAMTagAndValue(STRATIFIED_LIKELIHOOD_TAG_PREFIX + entry.getKey().charAt(0), composeLikelihoodTagValue(entry.getValue(), name)));
            }
        }
        return result;
    }

    private static String composeLikelihoodTagValue(final ReadLikelihoods<SVGenotypingContext.Allele> lk, final String name) {
        final int readIndex = lk.readIndex(0, name);
        if (readIndex < 0) {
            return VCFConstants.MISSING_VALUE_v4;
        } else {
            return IntStream.range(0, lk.numberOfAlleles())
             .mapToDouble(i -> lk.sampleMatrix(0).get(i, readIndex))
                    .mapToObj(d -> String.format("%.1f", d == 0.00 ? 0.0 : d * -10))
                    .collect(Collectors.joining(","));
        }
    }

    private static String composeLikelihoodCall(final ReadLikelihoods<SVGenotypingContext.Allele> lk, final String name, final double minimumPhredDiff) {
        final int readIndex = lk.readIndex(0, name);
        if (readIndex < 0) {
            return ".";
        } else {
            final double[] values = IntStream.range(0, lk.numberOfAlleles())
                    .mapToDouble(i -> lk.sampleMatrix(0).get(i, readIndex))
                    .toArray();
            double bestValue = values[0];
            double secondBestValue = Double.NEGATIVE_INFINITY;
            int bestIndex = 0;
            for (int i = 0; i < values.length; i++) {
                if (values[i] > bestValue) {
                    secondBestValue = bestValue;
                    bestValue = values[i];
                    bestIndex = i;
                } else if (values[i] > secondBestValue) {
                    secondBestValue = values[i];
                }
            }
            final double conf = (bestValue - secondBestValue) * 10;
            if (conf < minimumPhredDiff) {
                return Allele.NO_CALL_STRING;
            } else if (lk.alleles().get(bestIndex).isReference()) {
                return SVHaplotype.REF_HAPLOTYPE_NAME;
            } else {
                return SVHaplotype.ALT_HAPLOTYPE_NAME;
            }
        }
    }

    private static SAMRecord composeAlignedOutputTemplateSAMRecordForUnmapped(final SAMFileHeader header, final Template template, final SVGenotypingContext context, ReadLikelihoods<SVGenotypingContext.Allele> totalLikelihoods, Map<String, ReadLikelihoods<SVGenotypingContext.Allele>> stratifiedLikelihoods, double informativePhredDiff) {
        final SAMRecord sam = new SAMRecord(header);
        sam.setAlignmentStart(context.variant.getStart());
        sam.setReferenceName(context.variant.getContig());
        sam.setReadUnmappedFlag(false);
        final int length = template.fragments().stream().mapToInt(fragment -> fragment.length() + 20).sum() - 20;
        int nextBase = 0;
        final byte[] sequence = new byte[length];
        final List<AlignmentInterval> intervals = new ArrayList<>(10);
        for (final Template.Fragment fragment : template.fragments()) {
            final int offset = nextBase;
            final int tail = sequence.length - fragment.length() - nextBase;
            fragment.alignmentIntervals().stream()
                    .map(original -> new AlignmentInterval(original.referenceSpan,
                            original.startInAssembledContig + offset,
                            original.endInAssembledContig + offset,
                            addClippingToCigar(original.cigarAlong5to3DirectionOfContig,
                                    original.forwardStrand ? offset : tail,
                                    original.forwardStrand ? tail : offset, offset != 0 && original != fragment.alignmentIntervals().get(0)),
                            original.forwardStrand,
                            original.mapQual,
                            original.mismatches,
                            original.alnScore,
                            original.alnModType))
                    .forEach(intervals::add);
            fragment.copyBases(sequence, nextBase, 0, fragment.length());
            nextBase += fragment.length();
            if (nextBase < sequence.length) {
                Arrays.fill(sequence, nextBase, nextBase += 20, (byte) 'N');
            }
        }
        sam.setAttribute(SAMTag.SA.name(), AlignmentInterval.encode(intervals));
        sam.setReadBases(sequence);
        sam.setProperPairFlag(false);
        sam.setReadNegativeStrandFlag(false);
        sam.setMappingQuality(0);
        sam.setAttribute(SAMTag.RG.name(), TEMPLATE_READ_GROUP);
        final Collection<SAMRecord.SAMTagAndValue> likelihoodAnnotations = composeLikelihoodsAnnotations(template.name(), totalLikelihoods, stratifiedLikelihoods, informativePhredDiff);
        for (final SAMRecord.SAMTagAndValue annotation : likelihoodAnnotations) {
            sam.setAttribute(annotation.tag, annotation.value);
        }
        return sam;
    }

    private static Cigar addClippingToCigar(final Cigar original, final int left, final int right, final boolean hardIfUndefined) {
        if (left <= 0 && right <= 0) {
            return new Cigar(original.getCigarElements());
        } else {
            final List<CigarElement> newElements = new ArrayList<>(original.getCigarElements().size() + 2);
            if (left > 0) {
                if (original.getCigarElement(0).getOperator().isClipping()) {
                    newElements.add(new CigarElement(original.getCigarElement(0).getLength() + left,
                            original.getCigarElement(0).getOperator()));
                    newElements.addAll(original.getCigarElements().subList(1, original.numCigarElements()));
                } else {
                    newElements.add(new CigarElement(left, hardIfUndefined ? CigarOperator.H : CigarOperator.S));
                    newElements.addAll(original.getCigarElements());
                }
            } else {
                newElements.addAll(original.getCigarElements());
            }
            if (right > 0) {
                final CigarElement lastElement = newElements.get(newElements.size() - 1);
                if (lastElement.getOperator().isClipping()) {
                    newElements.set(newElements.size() - 1, new CigarElement(right + lastElement.getLength(), lastElement.getOperator()));
                } else {
                    newElements.add(new CigarElement(right, hardIfUndefined ? CigarOperator.H : CigarOperator.S));
                }
            }
            return new Cigar(newElements);
        }
    }

    private static List<SAMRecord> composeAlignedOutputHaplotypeSAMRecords(final SAMFileHeader header, final SVGenotypingContext context, final SVHaplotype haplotype, final RealignmentScoreParameters parameters) {
        final String variantUID = context.variant.getUniqueID();
        final AlignmentInterval primaryAlignment = haplotype.getReferenceAlignment().stream()
                .filter(ai -> !ai.cigarAlong5to3DirectionOfContig.containsOperator(CigarOperator.H))
                .min(Comparator.<AlignmentInterval>comparingInt(ai -> ai.mapQual == SAMRecord.UNKNOWN_MAPPING_QUALITY ? Integer.MAX_VALUE : -ai.mapQual)
                        .thenComparingInt(ai -> CigarUtils.countAlignedBases(ai.cigarAlong5to3DirectionOfContig))).orElse(null);
        return haplotype.getReferenceAlignment().stream()
                .map(ai -> {
                    final SAMRecord sam = ai.toSAMRecord(header, haplotype.getName(),
                            haplotype.getBases(),
                            ai != primaryAlignment,
                            ai != primaryAlignment ? SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue() : 0, null);
                    sam.setAttribute(VARIANT_CONTEXT_TAG, variantUID);
                    sam.setReadPairedFlag(false);
                    composeAndSetSATag(sam, haplotype);
                    if (haplotype.isContig()) {
                        annotateContigOutputSAMRecord(sam, (SVContig) haplotype);
                    } else {
                        annotateHaplotypeOutputSAMRecord(sam, haplotype, parameters);
                    }
                    return sam;
                })
                .collect(Collectors.toList());
    }

    private static void annotateContigOutputSAMRecord(final SAMRecord outputRecord, final SVContig contig) {
        outputRecord.setAttribute(SAMTag.RG.name(), CONTIG_READ_GROUP);
        final double refScore = contig.getReferenceHaplotypeScore();
        contig.getAlternativeHaplotypeScore();
        final double altScore = contig.getAlternativeHaplotypeScore();
        final double callQual = contig.getCallQuality();
        final String call = callQual == 0 ? Allele.NO_CALL_STRING : (refScore > altScore ? SVHaplotype.REF_HAPLOTYPE_NAME : SVHaplotype.ALT_HAPLOTYPE_NAME);
        outputRecord.setAttribute(HAPLOTYPE_CALL_TAG, call);
        outputRecord.setAttribute(HAPLOTYPE_QUAL_TAG, String.format("%.2f", callQual));
        outputRecord.setAttribute(REFERENCE_SCORE_TAG, contig.getReferenceScoreString());
        outputRecord.setAttribute(ALTERNATIVE_SCORE_TAG, contig.getAlternativeScoreString());
        outputRecord.setAttribute(REFERENCE_ALIGNMENT_TAG, AlignmentInterval.encode(contig.geReferenceHaplotypeAlignment()) + ';');
        outputRecord.setAttribute(ALTERNATIVE_ALIGNMENT_TAG, AlignmentInterval.encode(contig.getAlternativeHaplotypeAlignment()) + ';');
        outputRecord.setReadPairedFlag(false);
        outputRecord.setDuplicateReadFlag(false);
        outputRecord.setSecondOfPairFlag(false);
    }

    private static void annotateHaplotypeOutputSAMRecord(final SAMRecord outputRecord, final SVHaplotype haplotype, final RealignmentScoreParameters parameters) {
        outputRecord.setAttribute(SAMTag.RG.name(), HAPLOTYPE_READ_GROUP);
        final boolean isReference = haplotype.isReference();
        final double refScore = isReference ? 0 : parameters.interHaplotypePenalty;
        final double altScore = isReference ? parameters.interHaplotypePenalty : 0;
        final double callQual = parameters.interHaplotypePenalty;
        final String call = haplotype.getName();
        outputRecord.setAttribute(HAPLOTYPE_CALL_TAG, call);
        outputRecord.setAttribute(HAPLOTYPE_QUAL_TAG, String.format("%.2f", callQual));
        outputRecord.setAttribute(REFERENCE_SCORE_TAG, String.format("%.2f", refScore));
        outputRecord.setAttribute(ALTERNATIVE_SCORE_TAG, String.format("%.2f", altScore));
        outputRecord.setReadPairedFlag(false);
        outputRecord.setDuplicateReadFlag(false);
        outputRecord.setSecondOfPairFlag(false);
    }

    private static void composeAndSetSATag(SAMRecord outputRecord, SVHaplotype haplotype) {
        final List<AlignmentInterval> otherAlignments = haplotype.getReferenceAlignment().stream()
                .filter(ai -> ai.referenceSpan.getStart() != outputRecord.getAlignmentStart() ||
                              !Objects.equals(ai.referenceSpan.getContig(), outputRecord.getReferenceName()) ||
                              ai.referenceSpan.getEnd() != outputRecord.getAlignmentEnd() ||
                              !Objects.equals(ai.cigarAlong5to3DirectionOfContig, outputRecord.getCigar()) ||
                              ai.forwardStrand == outputRecord.getReadNegativeStrandFlag())
                .collect(Collectors.toList());
        if (!otherAlignments.isEmpty()) {
            outputRecord.setAttribute(SAMTag.SA.name(), AlignmentInterval.encode(otherAlignments) + ';');
        }
    }

    private static Genotype composeGenotype(final SVGenotypingContext context,
                                            final ReadLikelihoods<SVGenotypingContext.Allele> totalLikelihoods,
                                            final GenotypeLikelihoodCalculator genotypeCalculator,
                                            final ReadLikelihoods<SVGenotypingContext.Allele> splitsReadlikelihoods,
                                            final ReadLikelihoods<SVGenotypingContext.Allele> discordantOrientationLikelihoods,
                                            final ReadLikelihoods<SVGenotypingContext.Allele> insertSizeLikelihoods,
                                            final RealignmentScoreParameters penalties,
                                            final double informativePhredLikelihoodDifferenceThreshold,
                                            final boolean emitStratifiedAlleleDepths,
                                            final boolean emitStratifiedLikelihoods) {
        totalLikelihoods.removeUniformativeReads(0.0);
        totalLikelihoods.normalizeLikelihoods(-0.1 * penalties.maximumLikelihoodDifferencePerTemplate);
        final GenotypeLikelihoods likelihoods = genotypeCalculator.genotypeLikelihoods(totalLikelihoods.sampleMatrix(0));
        final int pl[] = likelihoods.getAsPLs();
        final int ad[] = calculateAD(totalLikelihoods, informativePhredLikelihoodDifferenceThreshold);
        final int gq = GATKVariantContextUtils.calculateGQFromPLs(pl);
        final List<Allele> genotypeAlleles = genotypeAlleleles(context, likelihoods, gq);
        final Map<String, Object> stratifiedAttributes = composeAdditionalGenotypeAttributes(emitStratifiedAlleleDepths, emitStratifiedLikelihoods, informativePhredLikelihoodDifferenceThreshold, genotypeCalculator, splitsReadlikelihoods, insertSizeLikelihoods, discordantOrientationLikelihoods);
        return new GenotypeBuilder().name(totalLikelihoods.samples().get(0))
                .PL(pl)
                .GQ(gq)
                .AD(ad)
                .DP((int) MathUtils.sum(ad))
                .attributes(stratifiedAttributes)
                .alleles(genotypeAlleles).make();
    }

    private static Map<String, Object> composeAdditionalGenotypeAttributes(final boolean emitStratifiedAlleleDepths, final boolean emitStratifiedLikelihoods, double informativePhredLikelihoodDifferenceThreshold, GenotypeLikelihoodCalculator genotypeCalculator, ReadLikelihoods<SVGenotypingContext.Allele> splitsReadlikelihoods, ReadLikelihoods<SVGenotypingContext.Allele> insertSizeLikelihoods, ReadLikelihoods<SVGenotypingContext.Allele> discordantOrientationLikelihoods) {
        if (!emitStratifiedAlleleDepths && !emitStratifiedLikelihoods) {
            return Collections.emptyMap();
        } else {
            splitsReadlikelihoods.removeUniformativeReads(0.0);
            discordantOrientationLikelihoods.removeUniformativeReads(0.0);
            insertSizeLikelihoods.removeUniformativeReads(0.0);
            final Map<String, Object> stratifiedAttributes = new LinkedHashMap<>((emitStratifiedAlleleDepths ? 3 : 0) + (emitStratifiedLikelihoods ? 3 : 0));

            if (emitStratifiedAlleleDepths) {
                stratifiedAttributes.put(GATKSVVCFConstants.TEMPLATE_MAPPING_ALLELE_DEPTH, calculateAD(splitsReadlikelihoods, informativePhredLikelihoodDifferenceThreshold));
                stratifiedAttributes.put(GATKSVVCFConstants.INSERT_SIZE_ALLELE_DEPTH, calculateAD(insertSizeLikelihoods, informativePhredLikelihoodDifferenceThreshold));
                stratifiedAttributes.put(GATKSVVCFConstants.DISCORDANT_PAIR_ORIENTATION_ALLELE_DEPTH, calculateAD(discordantOrientationLikelihoods, informativePhredLikelihoodDifferenceThreshold));
            }
            if (emitStratifiedLikelihoods) {
                stratifiedAttributes.put(GATKSVVCFConstants.TEMPLATE_MAPPING_LIKELIHOODS,
                        genotypeCalculator.genotypeLikelihoods(splitsReadlikelihoods.sampleMatrix(0)).getAsPLs());
                stratifiedAttributes.put(GATKSVVCFConstants.INSERT_SIZE_LIKELIHOODS,
                        genotypeCalculator.genotypeLikelihoods(insertSizeLikelihoods.sampleMatrix(0)).getAsPLs());
                stratifiedAttributes.put(GATKSVVCFConstants.DISCORDANT_PAIR_ORIENTATION_LIKELIHOODS,
                        genotypeCalculator.genotypeLikelihoods(discordantOrientationLikelihoods.sampleMatrix(0)).getAsPLs());
            }
            return stratifiedAttributes;
        }
    }

    private static List<Allele> genotypeAlleleles(final SVGenotypingContext context, final GenotypeLikelihoods likelihoods, final int gq) {
        final int bestGenotypeIndex = MathUtils.maxElementIndex(likelihoods.getAsVector());
        return gq == 0
                                    ? Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)
                                    : (bestGenotypeIndex == 0
                                    ? Collections.nCopies(2, context.refAllele)
                                    : ((bestGenotypeIndex == 1)
                                    ? Arrays.asList(context.refAllele, context.altAllele)
                                    : Collections.nCopies(2, context.altAllele)));
    }

    private static int[] calculateAD(final ReadLikelihoods<SVGenotypingContext.Allele> likelihoods, final double informativePhredLikelihoodDifferenceThreshold) {
        final int[] result = new int[likelihoods.numberOfAlleles()];
        final LikelihoodMatrix<SVGenotypingContext.Allele> matrix = likelihoods.sampleMatrix(0);
        final int numberOfTemplates = matrix.numberOfReads();
        final int numberOfAlleles = matrix.numberOfAlleles();
        final double log10Threshold = informativePhredLikelihoodDifferenceThreshold * -.1;
        for (int i = 0; i < numberOfTemplates; i++) {
            double bestLk = matrix.get(0, i);
            int bestLkIndex = 0;
            double secondBestLk = Double.NEGATIVE_INFINITY;
            for (int j = 1; j < numberOfAlleles; j++) {
                final double lk = matrix.get(j, i);
                if (lk > bestLk) {
                    secondBestLk = bestLk;
                    bestLk = lk;
                    bestLkIndex = j;
                } else if (lk > secondBestLk) {
                    secondBestLk = lk;
                }
            }
            if (bestLk - secondBestLk >= log10Threshold) {
                result[bestLkIndex]++;
            }
        }
        return result;
    }

    /**
     * Adjusts the expected relative frequency of alleles in the heterozygous genotype in the likelihood calculator.
     */
    private static void adjustLikelihoodCalculatorAlleleFrequencies(final SVGenotypingContext context, final InsertSizeDistribution insertSizeDistribution, final GenotypeLikelihoodCalculator genotypeCalculator) {
        if (context.variant.isInsertion()) {
            final int length = context.variant.getStructuralVariantLength();
            final double insertAverageSize = insertSizeDistribution.mean();
            genotypeCalculator.setRelativeAlleleFrequency(insertAverageSize, Math.min(2 * insertAverageSize, insertAverageSize + length));
        } else if (context.variant.isDeletion()) {
            final double insertAverageSize = insertSizeDistribution.mean();
            final int length = context.variant.getStructuralVariantLength();
            genotypeCalculator.setRelativeAlleleFrequency(Math.min(insertAverageSize * 2, insertAverageSize + length), insertAverageSize);
        } else {
            genotypeCalculator.resetRelativeAlleleFrequency(); //reset to default uniform.
        }
    }

    private static void setMissingAlignmentScores(final RealignmentScoreParameters penalties, final List<Template> templates,
                                                  final TemplateMappingTable scoreTable) {
        scoreTable.calculateBestMappingScores();
        for (int t = 0; t < templates.size(); t++) {
            final double worstFirstAlignmentScore = scoreTable.getWorstAlignmentScore(t, 0);
            final double firstMissingAlignmentScore = worstFirstAlignmentScore - 0.1 * penalties.unmappedFragmentPenalty;
            scoreTable.applyMissingAlignmentScore(t, 0, firstMissingAlignmentScore);
            final double worstSecondAlignmentScore = scoreTable.getWorstAlignmentScore(t, 1);
            final double secondMissingAlignmentScore = worstSecondAlignmentScore - 0.1 * penalties.unmappedFragmentPenalty;
            scoreTable.applyMissingAlignmentScore(t, 1, secondMissingAlignmentScore);
        }
    }

    /**
     * Often when we realign reads versus the reference haplotype, that is a small part of the full reference,
     * some of these reads map differently with full confidence (mq == 60). We shall consider these mappings not trust-worthy
     * and set the MQ to 0 so that they don't have effect on likelihoods in the end.
     * @param context the genotyping context.
     * @param scoreTable the mapping score table with the remap information.
     */
    private static void setMappingQualityToZeroForFragmentsThatMapDifferentlyOnFullReference(final SVGenotypingContext context, final TemplateMappingTable scoreTable) {
        for (int t = 0; t < context.numberOfTemplates; t++) {
            for (int f = 0; f < 2; f++) {
                final List<AlignmentInterval> newRefMapping = f == 0
                        ? scoreTable.getMappingInfo(context.refHaplotypeIndex, t).firstAlignmentIntervals
                        : scoreTable.getMappingInfo(context.refHaplotypeIndex, t).secondAlignmentIntervals;
                if (newRefMapping == null || newRefMapping.isEmpty()) {
                    continue;
                }
                final Template.Fragment fragment = context.templates.get(t).fragments().get(f);
                final List<AlignmentInterval> oldRefMapping = fragment.alignmentIntervals();
                if (oldRefMapping == null || oldRefMapping.isEmpty()) {
                    continue;
                }
                if (newRefMapping.size() != oldRefMapping.size()) {
                    context.templates.get(t).fragments().get(f).setMappingQuality(0);
                } else {
                    for (final AlignmentInterval newInterval : newRefMapping) {
                        if (oldRefMapping.stream().noneMatch(old ->
                                old.startInAssembledContig == newInterval.startInAssembledContig
                                        && old.endInAssembledContig == newInterval.endInAssembledContig
                                        && CigarUtils.equals(old.cigarAlong5to3DirectionOfContig, newInterval.cigarAlong5to3DirectionOfContig)
                                        && old.referenceSpan.getContig().equals(context.altHaplotype.getReferenceSpan().getContig())
                                        && old.referenceSpan.getStart() == context.refHaplotype.getReferenceSpan().getStart() + newInterval.referenceSpan.getStart() - 1)) {
                            context.templates.get(t).fragments().get(f).setMappingQuality(0);
                            break;
                        }
                    }
                }
            }
        }
    }

    private static TemplateMappingTable remapTemplatesOnHaplotypes(final SVGenotypingContext context, RealignmentScoreParameters realignmentScoreArguments) {
        final TemplateMappingTable scoreTable =
                new TemplateMappingTable(context.templates, context.haplotypes);
        final List<SVHaplotype> refAndAlthaplotypes = Arrays.asList(
                context.refHaplotype, context.altHaplotype
        );
        final List<SVContig> contigs = context.haplotypes.stream().filter(SVHaplotype::isContig)
                .map(SVContig.class::cast)
                .collect(Collectors.toList());
        final List<byte[]> contigSequences = contigs.stream().map(SVHaplotype::getBases).collect(Collectors.toList());
        final List<byte[]> templateSequences = context.templates.stream()
                .flatMap(template -> template.fragments().stream().map(Template.Fragment::bases))
                .collect(Collectors.toList());
        final int numberOfContigs = contigs.size();
        for (final SVHaplotype haplotype : refAndAlthaplotypes) {
            final List<byte[]> allSequences = new ArrayList<>(numberOfContigs + context.numberOfTemplates);
            allSequences.addAll(contigSequences);
            allSequences.addAll(templateSequences);
            final int haplotypeIndex = scoreTable.indexOf(haplotype);
                final List<List<AlignmentInterval>> intervals = alignSequenceAgainstHaplotype(context, haplotype, allSequences);
            annotateContigMappingOnHaplotypeInformation(realignmentScoreArguments, contigs, numberOfContigs, haplotype, intervals);
            final List<List<AlignmentInterval>> templateAlignmentIntervals = intervals.subList(numberOfContigs, intervals.size());
            annotateTemplateMappingInformation(context, realignmentScoreArguments, scoreTable, haplotype, haplotypeIndex, templateAlignmentIntervals);
        }
        for (final SVHaplotype contig : contigs) {
            final int haplotypeIndex = scoreTable.indexOf(contig);
            final List<List<AlignmentInterval>> intervals = alignSequenceAgainstHaplotype(context, contig, templateSequences);
            annotateTemplateMappingInformation(context, realignmentScoreArguments, scoreTable, contig, haplotypeIndex, intervals);
        }
        if (scoreTable.containsNulls()) {
            throw new IllegalStateException("bad");
        }
        return scoreTable;
    }

    private static void annotateContigMappingOnHaplotypeInformation(RealignmentScoreParameters realignmentScoreArguments, List<SVContig> contigs, int numberOfContigs, SVHaplotype haplotype, List<List<AlignmentInterval>> intervals) {
        for (int i = 0; i < numberOfContigs; i++) {
            final SVContig contig = contigs.get(i);
            final List<AlignmentInterval> alignment = intervals.get(i);
            final RealignmentScore score = RealignmentScore.calculate(realignmentScoreArguments, haplotype.getBases(), contig.getBases(), alignment);
            if (haplotype.isReference()) {
                contig.setReferenceHaplotypeAlignment(alignment, score);
            } else {
                contig.setAlternativeHaplotypeAlignment(alignment, score);
            }
        }
    }

    private static void annotateTemplateMappingInformation(SVGenotypingContext context, RealignmentScoreParameters realignmentScoreArguments, TemplateMappingTable scoreTable, SVHaplotype haplotype, int haplotypeIndex, List<List<AlignmentInterval>> intervals) {
        final int numberOfAlignments = intervals.size();
        for (int i = 0, templateIndex = 0; i < numberOfAlignments; templateIndex++) {
            final Template template = context.templates.get(templateIndex);
            final List<AlignmentInterval> firstIntervals = intervals.get(i++);
            final List<AlignmentInterval> secondIntervals = intervals.get(i++);
            final TemplateMapping mappingInformation = TemplateMapping.fromAlignments(realignmentScoreArguments, haplotype,
                    template.fragments().get(0).bases(), firstIntervals,
                    template.fragments().get(1).bases(), secondIntervals);
            scoreTable.setMapping(haplotypeIndex, templateIndex, mappingInformation);
        }
    }

    private static List<List<AlignmentInterval>> alignSequenceAgainstHaplotype(final SVGenotypingContext context,
                                                                               final SVHaplotype haplotype,
                                                                               final List<byte[]> sequences) {
        final IntFunction<String> haplotypeRefNameByIndex = i -> i == 0 ? haplotype.getName() : null;
        try (final TransientBwaMemIndex index = prepareIndex(context, haplotype);
             final BwaMemAligner aligner = new BwaMemAligner(index.get())) {
            return alignIntoIntervals(haplotypeRefNameByIndex, aligner, sequences);
        }
    }

    private static List<List<AlignmentInterval>> alignIntoIntervals(final IntFunction<String> haplotypeRefNameByIndex, BwaMemAligner aligner, final List<byte[]> sequences) {
        final int numberOfSequenes = sequences.size();
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(sequences);
        final List<List<AlignmentInterval>> intervals = new ArrayList<>(numberOfSequenes);
        for (int i = 0; i < numberOfSequenes; i++) {
            intervals.add(BwaMemAlignmentUtils.toAlignmentIntervals(alignments.get(i), haplotypeRefNameByIndex, sequences.get(i).length));
        }
        return intervals;
    }

    private static TransientBwaMemIndex prepareIndex(final SVGenotypingContext context, final SVHaplotype haplotype) {
        final String haplotypeName = haplotype.getName();
        final String variantID = context.variant.getUniqueID();
        final String imageName = haplotypeName.startsWith(variantID) ? haplotypeName : context.variant.getUniqueID() + "/" + haplotypeName;
        return new TransientBwaMemIndex(imageName, haplotype.getBases());
    }

    /**
     * Caculate the likelihoods associated with the way individual template read/fragments map to haplotypes and contigs.
     * <p>The only aspect that pair-information that is being used here is that both reads in a template are considered to be part
     * of the same haplotype and a mixture is not possible; insert-size or pair orientation are not relevant here</p>
     * @param context genotyping context.
     * @param realignmentScoreArguments realignment scoring parameters.
     * @param scoreTable the re-mapping score-table.
     * @return never {@code null}.
     */
    private static ReadLikelihoods<SVGenotypingContext.Allele> calculateSplitReadLikelihoods(final SVGenotypingContext context,
                                                                                             final RealignmentScoreParameters realignmentScoreArguments,
                                                                                             final TemplateMappingTable scoreTable) {

        final ReadLikelihoods<SVGenotypingContext.Allele> likelihoods = context.newLikelihoods();
        final ReadLikelihoods<SVGenotypingContext.Allele> likelihoodsFirst = context.newLikelihoods();
        final ReadLikelihoods<SVGenotypingContext.Allele> likelihoodsSecond = context.newLikelihoods();

        final LikelihoodMatrix<SVGenotypingContext.Allele> sampleLikelihoods = likelihoods.sampleMatrix(0);
        final LikelihoodMatrix<SVGenotypingContext.Allele> sampleLikelihoodsFirst = likelihoodsFirst.sampleMatrix(0);
        final LikelihoodMatrix<SVGenotypingContext.Allele> sampleLikelihoodsSecond = likelihoodsSecond.sampleMatrix(0);

        sampleLikelihoodsFirst.fill(Double.NEGATIVE_INFINITY);
        sampleLikelihoodsSecond.fill(Double.NEGATIVE_INFINITY);
        final List<SVContig> contigs = context.haplotypes.stream().filter(SVHaplotype::isNeitherReferenceNorAlternative).map(SVContig.class::cast)
                .collect(Collectors.toList());
        for (int t = 0; t < context.numberOfTemplates; t++) {
            sampleLikelihoodsFirst.set(context.refAlleleIndex, t,
                    scoreTable.getMappingInfo(context.refHaplotypeIndex, t).firstAlignmentScore.orElse(0));
            sampleLikelihoodsSecond.set(context.refAlleleIndex, t,
                    scoreTable.getMappingInfo(context.refHaplotypeIndex, t).secondAlignmentScore.orElse(0));
            sampleLikelihoodsFirst.set(context.altAlleleIndex, t,
                    scoreTable.getMappingInfo(context.altHaplotypeIndex, t).firstAlignmentScore.orElse(0));
            sampleLikelihoodsSecond.set(context.altAlleleIndex, t,
                    scoreTable.getMappingInfo(context.altHaplotypeIndex, t).secondAlignmentScore.orElse(0));
        }

        final Set<String> altContigNames = new HashSet<>(context.variant.getSupportingContigIds());

        for (final SVContig contig : contigs) {
            final int mappingInfoIndex = context.haplotypes.indexOf(contig);
            double haplotypeAltScore = RealignmentScore.calculate(realignmentScoreArguments, context.haplotypes.get(context.altHaplotypeIndex).getBases(), contig.getBases(), contig.getAlternativeHaplotypeAlignment()).getLog10Prob();
            double haplotypeRefScore = RealignmentScore.calculate(realignmentScoreArguments, context.haplotypes.get(context.refHaplotypeIndex).getBases(), contig.getBases(), contig.geReferenceHaplotypeAlignment()).getLog10Prob();
            if (altContigNames.contains(contig.getName())) {
                haplotypeAltScore = 0;
                haplotypeRefScore = -.1 * realignmentScoreArguments.interHaplotypePenalty;
            }
            final double maxMQ = contig.getCallQuality();
            final double base = Math.max(haplotypeAltScore, haplotypeRefScore);
            haplotypeAltScore -= base;
            haplotypeRefScore -= base;
            // we cap the difference of scores by the contig mapping quality.
            // so that contigs that could map in several places in the genome has less weight when assigning
            // genotypes.
            if (haplotypeAltScore < haplotypeRefScore) {
                haplotypeAltScore = Math.max(haplotypeAltScore, haplotypeRefScore - 0.1 * maxMQ);
            } else {
                haplotypeRefScore = Math.max(haplotypeRefScore, haplotypeAltScore - 0.1 * maxMQ);
            }
            if (haplotypeRefScore == haplotypeAltScore) {
                haplotypeAltScore = haplotypeRefScore = Double.NEGATIVE_INFINITY;
            }

            // for each template we apply the scores thru the haplotype/contig `c`. We reduce/marginalize the likelihood
            // to take the maximum across all haplotype/contigs for that template.
            for (int t = 0; t < context.numberOfTemplates; t++) {
                final boolean noAlignment = !scoreTable.getMappingInfo(mappingInfoIndex, t).firstAlignmentScore.isPresent()
                        && !scoreTable.getMappingInfo(mappingInfoIndex, t).secondAlignmentScore.isPresent();
                if (noAlignment) continue;
                final double firstMappingScore = scoreTable.getMappingInfo(mappingInfoIndex, t).firstAlignmentScore.orElse(Double.NEGATIVE_INFINITY);
                final double secondMappingScore = scoreTable.getMappingInfo(mappingInfoIndex, t).secondAlignmentScore.orElse(Double.NEGATIVE_INFINITY);
                sampleLikelihoodsFirst.set(context.refAlleleIndex, t,
                        Math.max(firstMappingScore + haplotypeRefScore, sampleLikelihoodsFirst.get(context.refAlleleIndex, t)));
                sampleLikelihoodsFirst.set(context.altAlleleIndex, t,
                        Math.max(firstMappingScore + haplotypeAltScore, sampleLikelihoodsFirst.get(context.altAlleleIndex, t)));
                sampleLikelihoodsSecond.set(context.refAlleleIndex, t,
                        Math.max(secondMappingScore + haplotypeRefScore, sampleLikelihoodsSecond.get(context.refAlleleIndex, t)));
                sampleLikelihoodsSecond.set(context.altAlleleIndex, t,
                        Math.max(secondMappingScore + haplotypeAltScore, sampleLikelihoodsSecond.get(context.altAlleleIndex, t)));
            }
        }
        // we cap the likelihood difference for each allele haplotypes for each fragment in each template by th
        // mapping quality of the fragment so that those reads that may map to other locations in the genome count less
        // towards genotyping.
        for (int t = 0; t < context.numberOfTemplates; t++) {
            for (int f = 0; f < 2; f++) {
                final int maxMq = context.templates.get(t).fragments().get(f).getMappingQuality();
                final LikelihoodMatrix<SVGenotypingContext.Allele> matrix = f == 0 ? sampleLikelihoodsFirst : sampleLikelihoodsSecond;
                final double base = Math.max(matrix.get(context.refAlleleIndex, t), matrix.get(context.altAlleleIndex, t));
                final int maxIndex = matrix.get(context.refAlleleIndex, t) == base ? context.refAlleleIndex : context.altAlleleIndex;
                final int minIndex = maxIndex == context.refAlleleIndex ? context.altAlleleIndex : context.refAlleleIndex;
                matrix.set(minIndex, t, Math.max(matrix.get(maxIndex, t) - 0.1 * maxMq, matrix.get(minIndex, t)));
            }
        }

        for (int j = 0; j < sampleLikelihoods.numberOfReads(); j++) {
            final boolean considerFirstFragment = scoreTable.getMappingInfo(context.refHaplotypeIndex, j).crossesBreakPoint(context.refBreakPoints)
                    || scoreTable.getMappingInfo(context.altHaplotypeIndex, j).crossesBreakPoint(context.altBreakPoints);
            final boolean considerSecondFragment = scoreTable.getMappingInfo(context.refHaplotypeIndex, j).crossesBreakPoint(context.refBreakPoints)
                    || scoreTable.getMappingInfo(context.altHaplotypeIndex, j).crossesBreakPoint(context.altBreakPoints);

            //TODO  the following adjustment seem to increase accuracy but the maths are a bit sketchy so for now I
            //include it.
            // The idea is to penalize poorly modeled reads... those that allele haplotype likelihood difference is not
            // much larger than the magnitude of the likelihood of the best allele likelihood.
            // so for example a template fragment has lks Phred 3000 3100, there is 100 diff in Lk but is very small
            // when compare to how unlikely the read is even with the best haplotype (3000) so is 100 even relevant?
            // The code simply substract 3000 to those 100 (min 0) so in the end such a read to be of any value for Lk
            // calculation the worst allele lk would need to be at least 6000.
            //
            // at first may seem a reasonable thing to do but it it unclear how to weight down these reads... just take
            // the best lk seems arbitrary.
            //
            // This might be resolved differently by ignoring short variant differences so that the Lk is totally
            // defined by real SV variation and so the 3000 becomes closer to 0. Or that in general the 100 is actual
            // real SV differences rather that differences due to chance in the distribution of short variants between
            // contigs.
            //if (sampleLikelihoodsFirst == sampleLikelihoodsSecond) debug00(sampleLikelihoodsFirst, sampleLikelihoodsSecond);
            for (int k = 0; k < 2; k++) {
                final LikelihoodMatrix<SVGenotypingContext.Allele> matrix = k == 0 ? sampleLikelihoodsFirst : sampleLikelihoodsSecond;
                final double base = Math.max(matrix.get(context.refAlleleIndex, j), matrix.get(context.altAlleleIndex, j));
                final int maxIndex = matrix.get(context.refAlleleIndex, j) == base ? context.refAlleleIndex : context.altAlleleIndex;
                final int minIndex = maxIndex == context.refAlleleIndex ? context.altAlleleIndex : context.refAlleleIndex;
                matrix.set(minIndex, j, Math.min(matrix.get(maxIndex, j), matrix.get(minIndex, j) - scoreTable.getBestAlignmentScore(j, k)));
            }
            for (int i = 0; i < sampleLikelihoods.numberOfAlleles(); i++) {
                sampleLikelihoods.set(i, j, (considerFirstFragment ? sampleLikelihoodsFirst.get(i, j) : 0) +
                        ((considerSecondFragment) ? sampleLikelihoodsSecond.get(i, j) : 0));

            }
        }
        return likelihoods;
    }

    private static ReadLikelihoods<SVGenotypingContext.Allele> calculateDiscordantOrientationLikelihoods(final SVGenotypingContext context,
                                                                                                         final RealignmentScoreParameters penalties, final TemplateMappingTable scoreTable) {
        final ReadLikelihoods<SVGenotypingContext.Allele> result = context.newLikelihoods();
        final LikelihoodMatrix<SVGenotypingContext.Allele> matrix = result.sampleMatrix(0);
        for (int t = 0; t < context.numberOfTemplates; t++) {
            final TemplateMapping referenceMappingInfo = scoreTable.getMappingInfo(context.refHaplotypeIndex, t);
            final TemplateMapping alternativeMappingInfo = scoreTable.getMappingInfo(context.altHaplotypeIndex, t);
            if (!referenceMappingInfo.pairOrientation.isDefined()) continue;
            if (!alternativeMappingInfo.pairOrientation.isDefined()) continue;
            if (referenceMappingInfo.pairOrientation.isProper() == alternativeMappingInfo.pairOrientation.isProper())
                continue;
            if (referenceMappingInfo.pairOrientation.isProper()) { // && referenceMappingInfo.crossesBreakPoint(context.refBreakPoints)) {
                matrix.set(context.altAlleleIndex, t, penalties.improperPairPenalty * -.1);
            } else { // this must be true: alternativeMappingInfo.pairOrientation.isProper()
                matrix.set(context.refAlleleIndex, t, penalties.improperPairPenalty * -.1);
            }
        }
        return result;
    }

    private static ReadLikelihoods<SVGenotypingContext.Allele> calculateInsertSizeLikelihoods(final SVGenotypingContext context,
                                                                                              final RealignmentScoreParameters penalties,
                                                                                              final TemplateMappingTable scoreTable,
                                                                                              final InsertSizeDistribution dist) {
        final ReadLikelihoods<SVGenotypingContext.Allele> result = context.newLikelihoods();
        final LikelihoodMatrix<SVGenotypingContext.Allele> matrix = result.sampleMatrix(0);
        final double inverseOfLog10 = 1.0 / MathUtil.LOG_10_MATH.getLog_of_base();
        for (int t = 0; t < context.numberOfTemplates; t++) {
            final TemplateMapping referenceMappingInfo = scoreTable.getMappingInfo(context.refHaplotypeIndex, t);
            final TemplateMapping alternativeMappingInfo = scoreTable.getMappingInfo(context.altHaplotypeIndex, t);
            if (!referenceMappingInfo.pairOrientation.isProper()) continue;
            if (!alternativeMappingInfo.pairOrientation.isProper()) continue;
            final boolean accrossBreakPointsOnRef = referenceMappingInfo.crossesBreakPoint(context.refBreakPoints);
            final boolean accrossBreakPointsOnAlt = alternativeMappingInfo.crossesBreakPoint(context.altBreakPoints);
            if (accrossBreakPointsOnAlt || accrossBreakPointsOnRef) {
                @SuppressWarnings("ConstantConditions")
                final double refInsertSizeLk = dist.logProbability(referenceMappingInfo.insertSize.getAsInt()) * inverseOfLog10;
                @SuppressWarnings("ConstantConditions")
                final double altInsertSizeLk = dist.logProbability(alternativeMappingInfo.insertSize.getAsInt()) * inverseOfLog10;
                final double best = Math.max(refInsertSizeLk, altInsertSizeLk);
                final double worst = best + penalties.improperPairPenalty * -.1;
                matrix.set(context.refAlleleIndex, t, Math.max(worst, refInsertSizeLk));
                matrix.set(context.altAlleleIndex, t, Math.max(worst, altInsertSizeLk));
            }
        }
        return result;
    }

    private static boolean structuralVariantAlleleIsSupported(final SVContext ctx) {
        switch (ctx.getStructuralVariantType()) {
            case INS:
            case DEL:
            case DUP:
                return true;
            //case INV:
                //return !ctx.hasAttribute(GATKSVVCFConstants.IMPRECISE);
            default:
                return false;
        }
    }

    private void tearDown(@SuppressWarnings("unused") final JavaSparkContext ctx) {
        variantsSource = null;
    }
}
