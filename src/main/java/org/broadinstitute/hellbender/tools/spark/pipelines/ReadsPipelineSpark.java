package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.storage.StorageLevel;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.AddContextDataToReadSpark;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.JoinStrategy;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.ApplyBQSRUniqueArgumentCollection;
import org.broadinstitute.hellbender.tools.HaplotypeCallerSpark;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.tools.spark.transforms.ApplyBQSRSparkFn;
import org.broadinstitute.hellbender.tools.spark.transforms.BaseRecalibratorSparkFn;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.SerializableOpticalDuplicatesFinder;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;

import java.util.Collection;
import java.util.List;

/**
 * ReadsPipelineSpark is our standard pipeline that takes unaligned or aligned reads and runs BWA (if specified), MarkDuplicates,
 * BQSR, and HaplotypeCaller. The final result is analysis-ready variants.
 *
 *
 * <h3>Examples</h3>
 * <pre>
 * gatk ReadsPipelineSpark \
 *   -I gs://my-gcs-bucket/aligned_reads.bam \
 *   -R gs://my-gcs-bucket/reference.fasta \
 *   --known-sites gs://my-gcs-bucket/sites_of_variation.vcf \
 *   -O gs://my-gcs-bucket/output.vcf \
 *   -- \
 *   --sparkRunner GCS \
 *   --cluster my-dataproc-cluster
 * </pre>
 * <p>
 * To additionally align reads with BWA-MEM:
 * </p>
 * <pre>
 * gatk ReadsPipelineSpark \
 *   -I gs://my-gcs-bucket/unaligned_reads.bam \
 *   -R gs://my-gcs-bucket/reference.fasta \
 *   --known-sites gs://my-gcs-bucket/sites_of_variation.vcf \
 *   --align
 *   -O gs://my-gcs-bucket/output.vcf \
 *   -- \
 *   --sparkRunner GCS \
 *   --cluster my-dataproc-cluster
 * </pre>
 */

@CommandLineProgramProperties(
        summary = ReadsPipelineSpark.USAGE_SUMMARY,
        oneLineSummary = ReadsPipelineSpark.USAGE_ONE_LINE_SUMMARY,
        programGroup = ShortVariantDiscoveryProgramGroup.class
)

@DocumentedFeature
@BetaFeature
public class ReadsPipelineSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    static final String USAGE_ONE_LINE_SUMMARY = "Takes unaligned or aligned reads and runs BWA (if specified), MarkDuplicates, BQSR, and HaplotypeCaller to generate a VCF file of variants";
    static final String USAGE_SUMMARY = "Takes unaligned or aligned reads and runs BWA (if specified), MarkDuplicates, BQSR, and HaplotypeCaller. The final result is analysis-ready variants.";

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public boolean requiresReference() { return true; }

    @Argument(doc = "whether to perform alignment using BWA-MEM", fullName = "align", optional = true)
    private boolean align;

    @Argument(doc = "the known variants", fullName = BaseRecalibrator.KNOWN_SITES_ARG_FULL_NAME, optional = false)
    protected List<String> baseRecalibrationKnownVariants;

    @Argument(doc = "the output vcf", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String output;

    @Argument(doc = "the output bam", fullName = "output-bam", optional = true)
    protected String outputBam;

    @Argument(doc = "the join strategy for reference bases and known variants", fullName = "join-strategy", optional = true)
    private JoinStrategy joinStrategy = JoinStrategy.BROADCAST;

    @ArgumentCollection
    protected MarkDuplicatesSparkArgumentCollection markDuplicatesSparkArgumentCollection = new MarkDuplicatesSparkArgumentCollection();

    @ArgumentCollection
    public final BwaArgumentCollection bwaArgs = new BwaArgumentCollection();

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc = "all the command line arguments for BQSR and its covariates")
    private final RecalibrationArgumentCollection bqsrArgs = new RecalibrationArgumentCollection();

    @ArgumentCollection
    public final HaplotypeCallerSpark.ShardingArgumentCollection shardingArgs = new HaplotypeCallerSpark.ShardingArgumentCollection();

    /**
     * command-line arguments to fine tune the apply BQSR step.
     */
    @ArgumentCollection
    public ApplyBQSRUniqueArgumentCollection applyBqsrArgs = new ApplyBQSRUniqueArgumentCollection();

    @ArgumentCollection
    public HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();

    @Override
    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return BaseRecalibrationEngine.BQSR_REFERENCE_WINDOW_FUNCTION;
    }

    @Override
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() { return HaplotypeCallerEngine.getStandardHaplotypeCallerAnnotationGroups();}

    @Override
    public Collection<Annotation> makeVariantAnnotations() {
        final boolean referenceConfidenceMode = hcArgs.emitReferenceConfidence != ReferenceConfidenceMode.NONE;
        final Collection<Annotation> annotations = super.makeVariantAnnotations();
        return referenceConfidenceMode? HaplotypeCallerEngine.filterReferenceConfidenceAnnotations(annotations): annotations;
    }

    @Override
    protected void validateSequenceDictionaries(){
        //don't validate unaligned reads because we don't require them to have a sequence dictionary
        if( !align ){
            super.validateSequenceDictionaries();
        }
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        if (joinStrategy == JoinStrategy.BROADCAST && ! getReference().isCompatibleWithSparkBroadcast()){
            throw new UserException.Require2BitReferenceForBroadcast();
        }

        final JavaRDD<GATKRead> alignedReads;
        final SAMFileHeader header;
        final BwaSparkEngine bwaEngine;
        if (align) {
            bwaEngine = new BwaSparkEngine(ctx, referenceArguments.getReferenceFileName(), bwaArgs.indexImageFile, getHeaderForReads(), getReferenceSequenceDictionary());
            if (bwaArgs.singleEndAlignment) {
                alignedReads = bwaEngine.alignUnpaired(getReads());
            } else {
                // filter reads after alignment in the case of paired reads since filtering does not know about pairs
                final ReadFilter filter = makeReadFilter(bwaEngine.getHeader());
                alignedReads = bwaEngine.alignPaired(getUnfilteredReads()).filter(filter::test);
            }
            header = bwaEngine.getHeader();
        } else {
            bwaEngine = null;
            alignedReads = getReads();
            header = getHeaderForReads();
        }

        final JavaRDD<GATKRead> markedReads = MarkDuplicatesSpark.mark(alignedReads, header, markDuplicatesSparkArgumentCollection.duplicatesScoringStrategy, new SerializableOpticalDuplicatesFinder(), getRecommendedNumReducers(), markDuplicatesSparkArgumentCollection.dontMarkUnmappedMates);

        // The markedReads have already had the WellformedReadFilter applied to them, which
        // is all the filtering that MarkDupes and ApplyBQSR want. BQSR itself wants additional
        // filtering performed, so we do that here.
        //NOTE: this doesn't honor enabled/disabled commandline filters
        final ReadFilter bqsrReadFilter = ReadFilter.fromList(BaseRecalibrator.getBQSRSpecificReadFilterList(), header);

        JavaRDD<GATKRead> markedFilteredReadsForBQSR = markedReads.filter(read -> bqsrReadFilter.test(read));

        if (joinStrategy.equals(JoinStrategy.OVERLAPS_PARTITIONER)) {
            // the overlaps partitioner requires that reads are coordinate-sorted
            final SAMFileHeader readsHeader = header.clone();
            readsHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            markedFilteredReadsForBQSR = SparkUtils.sortReadsAccordingToHeader(markedFilteredReadsForBQSR, readsHeader, numReducers);
        }

        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        JavaRDD<GATKVariant> bqsrKnownVariants = variantsSparkSource.getParallelVariants(baseRecalibrationKnownVariants, getIntervals());

        JavaPairRDD<GATKRead, ReadContextData> rddReadContext = AddContextDataToReadSpark.add(ctx, markedFilteredReadsForBQSR, getReference(), bqsrKnownVariants, baseRecalibrationKnownVariants, joinStrategy, header.getSequenceDictionary(), shardingArgs.readShardSize, shardingArgs.readShardPadding);
        final RecalibrationReport bqsrReport = BaseRecalibratorSparkFn.apply(rddReadContext, header, getReferenceSequenceDictionary(), bqsrArgs);

        final Broadcast<RecalibrationReport> reportBroadcast = ctx.broadcast(bqsrReport);
        final JavaRDD<GATKRead> finalReads = ApplyBQSRSparkFn.apply(markedReads, reportBroadcast, header, applyBqsrArgs.toApplyBQSRArgumentCollection(bqsrArgs.PRESERVE_QSCORES_LESS_THAN));

        if (outputBam != null) { // only write output of BQSR if output BAM is specified
            writeReads(ctx, outputBam, finalReads, header);
        }

        // Run Haplotype Caller
        final ReadFilter hcReadFilter = ReadFilter.fromList(HaplotypeCallerEngine.makeStandardHCReadFilters(), header);
        final JavaRDD<GATKRead> filteredReadsForHC = finalReads.filter(read -> hcReadFilter.test(read));
        filteredReadsForHC.persist(StorageLevel.DISK_ONLY()); // without caching, computations are run twice as a side effect of finding partition boundaries for sorting
        final List<SimpleInterval> intervals = hasUserSuppliedIntervals() ? getIntervals() : IntervalUtils.getAllIntervalsForReference(header.getSequenceDictionary());
        HaplotypeCallerSpark.callVariantsWithHaplotypeCallerAndWriteOutput(ctx, filteredReadsForHC, header, getReference(), intervals, hcArgs, shardingArgs, numReducers, output, makeVariantAnnotations());

        if (bwaEngine != null) {
            bwaEngine.close();
        }
    }
}
