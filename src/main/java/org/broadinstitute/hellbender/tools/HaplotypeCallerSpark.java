package org.broadinstitute.hellbender.tools;

import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.SparkFiles;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.AssemblyRegionEvaluator;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.*;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSink;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Supplier;

/**
 * ********************************************************************************
 * *   This tool DOES NOT match the output of HaplotypeCaller.                    *
 * *   It is still under development and should not be used for production work.  *
 * *   For evaluation only.                                                       *
 * *   Use the non-spark HaplotypeCaller if you care about the results.           *
 * ********************************************************************************
 *
 * Call germline SNPs and indels via local re-assembly of haplotypes.
 *
 * <p>This is an implementation of {@link HaplotypeCaller} using spark to distribute the computation.
 * It is still in an early stage of development and does not yet support all the options that the non-spark version does.
 * Specifically it does not support the --dbsnp, --comp, and --bam-output options.</p>
 *
 * <h3>Usage Example</h3>
 * <pre>
 * gatk HaplotypeCallerSpark \
 * -R Homo_sapiens_assembly38.fasta \
 * -I input.bam \
 * -O output.vcf.gz
 * </pre>
 *
 */
@CommandLineProgramProperties(summary = "HaplotypeCaller on Spark", oneLineSummary = "HaplotypeCaller on Spark", programGroup = ShortVariantDiscoveryProgramGroup.class)
@DocumentedFeature
@BetaFeature
public final class HaplotypeCallerSpark extends AssemblyRegionWalkerSpark {
    private static final long serialVersionUID = 1L;

    public static final int DEFAULT_READSHARD_SIZE = 5000;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Single file to which variants should be written")
    public String output;

    public static class HaplotypeCallerAssemblyRegionArgumentCollection extends AssemblyRegionArgumentCollection {
        private static final long serialVersionUID = 1L;

        @Override
        protected int defaultMinAssemblyRegionSize() { return HaplotypeCaller.DEFAULT_MIN_ASSEMBLY_REGION_SIZE; }

        @Override
        protected int defaultMaxAssemblyRegionSize() { return HaplotypeCaller.DEFAULT_MAX_ASSEMBLY_REGION_SIZE; }

        @Override
        protected int defaultAssemblyRegionPadding() { return HaplotypeCaller.DEFAULT_ASSEMBLY_REGION_PADDING; }

        @Override
        protected int defaultMaxReadsPerAlignmentStart() { return HaplotypeCaller.DEFAULT_MAX_READS_PER_ALIGNMENT; }

        @Override
        protected double defaultActiveProbThreshold() { return HaplotypeCaller.DEFAULT_ACTIVE_PROB_THRESHOLD; }

        @Override
        protected int defaultMaxProbPropagationDistance() { return HaplotypeCaller.DEFAULT_MAX_PROB_PROPAGATION_DISTANCE; }
    }

    @Override
    protected AssemblyRegionArgumentCollection getAssemblyRegionArgumentCollection() {
        return new HaplotypeCallerAssemblyRegionArgumentCollection();
    }

    @ArgumentCollection
    public HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();

    @Override
    protected boolean includeReadsWithDeletionsInIsActivePileups() {
        return true;
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
    protected void runTool(JavaSparkContext ctx) {
        //TODO remove me when https://github.com/broadinstitute/gatk/issues/4303 are fixed
        if (output.endsWith(IOUtil.BCF_FILE_EXTENSION) || output.endsWith(IOUtil.BCF_FILE_EXTENSION + ".gz")) {
            throw new UserException.UnimplementedFeature("It is currently not possible to write a BCF file on spark.  See https://github.com/broadinstitute/gatk/issues/4303 for more details .");
        }
        Utils.validateArg(hcArgs.dbsnp.dbsnp == null, "HaplotypeCallerSpark does not yet support -D or --dbsnp arguments" );
        Utils.validateArg(hcArgs.comps.isEmpty(), "HaplotypeCallerSpark does not yet support -comp or --comp arguments" );
        Utils.validateArg(hcArgs.bamOutputPath == null, "HaplotypeCallerSpark does not yet support -bamout or --bamOutput");

        Utils.validate(getHeaderForReads().getSortOrder() == SAMFileHeader.SortOrder.coordinate, "The reads must be coordinate sorted.");
        logger.info("********************************************************************************");
        logger.info("The output of this tool DOES NOT match the output of HaplotypeCaller. ");
        logger.info("It is under development and should not be used for production work. ");
        logger.info("For evaluation only.");
        logger.info("Use the non-spark HaplotypeCaller if you care about the results. ");
        logger.info("********************************************************************************");
        try {
            super.runTool(ctx);
        } catch (Exception e) {
            if (e.getCause() instanceof UserException) {
                throw (UserException) e.getCause();
            } else {
                throw e;
            }
        }
    }

    @Override
    protected void processAssemblyRegions(JavaRDD<AssemblyRegionWalkerContext> rdd, JavaSparkContext ctx) {
        processAssemblyRegions(rdd, ctx, getHeaderForReads(), referenceArguments.getReferenceFileName(), hcArgs, output, makeVariantAnnotations(), logger);
    }

    private static void processAssemblyRegions(
            JavaRDD<AssemblyRegionWalkerContext> rdd,
            final JavaSparkContext ctx,
            final SAMFileHeader header,
            final String reference,
            final HaplotypeCallerArgumentCollection hcArgs,
            final String output,
            final Collection<Annotation> annotations,
            final Logger logger) {

        final VariantAnnotatorEngine variantannotatorEngine = new VariantAnnotatorEngine(annotations,  hcArgs.dbsnp.dbsnp, hcArgs.comps, hcArgs.emitReferenceConfidence != ReferenceConfidenceMode.NONE);

        final Path referencePath = IOUtils.getPath(reference);
        final ReferenceSequenceFile driverReferenceSequenceFile = new CachingIndexedFastaSequenceFile(referencePath);
        final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgs, false, false, header, driverReferenceSequenceFile, variantannotatorEngine);
        final String referenceFileName = referencePath.getFileName().toString();
        final Broadcast<HaplotypeCallerArgumentCollection> hcArgsBroadcast = ctx.broadcast(hcArgs);
        final Broadcast<VariantAnnotatorEngine> annotatorEngineBroadcast = ctx.broadcast(variantannotatorEngine);

        final JavaRDD<VariantContext> variants = rdd.mapPartitions(assemblyFunction(header, referenceFileName, hcArgsBroadcast, annotatorEngineBroadcast));

        variants.cache(); // without caching, computations are run twice as a side effect of finding partition boundaries for sorting
        try {
            VariantsSparkSink.writeVariants(ctx, output, variants, hcEngine.makeVCFHeader(header.getSequenceDictionary(), new HashSet<>()),
                    hcArgs.emitReferenceConfidence == ReferenceConfidenceMode.GVCF, new ArrayList<Number>(hcArgs.GVCFGQBands), hcArgs.genotypeArgs.samplePloidy);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(output, "writing failed", e);
        }
    }

    private static FlatMapFunction<Iterator<AssemblyRegionWalkerContext>, VariantContext> assemblyFunction(final SAMFileHeader header,
                                                                                                           final String referenceFileName,
                                                                                                           final Broadcast<HaplotypeCallerArgumentCollection> hcArgsBroadcast,
                                                                                                           final Broadcast<VariantAnnotatorEngine> annotatorEngineBroadcast) {
        return (FlatMapFunction<Iterator<AssemblyRegionWalkerContext>, VariantContext>) contexts -> {
            // HaplotypeCallerEngine isn't serializable but is expensive to instantiate, so construct and reuse one for every partition
            final ReferenceSequenceFile taskReferenceSequenceFile = taskReferenceSequenceFile(referenceFileName);
            final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgsBroadcast.value(), false, false, header, taskReferenceSequenceFile, annotatorEngineBroadcast.getValue());
            Iterator<Iterator<VariantContext>> iterators = Utils.stream(contexts).map(context -> {
                AssemblyRegion region = context.getAssemblyRegion();
                FeatureContext featureContext = context.getFeatureContext();
                return hcEngine.callRegion(region, featureContext).iterator();
            }).iterator();

            return Iterators.concat(iterators);
        };
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return HaplotypeCallerEngine.makeStandardHCReadFilters();
    }

    @Override
    public AssemblyRegionEvaluator assemblyRegionEvaluator() {
        return null; // not used (see assemblyRegionEvaluatorSupplierBroadcast)
    }

    @Override
    protected Broadcast<Supplier<AssemblyRegionEvaluator>> assemblyRegionEvaluatorSupplierBroadcast(final JavaSparkContext ctx) {
        final Path referencePath = IOUtils.getPath(referenceArguments.getReferenceFileName());
        final String referenceFileName = referencePath.getFileName().toString();
        final String pathOnExecutor = SparkFiles.get(referenceFileName);
        final ReferenceSequenceFile taskReferenceSequenceFile = new CachingIndexedFastaSequenceFile(IOUtils.getPath(pathOnExecutor));
        final Collection<Annotation> annotations = makeVariantAnnotations();
        final VariantAnnotatorEngine annotatorEngine = new VariantAnnotatorEngine(annotations,  hcArgs.dbsnp.dbsnp, hcArgs.comps, hcArgs.emitReferenceConfidence != ReferenceConfidenceMode.NONE);
        return assemblyRegionEvaluatorSupplierBroadcastFunction(ctx, hcArgs, getHeaderForReads(), taskReferenceSequenceFile, annotatorEngine);
    }

    private static Broadcast<Supplier<AssemblyRegionEvaluator>> assemblyRegionEvaluatorSupplierBroadcast(
            final JavaSparkContext ctx,
            final HaplotypeCallerArgumentCollection hcArgs,
            final SAMFileHeader header,
            final String reference,
            final Collection<Annotation> annotations) {
        final Path referencePath = IOUtils.getPath(reference);
        final String referenceFileName = referencePath.getFileName().toString();
        final ReferenceSequenceFile taskReferenceSequenceFile = taskReferenceSequenceFile(referenceFileName);
        final VariantAnnotatorEngine annotatorEngine = new VariantAnnotatorEngine(annotations,  hcArgs.dbsnp.dbsnp, hcArgs.comps, hcArgs.emitReferenceConfidence != ReferenceConfidenceMode.NONE);
        return assemblyRegionEvaluatorSupplierBroadcastFunction(ctx, hcArgs, header, taskReferenceSequenceFile, annotatorEngine);
    }

    private static ReferenceSequenceFile taskReferenceSequenceFile(final String referenceFileName) {
        final String pathOnExecutor = SparkFiles.get(referenceFileName);
        return new CachingIndexedFastaSequenceFile(IOUtils.getPath(pathOnExecutor));
    }

    private static Broadcast<Supplier<AssemblyRegionEvaluator>> assemblyRegionEvaluatorSupplierBroadcastFunction(
            final JavaSparkContext ctx,
            final HaplotypeCallerArgumentCollection hcArgs,
            final SAMFileHeader header,
            final ReferenceSequenceFile taskReferenceSequenceFile,
            final VariantAnnotatorEngine annotatorEngine) {
        Supplier<AssemblyRegionEvaluator> supplier = new Supplier<AssemblyRegionEvaluator>() {
            @Override
            public AssemblyRegionEvaluator get() {
                return new HaplotypeCallerEngine(hcArgs, false, false, header, taskReferenceSequenceFile, annotatorEngine);
            }
        };
        return ctx.broadcast(supplier);
    }

    /**
     * Call Variants using HaplotypeCaller on Spark and write out a VCF file.
     *
     * This may be called from any spark pipeline in order to call variants from an RDD of GATKRead
     * @param ctx the spark context
     * @param reads the reads variants should be called from
     * @param header the header that goes with the reads
     * @param reference the full path to the reference file (must have been added via {@code SparkContext#addFile()})
     * @param intervalShards the interval shards to restrict calling to
     * @param hcArgs haplotype caller arguments
     * @param shardingArgs arguments to control how the assembly regions are sharded
     * @param output the output path for the VCF
     * @param logger
     * @param strict whether to use the strict implementation (slower) for finding assembly regions to match walker version
     */
    public static void callVariantsWithHaplotypeCallerAndWriteOutput(
            final JavaSparkContext ctx,
            final JavaRDD<GATKRead> reads,
            final SAMFileHeader header,
            final SAMSequenceDictionary sequenceDictionary,
            final String reference,
            final List<ShardBoundary> intervalShards,
            final HaplotypeCallerArgumentCollection hcArgs,
            final AssemblyRegionReadShardArgumentCollection shardingArgs,
            final AssemblyRegionArgumentCollection assemblyRegionArgs,
            final boolean includeReadsWithDeletionsInIsActivePileups,
            final String output,
            final Collection<Annotation> annotations,
            final Logger logger,
            final boolean strict) {

        final Path referencePath = IOUtils.getPath(reference);
        final String referenceFileName = referencePath.getFileName().toString();
        Broadcast<Supplier<AssemblyRegionEvaluator>> assemblyRegionEvaluatorSupplierBroadcast = assemblyRegionEvaluatorSupplierBroadcast(ctx, hcArgs, header, reference, annotations);
        JavaRDD<AssemblyRegionWalkerContext> assemblyRegions = strict ?
                FindAssemblyRegionsSpark.getAssemblyRegionsStrict(ctx, reads, header, sequenceDictionary, referenceFileName, null, intervalShards, assemblyRegionEvaluatorSupplierBroadcast, shardingArgs, assemblyRegionArgs, includeReadsWithDeletionsInIsActivePileups, false) :
                FindAssemblyRegionsSpark.getAssemblyRegionsFast(ctx, reads, header, sequenceDictionary, referenceFileName, null, intervalShards, assemblyRegionEvaluatorSupplierBroadcast, shardingArgs, assemblyRegionArgs, includeReadsWithDeletionsInIsActivePileups, false);
        processAssemblyRegions(assemblyRegions, ctx, header, reference, hcArgs, output, annotations, logger);
    }
}
