package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.storage.StorageLevel;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalLocator;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SerializableComparator;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.gcs.BamBucketIoUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Composes a bam file with assembled contigs conveniently annotated for the SV indel genotyping step.
 *
 * <p>
 *     Currently only SV indels are supported: {@code SVTYPE=INV} or {@code SVTYPE=DEL}.
 * </p>
 * <p>
 *     The main output alignment (-O|--output output.bam) is composed of unmapped records representing haplotypes and contigs with coordinates
 *     that match the variant they are relevant to. This is the most convenient format for
 *     the genotyping tool downstream. Optionally one can request to generate an additional re-aligned version where all records are
 *     aligned against the reference around the relevant variant (-alnO|--alignedOutput alnout.bam).
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
 *     overlap any of the break-points of that variant in the input assemblies file (argument -C|--contigs xxx.bam|sam).
 *     These are all part of the {@code "CTG"} read-group.
 * </p>
 * <p>All records are marked as unmapped, unpaired and pass vendor filters</p>
 * <p>
 *    Also these have a number of annotations indicating amongst other things what allele/haplotype they support:
 *     <dl>
 *         <dt>HP</dt><dd>the haplotype the support, {@code 'ref'} for reference, {@code 'alt'} for alternative and {@code '.'} for neither. (e.g. {@code "HP:Z:ref", "HP:Z:alt", "HP:Z:."}).</dd>
 *         <dt>HQ</dt><dd>the confidence in the support for the haplotype in {@code "HP"} (e.g. {@code "HQ:Z:10.0"}).
 *                        This value is equal to the difference between the score for the supported and the other haplotype</dd></dt>
 *         <dt>RS</dt><dd>the reference support score screen (e.g. {@code "RS:Z:-100,55,1,2,30,0"}).</dd>
 *         <dt>XS</dt><dd>the alternative allele score screen (e.g. {@code "XS:Z:-1241,134,2,1,5,1}).</dd>
 *         <dt>RA</dt><dd>alignment versus the "ref" haplotype for contigs and vs the reference for haplotypes</dd>
 *         <dt>XA</dt><dd>alignment versus the "alt" haplotype for contigs. It does not apply to haplotypes</dd>
 *         <dt>VC</dt><dd>coordinate of the targeted variant {@code chr:pos}</dd>
 *     </dl>
 * </p>
 * <p>
 *     {@code RS} and {@code XS} score annotations follows this format:
 * </p>
 *     <pre>Score,matches,mismatches,indels,indelLength,reversals</pre>
 * <p>
 *     Where:
 *     <dl>
 *         <dt>score</dt><dd>score the contig given the happlotype (ref. for {@code RS}, alt. for {@code XS})</dd>
 *         <dt>matches</dt><dd>number of base call matches</dd>
 *         <dt>mismatches</dt><dd>number of base call mismatches</dd>
 *         <dt>indels</dt><dd>number of insertion or deletions (i.e. gap-openings)</dd>
 *         <dt>indelLength</dt><dd>total length of insertion and deletions (i.e. gap-openings + gap-extensions)</dd>
 *         <dt>reversals</dt><dd>number of those strand switches between contiguous aligned intervals in the contig</dd>
 *     </dl>
 * </p>
 */
@CommandLineProgramProperties(
        summary = "composes contigs file for genotyping",
        oneLineSummary = "composes contigs file for genotyping",
        usageExample = "--variant ins_and_dels.vcf --contigs assemblies.bam --reference my-ref.fa " +
                       "--shardSize 10000 --paddingSize 300 --output genotyping-contigs.bam",
        programGroup = StructuralVariantDiscoveryProgramGroup.class )
@BetaFeature
public class ComposeStructuralVariantHaplotypesSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    public static final String CONTIGS_FILE_SHORT_NAME = "C";
    public static final String CONTIGS_FILE_FULL_NAME = "contigs";
    public static final String SHARD_SIZE_SHORT_NAME = "sz";
    public static final String SHARD_SIZE_FULL_NAME = "shardSize";
    public static final String PADDING_SIZE_SHORT_NAME = "pd";
    public static final String PADDING_SIZE_FULL_NAME = "paddingSize";
    public static final String ALIGNED_OUTPUT_SHORT_NAME = "alnOut";
    public static final String ALIGNED_OUTPUT_FULL_NAME = "alignedOutput";
    public static final String SAMPLE_SHORT_NAME = "S";
    public static final String SAMPLE_FULL_NAME = "sampleName";

    public static final int DEFAULT_SHARD_SIZE = 10_000;
    public static final int DEFAULT_PADDING_SIZE = 300; // ideally larger than a read length.
    public static final String DEFAULT_SAMPLE = "sample";

    public static final String HAPLOTYPE_READ_GROUP = "HAP";
    public static final String CONTIG_READ_GROUP = "CTG";
    public static final String HAPLOTYPE_CALL_TAG = "HP";
    public static final String HAPLOTYPE_QUAL_TAG = "HQ";
    public static final String REFERENCE_SCORE_TAG = "RS";
    public static final String REFERENCE_ALIGNMENT_TAG = "RA";
    public static final String ALTERNATIVE_SCORE_TAG = "XS";
    public static final String ALTERNATIVE_ALIGNMENT_TAG = "XA";
    public static final String VARIANT_CONTEXT_TAG = "VC";

    public static final String REF_CONTIG_NAME = "ref";
    public static final String ALT_CONTIG_NAME = "alt";

    @Argument(doc = "shard size",
              shortName = SHARD_SIZE_SHORT_NAME,
              fullName = SHARD_SIZE_FULL_NAME,
    optional = true)
    private int shardSize = DEFAULT_SHARD_SIZE;

    @Argument(doc = "sample name to use in output BAM files read-groups",
              shortName = SAMPLE_SHORT_NAME,
              fullName = SAMPLE_FULL_NAME,
              optional = true)
    private String sample = DEFAULT_SAMPLE;

    @Argument(doc ="padding size",
              shortName = PADDING_SIZE_SHORT_NAME,
              fullName = PADDING_SIZE_FULL_NAME,
              optional = true)
    private int paddingSize = DEFAULT_PADDING_SIZE;

    @Argument(doc = "input variant file",
              shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
              fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME)
    private String variantsFileName = null;

    @Argument(doc = "aligned contig file",
              fullName = CONTIGS_FILE_FULL_NAME,
              shortName = CONTIGS_FILE_SHORT_NAME
    )
    private String alignedContigsFileName = null;

    @Argument(doc = "output bam file with contigs per variant",
              fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputFileName = null;

    @Argument(doc = "aligned output bam file with haplotypes and contigs aligned against the reference",
            fullName = ALIGNED_OUTPUT_FULL_NAME,
            shortName = ALIGNED_OUTPUT_SHORT_NAME,
            optional = true
    )
    private String alignedOutputFileName = null;

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        Utils.nonNull(ctx);
        Utils.nonNull(alignedContigsFileName);

        final ReadsSparkSource alignedContigs = new ReadsSparkSource(ctx);
        final VariantsSparkSource variantsSource = new VariantsSparkSource(ctx);

        final List<SimpleInterval> wholeGenomeIntervals = IntervalUtils.getAllIntervalsForReference(getReferenceSequenceDictionary());

        final List<SimpleInterval> targetIntervals = hasIntervals()
                ? getIntervals()
                : wholeGenomeIntervals;


        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        final Broadcast<SAMSequenceDictionary> sequenceDictionaryBroadcast = ctx.broadcast(sequenceDictionary);
        final Broadcast<IntervalsSkipList<SimpleInterval>> shardIntervalsBroadcast = composeShardBroadcastWholeReference(ctx, sequenceDictionary);

        final int paddingSize = this.paddingSize;
        final JavaRDD<SVContext> variants = variantsSource
                .getParallelVariantContexts(variantsFileName, getIntervals())
                .filter(vc -> supportedVariant(vc,paddingSize))
                .map(SVContext::of);
        final JavaRDD<GATKRead> contigs = alignedContigs
                .getParallelReads(alignedContigsFileName, referenceArguments.getReferenceFileName(),
                        new TraversalParameters(targetIntervals, false));
        final JavaPairRDD<SVContext, Tuple2<List<String>, List<SimpleInterval>>> variantsAndRelevantContigIntervals =
                composeVariantAndRelevantContigsIntervals(ctx, shardIntervalsBroadcast, sequenceDictionaryBroadcast, contigs, variants);

        final Map<String, SVContext> variantsByUID =
                variants.mapToPair(v -> new Tuple2<>(v.getUniqueID(), v)).collectAsMap();
        final Map<String, List<String>> variantUIDByContigName = variantsAndRelevantContigIntervals.flatMapToPair(
                t -> t._2()._1().stream().map(readName -> new Tuple2<>(readName, (List<String>) new ArrayList<>(Collections.singleton(t._1().getUniqueID())))).iterator())
                .reduceByKey((a, b) -> {a.addAll(b); return a; }).collectAsMap();

        final SVIntervalLocator locator = SVIntervalLocator.of(sequenceDictionary);
        final SVIntervalTree<ShardBoundary> shardTree = wholeGenomeIntervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, shardSize, 0, sequenceDictionary).stream())
                .collect(locator.toTreeCollector(Function.identity()));

        final Broadcast<SVIntervalTree<ShardBoundary>> shardTreeBCast = ctx.broadcast(shardTree);
        final Broadcast<SVIntervalLocator> locatorBCast = ctx.broadcast(locator);

        final List<ShardBoundary> relevantShards = variantsAndRelevantContigIntervals
                .mapPartitions(it -> {
                    final SVIntervalTree<ShardBoundary> mapShardTree = shardTreeBCast.getValue();
                    final SVIntervalLocator mapLocator = locatorBCast.getValue();
                    final Set<ShardBoundary> bundaries = new LinkedHashSet<>();
                    while (it.hasNext()) {
                        final List<SimpleInterval> list = it.next()._2()._2();
                        for (final SimpleInterval interval : list) {
                            final Iterator<SVIntervalTree.Entry<ShardBoundary>> entryIt =
                                        mapShardTree.overlappers(mapLocator.toSVInterval(interval));
                            while (entryIt.hasNext()) {
                                bundaries.add(entryIt.next().getValue());
                            }
                        }
                    }
                    return bundaries.iterator();
                }).distinct().collect();
        final List<SimpleInterval> relevantShardIntervals = relevantShards.stream()
                .map(ShardBoundary::getInterval).collect(Collectors.toList());

        final Broadcast<Map<String, List<String>>> variantUIDByContigNameBCast = ctx.broadcast(variantUIDByContigName);
        final Broadcast<Map<String, SVContext>> variantsByUIDBCast = ctx.broadcast(variantsByUID);

        final JavaRDD<GATKRead> primaryContings = alignedContigs
                .getParallelReads(alignedContigsFileName, referenceArguments.getReferenceFileName(),
                        new TraversalParameters(relevantShardIntervals, false));
        final JavaPairRDD<SVContext, List<AlignedContig>> result =
                primaryContings.mapPartitionsToPair(it -> {
                    final Map<String, List<String>> mapVariantUIDByContigName = variantUIDByContigNameBCast.getValue();
                    final Map<String, SVContext> variantByUID = variantsByUIDBCast.getValue();
                    final Map<SVContext, List<GATKRead>> accumulator = new LinkedHashMap<>();
                    while (it.hasNext()) {
                        final GATKRead read = it.next();
                        final List<String> uids = mapVariantUIDByContigName.getOrDefault(read.getName(), Collections.emptyList());
                        for (final String uid : uids) {
                            final SVContext svc = variantByUID.get(uid);
                            if (svc != null) {
                                final List<GATKRead> alignedContig = accumulator.computeIfAbsent(svc, (key) -> new ArrayList<>());
                                alignedContig.add(read);
                            }
                        }
                    }
                    return accumulator.entrySet().stream()
                            .map(e -> new Tuple2<>(e.getKey().getUniqueID(), new Tuple2<>(e.getKey(), e.getValue())))
                            .iterator();
                })
                // JavaPairRDD<String,Tuple2<SVContext, List<GATKRead>>> -> <SVContent, List<GATKRead>> :
                .reduceByKey((t1, t2) -> {
                    t1._2().addAll(t2._2());
                    return t1;
                }).mapToPair(Tuple2::_2)
                // JavaPairRDD<SVContext, List<GATKRead>> -> <SVContext, List<AlignedContig>> :
                .mapValues(readList -> readList.stream()
                           .collect(Collectors.groupingBy(GATKRead::getName)).values().stream() //: Stream<List<GATKRead>>.
                                .map(ll -> AlignedContig.of(removeRepeatedReads(ll))) // AlignmentContig.of(Iterable<GATKRead>) : AlignedContig.
                                .collect(Collectors.toList()));

        processVariants(result, sequenceDictionary, sequenceDictionaryBroadcast, shardIntervalsBroadcast);

    }

    private static List<GATKRead> removeRepeatedReads(final List<GATKRead> original) {
        if (original.size() < 1) {
            return original;
        } else if (original.size() == 2 && sameRead(original.get(0), original.get(1))) {
            return Collections.singletonList(original.get(0));
        } else {

            final List<GATKRead> result = new ArrayList<>(original.size());
            original.stream().sorted(GATK_WITHIN_READ_RECORD_COMPARATOR).forEach(r -> {
                if (result.isEmpty() || !sameRead(result.get(result.size() - 1), r)) {
                    result.add(r);
                }
            });
            return result;
        }
    }


    private static boolean sameRead(final GATKRead r1, final GATKRead r2) {
        if (r1 == r2) {
            return true;
        } else if (!r1.getName().equals(r2.getName())) {
            return false;
        } else if (!r1.getContig().equals(r2.getContig())) {
            return false;
        } else if (r1.isUnmapped() != r2.isUnmapped()) {
            return false;
        } else if (r1.isPaired() != r2.isPaired()) {
            return false;
        } else if (Objects.equals(r1.getReadGroup(), r2.getReadGroup())) {
            return false;
        }
        if (!r1.isUnmapped()) {
            if (r1.getStart() != r2.getStart()) {
                return false;
            } else if (r1.getUnclippedStart() != r2.getUnclippedStart()) {
                return false;
            } else if (r1.getEnd() != r2.getEnd()) {
                return false;
            } else if (r1.getUnclippedEnd() != r2.getUnclippedEnd()) {
                return false;
            } else if (r1.getMappingQuality() != r2.getMappingQuality()) {
                return false;
            } else if (r1.isReverseStrand() != r2.isReverseStrand()) {
                return false;
            } else if (r1.isSecondaryAlignment() != r2.isSecondaryAlignment()) {
                return false;
            } else if (r1.isSupplementaryAlignment() != r2.isSupplementaryAlignment()) {
                return false;
            } else if (!CigarUtils.equals(r1.getCigar(), r2.getCigar())) {
                return false;
            }
        }
        if (r1.isPaired()) {
            if (r1.mateIsUnmapped() != r2.mateIsUnmapped()) {
                return false;
            } else if (!r1.mateIsUnmapped()) {
              if (r1.mateIsReverseStrand() != r2.mateIsReverseStrand()) {
                  return false;
              } else if (!r1.getMateContig().equals(r2.getMateContig())) {
                  return false;
              } else if (r1.getMateStart() != r2.getMateStart()) {
                  return false;
              }
            }
        }
        if (r1.getBases() != r2.getBases()) {
            return false;
        } else if (r1.getBases() != null && !Arrays.equals(r1.getBases(), r2.getBases())) {
            return false;
        }
        return true;
    }

    private static final Comparator<GATKRead> GATK_WITHIN_READ_RECORD_COMPARATOR =
            Comparator.comparing(GATKRead::getContig).thenComparingInt(GATKRead::getStart);

    private static boolean supportedVariant(final VariantContext vc, final int paddingSize) {
        final StructuralVariantType type = vc.getStructuralVariantType();

        if (type == null) {
            return false;
        } else if (vc.getAlternateAlleles().size() != 1) {
            return false;
        } else if (vc.hasAttribute(GATKSVVCFConstants.IMPRECISE)) {
            return false;
        } else if (type == StructuralVariantType.INS || type == StructuralVariantType.DEL) {
            return vc.hasAttribute(GATKSVVCFConstants.SVLEN); // for now we skip indels without SVLEN... there are some, perhaps a bug
        } else if (type == StructuralVariantType.DUP) {
            return true;
        } else if (type == StructuralVariantType.INV) {
            return true;
        } else {
            return false;
        }
    }

    private static final Pattern SA_SPLIT_PATTERN = Pattern.compile(";");

    private JavaPairRDD<SVContext, Tuple2<List<String>, List<SimpleInterval>>> composeVariantAndRelevantContigsIntervals(
            final JavaSparkContext ctx, final Broadcast<IntervalsSkipList<SimpleInterval>> shardIntervalsBroadcast,
                                        final Broadcast<SAMSequenceDictionary> sequenceDictionaryBroadcast,
            final JavaRDD<GATKRead> contigs, final JavaRDD<SVContext> variants) {

        final JavaPairRDD<SimpleInterval, List<Tuple2<SimpleInterval,GATKRead>>> contigsInShards =
                groupInShards(contigs, ComposeStructuralVariantHaplotypesSpark::readIntervalList, shardIntervalsBroadcast);
        final int paddingSize = this.paddingSize;

        final JavaPairRDD<SimpleInterval, List<Tuple2<SimpleInterval, SVContext>>> variantsInShards =
                groupInShards(variants, (v) -> v.getBreakPointIntervals(paddingSize, sequenceDictionaryBroadcast.getValue(), true), shardIntervalsBroadcast);

        final JavaPairRDD<SimpleInterval, Tuple2<List<Tuple2<SimpleInterval, GATKRead>>,
                List<Tuple2<SimpleInterval, SVContext>>>> contigAndVariantsInShards
                = contigsInShards.join(variantsInShards);

        final JavaPairRDD<SVContext, Tuple2<List<String>, List<SimpleInterval>>> contigsPerVariantInterval
                = contigAndVariantsInShards.flatMapToPair(t -> {
            final List<Tuple2<SimpleInterval, SVContext>> vars = t._2()._2();
            final List<Tuple2<SimpleInterval, GATKRead>> ctgs = t._2()._1();

            return vars.stream()
                .map(v -> {
                        final List<GATKRead> cs = ctgs.stream()
                                .filter(ctg -> v._1().overlaps(ctg._1()))
                                .map(Tuple2::_2)
                                .collect(Collectors.toList());
                        final List<SimpleInterval> its = cs.stream()
                                .flatMap(GATKREAD_TO_ALL_ALIGNMENT_STARTS)
                                .collect(Collectors.toList());
                        final List<String> names = cs.stream()
                                .map(GATKRead::getName)
                                .distinct()
                                .collect(Collectors.toList());
                        return new Tuple2<>(v._2(), new Tuple2<>( names, its)); })
                    .collect(Collectors.toList()).iterator();
                });

        return contigsPerVariantInterval;
    }

    private Broadcast<IntervalsSkipList<SimpleInterval>> composeShardBroadcastWholeReference(final JavaSparkContext ctx, SAMSequenceDictionary sequenceDictionary) {
        final List<SimpleInterval> intervals = IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        // use unpadded shards (padding is only needed for reference bases)
        final List<ShardBoundary> shardBoundaries = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, shardSize, 0, sequenceDictionary).stream())
                .collect(Collectors.toList());
        final IntervalsSkipList<SimpleInterval> shardIntervals = new IntervalsSkipList<>(shardBoundaries.stream()
                .map(ShardBoundary::getPaddedInterval)
                .collect(Collectors.toList()));

        return ctx.broadcast(shardIntervals);
    }

    private static final SerializableFunction<GATKRead, Stream<SimpleInterval>> GATKREAD_TO_ALL_ALIGNMENT_STARTS =
            ComposeStructuralVariantHaplotypesSpark::primaryAlignmentInterval;

    private static Stream<SimpleInterval> primaryAlignmentInterval(final GATKRead record) {
        final Stream<SimpleInterval> recordCoords = Stream.of(new SimpleInterval(record.getAssignedContig(), record.getAssignedStart(), record.getAssignedStart()));
        if (!record.hasAttribute(SAMTag.SA.name())) {
            return recordCoords;
        } else {
            final Stream<SimpleInterval> saCoords = SA_SPLIT_PATTERN.splitAsStream(
                    record.getAttributeAsString(SAMTag.SA.name()))
                    .filter(s -> !s.isEmpty())
                    .map(AlignmentInterval::new)
                    .map(ai -> ai.referenceSpan.getStartInterval());
            return Stream.concat(recordCoords, saCoords);
        }
    }

    private JavaPairRDD<SVContext, List<GATKRead>> composeOverlappingContigRecordsPerVariant(
            final JavaSparkContext ctx, final JavaRDD<GATKRead> contigs, final JavaRDD<SVContext> variants) {
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        final List<SimpleInterval> intervals = hasIntervals()
                ? getIntervals()
                : IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        // use unpadded shards (padding is only needed for reference bases)
        final List<ShardBoundary> shardBoundaries = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, shardSize, 0, sequenceDictionary).stream())
                .collect(Collectors.toList());
        final IntervalsSkipList<SimpleInterval> shardIntervals = new IntervalsSkipList<>(shardBoundaries.stream()
                .map(ShardBoundary::getPaddedInterval)
                .collect(Collectors.toList()));
        final Broadcast<SAMSequenceDictionary> dictionaryBroadcast = ctx.broadcast(sequenceDictionary);

        final Broadcast<IntervalsSkipList<SimpleInterval>> shardIntervalsBroadcast = ctx.broadcast(shardIntervals);

        final JavaPairRDD<SimpleInterval, List<Tuple2<SimpleInterval,GATKRead>>> contigsInShards =
            groupInShards(contigs, ComposeStructuralVariantHaplotypesSpark::readIntervalList, shardIntervalsBroadcast);
        final int paddingSize = this.paddingSize;

        final JavaPairRDD<SimpleInterval, List<Tuple2<SimpleInterval, SVContext>>> variantsInShards =
            groupInShards(variants, (v) -> v.getBreakPointIntervals(paddingSize, dictionaryBroadcast.getValue(), true), shardIntervalsBroadcast);

        final JavaPairRDD<SimpleInterval, Tuple2<List<Tuple2<SimpleInterval, GATKRead>>,
                    List<Tuple2<SimpleInterval, SVContext>>>> contigAndVariantsInShards
                = contigsInShards.join(variantsInShards);


        final JavaPairRDD<SVContext, List<GATKRead>> contigsPerVariantInterval
                = contigAndVariantsInShards.flatMapToPair(t -> {
                    final List<Tuple2<SimpleInterval, SVContext>> vars = t._2()._2();
                    final List<Tuple2<SimpleInterval, GATKRead>> ctgs = t._2()._1();
                    return vars.stream()
                            .map(v -> {
                                final List<GATKRead> cs = ctgs.stream()
                                        .filter(ctg -> v._1().overlaps(ctg._1()))
                                        .map(Tuple2::_2)
                                        .collect(Collectors.toList());

                                return new Tuple2<>(v._2(), cs);})
                            .collect(Collectors.toList()).iterator();
                });

        final Function<SVContext, String> variantId =
                (SerializableFunction<SVContext, String>) ComposeStructuralVariantHaplotypesSpark::variantId;
        final Function2<List<GATKRead>, List<GATKRead>, List<GATKRead>> readListMerger
                = (a, b) -> Stream.concat(a.stream(), b.stream()).collect(Collectors.toList());

        // Merge contig lists on the same variant-context coming from different intervals
        // into one.
        final JavaPairRDD<SVContext, List<GATKRead>> contigsPerVariant = contigsPerVariantInterval
                .mapToPair(t -> new Tuple2<>(variantId.apply(t._1()), t))
                .reduceByKey((a, b) -> new Tuple2<>(a._1(), readListMerger.call(a._2(), b._2())))
                .mapToPair(Tuple2::_2);
        return contigsPerVariant;
    }

    private static String variantId(final SVContext variant) {
        if (variant.getID() != null && !variant.getID().isEmpty()) {
            return variant.getID();
        } else {
            final int length = variant.getStructuralVariantLength();
            return "var_" + variant.getAlternateAllele(0).getDisplayString() + "_" + length;
        }
    }

    private static List<SimpleInterval> readIntervalList(final GATKRead contig) {
        if (contig.isUnmapped()) {
            return Collections.emptyList();
        } else {
            return Collections.singletonList(new SimpleInterval(contig.getContig(), contig.getStart(), contig.getEnd()));
        }
    }

    private <T> JavaPairRDD<SimpleInterval, List<Tuple2<SimpleInterval, T>>> groupInShards(
                final JavaRDD<T> elements,
                final org.apache.spark.api.java.function.Function<T, List<SimpleInterval>> intervalsOf,
                final Broadcast<IntervalsSkipList<SimpleInterval>> shards) {
        final PairFlatMapFunction<T, SimpleInterval, Tuple2<SimpleInterval, T>> flatMapIntervals =
                t -> intervalsOf.call(t).stream().map(i -> new Tuple2<>(i, new Tuple2<>(i,t))).iterator();

        return elements
                .flatMapToPair(flatMapIntervals)
                .flatMapToPair(t -> shards.getValue().getOverlapping(t._1()).stream().map(i -> new Tuple2<>(i, t._2())).iterator())
                .aggregateByKey(new ArrayList<Tuple2<SimpleInterval, T>>(10),
                        (l1, c) -> { l1.add(c); return l1;},
                        (l1, l2) -> {l1.addAll(l2); return l1;});
    }

    protected void processVariants(final JavaPairRDD<SVContext, List<AlignedContig>> variantsAndOverlappingUniqueContigs,
                                   final SAMSequenceDictionary dictionary,
                                   final Broadcast<SAMSequenceDictionary> broadcastDictionary,
                                   final Broadcast<IntervalsSkipList<SimpleInterval>> shardBroadcast) {

        final SAMFileHeader outputHeader = composeOutputHeader(dictionary);
        try (final SAMFileWriter outputWriter = BamBucketIoUtils.makeWriter(outputFileName, outputHeader, true);
             final SAMFileWriter alignedOutputWriter = alignedOutputFileName == null
                ? null : BamBucketIoUtils.makeWriter(alignedOutputFileName, outputHeader, true)) {

            final SerializableComparator<SimpleInterval> shardSorter = createLocatableSparkSorter(broadcastDictionary);

            final JavaPairRDD<SVContext, VariantHaplotypesAndContigsComposite> alignedContigs = variantsAndOverlappingUniqueContigs.mapToPair(tuple -> {
                        final SVContext vc = tuple._1();
                        final List<AlignedContig> contigs = tuple._2();
                        final List<SVHaplotype> haplotypes = vc.composeHaplotypesBasedOnReference(paddingSize, getReference(), contigs);
                        final SVHaplotype referenceHaplotype = haplotypes.get(0);
                        final SVHaplotype alternativeHaplotype = haplotypes.get(1);
                        final List<AlignedContig> referenceAlignedContigs = referenceHaplotype.alignContigs(contigs);
                        final List<AlignedContig> alternativeAlignedContigs = alternativeHaplotype.alignContigs(contigs);
                        return new Tuple2<>(vc, new VariantHaplotypesAndContigsComposite(referenceHaplotype, alternativeHaplotype, contigs, referenceAlignedContigs, alternativeAlignedContigs));
            });


            final JavaPairRDD<SVContext, VariantHaplotypesAndContigsComposite> cachedIfNeededAlignedContigs
                    = alignedOutputWriter == null ? alignedContigs : alignedContigs.persist(StorageLevel.MEMORY_AND_DISK());

            cachedIfNeededAlignedContigs
                    .mapToPair(tuple -> new Tuple2<>(shardBroadcast.getValue().getOverlapping(tuple._1().getStartPositionInterval()).get(0), tuple))
                    .combineByKey(v -> Stream.of(v).collect(Collectors.toList()),
                            (l, v) -> {l.add(v); return l;} , (l1, l2) -> { l1.addAll(l2); return l1; })
                    .sortByKey(shardSorter)
                    .mapValues(l -> { l.sort(Comparator.comparingInt(v -> v._1().getStart())); return l; })
                    .values()
                    .flatMap(List<Tuple2<SVContext, VariantHaplotypesAndContigsComposite>>::iterator)
                    .toLocalIterator().forEachRemaining(tuple -> {
                        final SVContext variant = tuple._1();
                        final VariantHaplotypesAndContigsComposite haplotypesAndContigs = tuple._2();
                        final int numberOfContigs = haplotypesAndContigs.numberOfContigs();
                        outputWriter.addAlignment(haplotypesAndContigs.composeOutputReferenceHaplotypeRecord(outputHeader, variant));
                        outputWriter.addAlignment(haplotypesAndContigs.composeOutputAlternativeHaplotypeRecord(outputHeader, variant));
                        for (int i = 0; i < numberOfContigs; i++) {
                            outputWriter.addAlignment(haplotypesAndContigs.composeOutputSAMRecord(outputHeader, variant, i));
                        }
                    });

            composeAlignedOutput(shardBroadcast, outputHeader, alignedOutputWriter, shardSorter, cachedIfNeededAlignedContigs);
        }
    }

    private void composeAlignedOutput(Broadcast<IntervalsSkipList<SimpleInterval>> shardBroadcast, SAMFileHeader outputHeader, SAMFileWriter alignedOutputWriter, SerializableComparator<SimpleInterval> shardSorter, JavaPairRDD<SVContext, VariantHaplotypesAndContigsComposite> cachedIfNeededAlignedContigs) {
        if (alignedOutputWriter != null) {
            cachedIfNeededAlignedContigs.flatMapToPair(tuple -> {
                final SVContext variant = tuple._1();
                final VariantHaplotypesAndContigsComposite haplotypesAndContigs = tuple._2();
                final List<SAMRecord> result = new ArrayList<>();

                result.addAll(haplotypesAndContigs.composeAlignedOutputReferenceHaplotypeSAMRecords(null, variant));
                result.addAll(haplotypesAndContigs.composeAlignedOutputAlternativeHaplotypeSAMRecords(null, variant));
                for (int i = 0; i < haplotypesAndContigs.numberOfContigs(); i++) {
                    result.addAll(haplotypesAndContigs.composeAlignedOutputSAMRecords(null, variant, i));
                }
                final IntervalsSkipList<SimpleInterval> shards = shardBroadcast.getValue();
                return result.stream()
                        .map(r -> new Tuple2<>(shards.getOverlapping(new SimpleInterval(r.getContig(), r.getStart(), r.getStart())).get(0), r)).iterator();
            })
            .combineByKey((v) -> Stream.of(v).collect(Collectors.toList()) ,
                            (l, v) -> { l.add(v); return l; }, (l1, l2) -> { l1.addAll(l2); return l1; } )
            .sortByKey(shardSorter)
            .mapValues(unsorted -> { unsorted.sort(Comparator.comparingInt(SAMRecord::getAlignmentStart)); return unsorted; })
            .values()
            .flatMap(List::iterator)
            .toLocalIterator().forEachRemaining(record -> {
                record.setHeader(outputHeader);
                alignedOutputWriter.addAlignment(record);
            });
        }
        ;
    }

    private <T extends Locatable> SerializableComparator<T> createLocatableSparkSorter(final Broadcast<SAMSequenceDictionary> broadcastDictionary) {
        return (a, b) -> {
            final String aContig = a.getContig();
            final String bContig = b.getContig();
            if (aContig.equals(bContig)) {
                return Integer.compare(a.getStart(), b.getStart());
            } else {
                final SAMSequenceDictionary myDictionary = broadcastDictionary.getValue();
                return Integer.compare(myDictionary.getSequenceIndex(aContig), myDictionary.getSequenceIndex(bContig));
            }
        };
    }

    private SAMFileHeader composeOutputHeader(SAMSequenceDictionary dictionary) {
        final SAMFileHeader outputHeader = new SAMFileHeader();
        final SAMProgramRecord programRecord = new SAMProgramRecord(getProgramName());
        programRecord.setCommandLine(getCommandLine());
        outputHeader.setSequenceDictionary(dictionary);
        outputHeader.addProgramRecord(programRecord);
        outputHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord contigsReadGroup = new SAMReadGroupRecord(CONTIG_READ_GROUP);
        contigsReadGroup.setSample(sample);
        final SAMReadGroupRecord haplotypesReadGroup = new SAMReadGroupRecord(HAPLOTYPE_READ_GROUP);
        haplotypesReadGroup.setSample(sample);
        outputHeader.addReadGroup(haplotypesReadGroup);
        outputHeader.addReadGroup(contigsReadGroup);
        return outputHeader;
    }

    /**
     * Class to represent and calculate the aligned contig score.
     */
    private static class AlignedContigScore {

        public static final int STRAND_SWITCH_COST = 60;
        public static final double MATCH_COST = 0.01;
        public static final int MISMATCH_COST = 30;
        public static final double GAP_OPEN_COST = 45;
        public static final double GAP_EXTEND_COST = 3;

        final int totalReversals;
        final int totalIndels;
        final int totalMatches;
        final int totalMismatches;
        final int totalIndelLength;

        public AlignedContigScore(final int reversals, final int indels, final int matches, final int mismatches, final int totalIndelLength) {
            this.totalReversals = reversals;
            this.totalIndels = indels;
            this.totalMatches = matches;
            this.totalMismatches = mismatches;
            this.totalIndelLength = totalIndelLength;
        }

        public double getValue() {
            return -(int) Math.round(totalMatches * MATCH_COST
                    + totalMismatches * MISMATCH_COST
                    + totalIndels * GAP_OPEN_COST
                    + (totalIndelLength - totalIndels) * GAP_EXTEND_COST
                    + totalReversals * STRAND_SWITCH_COST);
        }

        public String toString() {
            return  getValue() + ":" + Utils.join(",", totalMatches, totalMismatches,
                    totalIndels, totalIndelLength, totalReversals);
        }
    }

    /**
     * Represent the information derived of collecting contigs and haplotypes that are relevant to a variant context
     * and align them versus the reference and alternative haplotype.
     */
    @DefaultSerializer(VariantHaplotypesAndContigsComposite.Serializer.class)
    static class VariantHaplotypesAndContigsComposite implements Serializable {

        private static final long serialVersionUID = 1L;
        private static final SerializableComparator<AlignmentInterval> ALIGNMENT_INTERVAL_PRIMARY_SECONDARY_SORTER =
                (a, b) -> {
                    if (a.mapQual != b.mapQual) {
                        return Integer.compare(b.mapQual == SAMRecord.UNKNOWN_MAPPING_QUALITY ? -1 : b.mapQual,
                                a.mapQual == SAMRecord.UNKNOWN_MAPPING_QUALITY ? - 1: a.mapQual);
                    } else if (a.alnScore != b.alnScore) {
                        return Integer.compare(b.alnScore == AlignmentInterval.NO_AS? -1 : b.alnScore,
                                a.alnScore == AlignmentInterval.NO_AS? -1 : a.alnScore);
                    } else {
                        final int bAlignedBases = CigarUtils.countAlignedBases(b.cigarAlong5to3DirectionOfContig);
                        final int aAlignedBases = CigarUtils.countAlignedBases(a.cigarAlong5to3DirectionOfContig);
                        if (bAlignedBases != aAlignedBases) {
                            return Integer.compare(bAlignedBases, aAlignedBases);
                        } else if (a.mismatches != AlignmentInterval.NO_NM && b.mismatches != AlignmentInterval.NO_NM && a.mismatches != b.mismatches) {
                            return Integer.compare(b.mismatches, a.mismatches);
                        } else if (!a.referenceSpan.getContig().equals(b.referenceSpan.getContig())) {
                            return b.referenceSpan.getContig().compareTo(a.referenceSpan.getContig());
                        } else if (a.referenceSpan.getStart() != b.referenceSpan.getStart()) {
                            return Integer.compare(b.referenceSpan.getStart(), a.referenceSpan.getStart());
                        } else if (a.referenceSpan.getEnd() != b.referenceSpan.getEnd()) {
                            return Integer.compare(a.referenceSpan.getStart(), b.referenceSpan.getStart());
                        } else if (a.startInAssembledContig != b.startInAssembledContig) {
                            return Integer.compare(b.startInAssembledContig, a.startInAssembledContig);
                        } else if (a.endInAssembledContig != b.endInAssembledContig) {
                            return Integer.compare(a.endInAssembledContig, b.endInAssembledContig);
                        } else { // last resort, compare the cigar representation as text.... this probrably never happens though.
                            return a.cigarAlong5to3DirectionOfContig.toString().compareTo(b.cigarAlong5to3DirectionOfContig.toString());
                        }
                    }
                };
        private final SVHaplotype referenceHaplotype;
        private final SVHaplotype alternativeHaplotype;
        private final AlignedContig[] originalAlignments;
        private final AlignedContig[] referenceAlignments;
        private final AlignedContig[] alternativeAlignments;
        private final String[] referenceScoreTagValue;
        private final String[] alternativeScoreTagValue;
        private final String[] hpTagValue;
        private final double[] hpQualTagValue;

        VariantHaplotypesAndContigsComposite(final SVHaplotype referenceHaplotype, final SVHaplotype alternativeHaplotype,
                                             final List<AlignedContig> originalAlignment, final List<AlignedContig> referenceAlignedContigs,
                                             final List<AlignedContig> alternativeAlignedContigs) {
            final int numberOfContigs = referenceAlignedContigs.size();
            this.referenceHaplotype = referenceHaplotype;
            this.alternativeHaplotype = alternativeHaplotype;
            this.originalAlignments = new AlignedContig[numberOfContigs];
            this.referenceAlignments = new AlignedContig[numberOfContigs];
            this.alternativeAlignments = new AlignedContig[numberOfContigs];
            this.referenceScoreTagValue = new String[numberOfContigs];
            this.alternativeScoreTagValue = new String[numberOfContigs];
            this.hpTagValue = new String[numberOfContigs];
            this.hpQualTagValue = new double[numberOfContigs];
            for (int i = 0; i < numberOfContigs; i++) {
                this.originalAlignments[i] = originalAlignment.get(i);
                referenceAlignments[i] = referenceAlignedContigs.get(i);
                alternativeAlignments[i] = alternativeAlignedContigs.get(i);
                final AlignedContigScore referenceScore = calculateAlignedContigScore(referenceAlignments[i]
                        = referenceAlignedContigs.get(i));
                final AlignedContigScore alternativeScore = calculateAlignedContigScore(alternativeAlignments[i]
                        = alternativeAlignedContigs.get(i));
                referenceScoreTagValue[i] = referenceScore.toString();
                alternativeScoreTagValue[i] = alternativeScore.toString();
                hpTagValue[i] = calculateHPTag(referenceScore.getValue(), alternativeScore.getValue());
                hpQualTagValue[i] = calculateHPQualTag(referenceScore.getValue(), alternativeScore.getValue());
            }
        }

        public VariantHaplotypesAndContigsComposite(final Kryo kryo, final Input input) {
            final int numberOfContigs = input.readInt();
            this.originalAlignments = new AlignedContig[numberOfContigs];
            this.referenceAlignments = new AlignedContig[numberOfContigs];
            this.alternativeAlignments = new AlignedContig[numberOfContigs];
            this.referenceScoreTagValue = new String[numberOfContigs];
            this.alternativeScoreTagValue = new String[numberOfContigs];
            this.hpTagValue = new String[numberOfContigs];
            this.hpQualTagValue = new double[numberOfContigs];
            this.referenceHaplotype = (SVHaplotype) kryo.readClassAndObject(input);
            this.alternativeHaplotype = (SVHaplotype) kryo.readClassAndObject(input);
            for (int i = 0; i < numberOfContigs; i++) {
                this.originalAlignments[i] = kryo.readObject(input, AlignedContig.class);
                this.referenceAlignments[i] = kryo.readObject(input, AlignedContig.class);
                this.alternativeAlignments[i] = kryo.readObject(input, AlignedContig.class);
                this.referenceScoreTagValue[i] = input.readString();
                this.alternativeScoreTagValue[i] = input.readString();
                this.hpTagValue[i] = input.readString();
                this.hpQualTagValue[i] = input.readDouble();
            }
        }

        private int numberOfContigs() {
            return referenceAlignments.length;
        }

        private static AlignedContigScore calculateAlignedContigScore(final AlignedContig ctg) {
            final List<AlignmentInterval> intervals = ctg.alignmentIntervals.stream()
                    .sorted(Comparator.comparing(ai -> ai.startInAssembledContig))
                    .collect(Collectors.toList());
            int totalReversals = 0;
            int totalIndels = 0;
            int totalMatches = 0;
            int totalMismatches = 0;
            int totalIndelLength = 0;
            for (int i = 0; i < intervals.size(); i++) {
                final AlignmentInterval ai = intervals.get(i);
                if (i > 0) {
                    final AlignmentInterval prev = intervals.get(i - 1);
                    if (prev.forwardStrand != ai.forwardStrand) {
                        totalReversals++;
                    } else {
                        final AlignmentInterval left = ai.forwardStrand ? prev : ai;
                        final AlignmentInterval right = ai.forwardStrand ? ai : prev;
                        if (left.referenceSpan.getEnd() < right.referenceSpan.getStart()) {
                            totalIndels++;
                            totalIndelLength += right.referenceSpan.getStart() - left.referenceSpan.getEnd();
                        }
                        if (left.endInAssembledContig < right.startInAssembledContig) {
                            totalIndels++;
                            totalIndelLength += right.startInAssembledContig - left.endInAssembledContig - 1;
                        }
                    }
                }
                final int matches = ai.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                        .filter(ce -> ce.getOperator().isAlignment())
                        .mapToInt(CigarElement::getLength).sum();
                final int misMatches = ai.mismatches;
                final int indelCount = (int) ai.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                        .filter(ce -> ce.getOperator().isIndel())
                        .count();
                final int indelLengthSum = ai.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                        .filter(ce -> ce.getOperator().isIndel())
                        .mapToInt(CigarElement::getLength).sum();
                totalIndels += indelCount;
                totalMatches += matches;
                totalMismatches += misMatches;
                totalIndelLength += indelLengthSum;
            }
            if (intervals.isEmpty()) {
                totalIndelLength += ctg.contigSequence.length;
                totalIndels++;
            } else {
                if (intervals.get(0).startInAssembledContig > 1) {
                    totalIndelLength += intervals.get(0).startInAssembledContig - 1;
                    totalIndels++;
                }
                if (intervals.get(intervals.size() - 1).endInAssembledContig < ctg.contigSequence.length) {
                    totalIndelLength += ctg.contigSequence.length - intervals.get(intervals.size() - 1).endInAssembledContig;
                    totalIndels++;
                }
            }
            return new AlignedContigScore(totalReversals, totalIndels, totalMatches, totalMismatches, totalIndelLength);
        }

        private static String calculateHPTag(final double referenceScore, final double alternativeScore) {
            if (Double.isNaN(referenceScore) && Double.isNaN(alternativeScore)) {
                return VCFConstants.EMPTY_INFO_FIELD;
            } else if (Double.isNaN(referenceScore)) {
                return ComposeStructuralVariantHaplotypesSpark.ALT_CONTIG_NAME;
            } else if (Double.isNaN(alternativeScore)) {
                return ComposeStructuralVariantHaplotypesSpark.REF_CONTIG_NAME;
            } else {
                switch (Double.compare(alternativeScore, referenceScore)) {
                    case 0:
                        return VCFConstants.EMPTY_INFO_FIELD;
                    case -1:
                        return ComposeStructuralVariantHaplotypesSpark.REF_CONTIG_NAME;
                    default:
                        return ComposeStructuralVariantHaplotypesSpark.ALT_CONTIG_NAME;
                }
            }
        }

        private static double calculateHPQualTag(final double referenceScore, final double alternativeScore) {
            if (Double.isNaN(referenceScore) || Double.isNaN(alternativeScore)) {
                return Double.NaN;
            } else {
                switch (Double.compare(alternativeScore, referenceScore)) {
                    case 0:
                        return 0;
                    case -1:
                        return referenceScore - alternativeScore;
                    default:
                        return alternativeScore - referenceScore;
                }
            }
        }

        private SAMRecord composeOutputSAMRecord(final SAMFileHeader outputHeader, final SVContext vc, final int index) {

            final AlignedContig originalContig = originalAlignments[index];
            final SAMRecord outputRecord = new SAMRecord(outputHeader);
            outputRecord.setAttribute(SAMTag.RG.name(), CONTIG_READ_GROUP);
            outputRecord.setAttribute(HAPLOTYPE_CALL_TAG, hpTagValue[index]);
            outputRecord.setAttribute(HAPLOTYPE_QUAL_TAG, "" + hpQualTagValue[index]);
            outputRecord.setAttribute(REFERENCE_SCORE_TAG, referenceScoreTagValue[index]);
            outputRecord.setAttribute(ALTERNATIVE_SCORE_TAG, alternativeScoreTagValue[index]);
            outputRecord.setAttribute(SAMTag.SA.name(), composeSupplementaryLikeString(originalContig) + ';');
            outputRecord.setAttribute(REFERENCE_ALIGNMENT_TAG, composeSupplementaryLikeString(referenceAlignments[index]) + ';');
            outputRecord.setAttribute(ALTERNATIVE_ALIGNMENT_TAG, composeSupplementaryLikeString(alternativeAlignments[index]) + ';');
            outputRecord.setAttribute(VARIANT_CONTEXT_TAG, vc.getUniqueID());
            outputRecord.setReadName(originalContig.contigName);
            outputRecord.setReadPairedFlag(false);
            outputRecord.setDuplicateReadFlag(false);
            outputRecord.setSecondOfPairFlag(false);
            outputRecord.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
            outputRecord.setReadNegativeStrandFlag(false);
            final byte[] bases = originalContig.contigSequence;
            outputRecord.setReadBases(bases);
            outputRecord.setReferenceName(vc.getContig());
            outputRecord.setAlignmentStart(vc.getStart());
            outputRecord.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
            outputRecord.setReadUnmappedFlag(true);
            return outputRecord;
        }

        private static String composeSupplementaryLikeString(final AlignedContig referenceAlignment) {
            return composeSupplementaryLikeString(referenceAlignment.alignmentIntervals);
        }

        private static String composeSupplementaryLikeString(final AlignmentInterval interval) {
            return composeSupplementaryLikeString(Collections.singletonList(interval));
        }

        private static String composeSupplementaryLikeString(final List<AlignmentInterval> intervals) {
            if (intervals.isEmpty()) {
                return new Cigar().toString();
            } else {
                return intervals.stream()
                        .map(AlignmentInterval::toSumpplementaryAlignmentString)
                        .collect(Collectors.joining(";"));
            }
        }

        private List<SAMRecord> composeAlignedOutputReferenceHaplotypeSAMRecords(final SAMFileHeader header, final SVContext variant) {
            return composeAlignedOutputHaplotypeSAMRecords(header, variant, referenceHaplotype);
        }

        private List<SAMRecord> composeAlignedOutputAlternativeHaplotypeSAMRecords(final SAMFileHeader header, final SVContext variant) {
            return composeAlignedOutputHaplotypeSAMRecords(header, variant, alternativeHaplotype);
        }

        private List<SAMRecord> composeAlignedOutputHaplotypeSAMRecords(SAMFileHeader header, SVContext variant, SVHaplotype haplotype) {
            final List<AlignmentInterval> intervals = haplotype.getReferenceAlignmentIntervals()
                    .stream().sorted(ALIGNMENT_INTERVAL_PRIMARY_SECONDARY_SORTER).collect(Collectors.toList());
            final List<SAMRecord> result = new ArrayList<>(intervals.size());
            final String vcUID = variant.getUniqueID();
            for (int i = 0; i < intervals.size(); i++) {
                final AlignmentInterval interval = intervals.get(i);
                final SAMRecord record = intervals.get(i).toSAMRecord(header, haplotype.getName(), haplotype.getBases(), i != 0, SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue(), null);
                record.setReadNegativeStrandFlag(!interval.forwardStrand);
                record.setReadName(haplotype.getName());
                record.setAttribute(SAMTag.RG.name(), HAPLOTYPE_READ_GROUP);
                record.setAttribute(HAPLOTYPE_CALL_TAG, haplotype.getName());
                record.setAttribute(HAPLOTYPE_QUAL_TAG, "" + 60);
                record.setAttribute(VARIANT_CONTEXT_TAG, vcUID);
                record.setSupplementaryAlignmentFlag(i != 0);
                result.add(record);
            }
            final List<String> saTagValues = intervals.stream()
                    .map(AlignmentInterval::toSumpplementaryAlignmentString)
                    .collect(Collectors.toList());
            if (result.size() > 1) {
                for (int i = 0; i < result.size(); i++) {
                    final int idx = i;
                    result.get(i).setAttribute(SAMTag.SA.name(),
                            IntStream.range(0, result.size())
                                    .filter(ii -> ii != idx)
                                    .mapToObj(saTagValues::get)
                                    .collect(Collectors.joining(";")) + ";");
                }
            }
            return result;
        }

        private List<SAMRecord> composeAlignedOutputSAMRecords(final SAMFileHeader header,
                                                           final SVContext vc, final int index) {
            final AlignedContig alignment = originalAlignments[index];
            final List<SAMRecord> result = new ArrayList<>(alignment.alignmentIntervals.size());
            result.add(alignment.alignmentIntervals.get(0).toSAMRecord(header, alignment.contigName, alignment.contigSequence, false, 0, Collections.emptyList()));
            final String vcUID = vc.getUniqueID();
            for (int i = 1; i < alignment.alignmentIntervals.size(); i++) {
                result.add(alignment.alignmentIntervals.get(i).toSAMRecord(header, alignment.contigName, alignment.contigSequence, true, SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue(), Collections.emptyList()));
            }
            for (final SAMRecord record : result) {
                record.setReadName(alignment.contigName);
                record.setAttribute(SAMTag.RG.name(), CONTIG_READ_GROUP);
                record.setAttribute(HAPLOTYPE_CALL_TAG, hpTagValue[index]);
                record.setAttribute(HAPLOTYPE_QUAL_TAG, "" + hpQualTagValue[index]);
                record.setAttribute(REFERENCE_SCORE_TAG, referenceScoreTagValue[index]);
                record.setAttribute(ALTERNATIVE_SCORE_TAG, alternativeScoreTagValue[index]);
                record.setAttribute(VARIANT_CONTEXT_TAG, vcUID);
                record.setAttribute(REFERENCE_ALIGNMENT_TAG, composeSupplementaryLikeString(referenceAlignments[index]) + ';');
                record.setAttribute(ALTERNATIVE_ALIGNMENT_TAG, composeSupplementaryLikeString(alternativeAlignments[index]) + ';');
            }
            final List<String> saTagValues = alignment.alignmentIntervals.stream()
                    .map(AlignmentInterval::toSumpplementaryAlignmentString)
                    .collect(Collectors.toList());

            if (result.size() > 1) {
                for (int i = 0; i < result.size(); i++) {
                    final int idx = i;
                    result.get(i).setAttribute(SAMTag.SA.name(),
                            IntStream.range(0, result.size())
                                    .filter(ii -> ii != idx)
                                    .mapToObj(saTagValues::get)
                                    .collect(Collectors.joining(";")) + ";");
                }
            }
            return result;
        }

        private SAMRecord composeOutputReferenceHaplotypeRecord(final SAMFileHeader outputHeader, final SVContext variant) {
            return composeOutputReferenceHaplotypeRecord(outputHeader, referenceHaplotype, variant);
        }

        private SAMRecord composeOutputAlternativeHaplotypeRecord(final SAMFileHeader outputHeader, final SVContext variant) {
            return composeOutputReferenceHaplotypeRecord(outputHeader, alternativeHaplotype, variant);
        }

        private SAMRecord composeOutputReferenceHaplotypeRecord(final SAMFileHeader outputHeader, final SVHaplotype haplotype,
                                                                final SVContext variant) {
                final SAMRecord result = new SAMRecord(outputHeader);
            result.setReadName(haplotype.getName());
            result.setReadName(haplotype.getName());
            result.setReadUnmappedFlag(true);
            result.setReferenceName(variant.getContig());
            result.setAlignmentStart(variant.getStart());
            result.setAttribute(SAMTag.RG.name(), HAPLOTYPE_READ_GROUP);
            result.setAttribute(HAPLOTYPE_CALL_TAG, haplotype.getName());
            result.setAttribute(HAPLOTYPE_QUAL_TAG, "" + 60);
            result.setAttribute(REFERENCE_SCORE_TAG, "" + 0);
            result.setAttribute(ALTERNATIVE_SCORE_TAG, "" + 60);
            result.setAttribute(SAMTag.SA.name(), composeSupplementaryLikeString(haplotype.getReferenceAlignmentIntervals()) + ';');
            result.setAttribute(VARIANT_CONTEXT_TAG, variant.getUniqueID());
            result.setReadPairedFlag(false);
            result.setDuplicateReadFlag(false);
            result.setSecondOfPairFlag(false);
            result.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
            result.setReadNegativeStrandFlag(false);
            final byte[] bases = haplotype.getBases();
            result.setReadBases(bases);
            result.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
            return result;
        }

        public static class Serializer extends com.esotericsoftware.kryo.Serializer<VariantHaplotypesAndContigsComposite> {

            @Override
            public void write(final Kryo kryo, final Output output, final VariantHaplotypesAndContigsComposite object) {
                final int numberOfContigs = object.numberOfContigs();
                output.writeInt(numberOfContigs);
                kryo.writeClassAndObject(output, object.referenceHaplotype);
                kryo.writeClassAndObject(output, object.alternativeHaplotype);
                for (int i = 0; i < numberOfContigs; i++) {
                    kryo.writeObject(output, object.originalAlignments[i]);
                    kryo.writeObject(output, object.referenceAlignments[i]);
                    kryo.writeObject(output, object.alternativeAlignments[i]);
                    output.writeString(object.referenceScoreTagValue[i]);
                    output.writeString(object.alternativeScoreTagValue[i]);
                    output.writeString(object.hpTagValue[i]);
                    output.writeDouble(object.hpQualTagValue[i]);
                }
            }

            @Override
            public VariantHaplotypesAndContigsComposite read(final Kryo kryo, final Input input, Class<VariantHaplotypesAndContigsComposite> type) {
                return new VariantHaplotypesAndContigsComposite(kryo, input);
            }
        }
    }
}

