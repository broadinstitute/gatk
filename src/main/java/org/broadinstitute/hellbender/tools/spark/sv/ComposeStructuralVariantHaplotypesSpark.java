package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.gcs.BamBucketIoUtils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.report.GATKReportColumnFormat;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Composes a bam file with assembled contigs conveniently realigned and annotated for the SV indel genotyping.
 *
 * <p>
 *     Currently only SV indels are supported: {@code SVTYPE=INV} or {@code SVTYPE=DEL}.
 * </p>
 * <p>
 *     The resulting alignment file (BAM or SAM formatted base on the extension of the name provided)
 *     will be sorted by. Each supported variant position will be overlapped by a number of records:
 * </p>
 * <p>
 *     There are always two synthetic records representing the referene haplotype and the alternative allele haplotype
 *     for that variant. These receive names {@code var_id:ref} and {@code var_id:alt} respectively where {@code var_id}
 *     identifies the variant by its coordinates ({@code var_chrName_startPosition}).
 *     <br/>
 *     These records have as read group {@code {@value #HAPLOTYPE_READ_GROUP}}.
 * </p>
 * <p>
 *     Then for each assembly contig that overlapped that variant we have a number of records which depends on how these
 *     align against the reference (a single record, or several supplementary alignments). Their id is {@code var_id:assembly_id:contig_id}.
 *     <br/>
 *     These records have as read group {@code {@value #CONTIG_READ_GROUP}}.
 *     <br/>
 *     Also the have a number of annotations indicating what haplotype (ref or alt) the seem to support and with what
 *     confidence:
 *     <dl>
 *         <dt>HP</dt><dd>the haplotype the support, {@code 'ref'} for reference, {@code 'alt'} for alternative and {@code '.'} for neither. (e.g. {@code "HP:Z:ref", "HP:Z:alt", "HP:Z:."}).</dd>
 *         <dt>RS</dt><dd>the reference support score screen (e.g. {@code "RS:Z:-100,55,1,2,30,1"}).</dd>
 *         <dt>XS</dt><dd>the alternative allele score screen (e.g. {@code "XS:Z:-1241,134,2,1,5,1}).</dd>
 *         <dt>HQ</dt><dd>the confidence in the support for the haplotype in {@code "HP"} (e.g. {@code "HQ:Z:10.0"}).
 *                        This value is equal to the difference between the score for the supported and the other haplotype</dd></dt>
 *         <dt>VC</dt><dd>coordinate of the targeted variant {@code chr:pos}</dd>
 *     </dl>
 * </p>
 */
@CommandLineProgramProperties(
        summary = "composes contigs file for genotyping",
        oneLineSummary = "composes contigs file for genotyping",
        usageExample = "--variant ins_and_dels.vcf --contigs assemblies.bam --reference my-ref.fa " +
                       "--shardSize 10000 --paddingSize 300 --output genotyping-contigs.bam",
        programGroup = StructuralVariationSparkProgramGroup.class )
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

    public static final int DEFAULT_SHARD_SIZE = 10_000;
    public static final int DEFAULT_PADDING_SIZE = 50;
    private static final int FASTA_BASES_PER_LINE = 60;

    public static final String HAPLOTYPE_READ_GROUP = "HAP";
    public static final String CONTIG_READ_GROUP = "CTG";
    public static final String HAPLOTYPE_CALL_TAG = "HP";
    public static final String HAPLOTYPE_QUAL_TAG = "HQ";
    public static final String REFERENCE_SCORE_TAG = "RS";
    public static final String ALTERNATIVE_SCORE_TAG = "XS";
    public static final String VARIANT_CONTEXT_TAG = "VC";
    public static final String ALIGNMENT_SEQ_NAME = "seq";


    @Argument(doc = "shard size",
              shortName = SHARD_SIZE_SHORT_NAME,
              fullName = SHARD_SIZE_FULL_NAME,
    optional = true)
    private int shardSize = DEFAULT_SHARD_SIZE;

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

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        Utils.nonNull(ctx);
        Utils.nonNull(alignedContigsFileName);
        final ReadsSparkSource alignedContigs = new ReadsSparkSource(ctx);
        final VariantsSparkSource variantsSource = new VariantsSparkSource(ctx);

        final JavaRDD<GATKRead> contigs = alignedContigs.getParallelReads(alignedContigsFileName, referenceArguments.getReferenceFileName(), getIntervals());
        final JavaRDD<StructuralVariantContext> variants = variantsSource
                .getParallelVariantContexts(variantsFileName, getIntervals())
                .filter(ComposeStructuralVariantHaplotypesSpark::supportedVariant)
                .map(vc -> StructuralVariantContext.create(vc));

        final JavaPairRDD<StructuralVariantContext, List<GATKRead>> variantOverlappingContigs = composeOverlappingContigRecordsPerVariant(ctx, contigs, variants);
        processVariants(ctx, variantOverlappingContigs, getReferenceSequenceDictionary(), alignedContigs);

    }

    private static boolean supportedVariant(final VariantContext vc) {
        final StructuralVariantType type = vc.getStructuralVariantType();
        if (type == null) {
            return false;
        } else if (type == StructuralVariantType.INS || type == StructuralVariantType.DEL) {
            return vc.getAlternateAlleles().size() == 1;
        } else  {
            return false;
        }
    }

    private JavaPairRDD<StructuralVariantContext, List<GATKRead>> composeOverlappingContigRecordsPerVariant(final JavaSparkContext ctx, final JavaRDD<GATKRead> contigs, final JavaRDD<StructuralVariantContext> variants) {
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        final List<SimpleInterval> intervals = hasIntervals() ? getIntervals() : IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
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

        final JavaPairRDD<SimpleInterval, List<Tuple2<SimpleInterval, StructuralVariantContext>>> variantsInShards =
            groupInShards(variants, (v) -> v.getBreakPointIntervals(paddingSize, dictionaryBroadcast.getValue()), shardIntervalsBroadcast);

        final JavaPairRDD<SimpleInterval, Tuple2<List<Tuple2<SimpleInterval, GATKRead>>, List<Tuple2<SimpleInterval, StructuralVariantContext>>>> contigAndVariantsInShards =
                contigsInShards.join(variantsInShards);


        final JavaPairRDD<StructuralVariantContext, List<GATKRead>> contigsPerVariantInterval =
                contigAndVariantsInShards.flatMapToPair(t -> {
                    final List<Tuple2<SimpleInterval, StructuralVariantContext>> vars = t._2()._2();
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

        final Function<StructuralVariantContext, String> variantId = (Function<StructuralVariantContext, String> & Serializable) ComposeStructuralVariantHaplotypesSpark::variantId;
        final Function2<List<GATKRead>, List<GATKRead>, List<GATKRead>> readListMerger = (a, b) -> Stream.concat(a.stream(), b.stream()).collect(Collectors.toList());

        // Merge contig lists on the same variant-context coming from different intervals
        // into one.
        final JavaPairRDD<StructuralVariantContext, List<GATKRead>> contigsPerVariant = contigsPerVariantInterval
                .mapToPair(t -> new Tuple2<>(variantId.apply(t._1()), t))
                .reduceByKey((a, b) -> new Tuple2<>(a._1(), readListMerger.call(a._2(), b._2())))
                .mapToPair(Tuple2::_2);
        return contigsPerVariant;
    }

    private static String variantId(final StructuralVariantContext variant) {
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

    private <T> JavaPairRDD<SimpleInterval, List<Tuple2<SimpleInterval, T>>> groupInShards(final JavaRDD<T> elements, final org.apache.spark.api.java.function.Function<T, List<SimpleInterval>> intervalsOf,
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

    protected void processVariants(final JavaSparkContext ctx,
                                   final JavaPairRDD<StructuralVariantContext, List<GATKRead>> variantsAndOverlappingContigRecords,
                                   final SAMSequenceDictionary dictionary, final ReadsSparkSource readSource) {

        final SAMFileHeader outputHeader = new SAMFileHeader();
        final SAMProgramRecord programRecord = new SAMProgramRecord(getProgramName());
        programRecord.setCommandLine(getCommandLine());
        outputHeader.setSequenceDictionary(dictionary);
        outputHeader.addProgramRecord(programRecord);
        outputHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        outputHeader.addReadGroup(new SAMReadGroupRecord(HAPLOTYPE_READ_GROUP));
        outputHeader.addReadGroup(new SAMReadGroupRecord(CONTIG_READ_GROUP));
        final SAMFileWriter outputWriter = BamBucketIoUtils.makeWriter(outputFileName, outputHeader, false);

        final JavaPairRDD<StructuralVariantContext, List<AlignedContig>> variantsAndOverlappingUniqueContigs
                = variantsAndOverlappingContigRecords
                .mapValues(l -> l.stream().collect(Collectors.groupingBy(GATKRead::getName)))
                .mapValues(m -> m.values().stream()
                        .map(l -> l.stream()
                                .map(ComposeStructuralVariantHaplotypesSpark::convertToAlignedContig)
                                .reduce(ComposeStructuralVariantHaplotypesSpark::mergeAlignedContigs).orElseThrow(IllegalStateException::new))
                        .collect(Collectors.toList()));
        Utils.stream(variantsAndOverlappingUniqueContigs.toLocalIterator())
                .map(t -> resolvePendingContigs(t, readSource))
                .forEach(t -> {
                    final StructuralVariantContext vc = t._1();
                    final List<AlignedContig> contigs = t._2();
                    final int maxLength = contigs.stream()
                            .mapToInt(a -> a.contigSequence.length)
                            .max().orElse(paddingSize);
                    final Haplotype referenceHaplotype = vc.composeHaplotypeBasedOnReference(0, maxLength * 2, getReference());
                    //referenceHaplotype.setGenomeLocation(null);
                    final Haplotype alternativeHaplotype = vc.composeHaplotypeBasedOnReference(1, maxLength * 2, getReference());
                    //alternativeHaplotype.setGenomeLocation(null);
                    final String idPrefix = String.format("var_%s_%d", vc.getContig(), vc.getStart());
                    outputHaplotypesAsSAMRecords(outputHeader, outputWriter, referenceHaplotype, alternativeHaplotype, idPrefix);
                    final Map<String, AlignedContig> referenceAlignedContigs = alignContigsAgainstHaplotype(ctx, referenceHaplotype, contigs);
                    final Map<String, AlignedContig> alternativeAlignedContigs = alignContigsAgainstHaplotype(ctx, alternativeHaplotype, contigs);
                    contigs.forEach(contig -> {
                        final AlignedContig referenceAlignment = referenceAlignedContigs.get(contig.contigName);
                        final AlignedContig alternativeAlignment = alternativeAlignedContigs.get(contig.contigName);
                        final AlignedContigScore referenceScore = calculateAlignedContigScore(referenceAlignment);
                        final AlignedContigScore alternativeScore = calculateAlignedContigScore(alternativeAlignment);
                        final String hpTagValue = calculateHPTag(referenceScore.getValue(), alternativeScore.getValue());
                        final double hpQualTagValue = calculateHPQualTag(referenceScore.getValue(), alternativeScore.getValue());
                        if (!referenceAlignment.alignmentIntervals.isEmpty()) {
                            final List<SAMRecord> records = convertToSAMRecords(referenceAlignment, outputHeader, t._1().getContig(), referenceHaplotype.getGenomeLocation().getStart(), idPrefix, hpTagValue, hpQualTagValue, referenceScore, alternativeScore, t._1().getStart());
                            records.forEach(outputWriter::addAlignment);
                        } else {
                            final SAMRecord outputRecord = convertToUnmappedSAMRecord(outputHeader, t._1().getContig(), t._1().getStart(), idPrefix, contig, referenceScore, alternativeScore, hpTagValue, hpQualTagValue);
                            outputWriter.addAlignment(outputRecord);
                        }
                    });
                });
        outputWriter.close();
    }

    private SAMRecord convertToUnmappedSAMRecord(final SAMFileHeader outputHeader, final StructuralVariantContext vc,
                                                 final AlignedContig originalContig, final AlignedContigScore referenceScore,
                                                 final AlignedContigScore alternativeScore,
                                                 final AlignedContig referenceAlignment, final AlignedContig alternativeAlignment,
                                                 final String hpTagValue, double hpQualTagValue) {
        final SAMRecord outputRecord = new SAMRecord(outputHeader);
        outputRecord.setAttribute(SAMTag.RG.name(), CONTIG_READ_GROUP);
        outputRecord.setAttribute(HAPLOTYPE_CALL_TAG, hpTagValue);
        outputRecord.setAttribute(HAPLOTYPE_QUAL_TAG, "" + hpQualTagValue);
        outputRecord.setAttribute(REFERENCE_SCORE_TAG, "" + referenceScore);
        outputRecord.setAttribute(ALTERNATIVE_SCORE_TAG, "" + alternativeScore);
        outputRecord.setAttribute(VARIANT_CONTEXT_TAG, "" + referenceContig + ":" + start);
        outputRecord.setReadName(idPrefix + ":" + originalContig.contigName);
        outputRecord.setReadPairedFlag(false);
        outputRecord.setDuplicateReadFlag(false);
        outputRecord.setSecondOfPairFlag(false);
        outputRecord.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
        outputRecord.setReadNegativeStrandFlag(false);
        final byte[] bases = originalContig.contigSequence;
        outputRecord.setReadBases(bases);
        outputRecord.setReferenceName(referenceContig);
        outputRecord.setAlignmentStart(start);
        outputRecord.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
        outputRecord.setReadUnmappedFlag(true);
        return outputRecord;
    }

    private static List<SAMRecord> convertToSAMRecords(final AlignedContig alignment, final SAMFileHeader header, final String referenceContig, final int contigStart, final String idPrefix, final String hpTagValue, final double hpQualTagValue, final AlignedContigScore referenceScore, final AlignedContigScore alternativeScore, final int variantStart) {
        final List<SAMRecord> result = new ArrayList<>(alignment.alignmentIntervals.size());
        result.add(alignment.alignmentIntervals.get(0).toSAMRecord(header, alignment, false));
        for (int i = 1; i < alignment.alignmentIntervals.size(); i++) {
            result.add(alignment.alignmentIntervals.get(i).toSAMRecord(header, alignment, true));
        }
        for (final SAMRecord record : result) {
            record.setReadName(idPrefix + ":" + alignment.contigName);
            record.setReferenceName(referenceContig);
            record.setAlignmentStart(record.getAlignmentStart() + contigStart - 1);
            record.setAttribute(SAMTag.RG.name(), CONTIG_READ_GROUP);
            record.setAttribute(HAPLOTYPE_CALL_TAG, hpTagValue);
            record.setAttribute(HAPLOTYPE_QUAL_TAG, "" + hpQualTagValue);
            record.setAttribute(REFERENCE_SCORE_TAG, "" + referenceScore);
            record.setAttribute(ALTERNATIVE_SCORE_TAG, "" + alternativeScore);
            record.setAttribute(VARIANT_CONTEXT_TAG, "" + referenceContig + ":" + variantStart);
        }
        final List<String> saTagValues = result.stream()
                .map(record ->
                    Utils.join(",", record.getReferenceName(), record.getAlignmentStart(), record.getReadNegativeStrandFlag() ? "+" : "-", record.getCigarString(), record.getMappingQuality(), "" + record.getAttribute(SAMTag.NM.name()))
                )
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

    private AlignedContigScore calculateAlignedContigScore(final AlignedContig ctg) {
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

    private String calculateHPTag(final double referenceScore, final double alternativeScore) {
        if (Double.isNaN(referenceScore) && Double.isNaN(alternativeScore)) {
            return ".";
        } else if (Double.isNaN(referenceScore)) {
            return "alt";
        } else if (Double.isNaN(alternativeScore)) {
            return "ref";
        } else {
            switch (Double.compare(alternativeScore, referenceScore)) {
                case 0: return ".";
                case -1: return "ref";
                default: return "alt";
            }
        }
    }

    private double calculateHPQualTag(final double referenceScore, final double alternativeScore) {
        if (Double.isNaN(referenceScore) || Double.isNaN(alternativeScore)) {
            return Double.NaN;
        } else {
            switch (Double.compare(alternativeScore, referenceScore)) {
                case 0: return 0;
                case -1: return referenceScore - alternativeScore;
                default: return alternativeScore - referenceScore;
            }
        }
    }

    private final Map<String,AlignedContig> alignContigsAgainstHaplotype(final JavaSparkContext ctx, final Haplotype haplotype, final List<AlignedContig> contigs) {
        File fastaFile = null;
        File imageFile = null;
        final SAMFileHeader alignmentHeader = new SAMFileHeader();
        alignmentHeader.addSequence(new SAMSequenceRecord(ALIGNMENT_SEQ_NAME, haplotype.length()));
        try {

            fastaFile = createFastaFromHaplotype(haplotype);
            imageFile = createImageFileFromFasta(fastaFile);
            final Stream<AlignedContig> alignedContigSegments;
            if (contigs.size() > 10) { // just bother to sparkify it if there is a decent number of contigs.
                final Map<String, GATKRead> contigsByName = contigs.stream()
                        .collect(Collectors.toMap(contig -> contig.contigName, contig -> convertToUnmappedGATKRead(contig, alignmentHeader)));
                final BwaSparkEngine bwa = new BwaSparkEngine(ctx, imageFile.getPath(), alignmentHeader, alignmentHeader.getSequenceDictionary());
                alignedContigSegments = bwa.alignSingletons(ctx.parallelize(new ArrayList<>(contigsByName.values())))
                        .mapToPair(r -> new Tuple2<>(r.getName(), r))
                        .groupByKey()
                        .map(t -> Utils.stream(t._2()).collect(Collectors.toList()))
                        .collect()
                        .stream()
                        .map(l ->  l.stream()
                                    .map(ComposeStructuralVariantHaplotypesSpark::convertToAlignedContig)
                                    .reduce(ComposeStructuralVariantHaplotypesSpark::mergeAlignedContigs).get());
                bwa.close();
            } else {
                final BwaMemIndex index = new BwaMemIndex(imageFile.getPath());
                final BwaMemAligner aligner = new BwaMemAligner(index);
                final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(contigs, a -> a.contigSequence);
                final List<String> refNames = alignmentHeader.getSequenceDictionary().getSequences().stream().map(SAMSequenceRecord::getSequenceName).collect(Collectors.toList());
                alignedContigSegments = IntStream.range(0, contigs.size())
                        .mapToObj(i -> new Tuple2<>(contigs.get(i), alignments.get(i)))
                        .map(t -> new Tuple2<>(t._1(), t._2().stream()
                                .filter(bwa -> bwa.getRefId() >= 0) // remove unmapped reads.
                                .filter(bwa -> SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(bwa.getSamFlag())) // ignore secondary alignments.
                                .map(bma -> new AlignmentInterval(Utils.nonNull(bma), refNames, t._1().contigSequence.length))
                                .collect(Collectors.toList())))
                        .map(t -> new AlignedContig(t._1().contigName, t._1().contigSequence, t._2()));
            }
            return alignedContigSegments.collect(Collectors.toMap(a ->a.contigName, a -> a));
                   // .map(a -> {
                       // final Tuple2<List<AlignmentInterval>, Double> bestConfiguration =
                     //           PlaygroundExtract.pickAnyBestConfiguration(a, Collections.singleton("seq"));
                   //     return new AlignedContig(a.contigName, a.contigSequence, bestConfiguration._1(), bestConfiguration._2());
                   // })
                   // .collect(Collectors.toMap(a -> a.contigName, a -> a));
        } finally {
            Stream.of(fastaFile, imageFile)
                    .filter(Objects::nonNull)
                    .forEach(File::delete);
        }
    }

    private File createFastaFromHaplotype(final Haplotype haplotype) {
        final File result;
        try {
            result = File.createTempFile("gatk-sv-bwa-tmp-ref-", ".fasta");
        } catch (final IOException ex) {
            throw new GATKException("cound not optain a file location for a temporary haplotype reference fasta file", ex);
        }
        result.deleteOnExit();
        try (final PrintWriter fastaWriter = new PrintWriter(new FileWriter(result))) {

            fastaWriter.println(">" + ALIGNMENT_SEQ_NAME);
            final byte[] bases = haplotype.getBases();
            int nextIdx = 0;
            while (nextIdx < bases.length) {
                fastaWriter.println(new String(bases, nextIdx, Math.min(bases.length - nextIdx, FASTA_BASES_PER_LINE)));
                nextIdx += FASTA_BASES_PER_LINE;
            }
        } catch (final IOException ex) {
            throw new GATKException("could not write haplotype reference fasta file '" + result + "'", ex);
        }
        return result;
    }

    private File createImageFileFromFasta(final File fasta) {
        final File result = new File(fasta.getPath().replaceAll("\\.*$",".img"));
        try {
            BwaMemIndex.createIndexImageFromFastaFile(fasta.getPath(), result.getPath());
        } catch (final Throwable ex) {
            throw new GATKException("problems indexing fasta file '" + fasta + "' into '" + result + "'", ex);
        }
        result.deleteOnExit();
        return result;
    }

    private void outputHaplotypesAsSAMRecords(final SAMFileHeader outputHeader,
                                              final SAMFileWriter outputWriter,
                                              final Haplotype referenceHaplotype,
                                              final Haplotype alternativeHaplotype,
                                              final String idPrefix) {
        final Consumer<SAMRecord> haplotypeExtraSetup = r -> {
            //r.setReferenceName(t._1().getContig());
            //
            // r.setAlignmentStart(t._1().getStart());
            r.setAttribute(SAMTag.RG.name(), HAPLOTYPE_READ_GROUP);
            r.setMappingQuality(60);
        };
        outputWriter.addAlignment(referenceHaplotype.toSAMRecord(outputHeader, idPrefix + ":ref",
                haplotypeExtraSetup));
        outputWriter.addAlignment(alternativeHaplotype.toSAMRecord(outputHeader, idPrefix + ":alt",
                haplotypeExtraSetup));
    }

    private <V extends VariantContext> Tuple2<V, List<AlignedContig>> resolvePendingContigs(final Tuple2<V, List<AlignedContig>> vc, final ReadsSparkSource s) {
        logger.debug("VC " + vc._1().getContig() + ":" + vc._1().getStart() + "-" + vc._1().getEnd());
        final List<AlignedContig> updatedList = vc._2().stream()
                .map(a -> {
                    if (a.contigSequence[0] != 0 && a.contigSequence[a.contigSequence.length - 1] != 0) {
                        return a;
                    } else {
                        final String targetName = a.contigName;
                        final List<SimpleInterval> locations = a.alignmentIntervals.stream()
                                .map(ai -> ai.referenceSpan).collect(Collectors.toList());
                        final GATKRead candidate = s.getParallelReads(alignedContigsFileName, referenceArguments.getReferenceFileName(),
                                        new TraversalParameters(locations, false))
                                .filter(rr -> rr.getName().equals(targetName))
                                .filter(rr -> !rr.getCigar().containsOperator(CigarOperator.H))
                                .first();
                        final AlignedContig result = addAlignment(a, candidate);
                        if (result.contigSequence[0] != 0 && a.contigSequence[result.contigSequence.length - 1] != 0) {
                            logger.warn("Contig " + result.contigName + " " + readAlignmentString(candidate) + " "
                                    + candidate.getAttributeAsString("SA") + " gave-up!");
                            return null;
                        } else {
                            return result;
                        }
                    }})
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
        return new Tuple2<>(vc._1, updatedList);
    }

    private static String readAlignmentString(final GATKRead r) {
        return r.getContig() + "," + r.getStart() + "," + (r.isReverseStrand() ?  "-" : "+") + "," + r.getCigar() + ","
                + r.getMappingQuality() + "," + r.getAttributeAsString("NM");
    }

    private static AlignedContig convertToAlignedContig(final GATKRead read) {
        if (read.isUnmapped()) {
            return new AlignedContig(read.getName(), read.getBases(), Collections.emptyList());
        } else {
            final int contigLength = CigarUtils.readLength(read.getCigar());
            final byte[] contigBases = new byte[contigLength];
            final byte[] readBases = read.getBases();
            if (read.isReverseStrand()) {
                SequenceUtil.reverseComplement(readBases);
            }
            final int startPosition = 1 + (read.isReverseStrand() ? read.getRightHardClipLength() : read.getLeftHardClipLength());
            final int endPosition = contigLength - (read.isReverseStrand() ? read.getLeftHardClipLength() : read.getRightHardClipLength());
            for (int contigIdx = startPosition - 1, readIdx = 0; contigIdx < endPosition; contigIdx++, readIdx++) {
                if (contigBases[contigIdx] == 0) {
                    contigBases[contigIdx] = readBases[readIdx];
                } else if (contigBases[contigIdx] != readBases[readIdx]) {
                    throw new IllegalArgumentException("it seems that there is a base call conflict between overlapping alternative read alignments: " + read.getName() + " at " + (contigIdx + 1));
                }
            }
            return new AlignedContig(read.getName(), contigBases, Collections.singletonList(new AlignmentInterval(read)), false);
        }
    }

    private static GATKRead convertToUnmappedGATKRead(final AlignedContig contig, final SAMFileHeader header) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadBases(contig.contigSequence);
        record.setReadName(contig.contigName);
        record.setReadUnmappedFlag(true);
        return new SAMRecordToGATKReadAdapter(record);
    }

    private static AlignedContig addAlignment(final AlignedContig contig, final GATKRead read) {
        return mergeAlignedContigs(contig, convertToAlignedContig(read));
    }

    private static AlignedContig mergeAlignedContigs(final AlignedContig a, final AlignedContig b) {
        if (!a.contigName.equals(b.contigName)) {
            throw new IllegalArgumentException("trying to merge contigs with different names");
        }
        final byte[] aBases = a.contigSequence;
        final byte[] bBases = b.contigSequence;
        if (aBases.length != bBases.length) {
            throw new IllegalArgumentException("same contig cannot have different lengths!");
        }
        final byte[] cBases = new byte[aBases.length];
        for (int i = 0; i < aBases.length; i++) {
            final byte aBase = aBases[i];
            final byte bBase = bBases[i];
            if (aBase == bBase) {
                cBases[i] = aBase;
            } else if (aBase == 0) {
                cBases[i] = bBase;
            } else if (bBase == 0) {
                cBases[i] = aBase;
            } else {
                throw new IllegalStateException("conflict between read calls on the same contig " + a.contigName);
            }
        }
        return new AlignedContig(a.contigName, cBases,
                Stream.concat(a.alignmentIntervals.stream(), b.alignmentIntervals.stream())
                      .sorted(Comparator.comparing(ai -> ai.startInAssembledContig))
                      .collect(Collectors.toList()), false);
    }
}

