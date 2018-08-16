package org.broadinstitute.hellbender.engine.spark;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.ReferenceShard;
import org.broadinstitute.hellbender.engine.VariantShard;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.KV;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.testutils.FakeReferenceSource;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.MinimalVariant;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * ReadsPreprocessingPipelineTestData contains coordinated test data that can be used in the many transforms that
 * are a part of the ReadsPreprocessingPipeline.
 */
public class ReadsPreprocessingPipelineSparkTestData {
    private final List<KV<Integer, Integer>> readStartLength;

    private final List<GATKRead> reads;
    private final List<SimpleInterval> readIntervals;
    private final List<SimpleInterval> allIntervals;
    private final List<GATKVariant> variants;
    private final List<KV<GATKRead, ReferenceBases>> kvReadsRefBases;
    private final List<KV<GATKRead, GATKVariant>> kvReadVariant;
    private final List<KV<GATKRead, Iterable<GATKVariant>>> kvReadiVariant;
    private final List<KV<GATKRead, ReadContextData>> kvReadContextData;


    /**
     * ReadsPreprocessingPipelineTestData holds a bunch of connected data for testing classes that work with
     * reads, variants, references bases and pairing those types together.
     * @param clazz The class to be used to back the GATKRead, either Read.class, or SAMRecord.class.
     */
    public ReadsPreprocessingPipelineSparkTestData(Class<?> clazz) {
        readStartLength = Arrays.asList(KV.of(100, 50), KV.of(140, 100),
                KV.of(ReferenceShard.REFERENCE_SHARD_SIZE, 10),
                KV.of(3*ReferenceShard.REFERENCE_SHARD_SIZE - 1, 10));

        reads = Lists.newArrayList(
                makeRead("1", readStartLength.get(0), 1, clazz),
                makeRead("1", readStartLength.get(1), 2, clazz),
                makeRead("1", readStartLength.get(2), 3, clazz),
                makeRead("1", readStartLength.get(3), 4, clazz),
                makeRead("2", readStartLength.get(2), 5, clazz)
                );

        readIntervals = Lists.newArrayList(
                makeInterval("1", readStartLength.get(0)),
                makeInterval("1", readStartLength.get(1)),
                makeInterval("1", readStartLength.get(2)),
                makeInterval("1", readStartLength.get(3)),
                makeInterval("2", readStartLength.get(2))
                );

        // The first two reads are mapped onto the same reference shard. The ReferenceBases returned should
        // be from the start of the first read [rStartLength.get(0).getKey()] to the end
        // the second [rStartLength.get(1).getKey() + rStartLength.get(1).getValue()-1].
        SimpleInterval spannedReadInterval =
                new SimpleInterval("1", readStartLength.get(0).getKey(), readStartLength.get(1).getKey() + readStartLength.get(1).getValue()-1);

        allIntervals = Lists.newArrayList(readIntervals.iterator());
        allIntervals.add(spannedReadInterval);

        kvReadsRefBases = Arrays.asList(
                KV.of(reads.get(0), getBases("1", reads.get(0).getStart(), reads.get(0).getEnd())),
                KV.of(reads.get(1), getBases("1", reads.get(1).getStart(), reads.get(1).getEnd())),
                KV.of(reads.get(2), getBases("1", reads.get(2).getStart(), reads.get(2).getEnd())),
                KV.of(reads.get(3), getBases("1", reads.get(3).getStart(), reads.get(3).getEnd())),
                KV.of(reads.get(4), getBases("2", reads.get(4).getStart(), reads.get(4).getEnd()))
        );

        variants = Lists.newArrayList(
                new MinimalVariant(new SimpleInterval("1", 170, 180), true, false),
                new MinimalVariant(new SimpleInterval("1", 210, 220), false, true),
                new MinimalVariant(new SimpleInterval("1", ReferenceShard.REFERENCE_SHARD_SIZE,
                        ReferenceShard.REFERENCE_SHARD_SIZE), true, false),
                new MinimalVariant(new SimpleInterval("1", 3 * ReferenceShard.REFERENCE_SHARD_SIZE - 2,
                        3 * ReferenceShard.REFERENCE_SHARD_SIZE + 2), false, true),
                new MinimalVariant(new SimpleInterval("2", ReferenceShard.REFERENCE_SHARD_SIZE,
                        ReferenceShard.REFERENCE_SHARD_SIZE), false, true)
        );

        kvReadVariant = Arrays.asList(
                KV.of(reads.get(1), variants.get(0)),
                KV.of(reads.get(1), variants.get(1)),
                KV.of(reads.get(2), variants.get(2)),
                KV.of(reads.get(3), variants.get(3)),    // The read and variant span two variant shards, that's
                KV.of(reads.get(3), variants.get(3)),     // why there are two of them (2,3).
                KV.of(reads.get(4), variants.get(4))
        );
        final KV<GATKRead, GATKVariant> readNullVariant = KV.of(reads.get(0), null);

        Iterable<GATKVariant> variant10 = Lists.newArrayList(kvReadVariant.get(1).getValue(), kvReadVariant.get(0).getValue());
        Iterable<GATKVariant> variant2 = Lists.newArrayList(kvReadVariant.get(2).getValue());
        Iterable<GATKVariant> variant3 = Lists.newArrayList(kvReadVariant.get(3).getValue());
        Iterable<GATKVariant> variant4 = Lists.newArrayList(kvReadVariant.get(5).getValue());
        Iterable<GATKVariant> noVariant = Lists.newArrayList(new ArrayList<>());

        kvReadiVariant = Arrays.asList(
                KV.of(kvReadVariant.get(0).getKey(), variant10),
                KV.of(kvReadVariant.get(2).getKey(), variant2),
                KV.of(kvReadVariant.get(3).getKey(), variant3),
                KV.of(kvReadVariant.get(5).getKey(), variant4),
                KV.of(reads.get(0), noVariant)
        );

        kvReadContextData = Arrays.asList(
                KV.of(kvReadsRefBases.get(0).getKey(), new ReadContextData(kvReadsRefBases.get(0).getValue(), Lists.newArrayList())),
                KV.of(kvReadsRefBases.get(1).getKey(), new ReadContextData(kvReadsRefBases.get(1).getValue(), kvReadiVariant.get(0).getValue())),
                KV.of(kvReadsRefBases.get(2).getKey(), new ReadContextData(kvReadsRefBases.get(2).getValue(), kvReadiVariant.get(1).getValue())),
                KV.of(kvReadsRefBases.get(3).getKey(), new ReadContextData(kvReadsRefBases.get(3).getValue(), kvReadiVariant.get(2).getValue())),
                KV.of(kvReadsRefBases.get(4).getKey(), new ReadContextData(kvReadsRefBases.get(4).getValue(), kvReadiVariant.get(3).getValue()))
        );
    }

    /**
     * makeRead creates a read backed by either SAMRecord or Google model Read.
     * @param startLength the key is the start of the read, the value is the length.
     * @param i name, note that if a different i is used, then two otherwise identical reads are not equal.
     * @param clazz either Google model Read or SAMRecord
     * @return a new GAKTRead with either a Google model backed or SAMRecord backed read.
     */
    public static GATKRead makeRead(String contig, KV<Integer, Integer> startLength, int i, Class<?> clazz) {
        return makeRead(contig, startLength.getKey(), startLength.getValue(),i, clazz);
    }

    /**
     * makeRead creates a read backed by either SAMRecord or Google model Read.
     * @param start start position of the read
     * @param length length of the read
     * @param i name, note that if a different i is used, then two otherwise identical reads are not equal.
     * @param clazz either Google model Read or SAMRecord
     * @return a new GAKTRead with either a Google model backed or SAMRecord backed read.
     */
    public static GATKRead makeRead(String contig, int start, int length, int i, Class<?> clazz) {
        if (clazz == SAMRecord.class) {
            return ArtificialReadUtils.createHeaderlessSamBackedRead(Integer.toString(i), contig, start, length);
        } else {
            throw new GATKException("invalid GATKRead type");
        }
    }

    private SimpleInterval makeInterval(String contig, KV<Integer, Integer> startLength) {
        return new SimpleInterval(contig, startLength.getKey(), startLength.getKey() + startLength.getValue() - 1);
    }

    private ReferenceBases getBases(String contig, int start, int end) {
        return FakeReferenceSource.bases(new SimpleInterval(contig, start, end));
    }

    public List<SimpleInterval> getAllIntervals() {
        return allIntervals;
    }

    public List<GATKRead> getReads() {
        return reads;
    }

    public List<KV<GATKRead, ReferenceBases>> getKvReadsRefBases() {
        return kvReadsRefBases;
    }

    public List<GATKVariant> getVariants() {
        return variants;
    }

    public List<KV<GATKRead, ReadContextData>> getKvReadContextData() {
        return kvReadContextData;
    }

    @Test
    public static void verifyDivisibilityWithRefShard() {
        // We want the ratio between the two shard types to be an int so we can use them more easily for testing.
        Assert.assertEquals(Math.floorMod(ReferenceShard.REFERENCE_SHARD_SIZE, VariantShard.VARIANT_SHARDSIZE), 0);
    }

    public List<KV<GATKRead, Iterable<GATKVariant>>> getKvReadiVariant() {
        return kvReadiVariant;
    }
}
