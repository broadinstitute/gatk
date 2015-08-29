package org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates;

import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import com.google.common.base.Stopwatch;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.dataflow.DoFnWLog;
import org.broadinstitute.hellbender.utils.dataflow.lcollections.LocalCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;

import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

/**
 * The main transform for MarkDuplicates. Takes reads, returns them with some reads marked as duplicates.
 */
public final class MarkDuplicatesFromShardsDataflowTransform extends PTransform<PCollection<KV<String, Iterable<GATKRead>>>, PCollection<GATKRead>> {
    private static final long serialVersionUID = 1l;


    /**
     *  The header for all reads processed.
     */
    private final PCollectionView<SAMFileHeader> header;
    /**
     * Contains methods for finding optical duplicates.
     */
    private final PCollectionView<OpticalDuplicateFinder> finder;

    public MarkDuplicatesFromShardsDataflowTransform(final PCollectionView<SAMFileHeader> header, final PCollectionView<OpticalDuplicateFinder> finder) {
        this.header = header;
        this.finder = finder;
    }

    @Override
    public PCollection<GATKRead> apply(final PCollection<KV<String, Iterable<GATKRead>>> preads) {
        return preads.apply(ParDo.named("MarkDuplicatesFromShardsDataflow")
            .withSideInputs(header, finder)
            .of(new DoFnWLog<KV<String, Iterable<GATKRead>>, GATKRead>("MarkDuplicatesFromShardsDataflow") {
                private static final long serialVersionUID = 1l;
                private transient Stopwatch processing;
                private transient Stopwatch shipping;

                @Override
                public void startBundle(Context c) throws Exception {
                    super.startBundle(c);
                    processing = Stopwatch.createUnstarted();
                    shipping = Stopwatch.createUnstarted();
                }

                @Override
                public void finishBundle(Context c) throws Exception {
                    super.finishBundle(c);
                    log.info("MarkDuplicatesFromShardsDataflow bundle worked: "+processing.elapsed(TimeUnit.MILLISECONDS)+" ms. Shipped: "+shipping.elapsed(TimeUnit.MILLISECONDS)+" ms.");
                }

                // Convert from PCollection to LocalCollection, run the computation there, convert back to PCollection.
                @Override
                public void processElement(ProcessContext c) throws Exception {
                    processing.start();
                    final SAMFileHeader h = c.sideInput(header);
                    final OpticalDuplicateFinder f = c.sideInput(finder);
                    LocalCollection<GATKRead> answer = localApply(h, f, c.element().getValue());
                    processing.stop();
                    shipping.start();
                    for (GATKRead r : answer.iterable()) {
                        c.output(r);
                    }
                    shipping.stop();
                }
            }));
    }

    public static LocalCollection<GATKRead> localApply(final SAMFileHeader header, final OpticalDuplicateFinder finder, Iterable<GATKRead> readsIterable) {

        // read input, make our own copy.
        LocalCollection<GATKRead> reads = LocalCollection.of(readsIterable)
            .map(MarkDuplicatesUtils::markNonDuplicate);

        if (reads.size()==0) return reads;

        // we specify a minimum size to guarantee we always get 2 elements in the array, even if e.g. all the reads are primary.
        ArrayList<LocalCollection<GATKRead>> readsPartitioned = reads.partitionBy(MarkDuplicates::primaryNonPrimaryAlignment, 2);
        LocalCollection<GATKRead> fragments = readsPartitioned.get(MarkDuplicates.ReadsPartition.PRIMARY.ordinal());
        LocalCollection<GATKRead> fragmentsTransformed = transformFragments(header, fragments);

        LocalCollection<GATKRead> pairs = readsPartitioned.get(MarkDuplicates.ReadsPartition.PRIMARY.ordinal());
        LocalCollection<GATKRead> pairsTransformed = transformReads(header, finder, pairs);

        LocalCollection<GATKRead> not_primary = readsPartitioned.get(MarkDuplicates.ReadsPartition.NOT_PRIMARY.ordinal());

        LocalCollection<GATKRead> ret = LocalCollection.union(fragmentsTransformed, pairsTransformed, not_primary)
            .map(MarkDuplicatesFromShardsDataflowTransform::removeHeader);

        return ret;
    }


    private static LocalCollection<GATKRead> transformFragments(final SAMFileHeader header, LocalCollection<GATKRead> fragments) {
        return fragments
            // group by key (genomic location etc)
            .groupBy((read) -> MarkDuplicatesUtils.keyForFragment(header, read))
            // in each group, mark the duplicates.
            // false = do not make a copy. That's safe here because we know exactly how LocalCollection works.
            .mapTransform((reads)->MarkDuplicatesUtils.transformFragmentsFn(reads,false))
            // go back to a single list
            .flatten();
    }

    /**
     * (1) keyReadsByName: label each read with its read group and read name.
     * (2) GroupByKey: group together reads with the same group and name.
     * (3) keyPairedEndsWithAlignmentInfo:
     *   (a) Sort each group of reads (see GATKOrder below).
     *   (b) Pair consecutive reads into PairedEnds. In most cases there will only be two reads
     *       with the same name. TODO: explain why there might be more.
     *   (c) Label each read with alignment information: Library, reference index,
     *       stranded unclipped start and reverse strand.
     *   (d) Leftover reads are emitted, unmodified, as an unpaired end.
     * (4) GroupByKey: Group PairedEnds that share alignment information. These pairs
     *     are duplicates of each other.
     * (5) markDuplicatePairs:
     *   (a) For each group created by (4), sort the pairs by score and mark all but the
     *       highest scoring as duplicates.
     *   (b) Determine which duplicates are optical duplicates and increase the overall count.
     */
    static LocalCollection<GATKRead> transformReads(final SAMFileHeader header,
                                                    final OpticalDuplicateFinder finder, final LocalCollection<GATKRead> pairs) {

        return pairs
            // only keep the reads that have a mapped mate
            .filter(ReadUtils::readHasMappedMate)
            // group reads by genomic location etc.
            .groupBy((read) -> MarkDuplicatesUtils.keyForRead(header, read))
            // use the groups to build paired reads
            .mapTransform((reads) -> MarkDuplicatesUtils.pairTheReads(header, reads))
            // merge the groups together
            .flatten()
            // now group by pair key
            .groupBy((pair) -> MarkDuplicatesUtils.keyForPair(header, pair))
            // duplicates would have the same pair key, find them.
            // For each group, pass the list of pairs to markDuplicatePairedEndsFn and get a list of reads in return.
            // false = no copy. This is safe because no one else will see those pairs.
            .mapTransform((pairGroup) -> MarkDuplicatesUtils.markDuplicatePairedEndsFn(finder, pairGroup, false))
            // remove the groups, now we just have a list of reads.
            .flatten();
    }


    public static GATKRead removeHeader(GATKRead r) {
        if (r instanceof SAMRecordToGATKReadAdapter) {
            SAMRecordToGATKReadAdapter record = (SAMRecordToGATKReadAdapter)r;
            record.setHeader(null);
            return record;
        }
        return r;
    }

}