package org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates;

import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionList;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;

/**
 * The main transform for MarkDuplicates. Takes reads, returns reads with some reads marked as duplicates.
 */
public class MarkDuplicates extends PTransform<PCollection<GATKRead>, PCollection<GATKRead>> {
    private static final long serialVersionUID = 1l;

    //private vs non-primary alignments
    public enum ReadsPartition {
        PRIMARY, NOT_PRIMARY
    }

    /**
     * The header for all reads processed.
     */
    private final PCollectionView<SAMFileHeader> header;

    /**
     * Optical duplicate finder finds the physical location of the reads and determines if they are close enough to
     * be considered optical duplicates.
     */
    private final PCollectionView<OpticalDuplicateFinder> opticalDuplicateFinder;

    public MarkDuplicates(final PCollectionView<SAMFileHeader> header,
        final PCollectionView<OpticalDuplicateFinder> opticalDuplicateFinder) {
        this.header = header;
        this.opticalDuplicateFinder = opticalDuplicateFinder;
    }

    @Override
    public PCollection<GATKRead> apply(final PCollection<GATKRead> preads) {
        final PCollectionList<GATKRead> readsPartitioned = partitionReadsByPrimaryNonPrimaryAlignment(preads);

        final PCollection<GATKRead> fragments = readsPartitioned.get(ReadsPartition.PRIMARY.ordinal());
        final PCollection<GATKRead> fragmentsTransformed = MarkDuplicatesUtils.transformFragments(header, fragments);

        fragmentsTransformed.apply(Count.globally())
            .apply(DataflowUtils.convertToString())
            .apply(TextIO.Write.to("tmp-fragmentsTransformed"));

        final PCollection<GATKRead> pairs = readsPartitioned.get(ReadsPartition.PRIMARY.ordinal());
        final PCollection<GATKRead> pairsTransformed = MarkDuplicatesUtils.transformReads(header, opticalDuplicateFinder, pairs);

        pairsTransformed.apply(Count.globally())
            .apply(DataflowUtils.convertToString())
            .apply(TextIO.Write.to("tmp-pairsTransformed"));

        //no work on those
        final PCollection<GATKRead> not_primary = readsPartitioned.get(ReadsPartition.NOT_PRIMARY.ordinal());

        not_primary.apply(Count.globally())
            .apply(DataflowUtils.convertToString())
            .apply(TextIO.Write.to("tmp-not_primary"));

        return PCollectionList.of(fragmentsTransformed).and(pairsTransformed).and(not_primary).apply(Flatten.<GATKRead>pCollections());
    }

    public static int primaryNonPrimaryAlignment(GATKRead read) {
        return read.isSecondaryAlignment() || read.isSupplementaryAlignment() || read.isUnmapped() ? ReadsPartition.NOT_PRIMARY.ordinal() : ReadsPartition.PRIMARY.ordinal();
    }

    private static PCollectionList<GATKRead> partitionReadsByPrimaryNonPrimaryAlignment(PCollection<GATKRead> preads) {
        return preads.apply(Partition.of(2, new Partition.PartitionFn<GATKRead>() {
            private static final long serialVersionUID = 1l;

            @Override
            public int partitionFor(final GATKRead read, final int n) {
                return primaryNonPrimaryAlignment(read);
            }
        }));
    }


}
