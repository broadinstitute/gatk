package org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates;

import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionList;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;

/**
 * The main transform for MarkDuplicates. Takes reads, returns reads with some reads marked as duplicates.
 */
public class MarkDuplicates extends PTransform<PCollection<GATKRead>, PCollection<GATKRead>> {
    private static final long serialVersionUID = 1l;

    //private vs non-primary alignments
    private enum ReadsPartition {
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

        final PCollection<GATKRead> pairs = readsPartitioned.get(ReadsPartition.PRIMARY.ordinal());
        final PCollection<GATKRead> pairsTransformed = MarkDuplicatesUtils.transformReads(header, opticalDuplicateFinder, pairs);

        //no work on those
        final PCollection<GATKRead> not_primary = readsPartitioned.get(ReadsPartition.NOT_PRIMARY.ordinal());

        return PCollectionList.of(fragmentsTransformed).and(pairsTransformed).and(not_primary).apply(Flatten.<GATKRead>pCollections());
    }

    private static PCollectionList<GATKRead> partitionReadsByPrimaryNonPrimaryAlignment(PCollection<GATKRead> preads) {
        return preads.apply(Partition.of(2, new Partition.PartitionFn<GATKRead>() {
            private static final long serialVersionUID = 1l;

            @Override
            public int partitionFor(final GATKRead read, final int n) {
                return read.isSecondaryAlignment() || read.isSupplementaryAlignment() || read.isUnmapped() ? ReadsPartition.NOT_PRIMARY.ordinal() : ReadsPartition.PRIMARY.ordinal();
            }
        }));
    }


}
