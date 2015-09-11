package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;

public class MarkDuplicatesSpark  {

    public static JavaRDD<GATKRead> mark(final JavaRDD<GATKRead> reads, final SAMFileHeader header,
                                         final OpticalDuplicateFinder opticalDuplicateFinder) {

        JavaRDD<GATKRead> fragments = reads.filter(v1 -> !isNonPrimary(v1));
        JavaRDD<GATKRead> nonPrimaryReads = reads.filter(v1 -> isNonPrimary(v1));
        JavaRDD<GATKRead> pairsTransformed = MarkDuplicatesSparkUtils.transformReads(header, opticalDuplicateFinder, fragments);

        JavaRDD<GATKRead> fragmentsTransformed = MarkDuplicatesSparkUtils.transformFragments(header, fragments);
        return fragmentsTransformed.union(pairsTransformed).union(nonPrimaryReads);
    }

    private static boolean isNonPrimary(GATKRead read) {
        return read.isSecondaryAlignment() || read.isSupplementaryAlignment() || read.isUnmapped();
    }

}