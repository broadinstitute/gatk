package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.KolmogorovSmirnovCalculator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

@CommandLineProgramProperties(summary="Find distribution of fragment lengths.",
        oneLineSummary="Find distribution of fragment lengths.",
        programGroup = StructuralVariationSparkProgramGroup.class)
public final class FragmentLengthDistributionSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {
        final SAMFileHeader header = getHeaderForReads();
        if ( header.getSortOrder() != SAMFileHeader.SortOrder.coordinate ) {
            throw new GATKException("The reads must be coordinate sorted.");
        }
        final JavaRDD<GATKRead> reads = getUnfilteredReads();

        final int nReadGroups = header.getReadGroups().size();
        final List<ReadMetadata.FragmentLengthCounts> perPartitionStatistics =
                reads.mapPartitions(readItr ->
                        Collections.singletonList(new ReadMetadata.FragmentLengthCounts(readItr, nReadGroups)).iterator())
                        .collect();
        final ReadMetadata.FragmentLengthCounts combinedCounts =
                new ReadMetadata.FragmentLengthCounts(perPartitionStatistics, nReadGroups);
        final Map<String,long[]> countMap = combinedCounts.getReadGroupToFragmentSizeCountMap();
        final TreeMap<String,KolmogorovSmirnovCalculator> ksMap = new TreeMap<>();
        for ( final Map.Entry<String,long[]> entry : countMap.entrySet() ) {
            ksMap.put(entry.getKey(), new KolmogorovSmirnovCalculator(entry.getValue()));
        }
        final Broadcast<Map<String,KolmogorovSmirnovCalculator>> broadcastKSMap = ctx.broadcast(ksMap);
        try {
            final List<String> locs = reads.mapPartitions(readItr -> {
                final Map<String,Deque<GATKRead>> deQueMap = new HashMap<>(SVUtils.hashMapCapacity(nReadGroups));
                final List<String> evidence = new ArrayList<>();
                while (readItr.hasNext()) {
                    final GATKRead read = readItr.next();
                    int tLen = Math.abs(read.getFragmentLength());
                    if ( tLen != 0 && read.getMappingQuality() >= 60 &&
                            read.isPaired() &&
                            read.isReverseStrand() &&
                            !read.mateIsReverseStrand() &&
                            !read.isSecondaryAlignment() &&
                            !read.isSupplementaryAlignment() &&
                            !read.isDuplicate() ) {
                    }
                }
                return evidence.iterator();
            }).collect();
        }
        finally {
            broadcastKSMap.destroy();
        }
    }
}
