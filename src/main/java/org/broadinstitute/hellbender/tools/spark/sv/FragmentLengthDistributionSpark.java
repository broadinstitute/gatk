package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.KolmogorovSmirnovCalculator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

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
        final List<ReadMetadata.PartitionStatistics> perPartitionStatistics =
                reads.mapPartitions(readItr ->
                        Collections.singletonList(new ReadMetadata.PartitionStatistics(readItr, nReadGroups)).iterator())
                        .collect();
        final Map<String,long[]> countMap =
                new ReadMetadata.PartitionStatistics(perPartitionStatistics, nReadGroups).getReadGroupToFragmentSizeCountMap();
        final TreeMap<String,KolmogorovSmirnovCalculator> ksMap = new TreeMap<>();
        for ( final Map.Entry<String,long[]> entry : countMap.entrySet() ) {
            ksMap.put(entry.getKey(), new KolmogorovSmirnovCalculator(entry.getValue(),.05f));
        }
    }
}
