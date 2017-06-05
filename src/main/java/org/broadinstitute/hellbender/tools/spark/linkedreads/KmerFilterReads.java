package org.broadinstitute.hellbender.tools.spark.linkedreads;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

@CommandLineProgramProperties(summary = "Print reads from the input BAM", oneLineSummary = "PrintReads on Spark", programGroup = SparkProgramGroup.class)
public class KmerFilterReads extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "upstream kmers to limit reads by",
            fullName = "kmers",
            optional = false)
    public Set<String> kmers = new HashSet<>();


    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    public String output;

    @Override
    public boolean requiresReads() {
        return true;
    }


    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final UpstreamKmerReadFilter upstreamKmerReadFilter = new UpstreamKmerReadFilter(kmers, getReference());
        final JavaRDD<GATKRead> reads = getReads().filter(upstreamKmerReadFilter::test);
        writeReads(ctx, output, reads);
    }

    @Override
    public boolean requiresReference() {
        return true;
    }


}
