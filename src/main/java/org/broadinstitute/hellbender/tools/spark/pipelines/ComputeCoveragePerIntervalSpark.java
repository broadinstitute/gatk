package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ExomeReadCounts;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollections;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "Counts reads in the input BAM", oneLineSummary = "Counts reads in a BAM file", programGroup = SparkProgramGroup.class)
public final class ComputeCoveragePerIntervalSpark extends SparkCommandLineProgram {

    private static final long serialVersionUID = 1l;

    @ArgumentCollection
    public ReadInputArgumentCollection readArguments= new RequiredReadInputArgumentCollection();;

    @Argument(
            doc = "File containing the targets for analysis",
            shortName = ExomeReadCounts.EXOME_FILE_SHORT_NAME,
            fullName = ExomeReadCounts.EXOME_FILE_FULL_NAME,
            optional = false
    )
    public File targetIntervalFile;

    @ArgumentCollection
    public IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String out;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        if ( readArguments.getReadFilesNames().size() != 1 ) {
            throw new UserException("This tool only accepts a single bam/sam/cram as input");
        }
        final String bam = readArguments.getReadFilesNames().get(0);
        final ReadsSparkSource readSource = new ReadsSparkSource(ctx);

        final SAMFileHeader readsHeader = ReadsSparkSource.getHeader(ctx, bam);

        /**
         * If no intervals are provided, we want all reads, mapped and unmapped.
         */
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(readsHeader.getSequenceDictionary())
                : null;

        final TargetCollection<Locatable> result = getTargetCollection();

        final IntervalsSkipList<Locatable> isl = new IntervalsSkipList<>(result.targets());
        final Broadcast<IntervalsSkipList<Locatable>> islBroad = ctx.broadcast(isl);

        //NOTE: we hit some serialization error with lambdas in WellformedReadFilter
        //so for now we'll write some filters out explicitly.
//        final ReadFilter readFilter = ctx.broadcast(new WellformedReadFilter(readsHeader));
        final JavaRDD<GATKRead> rawReads = readSource.getParallelReads(bam, intervals);
        final JavaRDD<GATKRead> reads = rawReads.filter(read -> !read.isUnmapped() && read.getStart() <= read.getEnd() && !read.isDuplicate());
        final Map<Locatable, Long> byKey = reads.flatMap(read -> islBroad.getValue().getOverlapping(new SimpleInterval(read))).countByValue();

        final SortedMap<Locatable, Object> byKeySorted = new TreeMap<>(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
        byKeySorted.putAll(byKey);

        print(byKeySorted, System.out);

        if (out != null){
            final File file = new File(out);
            try(final OutputStream outputStream = BucketUtils.createFile(file.getPath(), null);
                final PrintStream ps = new PrintStream(outputStream)) {
                print(byKeySorted, ps);
            } catch(final IOException e){
                throw new UserException.CouldNotCreateOutputFile(file, e);
            }
        }
    }

    private void print(SortedMap<Locatable, Object> byKeySorted, PrintStream ps) {
        ps.println("CONTIG" + "\t" + "START" + "\t" + "END" + "\t" + "<ALL>");
        for (final Locatable loc : byKeySorted.keySet()){
            ps.println(loc.getContig() + "\t" + loc.getStart() + "\t" + loc.getEnd() + "\t" + byKeySorted.get(loc));
        }
    }

    private TargetCollection<Locatable> getTargetCollection() {
        //Target intervals
        Utils.regularReadableUserFile(targetIntervalFile);
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(targetIntervalFile);
        final Class<? extends Feature> featureType = codec.getFeatureType();
        final TargetCollection<BEDFeature> result;
        if (BEDFeature.class.isAssignableFrom(featureType)) {
            @SuppressWarnings("unchecked")
            final FeatureCodec<BEDFeature, ?> bedFeatureCodec = (FeatureCodec<BEDFeature, ?>) codec;
            result = TargetCollections.fromBEDFeatureFile(targetIntervalFile, bedFeatureCodec);
        } else {
            throw new UserException.BadInput(String.format("currently only BED formatted exome file are supported. '%s' does not seem to be a BED file",targetIntervalFile.getAbsolutePath()));
        }
        //NOTE: java does not allow conversion of List<BEDFeature> to List<Locatable>.
        //So we do this trick
        final List<Locatable> targets = result.targets().stream().collect(Collectors.toList());
        return new HashedListTargetCollection<>(targets);
    }
}
