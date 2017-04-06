package org.broadinstitute.hellbender.tools.spark.sv;

import avro.shaded.com.google.common.collect.Iterables;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.avro.test.Simple;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.storage.StorageLevel;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Tool to compose an alignment file with only the read evidence that was used to call a set of variants.
 *
 * @author &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(summary="Composes an aligment files with the reads that compose the evidence to call a set of variants. ",
        oneLineSummary="Extract reads that are evidence for called variants",
        programGroup = StructuralVariationSparkProgramGroup.class)
public final class EvidenceSamForVariantsSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    public static final String ASSEMBLY_DIRECTORY_SHORT_NAME = "adir";
    public static final String ASSEMBLY_DIRECTORY_FULL_NAME = "assemblyFastqDirectory";
    public static final String ASSEMBLY_FILE_NAME_FORMAT_SHORT_NAME = "aname";
    public static final String ASSEMBLY_FILE_NAME_FORMAT_FULL_NAME = "assemblyFileNameFormat";
    public static final String SHARD_SIZE_SHORT_NAME = "ssize";
    public static final String SHARD_SIZE_FULL_NAME = "shardSize";
    public static final String KEEP_INTERVALS_SHORT_NAME = "keepL";
    public static final String KEEP_INTERVALS_FULL_NAME = "keepInterval";
    public static final String READ_INTERVALS_SHORT_NAME = "rio";
    public static final String READ_INTERVALS_FULL_NAME = "readIntervalsOutput";

    private static final Pattern FASTQ_MAPPING_PATTERN = Pattern.compile("^mapping=(\\S+):(\\d+);(\\S+)$");

    @Argument(doc = "the input variant calls file name", shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME)
    private String variantSourceName = null;

    @Argument(doc = "the output bam file name", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private File output = null;

    @Argument(doc = "directory where the assembly fastq files are located",
              shortName = ASSEMBLY_DIRECTORY_SHORT_NAME,
              fullName = ASSEMBLY_DIRECTORY_FULL_NAME)
    private File fastqDir = null;

    @Argument(doc = "format of the assembly fastq file name given the assembly id as the only argument. This string follows the printf scape format (e.g. 'assembly_%d.fastq'",
              shortName = ASSEMBLY_FILE_NAME_FORMAT_SHORT_NAME,
              fullName = ASSEMBLY_FILE_NAME_FORMAT_FULL_NAME)
    private String assemblyFileNameFormat = "assembly%d.fastq";

    @Argument(doc = "maximum distance between locations to be considered part of the same shard.",
              shortName = SHARD_SIZE_SHORT_NAME,
              fullName = SHARD_SIZE_FULL_NAME,
              optional = true)
    private int shardWidth = 100;

    @Argument(doc = "intervals to keep aligned. If specified, reads that fall out of this intervals (first aligned position) will be declared as unmapped",
              shortName = KEEP_INTERVALS_SHORT_NAME,
              fullName = KEEP_INTERVALS_FULL_NAME,
              optional = true)
    private List<String> keepAlignmenedIntervals = new ArrayList<>();

    @Argument(doc = "output file to dump the intervals that contain the reads of interest",
              shortName = READ_INTERVALS_SHORT_NAME,
              fullName = READ_INTERVALS_FULL_NAME,
              optional = true)
    private File readIntervalsOutFile = null;

    private transient SAMFileGATKReadWriter outputWriter;
    private transient VariantsSparkSource variantsSource;

    @Override
    public boolean requiresIntervals() { return true; }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        Utils.nonNull(ctx);
        try {
            setup(ctx);
            final JavaRDD<VariantContext> variants = variantsSource.getParallelVariantContexts(variantSourceName,
                    intervalArgumentCollection.getIntervals(this.getReferenceSequenceDictionary()));
            processVariants(variants,ctx);
        } finally {
            tearDown();
        }
    }

    public void setup(JavaSparkContext ctx) {
        if (shardWidth < 1) {
            throw new UserException.BadInput("the shard size must be 1 or greather: " + shardWidth);
        }
        try {
            final SAMFileHeader readsHeader = getHeaderForReads();
            // use the reference's dictionary rather than the input reads in case the
            // first is a subset of the second.
            readsHeader.setSequenceDictionary(this.getReferenceSequenceDictionary());
            outputWriter = outputAlignmentArguments.createSAMWriter(output, true, readsHeader);
        } catch (final Throwable th) {
            throw new UserException.CouldNotCreateOutputFile(output, new Exception(th.getMessage(), th));
        }
        variantsSource = new VariantsSparkSource(ctx);
    }

    public void tearDown() {
        if (outputWriter != null) {
            outputWriter.close();
            outputWriter = null;
        }
        variantsSource = null;
    }

    private void processVariants(final JavaRDD<VariantContext> rdd, final JavaSparkContext ctx) {
        final List<SimpleInterval> toKeepList = new IntervalArgumentCollection() {

            private static final long serialVersionUID = 1L;

            @Override
            protected List<String> getIntervalStrings() {
                return EvidenceSamForVariantsSpark.this.keepAlignmenedIntervals;
            }

            @Override
            protected void addToIntervalStrings(final String newInterval) {
                 throw new UnsupportedOperationException();
            }
        }.getIntervals(getReferenceSequenceDictionary());

        final JavaRDD<Integer> assemblyIds = rdd.map(this::variantAssemblyIds).flatMap(List::iterator).distinct();
        final JavaRDD<File> assemblyFiles = assemblyIds.map(this::assemblyFileForId);
        final JavaRDD<ReadPairInfo> fastqReadPairs = assemblyFiles.flatMap(this::fastqReadRecordsForFile).cache();
        final JavaRDD<ReadPairLocationInfo> inputReadLocations = fastqReadPairs.flatMap(rp -> rp.locations().iterator()).cache();
        final JavaRDD<ReadPairLocationInfo> mappedReadLocations = inputReadLocations.filter(rp -> rp.location != null);
        final TreeSet<SimpleInterval> toKeepSet = new TreeSet<>(LocatableComparator.LEXICOGRAPHIC);
        toKeepSet.addAll(toKeepList);
        processMappedReadLocations(mappedReadLocations, ctx, toKeepSet);


        final JavaRDD<ReadPairInfo> unmappedReadPairs = fastqReadPairs.filter(ReadPairInfo::isUnmapped).cache();
        if (!unmappedReadPairs.isEmpty()) {
            logger.info("There are completely unmapped read pairs and so we need to process unmapped reads in the input bam");
            processUnmappedReadPairs(unmappedReadPairs);
        }
    }

    private void processMappedReadLocations(final JavaRDD<ReadPairLocationInfo> mappedReadLocations, final JavaSparkContext ctx, final TreeSet<SimpleInterval> toKeepSet) {
        final JavaPairRDD<String, ReadPairLocationInfo> contigAndReadPairInfo =
                mappedReadLocations.mapToPair(rpl -> new Tuple2<>(rpl.location.getContig(), rpl));
        final JavaPairRDD<String, SimpleInterval> contigAndInterval =
                contigAndReadPairInfo.mapToPair(rpl -> new Tuple2<>(rpl._1(), rpl._2().location));
        final JavaPairRDD<String, TreeSet<SimpleInterval>> contigAndIntervalSet =
                contigAndInterval.combineByKey(EvidenceSamForVariantsSpark::createCombiner,
                        (c, v) -> mergeValue(c, v, shardWidth),
                        (a, b) -> mergeCombiners(a, b, shardWidth));
        if (readIntervalsOutFile != null) {
            try (final PrintWriter writer = new PrintWriter(new FileWriter(readIntervalsOutFile))) {
                contigAndIntervalSet.sortByKey().toLocalIterator().forEachRemaining(pair -> {
                    for (final SimpleInterval interval : pair._2()) {
                        writer.println(interval);
                    }
                });
            } catch (final Exception ex) {
                throw new UserException.CouldNotCreateOutputFile(readIntervalsOutFile, ex.getMessage(), ex);
            }
        }
        final JavaPairRDD<String, TreeSet<SimpleInterval>> contigAndIntervalShardSet =
                contigAndIntervalSet.mapToPair(rpl -> new Tuple2<>(rpl._1(), IntervalUtils.cutToShards(rpl._2(), shardWidth, () -> new TreeSet<>(LocatableComparator.LEXICOGRAPHIC))));


        final Broadcast<Map<String, TreeSet<SimpleInterval>>> intervalShardsByContig =
                ctx.broadcast(contigAndIntervalShardSet.collectAsMap());
        final List<SimpleInterval> allIntervals = intervalShardsByContig.getValue().entrySet().stream()
                .flatMap(entry -> entry.getValue().stream())
                .collect(Collectors.toList());

        final JavaPairRDD<SimpleInterval, ReadPairLocationInfo> shardAndReadPairInfo =
                contigAndReadPairInfo.flatMapToPair(rpl -> flatMapReadPairLocationInfoToShard(rpl._2(), intervalShardsByContig.getValue().get(rpl._1())));
        final JavaPairRDD<SimpleInterval, String> shardAndReadName = shardAndReadPairInfo.mapValues(r -> {
            return r.name;
        });
        final JavaPairRDD<SimpleInterval, HashSet<String>> shardAndReadNameSet = shardAndReadName.groupByKey()
                .mapValues(names -> CollectionUtils.collect(names, HashSet::new)).cache();
        final GenomeLocParser locParser = new GenomeLocParser(getReferenceSequenceDictionary());
        final Broadcast<GenomeLocParser> dictBroadcast = ctx.broadcast(locParser);
        final Broadcast<TreeSet<SimpleInterval>> toKeepBroadcast = ctx.broadcast(toKeepSet);
        final Broadcast<Map<SimpleInterval, HashSet<String>>> shardAndReadNameSetBroadCast = ctx.broadcast(shardAndReadNameSet.collectAsMap());
        final JavaPairRDD<SimpleInterval, GATKRead> shardAndReads = getUnfilteredReads(contigAndIntervalSet.flatMap(pair -> pair._2().iterator()).collect())
                .mapToPair(read -> new Tuple2<>(read, readLocation(read)))
                .mapToPair(pair -> new Tuple2<>(findShard(intervalShardsByContig.getValue().get(pair._2().getContig()), pair._2()), pair._1()))
                .filter(pair -> pair._1() != null);


        final JavaPairRDD<SimpleInterval, GATKRead> shardsNamedReads = shardAndReads
                .filter(pair -> shardAndReadNameSetBroadCast.getValue().get(pair._1()).contains(pair._2().getName()))
                .mapToPair(pair -> new Tuple2<>(pair._1(), updateBasedOnKeepSet(pair._2(), toKeepBroadcast)));

        shardsNamedReads.values()
                        .filter(r -> !r.isUnmapped() || !r.mateIsUnmapped())
                        .sortBy(r -> new ReadSortingAttributes(dictBroadcast.getValue(), r), true, 100)
                        .toLocalIterator().forEachRemaining(outputWriter::addRead);

        shardsNamedReads.values()
                .filter(r -> r.isUnmapped() && r.mateIsUnmapped())
                .sortBy(r -> r.getName() + "/" + (r.isFirstOfPair() ? "1" : "2"), true, 100)
                .toLocalIterator().forEachRemaining(outputWriter::addRead);
    }

    private GATKRead updateBasedOnKeepSet(final GATKRead read, Broadcast<TreeSet<SimpleInterval>> toKeepBroadcast) {
        final TreeSet<SimpleInterval> toKeepSet = toKeepBroadcast.getValue();
        final SimpleInterval readInterval = read.isUnmapped() ? null : new SimpleInterval(read.getContig(), read.getStart(), read.getStart());
        final SimpleInterval mateInterval = read.mateIsUnmapped() ? null : new SimpleInterval(read.getMateContig(), read.getMateStart(), read.getMateStart());
        final SimpleInterval newReadInterval = readInterval == null ? null : findShard(toKeepSet, readInterval) == null ? null : readInterval;
        final SimpleInterval newMateInterval = mateInterval == null ? null : findShard(toKeepSet, mateInterval) == null ? null : mateInterval;

        if (Objects.equals(readInterval, newReadInterval) && Objects.equals(mateInterval, newMateInterval)) {
            return read;
        } else {
            final GATKRead result = read.copy();
            if (!Objects.equals(mateInterval, newMateInterval )) {
                result.unmapMate();
                if (readInterval == null) {
                    result.resetPosition();
                }
            }
            if (!Objects.equals(readInterval, newReadInterval )) {
                result.unmap();
                if (newMateInterval == null) {
                    result.unmapMate();
                }
            }

            return result;
        }
    }

    private SimpleInterval findShard(final TreeSet<SimpleInterval> target, final SimpleInterval query) {
        if (target == null) {
            return null;
        } else if (query == null) {
            return null;
        } else {
            final NavigableSet<SimpleInterval> tailSet = target.tailSet(query, true);
            final NavigableSet<SimpleInterval> searchSet = tailSet.isEmpty() ? target.descendingSet() :
                    target.headSet(tailSet.first(), true).descendingSet();
            if (searchSet.isEmpty()) {
                return null;
            } else {
                final Comparator<? super SimpleInterval> comparator = Objects.requireNonNull(target.comparator());
                for (final SimpleInterval t : searchSet) {
                    if (t.overlaps(query)) {
                        return t;
                    } else if (comparator.compare(t,query) < 0) {
                        break;
                    }
                }
                return null;
            }
        }
    }

    private Iterator<Tuple2<SimpleInterval, ReadPairLocationInfo>> flatMapReadPairLocationInfoToShard(final ReadPairLocationInfo readPairLocationInfo, final TreeSet<SimpleInterval> shards) {
        final NavigableSet<SimpleInterval> headSet = shards.headSet(readPairLocationInfo.location, true);
        final NavigableSet<SimpleInterval> tailSet;
        final List<Tuple2<SimpleInterval, ReadPairLocationInfo>> result = new ArrayList<>(4);
        if (headSet.isEmpty()) {
            tailSet = shards.tailSet(readPairLocationInfo.location, true);
        } else {
            final SimpleInterval previous = headSet.last();
            tailSet = (previous.getEnd() >= readPairLocationInfo.location.getStart())
                    ? shards.tailSet(previous, true)
                    : shards.tailSet(readPairLocationInfo.location, true);
        }
        for (final SimpleInterval shard : tailSet) {
            if (shard.getStart() > readPairLocationInfo.location.getEnd()) {
                break;
            } else {
                result.add(new Tuple2<>(shard, readPairLocationInfo));
            }
        }
        return result.iterator();
    }

    private static TreeSet<SimpleInterval> createCombiner(final SimpleInterval interval) {
        final TreeSet<SimpleInterval> result = new TreeSet<>(LocatableComparator.LEXICOGRAPHIC);
        result.add(interval);
        return result;
    }

    private static TreeSet<SimpleInterval> mergeValue(final TreeSet<SimpleInterval> accu, final SimpleInterval interval, final int shardSize) {
        SimpleInterval adding = interval;
        final NavigableSet<SimpleInterval> headSet = accu.headSet(interval, true);
        final NavigableSet<SimpleInterval> updateSet;
        if (headSet.isEmpty())
            updateSet = accu;
        else {
            final Iterator<SimpleInterval>  headIt = headSet.descendingIterator();
            SimpleInterval last = interval;
            while (headIt.hasNext()) {
                final SimpleInterval next = headIt.next();
                if (next.getEnd() < interval.getStart() - shardSize) {
                    break;
                } else {
                    last = next;
                }
            }
            updateSet = accu.tailSet(last, true);
        }
        final Iterator<SimpleInterval> updateIt = updateSet.iterator();
        while (updateIt.hasNext()) {
            final SimpleInterval next = updateIt.next();
            if (next.getStart() < adding.getEnd() + shardSize) {
                updateIt.remove();
                adding = new SimpleInterval(adding.getContig(),
                        Math.min(adding.getStart(), next.getStart()),
                        Math.max(adding.getEnd(),   next.getEnd()));
            } else {
                break;
            }
        }
        updateSet.add(adding);
        return accu;
    }

    private static TreeSet<SimpleInterval> mergeCombiners(final TreeSet<SimpleInterval> accu, final TreeSet<SimpleInterval> addOn, final int shardSize) {
        if (addOn.isEmpty())
            return accu;
        final ArrayList<SimpleInterval> intervals = CollectionUtils.heapSort(Iterables.concat(accu, addOn), IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR::compare, ArrayList::new);
        final ArrayList<SimpleInterval> mergedIntervals = IntervalUtils.mergeLocatables(intervals, (a, b) -> IntervalUtils.distance(a, b) <= shardSize, SimpleInterval::new, ArrayList::new);
        accu.clear();
        accu.addAll(mergedIntervals);
        return accu;
    }

    private void processUnmappedReadPairs(final JavaRDD<ReadPairInfo> unmappedReadPairs) {

        final Set<String> names = unmappedReadPairs.map(ReadPairInfo::name).collect().stream().collect(Collectors.toSet());
        final ReadsDataSource dataSource = new ReadsDataSource(readArguments.getReadPaths());
        dataSource.setTraversalBounds(Collections.emptyList(), true);

        final Stream<GATKRead> unmapped = StreamSupport.stream(Spliterators.spliteratorUnknownSize(dataSource.queryUnmapped(), Spliterator.DISTINCT | Spliterator.IMMUTABLE), false);
        unmapped.filter(r -> names.contains(r.getName()))
                .forEach(outputWriter::addRead);
    }

    @Override
    public void onShutdown() {
        if (outputWriter == null)
            return;
        try {
            outputWriter.close();
        } catch (final Throwable th) {
            throw new UserException.CouldNotCreateOutputFile(output, new Exception(th.getMessage(), th));
        }
    }

    private List<Integer> variantAssemblyIds(final VariantContext vc) {
        if (!vc.hasAttribute(GATKSVVCFHeaderLines.ASSEMBLY_IDS)) {
            return Collections.emptyList();
        } else {
            return vc.getAttributeAsIntList(GATKSVVCFHeaderLines.ASSEMBLY_IDS, -1);
        }
    }

    private File assemblyFileForId(final int id) {
        final File candidate = new File(fastqDir, String.format(assemblyFileNameFormat, id));
        if (!candidate.isFile())
            throw new UserException.CouldNotReadInputFile(candidate, "cannot find the assembly file for id " + id);
        else
            return candidate;
    }

    private Iterator<ReadPairInfo> fastqReadRecordsForFile(final File fastq) {
        final List<ReadPairInfo> result = new ArrayList<>(1000);
        try (final FastqReader reader = new FastqReader(fastq)) {
            while (reader.hasNext()) {
            final FastqRecord read1 = reader.next();
            if (!reader.hasNext())
                throw new IllegalArgumentException("premature end of file");

            final FastqRecord read2 = reader.next();

            final String[] read1Parts = read1.getReadHeader().split("\\s+");
            final String[] read2Parts = read2.getReadHeader().split("\\s+");
            if (!read1Parts[0].endsWith("/1"))
                throw new IllegalArgumentException("first read in pair name does not finish with /1:" + read1Parts[0]);
            if (!read2Parts[0].endsWith("/2"))
                throw new IllegalArgumentException("second read in pair name does not finish with /2:" + read2Parts[0]);
            final String name = read1Parts[0].substring(0, read1Parts[0].length() - 2);
            if (!read2Parts[0].substring(0, read2Parts[0].length() - 2).equals(name))
                throw new IllegalArgumentException("first and second read in pair names do not match: " + read1Parts[0] + " and " + read2Parts[0]);
            final List<SimpleInterval> intervals = Stream.concat(Stream.of(read1Parts), Stream.of(read2Parts))
                    .filter(s -> s.startsWith("mapping="))
                    .filter(s -> !s.contains("unmapped"))
                    .filter(s -> s.matches(FASTQ_MAPPING_PATTERN.pattern()))
                    .map(s -> {
                        final Matcher matcher = FASTQ_MAPPING_PATTERN.matcher(s);
                        if (!matcher.matches()) {
                            throw new IllegalArgumentException("unexpected lack of match");
                        }
                        final String contig = matcher.group(1);
                        final int start = Integer.parseInt(matcher.group(2));
                        return new SimpleInterval(contig, start, start);
                    })
                    .collect(Collectors.toList());
            result.add(new ReadPairInfo(name, intervals));
            }
            return result.iterator();
        } catch (final Throwable th) {
            throw new UserException.CouldNotReadInputFile(fastq, new RuntimeException(th.getMessage(), th));
        }
    }

    private static class ReadPairLocationInfo implements Serializable {

        private static final long serialVersionUID = 1L;

        private final String name;
        private final SimpleInterval location;
        private final int occurrences;

        public ReadPairLocationInfo(final String name, final SimpleInterval location, final int occurrences) {
            this.name = name;
            this.location = location;
            this.occurrences = occurrences;
        }
    }

    private static class ReadPairInfo implements Serializable {

        private static final long serialVersionUID = 1L;
        private final String name;
        private final List<SimpleInterval> intervals;

        public ReadPairInfo(final String name, final List<SimpleInterval> intervals) {
            this.name = name;
            this.intervals = intervals;
        }

        public String name() {
            return name;
        }

        public boolean isUnmapped() {
            return intervals.isEmpty();
        }

        public List<ReadPairLocationInfo> locations() {
            if (isUnmapped()) {
                return Collections.singletonList(new ReadPairLocationInfo(name, null, 2));
            } else if (intervals.size() == 1) {
                return Collections.singletonList(new ReadPairLocationInfo(name, intervals.get(0), 2));
            } else {
                return Arrays.asList(new ReadPairLocationInfo(name, intervals.get(0), 1),
                                     new ReadPairLocationInfo(name, intervals.get(1), 1));
            }
        }
    }

        private SimpleInterval readLocation(final GATKRead read) {
            if (read.isUnmapped()) {
                if (read.mateIsUnmapped()) {
                    return null;
                } else {
                    return new SimpleInterval(read.getMateContig(), read.getMateStart(), read.getMateStart());
                }
            } else {
                return new SimpleInterval(read.getContig(), read.getStart(), read.getStart());
            }
        }

        private int readStartLocation(final GATKRead read) {
            if (read.isUnmapped()) {
                if (read.mateIsUnmapped()) {
                    return SAMRecord.NO_ALIGNMENT_START;
                } else {
                    return read.getMateStart();
                }
            } else {
                return read.getStart();
            }
        }


private static class ReadSortingAttributes implements Comparable<ReadSortingAttributes>, Serializable {

        private static final long serialVersionUID = 1L;

        private final GenomeLocParser dict;
        private final String contig;
        private final int start;
        private final String name;
        private final boolean isMapped;
        private final boolean hasCoordinates;
        private final boolean firstInPair;

        public SimpleInterval toInterval() {
            if (hasCoordinates)
                return new SimpleInterval(contig, start, start);
            else
                return null;
        }

        public ReadSortingAttributes(final GenomeLocParser dict, final GATKRead read) {
            this.dict = dict;
            this.name = read.getName();
            this.firstInPair = read.isFirstOfPair();
            this.isMapped = !read.isUnmapped();
            if (this.isMapped) {
                this.contig = read.getContig();
                this.start = read.getStart();
                this.hasCoordinates = true;
            } else if (!read.mateIsUnmapped()) {
                this.contig = read.getMateContig();
                this.start = read.getMateStart();
                this.hasCoordinates = true;
            } else {
                this.contig = null;
                this.start = 0;
                this.hasCoordinates = false;
            }
        }

       @Override
        public int compareTo(final ReadSortingAttributes other) {
            if (other.hasCoordinates == this.hasCoordinates) {
                if (this.hasCoordinates) {
                    final int otherContigIndex = dict.getContigIndex(other.contig);
                    final int thisContigIndex = dict.getContigIndex(this.contig);
                    if (otherContigIndex != thisContigIndex) {
                        return otherContigIndex < thisContigIndex ? 1 : -1;
                    } else if (this.start == other.start) {
                        return compareNameAndMapped(other);
                    } else {
                        return this.start < other.start ? -1 : 1;
                    }
                } else {
                    return compareNameAndMapped(other);
                }
            } else {
                return this.hasCoordinates ? -1 : 1;
            }
        }

        private int compareNameAndMapped(ReadSortingAttributes other) {
            if (this.name.equals(other.name)) {
                if (this.isMapped == other.isMapped) {
                    return this.firstInPair ? -1 : 1;
                } else {
                    return this.isMapped ? -1 : 1;
                }
            } else {
                return this.name.compareTo(other.name);
            }
        }
    }
}