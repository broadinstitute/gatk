package org.broadinstitute.hellbender.tools.genome;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.SampleCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetCoverageUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Spark WGS Coverage",
        oneLineSummary = "Spark WGS Coverage.",
        programGroup = CopyNumberProgramGroup.class)
public class SparkGenomeReadCounts extends GATKSparkTool {
    private static final long serialVersionUID = 1l;

    protected static final String BINSIZE_SHORT_NAME = "bins";
    protected static final String BINSIZE_LONG_NAME = "binsize";

    @Argument(doc = "The size of bins in bases to collect coverage into",
            fullName = BINSIZE_LONG_NAME,
            shortName = BINSIZE_SHORT_NAME,
            optional = true)
    protected static int binsize = 10000;

    protected static final String OUTPUT_FILE_SHORT_NAME = "o";
    protected static final String OUTPUT_FILE_LONG_NAME = "outputFile";

    @Argument(doc = "Output tsv file for the coverage counts",
            fullName = OUTPUT_FILE_LONG_NAME,
            shortName = OUTPUT_FILE_SHORT_NAME,
            optional = false
    )
    protected File outputFile;

    private static final Logger logger = LogManager.getLogger(SparkGenomeReadCounts.class);

    public void collectReads(final JavaSparkContext ctx) {
        if ( readArguments.getReadFilesNames().size() != 1 ) {
            throw new UserException("This tool only accepts a single bam/sam/cram as input");
        }
        final String bam = readArguments.getReadFilesNames().get(0);
        final ReadsSparkSource readSource = new ReadsSparkSource(ctx);

        final SampleCollection sampleCollection = new SampleCollection(getHeaderForReads());
        if(sampleCollection.sampleCount()>1){
            throw new UserException.BadInput("We do not support bams with more than one sample.");
        }
        final String sampleName = sampleCollection.sampleIds().get(0);
        final String[] comments = {"##fileFormat  = tsv",
                "##commandLine = " + getCommandLine(),
                String.format("##title = Coverage counts in %d base bins for WGS", binsize)};

        final ReadFilter filter = makeReadFilter();
//      TODO: Add interval parameter for getParallelReads to subset to non-X/Y/MT/etc. regions - issue filed #287

        logger.info("Starting Spark coverage collection...");
        final long coverageCollectionStartTime = System.currentTimeMillis();
        final JavaRDD<GATKRead> rawReads = readSource.getParallelReads(bam, referenceArguments.getReferenceFileName(), (int) this.bamPartitionSplitSize);
        final JavaRDD<GATKRead> reads = rawReads.filter(read -> filter.test(read));
        final SAMSequenceDictionary sequenceDictionary = getReferenceSequenceDictionary();
        final Map<SimpleInterval, Long> byKey = reads.map(read -> SparkGenomeReadCounts.createKey(read, sequenceDictionary)).countByValue();
        final long coverageCollectionEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished the spark coverage collection. Elapse of %d seconds",
                (coverageCollectionEndTime - coverageCollectionStartTime) / 1000));

        final List<SimpleInterval> fullGenomeBins = createFullGenomeBins(binsize, sequenceDictionary);
        List<Target> fullGenomeTargetCollection = createTargetCollectionFromSimpleInterval(fullGenomeBins);
        TargetCoverageUtils.writeTargetsAsBed(new File(outputFile.getAbsolutePath() + ".bed"), fullGenomeTargetCollection);

        final Set<SimpleInterval> missingGenomeBins = Sets.difference( new HashSet<>(fullGenomeBins),byKey.keySet());
        final Map<SimpleInterval, Long> byKeyMutable = new HashMap<>();
        byKeyMutable.putAll(byKey);
        missingGenomeBins.stream().forEach(m -> byKeyMutable.put(m, 0l));

        final SortedMap<SimpleInterval, Long> byKeySorted = new TreeMap<>(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
        byKeySorted.putAll(byKeyMutable);

        TargetCoverageUtils.writeTargetsWithCoverageFromSimpleInterval(outputFile, sampleName, byKeySorted, comments);
    }

    private static List<SimpleInterval> createFullGenomeBins(final int binsize, final SAMSequenceDictionary sequenceDictionary){
        final List<SimpleInterval> sequenceIntervals = IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        return IntervalUtils.cutToShards(sequenceIntervals, binsize);
    }

    private static SimpleInterval createKey(final GATKRead read, final SAMSequenceDictionary sequenceDictionary){
        final String contig = read.getContig();
        final int newStart = (read.getStart()/binsize)*binsize+1;
        final int newEnd = Math.min(newStart + binsize - 1, sequenceDictionary.getSequence(contig).getSequenceLength());

        return new SimpleInterval(contig, newStart, newEnd);
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        collectReads(ctx);
    }

    /**
     * Filter reads based on marked parameters
     * @return never {@code null}
     */
    public ReadFilter makeReadFilter() {
        return new WellformedReadFilter(getHeaderForReads())
                .and(ReadFilterLibrary.MAPPED)
                .and(ReadFilterLibrary.NOT_DUPLICATE)
                .and(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT)
                .and(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
    }

    private List<Target> createTargetCollectionFromSimpleInterval(final List<SimpleInterval> intervalList){
        return intervalList.stream().map(s->new Target(TargetCoverageUtils.createDummyTargetName(s),
                new SimpleInterval(s))).collect(Collectors.toList());
    }
}

