package org.broadinstitute.hellbender.tools.genome;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
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
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Calculate coverage on a WGS bam file using Spark.  This creates a set of pseudo-targets that span" +
                "the entire genome.  Use the 'binsize' parameter to specify the size of each interval.  By default, any " +
                "contigs X, Y, M, and MT are excluded.",
        oneLineSummary = "Calculate coverage on a WGS bam file using Spark",
        programGroup = CopyNumberProgramGroup.class)
public class SparkGenomeReadCounts extends GATKSparkTool {
    private static final long serialVersionUID = 1l;

    private static final Set<String> NONAUTOSOMALCONTIGS = new HashSet<>(Arrays.asList("X", "Y", "MT", "M", "x", "y",
            "m", "chrX", "chrY", "chrMT", "chrM", "chrm"));

    protected static final String DROP_NON_AUTOSOMES_SHORT_NAME = "noxy";
    protected static final String DROP_NON_AUTOSOMES_LONG_NAME = "noXYMT";

    @Argument(doc = "Remove X, Y and MT regions",
            fullName = DROP_NON_AUTOSOMES_LONG_NAME,
            shortName = DROP_NON_AUTOSOMES_SHORT_NAME,
            optional = true
    )
    protected boolean dropNonAutosomes = true;

    private static final Logger logger = LogManager.getLogger(SparkGenomeReadCounts.class);
    protected static final String BINSIZE_SHORT_NAME = "bins";
    protected static final String BINSIZE_LONG_NAME = "binsize";

    @Argument(doc = "The size of bins in bases to collect coverage into",
            fullName = BINSIZE_LONG_NAME,
            shortName = BINSIZE_SHORT_NAME,
            optional = true)
    protected int binsize = 10000;

    protected static final String OUTPUT_FILE_SHORT_NAME = "o";
    protected static final String OUTPUT_FILE_LONG_NAME = "outputFile";

    public static final String RAW_COV_OUTPUT_EXTENSION = ".raw_cov";

    @Argument(doc = "Output tsv file for the proportional coverage counts.  Raw coverage counts will also be written with extension '" + RAW_COV_OUTPUT_EXTENSION + "'",
            fullName = OUTPUT_FILE_LONG_NAME,
            shortName = OUTPUT_FILE_SHORT_NAME,
            optional = false
    )
    protected File outputFile;

    private void collectReads(final JavaSparkContext ctx) {
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
        final String[] commentsForRawCoverage = {"##fileFormat  = tsv",
                "##commandLine = " + getCommandLine(),
                String.format("##title = Coverage counts in %d base bins for WGS", binsize)};

        final ReadFilter filter = makeGenomeReadFilter();
        final SAMSequenceDictionary originalSequenceDictionary = getReferenceSequenceDictionary();
        SAMSequenceDictionary modifiedSequenceDictionary = originalSequenceDictionary;
        if (dropNonAutosomes){
            modifiedSequenceDictionary = dropContigsFromSequence(originalSequenceDictionary, NONAUTOSOMALCONTIGS);
        }
        final SAMSequenceDictionary sequenceDictionary = modifiedSequenceDictionary;

        logger.info("Starting Spark coverage collection...");
        final long coverageCollectionStartTime = System.currentTimeMillis();
        final JavaRDD<GATKRead> rawReads = readSource.getParallelReads(bam, referenceArguments.getReferenceFileName(), (int) this.bamPartitionSplitSize);
        final JavaRDD<GATKRead> reads = rawReads.filter(read -> filter.test(read));
        final JavaRDD<SimpleInterval> readIntervals = reads
                .filter(read -> sequenceDictionary.getSequence(read.getContig()) != null)
                .map(read -> SparkGenomeReadCounts.createKey(read, sequenceDictionary, binsize));
        final Map<SimpleInterval, Long> byKey = readIntervals.countByValue();
        final Set<SimpleInterval> readIntervalKeySet = byKey.keySet();
        final long totalReads = byKey.values().stream().mapToLong(v -> v).sum();
        final long coverageCollectionEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished the spark coverage collection with %d targets and %d reads. Elapse of %d seconds",
                readIntervalKeySet.size(), totalReads, (coverageCollectionEndTime - coverageCollectionStartTime) / 1000));

        final String[] commentsForProportionalCoverage = {commentsForRawCoverage[0], commentsForRawCoverage[1],
                String.format("##title = Proportional coverage counts in %d base bins for WGS (total reads: %d)",
                        binsize, totalReads)};

        logger.info("Creating full genome bins...");
        final long createGenomeBinsStartTime = System.currentTimeMillis();
        final List<SimpleInterval> fullGenomeBins = createFullGenomeBins(binsize, sequenceDictionary);
        List<Target> fullGenomeTargetCollection = createTargetCollectionFromSimpleInterval(fullGenomeBins);
        TargetCoverageUtils.writeTargetsAsBed(new File(outputFile.getAbsolutePath() + ".bed"), fullGenomeTargetCollection);
        final long createGenomeBinsEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished creating genome bins. Elapse of %d seconds",
                (createGenomeBinsEndTime - createGenomeBinsStartTime) / 1000));

        logger.info("Creating missing genome bins...");
        final long createMissingGenomeBinsStartTime = System.currentTimeMillis();
        logger.info("Creating missing genome bins: Creating a mutable mapping...");
        final Map<SimpleInterval, Long> byKeyMutable = new HashMap<>();
        byKeyMutable.putAll(byKey);

        logger.info("Creating missing genome bins: Populating mutable mapping with zero counts for empty regions...");
        fullGenomeBins.stream().forEach(b -> byKeyMutable.putIfAbsent(b, 0l));

        final long createMissingGenomeBinsEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished creating missing genome bins. Elapse of %d seconds",
                (createMissingGenomeBinsEndTime - createMissingGenomeBinsStartTime) / 1000));

        logger.info("Creating final map...");
        final long createFinalMapStartTime = System.currentTimeMillis();
        final SortedMap<SimpleInterval, Long> byKeySorted = new TreeMap<>(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
        byKeySorted.putAll(byKeyMutable);
        final long createFinalMapEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished creating final map. Elapse of %d seconds",
                (createFinalMapEndTime - createFinalMapStartTime) / 1000));

        logger.info("Creating proportional coverage... ");
        final long pCovFileStartTime = System.currentTimeMillis();
        final SortedMap<SimpleInterval, Double> byKeyProportionalSorted = new TreeMap<>(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
        byKeySorted.entrySet().stream().forEach(e -> byKeyProportionalSorted.put(e.getKey(), (double) e.getValue() / totalReads));
        final long pCovFileEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished creating proportional coverage map. Elapse of %d seconds",
                (pCovFileEndTime - pCovFileStartTime) / 1000));

        logger.info("Writing raw coverage file ...");
        final long writingCovFileStartTime = System.currentTimeMillis();
        TargetCoverageUtils.writeTargetsWithCoverageFromSimpleInterval(
                new File(outputFile.getAbsolutePath() + RAW_COV_OUTPUT_EXTENSION),
                sampleName, byKeySorted, commentsForRawCoverage,
                TargetCoverageUtils.determineRawCoverageSingleSampleColumnCollection(sampleName));
        final long writingCovFileEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished writing coverage file. Elapse of %d seconds",
                (writingCovFileEndTime - writingCovFileStartTime) / 1000));

        logger.info("Writing proportional coverage file ...");
        final long writingPCovFileStartTime = System.currentTimeMillis();
        TargetCoverageUtils.writeTargetsWithCoverageFromSimpleInterval(outputFile, sampleName, byKeyProportionalSorted,
                commentsForProportionalCoverage, TargetCoverageUtils.determineRawCoverageSingleSampleColumnCollection(sampleName));
        final long writingPCovFileEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished writing proportional coverage file. Elapse of %d seconds",
                (writingPCovFileEndTime - writingPCovFileStartTime) / 1000));
    }

    protected static SAMSequenceDictionary dropContigsFromSequence(final SAMSequenceDictionary originalSequenceDictionary, final Set<String> nonAutosomalContigs) {
        final List<SAMSequenceRecord> initialSequences = new ArrayList<>(originalSequenceDictionary.getSequences());
        final List<SAMSequenceRecord> finalSequences = initialSequences.stream().filter(s -> !nonAutosomalContigs.contains(s.getSequenceName())).collect(Collectors.toList());
        return new SAMSequenceDictionary(finalSequences);
    }

    private static List<SimpleInterval> createFullGenomeBins(final int binsize, final SAMSequenceDictionary sequenceDictionary){
        final List<SimpleInterval> sequenceIntervals = IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        return IntervalUtils.cutToShards(sequenceIntervals, binsize);
    }

    private static SimpleInterval createKey(final GATKRead read, final SAMSequenceDictionary sequenceDictionary, final int binsize){
        final String contig = read.getContig();
        final int newStart = (read.getStart()/binsize)*binsize+1;
        final int newEnd = Math.min(newStart + binsize - 1, sequenceDictionary.getSequence(contig).getSequenceLength());

        return new SimpleInterval(contig, newStart, newEnd);
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        Utils.regularReadableUserFile(outputFile);
        collectReads(ctx);
    }

    /**
     * Filter reads based on marked parameters
     * @return never {@code null}
     */
    private ReadFilter makeGenomeReadFilter() {
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

