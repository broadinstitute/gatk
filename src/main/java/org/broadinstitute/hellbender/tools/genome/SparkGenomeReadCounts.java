package org.broadinstitute.hellbender.tools.genome;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.SampleCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetWriter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Calculates read coverage on whole genome sequencing (WGS) alignments using Spark.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 *
 * <h3>Examples</h3>
 *
 * <p>
 *     The command encompasses empirically determined parameters for TCGA project data.
 *     You may obtain better results with different parameters.
 * </p>
 *
 * <p>For whole genome sequencing (WGS) data: </p>
 *
 * <pre>
 * java -Xmx4g -jar ${gatk_jar} SparkGenomeReadCounts \
 *   --input bam.bam \
 *   --reference ref_fasta.fa \
 *   --disableReadFilter NotDuplicateReadFilter \
 *   --outputFile base_filename.coverage.tsv
 * </pre>
 *
 * <p>
 *     The reference is optional if --disableSequenceDictionaryValidation is set to true.
 * </p>
 *
 * <p>
 *     For whole exome sequencing (WGS) data, use CalculateTargetCoverage instead.
 * </p>
 *
 */
@CommandLineProgramProperties(
        summary = "Calculate coverage on a WGS bam file using Spark.  This creates a set of pseudo-targets that span " +
                "the entire genome.  Use the 'binsize' parameter to specify the size of each interval.  By default, any " +
                "contigs X, Y, M, and MT are excluded.\n" +
        "Please see the " + SparkGenomeReadCounts.DROP_NON_AUTOSOMES_LONG_NAME + " option if using this tool on a non-human genome.",
        oneLineSummary = "Calculate coverage on a WGS bam file using Spark",
        programGroup = CopyNumberProgramGroup.class)
@DocumentedFeature
public class SparkGenomeReadCounts extends GATKSparkTool {
    private static final long serialVersionUID = 1l;

    private static final Set<String> NONAUTOSOMALCONTIGS = new HashSet<>(Arrays.asList("X", "Y", "MT", "M", "x", "y",
            "m", "chrX", "chrY", "chrMT", "chrM", "chrm"));

    protected static final String DROP_NON_AUTOSOMES_SHORT_NAME = "keepxy";
    protected static final String DROP_NON_AUTOSOMES_LONG_NAME = "keepXYMT";

    @Argument(doc = "Keep X, Y, GL*, NC_*, and MT regions.  If this option is not specified, these regions will be dropped, regardless of intervals specified.  Use -L (or -XL) and enable this option for exact specification of intervals.  This option may be removed in the future.",
            fullName = DROP_NON_AUTOSOMES_LONG_NAME,
            shortName = DROP_NON_AUTOSOMES_SHORT_NAME,
            optional = true
    )
    protected boolean keepNonAutosomes = false;

    private static final Logger logger = LogManager.getLogger(SparkGenomeReadCounts.class);
    protected static final String BINSIZE_SHORT_NAME = "bins";
    protected static final String BINSIZE_LONG_NAME = "binsize";

    @Argument(doc = "The size of bins for each interval specified.  E.g. chr2:100-200 --> chr2:100-150, chr2:151-200 if binsize = 50.",
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

    /**
     * Determine the intervals to consider for coverage collection.  Honors the keepAutosome parameter.
     *
     * <p>Developer's note:  This CLI will always set the attribute, intervals, to a non-null value.</p>
     *
     * @param rawIntervals Specified by the user.  If null, converts to SimpleIntervals specifying the entire
     *                     reference genome.  If keepNonAutosomes is NOT specified, it will prune these intervals (or the
     *                     ones specified by the user), to remove the contigs that are listed in the docs for
     *                     {@link SparkGenomeReadCounts#keepNonAutosomes}
     * @return Never {@code null}  Specified list of intervals.  These will be treated as if user had specified on the
     *  CLI.
     */
    @Override
    protected List<SimpleInterval> editIntervals(final List<SimpleInterval> rawIntervals) {
        List<SimpleInterval> modifiedIntervals = rawIntervals;
        if (rawIntervals == null) {
            modifiedIntervals = IntervalUtils.getAllIntervalsForReference(getReferenceSequenceDictionary());
        }

        if (keepNonAutosomes) {
            return modifiedIntervals;
        }

        // Enforce the elimination of certain contigs when proper option is set.
        logger.info("Dropping non-autosomes, as requested...");
        return modifiedIntervals.stream()
                .filter(s -> !NONAUTOSOMALCONTIGS.contains(s.getContig()))
                .filter(s -> !(s.getContig().startsWith("GL")) && !(s.getContig().startsWith("NC_")))
                .collect(Collectors.toList());
    }

    private void collectReads() {
        if ( readArguments.getReadFilesNames().size() != 1 ) {
            throw new UserException("This tool only accepts a single bam/sam/cram as input");
        }

        final SampleCollection sampleCollection = new SampleCollection(getHeaderForReads());
        if(sampleCollection.sampleCount()>1){
            throw new UserException.BadInput("We do not support bams with more than one sample.");
        }
        final String sampleName = sampleCollection.sampleIds().get(0);
        final String[] commentsForRawCoverage = {"##fileFormat  = tsv",
                "##commandLine = " + getCommandLine(),
                String.format("##title = Coverage counts in %d base bins for WGS", binsize)};

        final ReadFilter filter = makeGenomeReadFilter();
        final SAMSequenceDictionary sequenceDictionary = getReferenceSequenceDictionary();

        logger.info("Starting Spark coverage collection...");
        final long coverageCollectionStartTime = System.currentTimeMillis();
        final JavaRDD<GATKRead> rawReads = getReads();
        final JavaRDD<GATKRead> reads = rawReads.filter(read -> filter.test(read));

        //Note: using a field inside a closure will pull in the whole enclosing object to serialization
        // (which leads to bad performance and can blow up if some objects in the fields are not
        // Serializable - closures always use java Serializable and not Kryo)
        //Solution here is to use a temp variable for binsize because it's just an int.
        final int binsize_tmp = binsize;
        final JavaRDD<SimpleInterval> readIntervals = reads
                .filter(read -> sequenceDictionary.getSequence(read.getContig()) != null)
                .map(read -> SparkGenomeReadCounts.createKey(read, sequenceDictionary, binsize_tmp));
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
        final List<SimpleInterval> fullGenomeBins = createFullGenomeBins(binsize);
        List<Target> fullGenomeTargetCollection = createTargetListFromSimpleInterval(fullGenomeBins);
        TargetWriter.writeTargetsToFile(new File(outputFile.getAbsolutePath() + ".targets.tsv"), fullGenomeTargetCollection);
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
        ReadCountCollectionUtils.writeReadCountsFromSimpleInterval(new File(outputFile.getAbsolutePath() + RAW_COV_OUTPUT_EXTENSION), sampleName, byKeySorted, commentsForRawCoverage);
        final long writingCovFileEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished writing coverage file. Elapse of %d seconds",
                (writingCovFileEndTime - writingCovFileStartTime) / 1000));

        logger.info("Writing proportional coverage file ...");
        final long writingPCovFileStartTime = System.currentTimeMillis();
        ReadCountCollectionUtils.writeReadCountsFromSimpleInterval(outputFile, sampleName, byKeyProportionalSorted,
                commentsForProportionalCoverage);
        final long writingPCovFileEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished writing proportional coverage file. Elapse of %d seconds",
                (writingPCovFileEndTime - writingPCovFileStartTime) / 1000));
    }

    private List<SimpleInterval> createFullGenomeBins(final int binsize){
        return IntervalUtils.cutToShards(getIntervals(), binsize);
    }

    private static SimpleInterval createKey(final GATKRead read, final SAMSequenceDictionary sequenceDictionary, final int binsize){
        final String contig = read.getContig();
        final int newStart = (read.getStart()/binsize)*binsize+1;
        final int newEnd = Math.min(newStart + binsize - 1, sequenceDictionary.getSequence(contig).getSequenceLength());

        return new SimpleInterval(contig, newStart, newEnd);
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        collectReads();
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

    private List<Target> createTargetListFromSimpleInterval(final List<SimpleInterval> intervalList){
        return intervalList.stream().map(Target::new).collect(Collectors.toList());
    }
}

