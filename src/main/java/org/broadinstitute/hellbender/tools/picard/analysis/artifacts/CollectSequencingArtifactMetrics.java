package org.broadinstitute.hellbender.tools.picard.analysis.artifacts;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalListReferenceSequenceMask;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.picard.analysis.SinglePassSamProgram;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.SequencingArtifactMetrics.BaitBiasDetailMetrics;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.SequencingArtifactMetrics.BaitBiasSummaryMetrics;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.SequencingArtifactMetrics.PreAdapterDetailMetrics;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.SequencingArtifactMetrics.PreAdapterSummaryMetrics;
import org.broadinstitute.hellbender.utils.variant.DbSnpBitSetUtil;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.CodeUtil.getOrElse;

/**
 * Quantify substitution errors caused by mismatched base pairings during various
 * stages of sample / library prep.
 *
 * We measure two distinct error types - artifacts that are introduced before
 * the addition of the read1/read2 adapters ("pre adapter") and those that are
 * introduced after target selection ("bait bias"). For each of these, we provide
 * summary metrics as well as detail metrics broken down by reference context
 * (the ref bases surrounding the substitution event).
 *
 * For a deeper explanation, see Costello et al. 2013:
 * http://www.ncbi.nlm.nih.gov/pubmed/23303777
 *
 * @author mattsooknah
 *
 */
@CommandLineProgramProperties(
        summary = CollectSequencingArtifactMetrics.USAGE,
        oneLineSummary = CollectSequencingArtifactMetrics.USAGE,
        programGroup = QCProgramGroup.class
)
public final class CollectSequencingArtifactMetrics extends SinglePassSamProgram {
    static final String USAGE = "Produces metrics to quantify single-base sequencing artifacts from a SAM/BAM file";

    @Argument(doc = "An optional list of intervals to restrict analysis to.", optional = true)
    public File INTERVALS;

    @Argument(doc = "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis.", optional = true)
    public File DB_SNP;

    @Argument(shortName = "Q", doc = "The minimum base quality score for a base to be included in analysis.")
    public int MINIMUM_QUALITY_SCORE = 20;

    @Argument(shortName = "MQ", doc = "The minimum mapping quality score for a base to be included in analysis.")
    public int MINIMUM_MAPPING_QUALITY = 30;

    @Argument(shortName = "MIN_INS", doc = "The minimum insert size for a read to be included in analysis.")
    public int MINIMUM_INSERT_SIZE = 60;

    @Argument(shortName = "MAX_INS", doc = "The maximum insert size for a read to be included in analysis. Set to 0 to have no maximum.")
    public int MAXIMUM_INSERT_SIZE = 600;

    @Argument(shortName = "UNPAIRED", doc = "Include unpaired reads. If set to true then all paired reads will be included as well - " +
            "MINIMUM_INSERT_SIZE and MAXIMUM_INSERT_SIZE will be ignored.")
    public boolean INCLUDE_UNPAIRED = false;

    @Argument(shortName = "TANDEM", doc = "Set to true if mate pairs are being sequenced from the same strand, " +
            "i.e. they're expected to face the same direction.")
    public boolean TANDEM_READS = false;

    @Argument(doc = "When available, use original quality scores for filtering.")
    public boolean USE_OQ = true;

    @Argument(doc = "The number of context bases to include on each side of the assayed base.")
    public int CONTEXT_SIZE = 1;

    @Argument(doc = "If specified, only print results for these contexts in the detail metrics output. " +
                  "However, the summary metrics output will still take all contexts into consideration.", optional = true)
    public Set<String> CONTEXTS_TO_PRINT = new HashSet<>();

    private static final String UNKNOWN_LIBRARY = "UnknownLibrary";
    private static final String UNKNOWN_SAMPLE = "UnknownSample";

    private File preAdapterSummaryOut;
    private File preAdapterDetailsOut;
    private File baitBiasSummaryOut;
    private File baitBiasDetailsOut;

    private IntervalListReferenceSequenceMask intervalMask;
    private DbSnpBitSetUtil dbSnpMask;
    private SamRecordFilter recordFilter;

    private final Set<String> samples = new HashSet<>();
    private final Set<String> libraries = new HashSet<>();
    private final Map<String, ArtifactCounter> artifactCounters = new HashMap<>();

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> messages = new ArrayList<>();

        final int contextFullLength = 2 * CONTEXT_SIZE + 1;
        if (CONTEXT_SIZE < 0) messages.add("CONTEXT_SIZE cannot be negative");
        for (final String context : CONTEXTS_TO_PRINT) {
            if (context.length() != contextFullLength) {
                messages.add("Context " + context + " is not the length implied by CONTEXT_SIZE: " + contextFullLength);
            }
        }

        if (MINIMUM_INSERT_SIZE < 0) messages.add("MINIMUM_INSERT_SIZE cannot be negative");
        if (MAXIMUM_INSERT_SIZE < 0) messages.add("MAXIMUM_INSERT_SIZE cannot be negative");
        if (MAXIMUM_INSERT_SIZE > 0 && MAXIMUM_INSERT_SIZE < MINIMUM_INSERT_SIZE) {
            messages.add("MAXIMUM_INSERT_SIZE cannot be less than MINIMUM_INSERT_SIZE unless set to 0");
        }

        return messages.isEmpty() ? null : messages.toArray(new String[messages.size()]);
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        preAdapterSummaryOut = new File(OUTPUT + SequencingArtifactMetrics.PRE_ADAPTER_SUMMARY_EXT);
        preAdapterDetailsOut = new File(OUTPUT + SequencingArtifactMetrics.PRE_ADAPTER_DETAILS_EXT);
        baitBiasSummaryOut = new File(OUTPUT + SequencingArtifactMetrics.BAIT_BIAS_SUMMARY_EXT);
        baitBiasDetailsOut = new File(OUTPUT + SequencingArtifactMetrics.BAIT_BIAS_DETAILS_EXT);

        IOUtil.assertFileIsWritable(preAdapterSummaryOut);
        IOUtil.assertFileIsWritable(preAdapterDetailsOut);
        IOUtil.assertFileIsWritable(baitBiasSummaryOut);
        IOUtil.assertFileIsWritable(baitBiasDetailsOut);

        for (final SAMReadGroupRecord rec : header.getReadGroups()) {
            samples.add(getOrElse(rec.getSample(), UNKNOWN_SAMPLE));
            libraries.add(getOrElse(rec.getLibrary(), UNKNOWN_LIBRARY));
        }

        if (INTERVALS != null) {
            IOUtil.assertFileIsReadable(INTERVALS);
            intervalMask = new IntervalListReferenceSequenceMask(IntervalList.fromFile(INTERVALS).uniqued());
        }

        if (DB_SNP != null) {
            IOUtil.assertFileIsReadable(DB_SNP);
            dbSnpMask = new DbSnpBitSetUtil(DB_SNP, header.getSequenceDictionary());
        }

        // set record-level filters
        final List<SamRecordFilter> filters = new ArrayList<>();
        filters.add(new FailsVendorReadQualityFilter());
        filters.add(new NotPrimaryAlignmentFilter());
        filters.add(new DuplicateReadFilter());
        filters.add(new AlignedFilter(true)); // discard unmapped reads
        filters.add(new MappingQualityFilter(MINIMUM_MAPPING_QUALITY));
        if (!INCLUDE_UNPAIRED) {
            final int effectiveMaxInsertSize = (MAXIMUM_INSERT_SIZE == 0) ? Integer.MAX_VALUE : MAXIMUM_INSERT_SIZE;
            filters.add(new InsertSizeFilter(MINIMUM_INSERT_SIZE, effectiveMaxInsertSize));
        }
        recordFilter = new AggregateFilter(filters);

        // set up the artifact counters
        final String sampleAlias = StringUtil.join(",", new ArrayList<>(samples));
        for (final String library : libraries) {
            artifactCounters.put(library, new ArtifactCounter(sampleAlias, library, CONTEXT_SIZE, TANDEM_READS));
        }
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        // see if the whole read should be skipped
        if (recordFilter.filterOut(rec)) return;

        // check read group + library
        final String library = (rec.getReadGroup() == null) ? UNKNOWN_LIBRARY : getOrElse(rec.getReadGroup().getLibrary(), UNKNOWN_LIBRARY);
        if (!libraries.contains(library)) {
            // should never happen if SAM is valid
            throw new UserException("Record contains library that is missing from header: " + library);
        }

        // iterate over aligned positions
        for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
            for (int offset = 0; offset < block.getLength(); offset++) {
                // remember, these are 1-based!
                final int readPos = block.getReadStart() + offset;
                final int refPos = block.getReferenceStart() + offset;

                /**
                 * Skip regions outside of intervals.
                 *
                 * NB: IntervalListReferenceSequenceMask.get() has side-effects which assume
                 * that successive ReferenceSequence's passed to this method will be in-order
                 * (e.g. it will break if you call acceptRead() with chr1, then chr2, then chr1
                 * again). So this only works if the underlying iteration is coordinate-sorted.
                 */
                if (intervalMask != null && !intervalMask.get(ref.getContigIndex(), refPos)) continue;

                // skip dbSNP sites
                if (dbSnpMask != null && dbSnpMask.isDbSnpSite(ref.getName(), refPos)) continue;

                // skip the ends of the reference
                final int contextStartIndex = refPos - CONTEXT_SIZE - 1;
                final int contextFullLength = 2 * CONTEXT_SIZE + 1;
                if (contextStartIndex < 0 || contextStartIndex + contextFullLength > ref.length()) continue;

                // skip contexts with N bases
                final String context = StringUtil.bytesToString(ref.getBases(), contextStartIndex, contextFullLength).toUpperCase();
                if (context.contains("N")) continue;

                // skip low BQ sites
                if (failsBaseQualityCutoff(readPos, rec)) continue;

                // skip N bases in read
                final char readBase = Character.toUpperCase((char) rec.getReadBases()[readPos - 1]);
                if (readBase == 'N') continue;

                // count the base!
                artifactCounters.get(library).countRecord(context, readBase, rec);
            }
        }
    }

    @Override
    protected void finish() {
        final MetricsFile<PreAdapterSummaryMetrics, Integer> preAdapterSummaryMetricsFile = getMetricsFile();
        final MetricsFile<PreAdapterDetailMetrics, Integer> preAdapterDetailMetricsFile = getMetricsFile();
        final MetricsFile<BaitBiasSummaryMetrics, Integer> baitBiasSummaryMetricsFile = getMetricsFile();
        final MetricsFile<BaitBiasDetailMetrics, Integer> baitBiasDetailMetricsFile = getMetricsFile();

        for (final ArtifactCounter counter : artifactCounters.values()) {
            // build metrics
            counter.finish();

            // write metrics
            preAdapterSummaryMetricsFile.addAllMetrics(counter.getPreAdapterSummaryMetrics());
            baitBiasSummaryMetricsFile.addAllMetrics(counter.getBaitBiasSummaryMetrics());

            for (final PreAdapterDetailMetrics preAdapterDetailMetrics : counter.getPreAdapterDetailMetrics()) {
                if (CONTEXTS_TO_PRINT.size() == 0 || CONTEXTS_TO_PRINT.contains(preAdapterDetailMetrics.CONTEXT)) {
                    preAdapterDetailMetricsFile.addMetric(preAdapterDetailMetrics);
                }
            }
            for (final BaitBiasDetailMetrics baitBiasDetailMetrics : counter.getBaitBiasDetailMetrics()) {
                if (CONTEXTS_TO_PRINT.size() == 0 || CONTEXTS_TO_PRINT.contains(baitBiasDetailMetrics.CONTEXT)) {
                    baitBiasDetailMetricsFile.addMetric(baitBiasDetailMetrics);
                }
            }

        }

        preAdapterDetailMetricsFile.write(preAdapterDetailsOut);
        preAdapterSummaryMetricsFile.write(preAdapterSummaryOut);
        baitBiasDetailMetricsFile.write(baitBiasDetailsOut);
        baitBiasSummaryMetricsFile.write(baitBiasSummaryOut);
    }

    @Override
    protected boolean usesNoRefReads() { return false; }

    /**
     * Check if this read base fails the base quality cutoff.
     */
    private boolean failsBaseQualityCutoff(final int oneIndexedPos, final SAMRecord rec) {
        final byte qual;
        if (USE_OQ && rec.getOriginalBaseQualities() != null) {
            qual = rec.getOriginalBaseQualities()[oneIndexedPos - 1];
        } else {
            qual = rec.getBaseQualities()[oneIndexedPos - 1];
        }
        return (qual < MINIMUM_QUALITY_SCORE);
    }
}
