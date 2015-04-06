package org.broadinstitute.hellbender.tools.picard.illumina;

import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.IlluminaProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.illumina.AdapterMarker;
import org.broadinstitute.hellbender.utils.illumina.AdapterPair;
import org.broadinstitute.hellbender.utils.illumina.ClippingUtil;
import org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.IlluminaAdapterPair;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.ReservedTagConstants.XT;
import static htsjdk.samtools.SAMFileHeader.SortOrder;
import static htsjdk.samtools.SAMFileHeader.SortOrder.queryname;
import static htsjdk.samtools.SamReaderFactory.makeDefault;
import static htsjdk.samtools.util.CloserUtil.close;
import static htsjdk.samtools.util.CollectionUtil.makeList;
import static htsjdk.samtools.util.IOUtil.assertFileIsReadable;
import static htsjdk.samtools.util.IOUtil.assertFileIsWritable;
import static htsjdk.samtools.util.Log.getInstance;
import static htsjdk.samtools.util.SequenceUtil.reverseComplement;
import static htsjdk.samtools.util.StringUtil.stringToBytes;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.INPUT_SHORT_NAME;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.OUTPUT_SHORT_NAME;
import static org.broadinstitute.hellbender.utils.illumina.AdapterMarker.DEFAULT_ADAPTER_LENGTH;
import static org.broadinstitute.hellbender.utils.illumina.AdapterMarker.DEFAULT_NUM_ADAPTERS_TO_KEEP;
import static org.broadinstitute.hellbender.utils.illumina.AdapterMarker.DEFAULT_PRUNE_ADAPTER_LIST_AFTER_THIS_MANY_ADAPTERS_SEEN;
import static org.broadinstitute.hellbender.utils.illumina.ClippingUtil.MAX_ERROR_RATE;
import static org.broadinstitute.hellbender.utils.illumina.ClippingUtil.MAX_PE_ERROR_RATE;
import static org.broadinstitute.hellbender.utils.illumina.ClippingUtil.MIN_MATCH_BASES;
import static org.broadinstitute.hellbender.utils.illumina.ClippingUtil.MIN_MATCH_PE_BASES;
import static org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.IlluminaAdapterPair.DUAL_INDEXED;
import static org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.IlluminaAdapterPair.INDEXED;
import static org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.IlluminaAdapterPair.PAIRED_END;

/**
 * Command line program to mark the location of adapter sequences.
 * This also outputs a Histogram of metrics describing the clipped bases
 *
 * @author Tim Fennell (adapted by mborkan@broadinstitute.org)
 */
@CommandLineProgramProperties(
        usage = "Reads a SAM or BAM file and rewrites it with new adapter-trimming tags.\n" +
                "Clear any existing adapter-trimming tags (XT:i:).\n" +
                "Only works for unaligned files in query-name order.\n" +
                "Note: This is a utility program and will not be run in the pipeline.\n",
        usageShort = "Reads a SAM or BAM file and rewrites it with new adapter-trimming tags",
        programGroup = IlluminaProgramGroup.class
)
public class MarkIlluminaAdapters extends PicardCommandLineProgram {

    // The following attributes define the command-line arguments

    @Argument(shortName = INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "If output is not specified, just the metrics are generated",
            shortName = OUTPUT_SHORT_NAME, optional = true)
    public File OUTPUT;

    @Argument(doc = "Histogram showing counts of bases_clipped in how many reads", shortName = "M")
    public File METRICS;

    @Argument(doc = "The minimum number of bases to match over when clipping single-end reads.")
    public int MIN_MATCH_BASES_SE = MIN_MATCH_BASES;

    @Argument(doc = "The minimum number of bases to match over (per-read) when clipping paired-end reads.")
    public int MIN_MATCH_BASES_PE = MIN_MATCH_PE_BASES;

    @Argument(doc = "The maximum mismatch error rate to tolerate when clipping single-end reads.")
    public double MAX_ERROR_RATE_SE = MAX_ERROR_RATE;

    @Argument(doc = "The maximum mismatch error rate to tolerate when clipping paired-end reads.")
    public double MAX_ERROR_RATE_PE = MAX_PE_ERROR_RATE;

    @Argument(doc = "Which adapters sequences to attempt to identify and clip.")
    public List<IlluminaAdapterPair> ADAPTERS =
            makeList(INDEXED,
                    DUAL_INDEXED,
                    PAIRED_END
            );

    @Argument(doc = "For specifying adapters other than standard Illumina", optional = true)
    public String FIVE_PRIME_ADAPTER;
    @Argument(doc = "For specifying adapters other than standard Illumina", optional = true)
    public String THREE_PRIME_ADAPTER;

    @Argument(doc = "Adapters are truncated to this length to speed adapter matching.  Set to a large number to effectively disable truncation.")
    public int ADAPTER_TRUNCATION_LENGTH = DEFAULT_ADAPTER_LENGTH;

    @Argument(doc = "If looking for multiple adapter sequences, then after having seen this many adapters, shorten the list of sequences. " +
            "Keep the adapters that were found most frequently in the input so far. " +
            "Set to -1 if the input has a heterogeneous mix of adapters so shortening is undesirable.",
            shortName = "APT")
    public int PRUNE_ADAPTER_LIST_AFTER_THIS_MANY_ADAPTERS_SEEN = DEFAULT_PRUNE_ADAPTER_LIST_AFTER_THIS_MANY_ADAPTERS_SEEN;

    @Argument(doc = "If pruning the adapter list, keep only this many adapter sequences when pruning the list (plus any adapters that " +
            "were tied with the adapters being kept).")
    public int NUM_ADAPTERS_TO_KEEP = DEFAULT_NUM_ADAPTERS_TO_KEEP;

    private static final Log log = getInstance(MarkIlluminaAdapters.class);

    @Override
    protected String[] customCommandLineValidation() {
        if ((FIVE_PRIME_ADAPTER != null && THREE_PRIME_ADAPTER == null) || (THREE_PRIME_ADAPTER != null && FIVE_PRIME_ADAPTER == null)) {
            return new String[]{"Either both or neither of THREE_PRIME_ADAPTER and FIVE_PRIME_ADAPTER must be set."};
        } else {
            return null;
        }
    }

    @Override
    protected Object doWork() {
        assertFileIsReadable(INPUT);
        assertFileIsWritable(METRICS);

        final SamReader in = makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SortOrder order = in.getFileHeader().getSortOrder();
        SAMFileWriter out = null;
        if (OUTPUT != null) {
            assertFileIsWritable(OUTPUT);
            out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(), true, OUTPUT);
        }

        final Histogram<Integer> histo = new Histogram<Integer>("clipped_bases", "read_count");

        // Combine any adapters and custom adapter pairs from the command line into an array for use in clipping
        final AdapterPair[] adapters;
        {
            final List<AdapterPair> tmp = new ArrayList<AdapterPair>();
            tmp.addAll(ADAPTERS);
            if (FIVE_PRIME_ADAPTER != null && THREE_PRIME_ADAPTER != null) {
                tmp.add(new CustomAdapterPair(FIVE_PRIME_ADAPTER, THREE_PRIME_ADAPTER));
            }
            adapters = tmp.toArray(new AdapterPair[tmp.size()]);
        }

        ////////////////////////////////////////////////////////////////////////
        // Main loop that consumes reads, clips them and writes them to the output
        ////////////////////////////////////////////////////////////////////////
        final ProgressLogger progress = new ProgressLogger(log, 1000000, "Read");
        final SAMRecordIterator iterator = in.iterator();

        final AdapterMarker adapterMarker = new AdapterMarker(ADAPTER_TRUNCATION_LENGTH, adapters).
                setMaxPairErrorRate(MAX_ERROR_RATE_PE).setMinPairMatchBases(MIN_MATCH_BASES_PE).
                setMaxSingleEndErrorRate(MAX_ERROR_RATE_SE).setMinSingleEndMatchBases(MIN_MATCH_BASES_SE).
                setNumAdaptersToKeep(NUM_ADAPTERS_TO_KEEP).
                setThresholdForSelectingAdaptersToKeep(PRUNE_ADAPTER_LIST_AFTER_THIS_MANY_ADAPTERS_SEEN);

        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            final SAMRecord rec2 = rec.getReadPairedFlag() && iterator.hasNext() ? iterator.next() : null;
            rec.setAttribute(XT, null);

            // Do the clipping one way for PE and another for SE reads
            if (rec.getReadPairedFlag()) {
                // Assert that the input file is in query name order only if we see some PE reads
                if (order != queryname) {
                    throw new UserException("Input BAM file must be sorted by queryname");
                }

                if (rec2 == null) throw new UserException("Missing mate pair for paired read: " + rec.getReadName());
                rec2.setAttribute(XT, null);

                // Assert that we did in fact just get two mate pairs
                if (!rec.getReadName().equals(rec2.getReadName())) {
                    throw new UserException("Adjacent reads expected to be mate-pairs have different names: " +
                            rec.getReadName() + ", " + rec2.getReadName());
                }

                // establish which of pair is first and which second
                final SAMRecord first, second;

                if (rec.getFirstOfPairFlag() && rec2.getSecondOfPairFlag()) {
                    first = rec;
                    second = rec2;
                } else if (rec.getSecondOfPairFlag() && rec2.getFirstOfPairFlag()) {
                    first = rec2;
                    second = rec;
                } else {
                    throw new UserException("Two reads with same name but not correctly marked as 1st/2nd of pair: " + rec.getReadName());
                }

                adapterMarker.adapterTrimIlluminaPairedReads(first, second);
            } else {
                adapterMarker.adapterTrimIlluminaSingleRead(rec);
            }

            // Then output the records, update progress and metrics
            for (final SAMRecord r : new SAMRecord[]{rec, rec2}) {
                if (r != null) {
                    progress.record(r);
                    if (out != null) out.addAlignment(r);

                    final Integer clip = rec.getIntegerAttribute(XT);
                    if (clip != null) histo.increment(rec.getReadLength() - clip + 1);
                }
            }
        }

        if (out != null) out.close();

        // Lastly output the metrics to file
        final MetricsFile<?, Integer> metricsFile = getMetricsFile();
        metricsFile.setHistogram(histo);
        metricsFile.write(METRICS);

        close(in);
        return 0;
    }

    private class CustomAdapterPair implements AdapterPair {

        final String fivePrime, threePrime, fivePrimeReadOrder;
        final byte[] fivePrimeBytes, threePrimeBytes, fivePrimeReadOrderBytes;

        private CustomAdapterPair(final String fivePrime, final String threePrime) {
            this.threePrime = threePrime;
            this.threePrimeBytes = stringToBytes(threePrime);

            this.fivePrime = fivePrime;
            this.fivePrimeReadOrder = reverseComplement(fivePrime);
            this.fivePrimeBytes = stringToBytes(fivePrime);
            this.fivePrimeReadOrderBytes = stringToBytes(fivePrimeReadOrder);
        }

        public String get3PrimeAdapter() {
            return threePrime;
        }

        public String get5PrimeAdapter() {
            return fivePrime;
        }

        public String get3PrimeAdapterInReadOrder() {
            return threePrime;
        }

        public String get5PrimeAdapterInReadOrder() {
            return fivePrimeReadOrder;
        }

        public byte[] get3PrimeAdapterBytes() {
            return threePrimeBytes;
        }

        public byte[] get5PrimeAdapterBytes() {
            return fivePrimeBytes;
        }

        public byte[] get3PrimeAdapterBytesInReadOrder() {
            return threePrimeBytes;
        }

        public byte[] get5PrimeAdapterBytesInReadOrder() {
            return fivePrimeReadOrderBytes;
        }

        public String getName() {
            return "Custom adapter pair";
        }
    }
}
