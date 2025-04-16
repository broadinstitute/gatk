package org.broadinstitute.hellbender.tools;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Multiset;
import com.google.common.collect.HashMultiset;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.WorkflowOutput;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;

import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.MetadataUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.CachedOverlapDetector;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Compare two aligners. This tool takes in two BAM files, each of which is the output of a different aligner run
 * on the same set of reads. It compares the alignments in the two files and produces a summary of the differences.
 *
 */
@CommandLineProgramProperties(
        summary = CompareTwoAligners.USAGE_SUMMARY + CompareTwoAligners.USAGE_DETAILS,
        oneLineSummary = CompareTwoAligners.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class)
@DocumentedFeature
public class CompareTwoAligners extends CommandLineProgram {
    private static final String TOTAL_COUNTS_FIRST_BAM_FILE_NAME = "totalCountsFirstBamFile.tsv";
    private static final String TRUE_POSITIVE_COUNTS_FIRST_BAM_FILE_NAME = "truePositiveCountsFirstBamFile.tsv";
    private static final String TOTAL_COUNTS_SECOND_BAM_FILE_NAME = "totalCountsSecondBamFile.tsv";
    private static final String TRUE_POSITIVE_COUNTS_SECOND_BAM_FILE_NAME = "truePositiveCountsSecondBamFile.tsv";
    private static final String BINNED_CONCORDANCE_FILE_NAME = "binnedConcordance.tsv";

    static final String USAGE_SUMMARY = "Compare two input SAM/BAM/CRAM files.  ";
    static final String USAGE_DETAILS = "This tool compares the alignments in two BAM files, each of which is the output of a different aligner run on the same set of reads. " +
            "It produces a summary of the differences between the two sets of alignments." +
            "<h3>Usage example:</h3>" +
            "<pre>" +
            "gatk CompareTwoAligners \\<br />" +
            "     -firstBam first.bam \\<br />" +
            "     -secondBam second.bam \\<br />" +
            "     -O comparison.tsv" +
            "</pre>" +
            "\n";

    @Argument(doc = "The first BAM file to compare.")
    public File firstBam;

    @Argument(doc = "The second BAM file to compare.")
    public File secondBam;

    @Argument(doc = "The length of the bins to use for analyses.")
    public int binLength = 1000;

    @Argument(doc = "The length of the bins to use for for heat map.")
    public int binLengthHeatMap = 100000;

    @Argument(doc = "Require exact alignment position match for binned concordance.")
    public boolean requireExactAlignmentMatch = false;

    @Argument(
            doc = "Output directory to write comparison results to.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    @WorkflowOutput
    private File outputDir = null;

    // Typically in a GATK or Picard tool you have a logger, e.g.:
    private static final Logger logger = LogManager.getLogger();

    private static final ProgressLogger progressLogger = new ProgressLogger(logger, 1_000_000, "Processed", "records");

    private List<SimpleInterval> intervals;

    private EquallySpacedSearchableIntervalCollection searchableIntervalCollectionBinnedConcordance;

    private EquallySpacedSearchableIntervalCollection searchableIntervalCollectionSimpleConcordance;

    private SampleLocatableMetadata metadataFirstBam;
    private Multiset<SimpleInterval> intervalMultisetTotalCountFirstBam;
    private Multiset<SimpleInterval> intervalMultisetTruePositiveCountFirstBam;

    private SampleLocatableMetadata metadataSecondBam;
    private Multiset<SimpleInterval> intervalMultisetTotalCountSecondBam;
    private Multiset<SimpleInterval> intervalMultisetTruePositiveCountSecondBam;

    @ArgumentCollection
    final IntervalArgumentCollection intervalArgumentCollection = new RequiredIntervalArgumentCollection();

    @ArgumentCollection
    final ReferenceInputArgumentCollection reference = new RequiredReferenceInputArgumentCollection();

    // This class is used to store and read in a set of records with the same query name from the BAM file, i.e.
    // reads that come from the same template
    private static class TemplateReadCollection {
        private final String queryName;
        private final List<SAMRecord> records;
        private final SAMRecord firstInPairRecord;
        private final SAMRecord secondInPairRecord;

        /**
         * Note that this constructor changes the state of the iterator passed to it. Should only be used within this class.
         * @param iterator PeekableIterator of SAMRecords
         */
        public TemplateReadCollection(final PeekableIterator<SAMRecord> iterator) {
            if (!iterator.hasNext()) {
                // there must be more records in this iterator, throw exception
                throw new IllegalArgumentException("No more records in the iterator");
            }

            final SAMRecord firstRecord = iterator.next();
            this.queryName = firstRecord.getReadName();
            this.records = new ArrayList<>();
            this.records.add(firstRecord);
            progressLogger.record(firstRecord);
            while (iterator.hasNext()) {
                final SAMRecord record = iterator.peek();
                if (record.getReadName().equals(queryName)) {
                    records.add(iterator.next());
                    progressLogger.record(record);
                } else {
                    break;
                }
            }

            SAMRecord localFirstInPairRecord = null;
            SAMRecord localSecondInPairRecord = null;

            for (final SAMRecord record: records) {
                if (record.isSecondaryOrSupplementary()) {
                    continue;
                }

                if (record.getFirstOfPairFlag()) {
                    if (localFirstInPairRecord != null) {
                        throw new IllegalStateException("Multiple first in pair records found for query name: " + queryName);
                    }
                    localFirstInPairRecord = record;
                } else if (record.getSecondOfPairFlag()) {
                    if (localSecondInPairRecord != null) {
                        throw new IllegalStateException("Multiple second in pair records found for query name: " + queryName);
                    }
                    localSecondInPairRecord = record;
                }
            }

            if (localFirstInPairRecord == null || localSecondInPairRecord == null) {
                throw new IllegalStateException("Missing first or second in pair record for query name: " + queryName);
            }
            this.firstInPairRecord = localFirstInPairRecord;
            this.secondInPairRecord = localSecondInPairRecord;
        }

        // Getters
        public String getQueryName() {
            return queryName;
        }

        public List<SAMRecord> getRecords() {
            return records;
        }

        public SAMRecord getFirstInPairRecord() {
            return firstInPairRecord;
        }

        public SAMRecord getSecondInPairRecord() {
            return secondInPairRecord;
        }
    }

    private IntervalPair getIntervalPair(final SAMRecord firstRecord, final SAMRecord secondRecord) {
        SimpleInterval interval1;
        SimpleInterval interval2;

        if (requireExactAlignmentMatch && !firstRecord.getReadUnmappedFlag() && !secondRecord.getReadUnmappedFlag()) {
            if (!isSameAlignment(firstRecord, secondRecord)) {
                // If the alignments are not the same, skip this pair
                return new IntervalPair(null, null);
            }
        }

        if (firstRecord.getReadUnmappedFlag()) {
            // If either read is unmapped, skip this pair
            interval1 = null;
        } else {
            interval1 = searchableIntervalCollectionBinnedConcordance.getBinForPosition(firstRecord.getContig(), firstRecord.getStart());
        }
        if (secondRecord.getReadUnmappedFlag()) {
            // If either read is unmapped, skip this pair
            interval2 = null;
        } else {
            interval2 = searchableIntervalCollectionBinnedConcordance.getBinForPosition(secondRecord.getContig(), secondRecord.getStart());
        }

        return new IntervalPair(interval1, interval2);
    }

private boolean isSameAlignment(final SAMRecord firstRecord, final SAMRecord secondRecord) {
        // Check if the alignment start positions are the same
        if (!firstRecord.getContig().equals(secondRecord.getContig())) {
            return false;
        }

        if (firstRecord.getStart() != secondRecord.getStart()) {
            return false;
        }

        if (firstRecord.getEnd() != secondRecord.getEnd()) {
            return false;
        }

        return true;
    }

    private static class EquallySpacedSearchableIntervalCollection {

        private final Map<String, List<SimpleInterval>> contigToBins;
        private final List<SimpleInterval> bins; //we duplicate the bins for convenience
        private final int binLength;

        public EquallySpacedSearchableIntervalCollection(final SAMSequenceDictionary sequenceDictionary, final int binLength) {
            final List<SimpleInterval> referenceIntervals = IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
            this.contigToBins = new HashMap<>();
            this.binLength = binLength;
            this.bins = new ArrayList<>();
            for (final SimpleInterval interval : referenceIntervals) {
                for (int binStart = interval.getStart(); binStart <= interval.getEnd(); binStart += binLength) {
                    final int binEnd = FastMath.min(binStart + binLength - 1, interval.getEnd());
                    // add the bin to the appropriate contig list or create new list if it doesn't exist with the interval
                    final SimpleInterval bin = new SimpleInterval(interval.getContig(), binStart, binEnd);
                    contigToBins.computeIfAbsent(interval.getContig(), k -> new ArrayList<>()).add(bin);
                    this.bins.add(bin);
                }
            }
        }

        public List<SimpleInterval> getBinsForContig(final String contig) {
            // Check if contigToBins contains the specified contig, if not throw exception
            Utils.nonNull(contigToBins.get(contig),
                    "Contig " + contig + " not found in the bins collection. " +
                            "Please ensure that the reference sequence dictionary contains this contig.");
            return contigToBins.get(contig);
        }

        public List<SimpleInterval> getBins() {
            return bins;
        }

        public SimpleInterval getBinForPosition(final String contig, final int position) {
            Utils.nonNull(contig, "Contig cannot be null.");
            Utils.validateArg(position > 0, "Position must be greater than 0.");
            // Find the corresponding interval taking advantage of the fact that the bins are equally spaced
            List<SimpleInterval> bins = getBinsForContig(contig);
            final int contigEnd = bins.get(bins.size() - 1).getEnd();
            Utils.validateArg(position <= contigEnd,
                    "Position " + position + " exceeds the end of the contig " + contig + " which is at " + contigEnd + ".");
            return bins.get(
                    FastMath.min(
                            bins.size() - 1, // ensure we don't go out of bounds
                            FastMath.max(0, (position - 1) / binLength) // find the index of the bin for the position
                    )
            );
        }
    }

    private static class IntervalPair {
        private final SimpleInterval interval1;
        private final SimpleInterval interval2;

        public IntervalPair(SimpleInterval interval1, SimpleInterval interval2) {
            this.interval1 = interval1;
            this.interval2 = interval2;
        }

        public SimpleInterval getLeft() {
            return interval1;
        }

        public SimpleInterval getRight() {
            return interval2;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (obj == null || getClass() != obj.getClass()) return false;
            IntervalPair that = (IntervalPair) obj;
            return Objects.equals(interval1, that.interval1) && Objects.equals(interval2, that.interval2);
        }

        @Override
        public int hashCode() {
            return Objects.hash(interval1, interval2);
        }
    }

    private void updateIntervalCounts(
            final SAMRecord readA,
            final SAMRecord readB,
            final Multiset<SimpleInterval> totalCountA,
            final Multiset<SimpleInterval> truePositiveCountA,
            final Multiset<SimpleInterval> totalCountB,
            final Multiset<SimpleInterval> truePositiveCountB,
            final EquallySpacedSearchableIntervalCollection binProvider
    ) {
        // Update for the first BAM
        if (!readA.getReadUnmappedFlag()) {
            final SimpleInterval bin = binProvider.getBinForPosition(readA.getContig(), readA.getStart());
            totalCountA.add(bin);

            // If second read is mapped with the same alignment, count it as a true positive
            if (!readB.getReadUnmappedFlag() && isSameAlignment(readA, readB)) {
                truePositiveCountA.add(bin);
            }
        }

        // Update for the second BAM
        if (!readB.getReadUnmappedFlag()) {
            final SimpleInterval bin = binProvider.getBinForPosition(readB.getContig(), readB.getStart());
            totalCountB.add(bin);

            // If first read is unmapped but same alignment, count it as a true positive
            if (!readA.getReadUnmappedFlag() && isSameAlignment(readA, readB)) {
                truePositiveCountB.add(bin);
            }
        }
    }

    /**
     * Creates a SimpleCountCollection given metadata, a list of intervals, and a multiset of counts.
     */
    private SimpleCountCollection createSimpleCountCollection(
            final SampleLocatableMetadata metadata,
            final List<SimpleInterval> intervals,
            final Multiset<SimpleInterval> multiset
    ) {
        return new SimpleCountCollection(
                metadata,
                ImmutableList.copyOf(
                        intervals.stream()
                                .map(i -> new SimpleCount(i, multiset.count(i)))
                                .iterator()
                )
        );
    }

    @Override
    protected Object doWork() {
        final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().referenceSequence(reference.getReferencePath());

        final ReferenceDataSource referenceDataSource = ReferenceDataSource.of(reference.getReferencePath());
        intervals = intervalArgumentCollection.getIntervals(referenceDataSource.getSequenceDictionary());
        //intervalCachedOverlapDetector = new CachedOverlapDetector<>(intervals);
        searchableIntervalCollectionBinnedConcordance = new EquallySpacedSearchableIntervalCollection(
                referenceDataSource.getSequenceDictionary(),
                binLengthHeatMap
        );

        searchableIntervalCollectionSimpleConcordance = new EquallySpacedSearchableIntervalCollection(
                referenceDataSource.getSequenceDictionary(),
                binLength
        );

        final Map<IntervalPair, Integer> intervalPairCounts = new HashMap<>();

        try (final SamReader firstReader = samReaderFactory.open(firstBam);
             final SamReader secondReader = samReaderFactory.open(secondBam))
        {
            final PeekableIterator<SAMRecord> firstPeekableIterator = new PeekableIterator<>(firstReader.iterator());
            final PeekableIterator<SAMRecord> secondPeekableIterator = new PeekableIterator<>(secondReader.iterator());
            metadataFirstBam = MetadataUtils.fromHeader(firstReader.getFileHeader(), Metadata.Type.SAMPLE_LOCATABLE);
            intervalMultisetTotalCountFirstBam = HashMultiset.create(searchableIntervalCollectionSimpleConcordance.getBins().size());
            intervalMultisetTruePositiveCountFirstBam = HashMultiset.create(searchableIntervalCollectionSimpleConcordance.getBins().size());
            metadataSecondBam = MetadataUtils.fromHeader(secondReader.getFileHeader(), Metadata.Type.SAMPLE_LOCATABLE);
            intervalMultisetTotalCountSecondBam = HashMultiset.create(searchableIntervalCollectionSimpleConcordance.getBins().size());
            intervalMultisetTruePositiveCountSecondBam = HashMultiset.create(searchableIntervalCollectionSimpleConcordance.getBins().size());
            // Read in set of records with matching query names, i.e. keep reading in records until the query names changes
            while (firstPeekableIterator.hasNext() && secondPeekableIterator.hasNext()) {
                final TemplateReadCollection firstTemplateReads = new TemplateReadCollection(firstPeekableIterator);
                final TemplateReadCollection secondTemplateReads = new TemplateReadCollection(secondPeekableIterator);

                // Query names must match
                if (!firstTemplateReads.getQueryName().equals(secondTemplateReads.getQueryName())) {
                    throw new IllegalStateException("Query names do not match: " + firstTemplateReads.getQueryName() + " vs " + secondTemplateReads.getQueryName());
                }

                final IntervalPair firstInPairResult = getIntervalPair(firstTemplateReads.getFirstInPairRecord(), secondTemplateReads.getFirstInPairRecord());
                intervalPairCounts.put(firstInPairResult, intervalPairCounts.getOrDefault(firstInPairResult, 0) + 1);
                final IntervalPair secondInPairResult = getIntervalPair(firstTemplateReads.getSecondInPairRecord(), secondTemplateReads.getSecondInPairRecord());
                intervalPairCounts.put(secondInPairResult, intervalPairCounts.getOrDefault(secondInPairResult, 0) + 1);


                // Update interval counts for both BAMs (first/second in pair)
                updateIntervalCounts(
                        firstTemplateReads.getFirstInPairRecord(),
                        secondTemplateReads.getFirstInPairRecord(),
                        intervalMultisetTotalCountFirstBam,
                        intervalMultisetTruePositiveCountFirstBam,
                        intervalMultisetTotalCountSecondBam,
                        intervalMultisetTruePositiveCountSecondBam,
                        searchableIntervalCollectionSimpleConcordance
                );
                updateIntervalCounts(
                        firstTemplateReads.getSecondInPairRecord(),
                        secondTemplateReads.getSecondInPairRecord(),
                        intervalMultisetTotalCountFirstBam,
                        intervalMultisetTruePositiveCountFirstBam,
                        intervalMultisetTotalCountSecondBam,
                        intervalMultisetTruePositiveCountSecondBam,
                        searchableIntervalCollectionSimpleConcordance
                );

            }
            // Check that both files had the same number of records
            if (firstPeekableIterator.hasNext() || secondPeekableIterator.hasNext()) {
                throw new IllegalStateException("The two BAM files do not contain the same number of reads.");
            }

        } catch (IOException e) {
            throw new RuntimeIOException("Error opening file", e);
        }

        // Create a multiset of intervals and their counts
        final SimpleCountCollection totalReadCountsFirstBam =
                createSimpleCountCollection(metadataFirstBam, searchableIntervalCollectionSimpleConcordance.getBins(), intervalMultisetTotalCountFirstBam);
        final SimpleCountCollection truePositiveReadCountsFirstBam =
                createSimpleCountCollection(metadataFirstBam, searchableIntervalCollectionSimpleConcordance.getBins(), intervalMultisetTruePositiveCountFirstBam);
        final SimpleCountCollection totalReadCountsSecondBam =
                createSimpleCountCollection(metadataSecondBam, searchableIntervalCollectionSimpleConcordance.getBins(), intervalMultisetTotalCountSecondBam);
        final SimpleCountCollection truePositiveReadCountsSecondBam =
                createSimpleCountCollection(metadataSecondBam, searchableIntervalCollectionSimpleConcordance.getBins(), intervalMultisetTruePositiveCountSecondBam);

        // Write the results to the output directory
        Utils.nonNull(outputDir, "Output directory must be specified.");
        totalReadCountsFirstBam.write(new File(outputDir, TOTAL_COUNTS_FIRST_BAM_FILE_NAME));
        truePositiveReadCountsFirstBam.write(new File(outputDir, TRUE_POSITIVE_COUNTS_FIRST_BAM_FILE_NAME));
        totalReadCountsSecondBam.write(new File(outputDir, TOTAL_COUNTS_SECOND_BAM_FILE_NAME));
        truePositiveReadCountsSecondBam.write(new File(outputDir, TRUE_POSITIVE_COUNTS_SECOND_BAM_FILE_NAME));
        try (final FileWriter writer = new FileWriter(new File(outputDir, BINNED_CONCORDANCE_FILE_NAME))) {
            writer.write("Interval1\tInterval2\tCount\n");
            for (Map.Entry<IntervalPair, Integer> entry : intervalPairCounts.entrySet()) {
                IntervalPair key = entry.getKey();
                writer.write(String.format("%s\t%s\t%d%n", key.getLeft(), key.getRight(), entry.getValue()));
            }
        } catch (IOException e) {
            throw new RuntimeIOException("Error writing to output file", e);
        }
        return 0;
    }
}