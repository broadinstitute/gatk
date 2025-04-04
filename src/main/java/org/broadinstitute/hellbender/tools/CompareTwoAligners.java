package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
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

    @Argument(
            doc = "Output file to write comparison results to.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    @WorkflowOutput
    private File outputComparisonFile = null;

    // Typically in a GATK or Picard tool you have a logger, e.g.:
    private static final Logger logger = LogManager.getLogger();

    // Create a ProgressLogger that reports every 1,000,000 reads (you can pick any interval).
    private static final ProgressLogger progressLogger = new ProgressLogger(logger, 1_000_000, "Processed", "records");


    private List<SimpleInterval> intervals;

    private CachedOverlapDetector<SimpleInterval> intervalCachedOverlapDetector;

    private EquallySpacedSearchableIntervalCollection searchableIntervalCollection;

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

    private Pair<SimpleInterval, SimpleInterval> getIntervalPair(final SAMRecord firstRecord, final SAMRecord secondRecord) {
        SimpleInterval interval1;
        SimpleInterval interval2;

        if (firstRecord.getReadUnmappedFlag()) {
            // If either read is unmapped, skip this pair
            interval1 = null;
        } else {
            interval1 = searchableIntervalCollection.getBinForPosition(firstRecord.getContig(), firstRecord.getStart());
        }
        if (secondRecord.getReadUnmappedFlag()) {
            // If either read is unmapped, skip this pair
            interval2 = null;
        } else {
            interval2 = searchableIntervalCollection.getBinForPosition(secondRecord.getContig(), secondRecord.getStart());
        }

        return new ImmutablePair<>(interval1, interval2);
    }

    private static class EquallySpacedSearchableIntervalCollection {

        private final Map<String, List<SimpleInterval>> contigToBins;
        private final int binLength;

        public EquallySpacedSearchableIntervalCollection(final SAMSequenceDictionary sequenceDictionary, final int binLength) {
            final List<SimpleInterval> referenceIntervals = IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
            this.contigToBins = new HashMap<>();
            this.binLength = binLength;
            for (final SimpleInterval interval : referenceIntervals) {
                for (int binStart = interval.getStart(); binStart <= interval.getEnd(); binStart += binLength) {
                    final int binEnd = FastMath.min(binStart + binLength - 1, interval.getEnd());
                    // add the bin to the appropriate contig list or create new list if it doesn't exist with the interval
                    contigToBins.computeIfAbsent(interval.getContig(), k -> new ArrayList<>())
                            .add(new SimpleInterval(interval.getContig(), binStart, binEnd));
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

    @Override
    protected Object doWork() {
        final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().referenceSequence(reference.getReferencePath());

        final ReferenceDataSource referenceDataSource = ReferenceDataSource.of(reference.getReferencePath());
        intervals = intervalArgumentCollection.getIntervals(referenceDataSource.getSequenceDictionary());
        intervalCachedOverlapDetector = new CachedOverlapDetector<>(intervals);
        searchableIntervalCollection = new EquallySpacedSearchableIntervalCollection(
                referenceDataSource.getSequenceDictionary(),
                10000
        );

        final Map<Pair<SimpleInterval, SimpleInterval>, Integer> intervalPairCounts = new HashMap<>();
        try (final SamReader firstReader = samReaderFactory.open(firstBam);
             final SamReader secondReader = samReaderFactory.open(secondBam))
        {
            final PeekableIterator<SAMRecord> firstPeekableIterator = new PeekableIterator<>(firstReader.iterator());
            final PeekableIterator<SAMRecord> secondPeekableIterator = new PeekableIterator<>(secondReader.iterator());

            // Read in set of records with matching query names, i.e. keep reading in records until the query names changes

            while (firstPeekableIterator.hasNext() && secondPeekableIterator.hasNext()) {
                final TemplateReadCollection firstTemplateReads = new TemplateReadCollection(firstPeekableIterator);
                final TemplateReadCollection secondTemplateReads = new TemplateReadCollection(secondPeekableIterator);

                // Query names must match
                if (!firstTemplateReads.getQueryName().equals(secondTemplateReads.getQueryName())) {
                    throw new IllegalStateException("Query names do not match: " + firstTemplateReads.getQueryName() + " vs " + secondTemplateReads.getQueryName());
                }

                Pair<SimpleInterval, SimpleInterval> firstInPairResult = getIntervalPair(firstTemplateReads.getFirstInPairRecord(), secondTemplateReads.getFirstInPairRecord());
                intervalPairCounts.put(firstInPairResult, intervalPairCounts.getOrDefault(firstInPairResult, 0) + 1);
                Pair<SimpleInterval, SimpleInterval> secondInPairResult = getIntervalPair(firstTemplateReads.getSecondInPairRecord(), secondTemplateReads.getSecondInPairRecord());
                intervalPairCounts.put(secondInPairResult, intervalPairCounts.getOrDefault(secondInPairResult, 0) + 1);
                // Compare alignments start position and CIGAR strings.
//                if (firstRecord.getAlignmentStart() != secondRecord.getAlignmentStart() ||
//                        !firstRecord.getCigarString().equals(secondRecord.getCigarString())) {
//                    System.out.println("Alignments do not match for read: " + firstRecord.getReadName());
//                }
            }
            // Check that both files had the same number of records
            if (firstPeekableIterator.hasNext() || secondPeekableIterator.hasNext()) {
                throw new IllegalStateException("The two BAM files do not contain the same number of reads.");
            }

        } catch (IOException e) {
            throw new RuntimeIOException("Error opening file", e);
        }
        // Write the results to the output file
        Utils.nonNull(outputComparisonFile, "Output file must be specified.");
        try (final FileWriter writer = new FileWriter(outputComparisonFile)) {
            writer.write("Interval1\tInterval2\tCount\n");
            for (Map.Entry<Pair<SimpleInterval, SimpleInterval>, Integer> entry : intervalPairCounts.entrySet()) {
                Pair<SimpleInterval, SimpleInterval> key = entry.getKey();
                writer.write(String.format("%s\t%s\t%d%n", key.getLeft(), key.getRight(), entry.getValue()));
            }
        } catch (IOException e) {
            throw new RuntimeIOException("Error writing to output file", e);
        }
        return 0;
    }
}