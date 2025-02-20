package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.*;
import htsjdk.samtools.util.RuntimeIOException;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
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
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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

    private List<SimpleInterval> intervals;

    private CachedOverlapDetector<SimpleInterval> intervalCachedOverlapDetector;

    @ArgumentCollection
    final IntervalArgumentCollection intervalArgumentCollection = new RequiredIntervalArgumentCollection();

    @ArgumentCollection
    final ReferenceInputArgumentCollection reference = new RequiredReferenceInputArgumentCollection();

    @Override
    protected Object doWork() {
        final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().referenceSequence(reference.getReferencePath());

        final ReferenceDataSource referenceDataSource = ReferenceDataSource.of(reference.getReferencePath());
        intervals = intervalArgumentCollection.getIntervals(referenceDataSource.getSequenceDictionary());
        intervalCachedOverlapDetector = new CachedOverlapDetector<>(intervals);

        final Map<Pair<SimpleInterval, SimpleInterval>, Integer> intervalPairCounts = new HashMap<>();
        try (final SamReader firstReader = samReaderFactory.open(firstBam);
             final SamReader secondReader = samReaderFactory.open(secondBam))
        {
            final SAMRecordIterator firstIterator = firstReader.iterator();
            final SAMRecordIterator secondIterator = secondReader.iterator();

            while (firstIterator.hasNext() && secondIterator.hasNext()) {
                final SAMRecord firstRecord = firstIterator.next();
                final SAMRecord secondRecord = secondIterator.next();

                // Query names must match
                if (!firstRecord.getReadName().equals(secondRecord.getReadName())) {
                    throw new IllegalStateException("Query names do not match: " + firstRecord.getReadName() + " vs " + secondRecord.getReadName());
                }

                SimpleInterval interval1 = intervalCachedOverlapDetector.getOverlap(new SimpleInterval(firstRecord.getContig(), firstRecord.getStart(), firstRecord.getStart()));
                SimpleInterval interval2 = intervalCachedOverlapDetector.getOverlap(new SimpleInterval(secondRecord.getContig(), secondRecord.getStart(), secondRecord.getStart()));

                Pair<SimpleInterval, SimpleInterval> key = new ImmutablePair<>(interval1, interval2);
                intervalPairCounts.put(key, intervalPairCounts.getOrDefault(key, 0) + 1);
                // Compare alignments start position and CIGAR strings.
                if (firstRecord.getAlignmentStart() != secondRecord.getAlignmentStart() ||
                        !firstRecord.getCigarString().equals(secondRecord.getCigarString())) {
                    System.out.println("Alignments do not match for read: " + firstRecord.getReadName());
                }
            }
            // Check that both files had the same number of records
            if (firstIterator.hasNext() || secondIterator.hasNext()) {
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