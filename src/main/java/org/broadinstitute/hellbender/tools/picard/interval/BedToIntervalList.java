package org.broadinstitute.hellbender.tools.picard.interval;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.IntervalProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;
import java.io.IOException;

/**
 * @author nhomer
 */
@CommandLineProgramProperties(
        summary = "Converts a BED file to an Picard Interval List",
        oneLineSummary = "Converts a BED file to an Picard Interval List",
        programGroup = IntervalProgramGroup.class
)
public final class BedToIntervalList extends PicardCommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "The input BED file")
    public File INPUT;

    @Argument(shortName = StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, doc = "The sequence dictionary")
    public File SEQUENCE_DICTIONARY;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "The output Picard Interval List")
    public File OUTPUT;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
        IOUtil.assertFileIsWritable(OUTPUT);
        try {
            final SAMFileHeader header = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).getFileHeader(SEQUENCE_DICTIONARY);
            final IntervalList intervalList = new IntervalList(header);

            /**
             * NB: BED is zero-based, but a BEDCodec by default (since it is returns tribble Features) has an offset of one,
             * so it returns 1-based starts.  Ugh.  Set to zero.
             */
            final FeatureReader<BEDFeature> bedReader = AbstractFeatureReader.getFeatureReader(INPUT.getAbsolutePath(), new BEDCodec(BEDCodec.StartOffset.ZERO), false);
            final CloseableTribbleIterator<BEDFeature> iterator = bedReader.iterator();
            final ProgressLogger progressLogger = new ProgressLogger(logger, (int) 1e6);

            while (iterator.hasNext()) {
                final BEDFeature bedFeature = iterator.next();
                final String sequenceName = bedFeature.getContig();
                /**
                 * NB: BED is zero-based, so we need to add one here to make it one-based.  Please observe we set the start
                 * offset to zero when creating the BEDCodec.
                 */
                final int start = bedFeature.getStart() + 1;
                /**
                 * NB: BED is 0-based OPEN (which, for the end is equivalent to 1-based closed).
                 */
                final int end = bedFeature.getEnd();
                // NB: do not use an empty name within an interval
                String name = bedFeature.getName();
                if (name.isEmpty()) name = null;

                final SAMSequenceRecord sequenceRecord = header.getSequenceDictionary().getSequence(sequenceName);

                // Do some validation
                if (null == sequenceRecord) {
                    throw new GATKException(String.format("Sequence '%s' was not found in the sequence dictionary", sequenceName));
                } else if (start < 1) {
                    throw new GATKException(String.format("Start on sequence '%s' was less than one: %d", sequenceName, start));
                } else if (sequenceRecord.getSequenceLength() < start) {
                    throw new GATKException(String.format("Start on sequence '%s' was past the end: %d < %d", sequenceName, sequenceRecord.getSequenceLength(), start));
                } else if (end < 1) {
                    throw new GATKException(String.format("End on sequence '%s' was less than one: %d", sequenceName, end));
                } else if (sequenceRecord.getSequenceLength() < end) {
                    throw new GATKException(String.format("End on sequence '%s' was past the end: %d < %d", sequenceName, sequenceRecord.getSequenceLength(), end));
                } else if (end < start - 1) {
                    throw new GATKException(String.format("On sequence '%s', end < start-1: %d <= %d", sequenceName, end, start));
                }

                final Interval interval = new Interval(sequenceName, start, end, bedFeature.getStrand() == Strand.POSITIVE, name);
                intervalList.add(interval);

                progressLogger.record(sequenceName, start);
            }
            CloserUtil.close(bedReader);

            // Sort and write the output
            intervalList.uniqued().write(OUTPUT);

        } catch (final IOException e) {
            throw new RuntimeException(e);
        }

        return null;
    }
}
