package org.broadinstitute.hellbender.tools.picard.fastq;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.MiscProgramGroup;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

/**
 * Simple command line program that allows sub-sequences represented by an interval
 * list to be extracted from a reference sequence file.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = "Extracts one or more intervals described in an interval_list file " +
                "from a given reference sequence and writes them out in FASTA format. Requires a fasta index " +
                "file to be present.",
        usageShort = "Extracts intervals from a reference sequence, writing them to a FASTA file",
        programGroup = MiscProgramGroup.class
)
public class ExtractSequences extends CommandLineProgram {

    @Argument(doc="Interval list describing intervals to be extracted from the reference sequence.")
    public File INTERVAL_LIST;

    @Argument(shortName=StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc="Reference sequence file.")
    public File REFERENCE_SEQUENCE;

    @Argument(shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Output fasta file.")
    public File OUTPUT;

    @Argument(doc="Maximum line length for sequence data.")
    public int LINE_LENGTH = 80;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INTERVAL_LIST);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IOUtil.assertFileIsWritable(OUTPUT);

        final IntervalList intervals = IntervalList.fromFile(INTERVAL_LIST);
        final ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
        SequenceUtil.assertSequenceDictionariesEqual(intervals.getHeader().getSequenceDictionary(), ref.getSequenceDictionary());

        final BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT);

        for (final Interval interval : intervals) {
            final ReferenceSequence seq = ref.getSubsequenceAt(interval.getContig(), interval.getStart(), interval.getEnd());
            final byte[] bases = seq.getBases();
            if (interval.isNegativeStrand()) SequenceUtil.reverseComplement(bases);

            try {
                out.write(">");
                out.write(interval.getName());
                out.write("\n");

                for (int i=0; i<bases.length; ++i) {
                    if (i > 0 && i % LINE_LENGTH == 0) out.write("\n");
                    out.write(bases[i]);
                }

                out.write("\n");
            }
            catch (IOException ioe) {
                throw new RuntimeIOException("Error writing to file " + OUTPUT.getAbsolutePath(), ioe);
            }
        }

        CloserUtil.close(out);
        return null;
    }
}
