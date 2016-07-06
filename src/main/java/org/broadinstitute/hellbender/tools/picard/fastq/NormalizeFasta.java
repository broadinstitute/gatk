package org.broadinstitute.hellbender.tools.picard.fastq;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.FastaProgramGroup;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

/**
 * Little program to "normalize" a fasta file to ensure that all line of sequence are the
 * same length, and are a reasonable length!
 */
@CommandLineProgramProperties(
        summary = "Takes any file that conforms to the fasta format and " +
                "normalizes it so that all lines of sequence except the last line per named sequence " +
                "are of the same length.",
        oneLineSummary = "Normalizes lines of sequence in a fasta file to be of the same length",
        programGroup = FastaProgramGroup.class
)
public final class NormalizeFasta extends CommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.INPUT_SHORT_NAME, doc="The input fasta file to normalize.")
    public File INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="The output fasta file to write.")
    public File OUTPUT;

    @Argument(doc="The line length to be used for the output fasta file.")
    public int LINE_LENGTH=100;

    @Argument(doc="Truncate sequence names at first whitespace.")
    public boolean TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE=false;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        Utils.validateArg(!INPUT.getAbsoluteFile().equals(OUTPUT.getAbsoluteFile()), "Input and output cannot be the same file.");

        try (final ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(INPUT, TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE);
             final BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT)) {

            ReferenceSequence seq = null;
            while ((seq = ref.nextSequence()) != null) {
                final String name  = seq.getName();
                final byte[] bases = seq.getBases();

                try {
                    out.write(">");
                    out.write(name);
                    out.newLine();

                    if (bases.length == 0) {
                        logger.warn("Sequence " + name + " contains 0 bases.");
                    }
                    else {
                        for (int i=0; i<bases.length; ++i) {
                            if (i > 0 && i % LINE_LENGTH == 0) out.write("\n");
                            out.write(bases[i]);
                        }

                        out.write("\n");
                    }
                }
                catch (IOException ioe) {
                    throw new RuntimeIOException("Error writing to file " + OUTPUT.getAbsolutePath(), ioe);

                }
            }

        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }

        return null;
    }
}
