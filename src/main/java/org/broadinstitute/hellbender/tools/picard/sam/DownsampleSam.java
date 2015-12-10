package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 * Class to randomly downsample a BAM file while respecting that we should either get rid
 * of both ends of a pair or neither end of the pair!
 */
@CommandLineProgramProperties(
        summary = "Randomly down-sample a SAM/BAM/CRAM file to retain " +
                "a random subset of the reads. Mate-pairs are either both kept or both discarded. Reads marked as not primary " +
                "alignments are all discarded. Each read is given a probability P of being retained - results with the exact " +
                "same input in the same order and with the same value for RANDOM_SEED will produce the same results.",
        oneLineSummary = "Down-sample a SAM/BAM file to retain a random subset of the reads",
        programGroup = ReadProgramGroup.class
)
public final class DownsampleSam extends PicardCommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, 
	    doc = "The input SAM/BAM/CRAM file to downsample.")
    public File INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, 
	    doc = "The output, downsampled, SAM/BAM/CRAM file to write.")
    public File OUTPUT;

    @Argument(fullName = "random_seed",
              shortName = "rand",
              doc = "Random seed to use if reproducibilty is desired.  " +
            "Setting to null will cause multiple invocations to produce different results.")
    public Long RANDOM_SEED = 1L;

    @Argument(shortName = "P", doc = "The probability of keeping any individual read, between 0 and 1.")
    public double PROBABILITY = 1;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final Random r = RANDOM_SEED == null ? new Random() : new Random(RANDOM_SEED);
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        long total = 0;
        long kept = 0;

        try (final SAMFileWriter out = createSAMWriter(OUTPUT, REFERENCE_SEQUENCE, in.getFileHeader(), true)) {
            final Map<String, Boolean> decisions = new HashMap<>();

            final ProgressLogger progress = new ProgressLogger(logger, (int) 1e7, "Read");

            for (final SAMRecord rec : in) {
                if (rec.isSecondaryOrSupplementary()) continue;
                ++total;

                final String key = rec.getReadName();
                final Boolean previous = decisions.remove(key);
                final boolean keeper;

                if (previous == null) {
                    keeper = r.nextDouble() <= PROBABILITY;
                    if (rec.getReadPairedFlag()) decisions.put(key, keeper);
                } else {
                    keeper = previous;
                }

                if (keeper) {
                    out.addAlignment(rec);
                    ++kept;
                }

                progress.record(rec);
            }
        }
        finally {
            CloserUtil.close(in);
        }
        logger.info("Finished! Kept " + kept + " out of " + total + " reads.");

        return null;
    }
}
