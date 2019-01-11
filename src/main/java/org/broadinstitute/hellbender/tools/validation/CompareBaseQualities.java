package org.broadinstitute.hellbender.tools.validation;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SecondaryOrSupplementarySkippingIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.PositionalArguments;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

import static org.broadinstitute.hellbender.transformers.BQSRReadTransformer.constructStaticQuantizedMapping;

/**
 * Compares the base qualities of two SAM/BAM/CRAM files. The reads in the two files must have exactly the same
 * names and appear in the same order. The two files are each specified without a flag.
 * The tool summarizes the results in a text file.
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk CompareBaseQualities \
 *   input_1.bam \
 *   input_2.bam \
 *   -O diff.txt
 * </pre>
 *
 * <p><b>An example result with identical base qualities between two inputs</b></p>
 * <pre>
 * -----------CompareMatrix summary------------
 * all 10000 quality scores are the same
 *
 * ---------CompareMatrix full matrix (non-zero entries) ----------
 * QRead1	QRead2	diff	count
 * 40	40	0	10000
 * -----------CompareMatrix-binned summary------------
 * all 10000 quality scores are the same
 *
 * ---------CompareMatrix-binned full matrix (non-zero entries) ----------
 * QRead1	QRead2	diff	count
 * 40	40	0	10000
 * </pre>
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Compares the base qualities of two SAM/BAM/CRAM files. The reads in the two files must have " +
                "exactly the same names and appear in the same order.",
        oneLineSummary = "Compares the base qualities of two SAM/BAM/CRAM files",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
public final class CompareBaseQualities extends PicardCommandLineProgram {
    @PositionalArguments(minElements = 2, maxElements = 2)
    public List<File> samFiles;

    public static final String THROW_ON_DIFF_LONG_NAME = "throw-on-diff";

    public static final String STATIC_QUANTIZED_QUALS_LONG_NAME = "static-quantized-quals";

    public static final String ROUND_DOWN_QUANTIZED_LONG_NAME = "round-down-quantized";

    @Argument(
            doc = "Summary output file.",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true
    )
    protected String outputFilename = null;

    @Argument(
            doc = "Throw exception on difference.",
            fullName = THROW_ON_DIFF_LONG_NAME,
            optional = true
    )
    protected boolean throwOnDiff = false;

    /**
     * Return value is 0 if the two files have identical base qualities and non-zero otherwise.
     * Use static quantized quality scores to a given number of levels.
     */
    @Advanced
    @Argument(
            doc = "Use static quantized quality scores to a given number of levels " +
                    "(with --" + StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME+ ")",
            fullName = STATIC_QUANTIZED_QUALS_LONG_NAME,
            optional = true
    )
    public List<Integer> staticQuantizationQuals = new ArrayList<>();

    /**
     * Round down quantized only works with the {@link #STATIC_QUANTIZED_QUALS_LONG_NAME} option.  If enabled,
     * rounding is done in probability space to the nearest bin.  Otherwise, the value is rounded to the
     * nearest bin that is smaller than the current bin.
     *
     * Note: this option only works when {@link #STATIC_QUANTIZED_QUALS_LONG_NAME} is set.
     */
    @Advanced
    @Argument(
            doc = "Round down quality scores to nearest quantized value.",
            fullName = ROUND_DOWN_QUANTIZED_LONG_NAME,
            optional = true
    )
    public boolean roundDown = false;

    private byte[] staticQuantizedMapping;

    @Override
    protected Object doWork() {
        if (roundDown && (staticQuantizationQuals == null || staticQuantizationQuals.isEmpty())){
            throw new CommandLineException.BadArgumentValue(
                    ROUND_DOWN_QUANTIZED_LONG_NAME, "true",
                    "This option can only be used if " + STATIC_QUANTIZED_QUALS_LONG_NAME + " is also used.");
        }
        staticQuantizedMapping = constructStaticQuantizedMapping(staticQuantizationQuals, roundDown);

        IOUtil.assertFileIsReadable(samFiles.get(0));
        IOUtil.assertFileIsReadable(samFiles.get(1));

        final SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(VALIDATION_STRINGENCY);
        final SamReader reader1 = factory.referenceSequence(REFERENCE_SEQUENCE).open(samFiles.get(0));
        final SamReader reader2 = factory.referenceSequence(REFERENCE_SEQUENCE).open(samFiles.get(1));

        final SecondaryOrSupplementarySkippingIterator it1 =
                new SecondaryOrSupplementarySkippingIterator(reader1.iterator());
        final SecondaryOrSupplementarySkippingIterator it2 =
                new SecondaryOrSupplementarySkippingIterator(reader2.iterator());

        final CompareMatrix finalMatrix = new CompareMatrix(staticQuantizedMapping);

        final ProgressMeter progressMeter = new ProgressMeter(1.0);
        progressMeter.start();

        while (it1.hasCurrent() && it2.hasCurrent()) {
            final SAMRecord read1= it1.getCurrent();
            final SAMRecord read2= it2.getCurrent();

            progressMeter.update(read1);

            if (!read1.getReadName().equals(read2.getReadName())){
                throw new UserException.BadInput("files do not have the same exact order of reads:" +
                        read1.getReadName() + " vs " + read2.getReadName());
            }

            finalMatrix.add(read1.getBaseQualities(), read2.getBaseQualities());
            it1.advance();
            it2.advance();
        }
        progressMeter.stop();

        if (it1.hasCurrent() || it2.hasCurrent()){
            throw new UserException.BadInput("files do not have the same exact number of reads");
        }
        CloserUtil.close(reader1);
        CloserUtil.close(reader2);

        finalMatrix.printOutResults(outputFilename);

        if (throwOnDiff && finalMatrix.hasNonDiagonalElements()) {
            throw new UserException("Quality scores from the two BAMs do not match");
        }

        return finalMatrix.hasNonDiagonalElements() ? 1 : 0;
    }
}
