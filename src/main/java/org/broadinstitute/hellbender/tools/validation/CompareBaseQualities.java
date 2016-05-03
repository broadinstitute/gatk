package org.broadinstitute.hellbender.tools.validation;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SecondaryOrSupplementarySkippingIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

import static org.broadinstitute.hellbender.transformers.BQSRReadTransformer.constructStaticQuantizedMapping;

@CommandLineProgramProperties(
        summary = "USAGE: CompareBaseQualities <SAMFile1> <SAMFile2>\n" +
                "Compares the base qualities of two input SAM/BAM/CRAM files. The files must be sorted exactly the same name.",
        oneLineSummary = "Compares base qualities of two input SAM/BAM/CRAM files",
        programGroup = ReadProgramGroup.class
)
public final class CompareBaseQualities extends PicardCommandLineProgram {
    @PositionalArguments(minElements = 2, maxElements = 2)
    public List<File> samFiles;

    @Argument(doc="summary output file", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = true)
    protected String outputFilename = null;

    @Argument(doc="throw error on diff", shortName = "cd", fullName = "throwOnDiff", optional = true)
    protected boolean throwOnDiff = false;

    /**
     * Return value is 0 if the two files have identical base qualities and non-zero otherwise.
     * Use static quantized quality scores to a given number of levels.
     */
    @Advanced
    @Argument(fullName="static_quantized_quals", shortName = "SQQ", doc = "Use static quantized quality scores to a given number of levels (with -"+ StandardArgumentDefinitions.BQSR_TABLE_SHORT_NAME+ ")", optional=true)
    public List<Integer> staticQuantizationQuals = new ArrayList<>();

    /**
     * Round down quantized only works with the static_quantized_quals option.  When roundDown = false, rounding is done in
     * probability space to the nearest bin.  When roundDown = true, the value is rounded to the nearest bin
     * that is smaller than the current bin.
     * Note: it only works when static_quantized_quals is set.
     */
    @Advanced
    @Argument(fullName="round_down_quantized", shortName = "RDQ", doc = "Round quals down to nearest quantized qual", optional=true)
    public boolean roundDown = false;

    private byte[] staticQuantizedMapping;

    @Override
    protected Object doWork() {
        if (roundDown && (staticQuantizationQuals == null || staticQuantizationQuals.isEmpty())){
            throw new UserException.BadArgumentValue("round_down_quantized", "true", "This option can only be used if static_quantized_quals is also used");
        }
        staticQuantizedMapping = constructStaticQuantizedMapping(staticQuantizationQuals, roundDown);

        IOUtil.assertFileIsReadable(samFiles.get(0));
        IOUtil.assertFileIsReadable(samFiles.get(1));

        final SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(VALIDATION_STRINGENCY);
        final SamReader reader1 = factory.referenceSequence(REFERENCE_SEQUENCE).open(samFiles.get(0));
        final SamReader reader2 = factory.referenceSequence(REFERENCE_SEQUENCE).open(samFiles.get(1));

        final SecondaryOrSupplementarySkippingIterator it1 = new SecondaryOrSupplementarySkippingIterator(reader1.iterator());
        final SecondaryOrSupplementarySkippingIterator it2 = new SecondaryOrSupplementarySkippingIterator(reader2.iterator());

        final CompareMatrix finalMatrix = new CompareMatrix(staticQuantizedMapping);

        final ProgressMeter progressMeter = new ProgressMeter(1.0);
        progressMeter.start();

        while (it1.hasCurrent() && it2.hasCurrent()) {
            final SAMRecord read1= it1.getCurrent();
            final SAMRecord read2= it2.getCurrent();

            progressMeter.update(read1);

            if (!read1.getReadName().equals(read2.getReadName())){
                throw new UserException.BadInput("files do not have the same exact order of reads:" + read1.getReadName() + " vs " + read2.getReadName());
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
