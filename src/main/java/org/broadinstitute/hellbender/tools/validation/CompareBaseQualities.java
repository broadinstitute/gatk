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

import java.io.File;
import java.util.List;

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

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(samFiles.get(0));
        IOUtil.assertFileIsReadable(samFiles.get(1));

        final SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(VALIDATION_STRINGENCY);
        final SamReader reader1 = factory.referenceSequence(REFERENCE_SEQUENCE).open(samFiles.get(0));
        final SamReader reader2 = factory.referenceSequence(REFERENCE_SEQUENCE).open(samFiles.get(1));

        final SecondaryOrSupplementarySkippingIterator it1 = new SecondaryOrSupplementarySkippingIterator(reader1.iterator());
        final SecondaryOrSupplementarySkippingIterator it2 = new SecondaryOrSupplementarySkippingIterator(reader2.iterator());

        final CompareMatrix finalMatrix = new CompareMatrix();

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

        finalMatrix.printOutput(outputFilename);

        if (throwOnDiff && finalMatrix.hasNonDiagonalElements()) {
            throw new UserException("Quality scores from the two BAMs do not match");
        }

        return 0;
    }

}
