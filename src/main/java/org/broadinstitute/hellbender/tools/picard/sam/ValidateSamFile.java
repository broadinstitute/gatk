package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * Command line program wrapping SamFileValidator.
 *
 * @author Doug Voet
 */
@CommandLineProgramProperties(
        summary = "Read a SAM/BAM/CRAM file and report on its validity.",
        oneLineSummary = "Validates a SAM/BAM/CRAM file",
        programGroup = ReadProgramGroup.class
)
public final class ValidateSamFile extends PicardCommandLineProgram {

    public enum Mode {VERBOSE, SUMMARY}

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "Input SAM/BAM/CRAM file")
    public File INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file or standard out if missing",
            optional = true)
    public File OUTPUT;

    @Argument(shortName = "M",
            doc = "Mode of output")
    public Mode MODE = Mode.VERBOSE;

    @Argument(doc = "List of validation error types to ignore.", optional = true)
    public List<SAMValidationError.Type> IGNORE = new ArrayList<>();

    @Argument(shortName = "MO",
            doc = "The maximum number of lines output in verbose mode")
    public Integer MAX_OUTPUT = 100;

    @Argument(doc = "If true, only report errors and ignore warnings.")
    public boolean IGNORE_WARNINGS = false;

    @Argument(doc = "If true and input is a BAM file with an index file, also validates the index.")
    public boolean VALIDATE_INDEX = true;

    @Argument(shortName = "BISULFITE",
            doc = "Whether the SAM/BAM/CRAM file consists of bisulfite sequenced reads. " +
                    "If so, C->T is not counted as an error in computing the value of the NM tag.")
    public boolean IS_BISULFITE_SEQUENCED = false;

    @Argument(doc = "Relevant for a coordinate-sorted file containing read pairs only. " +
            "Maximum number of file handles to keep open when spilling mate info to disk. " +
            "Set this number a little lower than the per-process maximum number of file that may be open. " +
            "This number can be found by executing the 'ulimit -n' command on a Unix system.")
    public int MAX_OPEN_TEMP_FILES = 8000;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        ReferenceSequenceFile reference = null;
        if (REFERENCE_SEQUENCE != null) {
            IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
            reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);

        }
        final PrintWriter out;
        if (OUTPUT != null) {
            IOUtil.assertFileIsWritable(OUTPUT);
            try {
                out = new PrintWriter(OUTPUT);
            } catch (FileNotFoundException e) {
                // we already asserted this so we should not get here
                throw new GATKException("Unexpected exception", e);
            }
        } else {
            out = new PrintWriter(System.out);
        }

        final SamReaderFactory factory = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE)
                .validationStringency(ValidationStringency.SILENT)
                .enable(SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS);
        final SamReader samReader = factory.open(INPUT);

        if (samReader.type() != SamReader.Type.BAM_TYPE) VALIDATE_INDEX = false;

        factory.setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, VALIDATE_INDEX);
        factory.reapplyOptions(samReader);

        final SamFileValidator validator = new SamFileValidator(out, MAX_OPEN_TEMP_FILES);
        validator.setErrorsToIgnore(IGNORE);

        if (IGNORE_WARNINGS) {
            validator.setIgnoreWarnings(IGNORE_WARNINGS);
        }
        if (MODE == Mode.SUMMARY) {
            validator.setVerbose(false, 0);
        } else {
            validator.setVerbose(true, MAX_OUTPUT);
        }
        if (IS_BISULFITE_SEQUENCED) {
            validator.setBisulfiteSequenced(IS_BISULFITE_SEQUENCED);
        }
        if (VALIDATE_INDEX) {
            validator.setIndexValidationStringency(BamIndexValidator.IndexValidationStringency.EXHAUSTIVE);
        }
        if (IOUtil.isRegularPath(INPUT)) {
            // Do not check termination if reading from a stream
            validator.validateBamFileTermination(INPUT);
        }

        boolean isValid = false;

        switch (MODE) {
            case SUMMARY:
                isValid = validator.validateSamFileSummary(samReader, reference);
                break;
            case VERBOSE:
                isValid = validator.validateSamFileVerbose(samReader, reference);
                break;
        }
        out.flush();

        return isValid;
    }
}
