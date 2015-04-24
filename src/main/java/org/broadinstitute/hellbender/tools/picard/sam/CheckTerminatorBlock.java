package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.IOException;

/**
 * Simple class to check the terminator block of a SAM file.
 */
@CommandLineProgramProperties(
        usage = CheckTerminatorBlock.USAGE,
        usageShort = CheckTerminatorBlock.USAGE,
        programGroup = ReadProgramGroup.class
)
public final class CheckTerminatorBlock extends PicardCommandLineProgram {
    static final String USAGE = "Returns true if the gzip file's (e.g., BAM) last block is well-formed; false otherwise";

    @Argument(shortName= StandardArgumentDefinitions.INPUT_SHORT_NAME, doc="The block compressed file to check.")
    public File INPUT;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        try {
            final BlockCompressedInputStream.FileTermination term = BlockCompressedInputStream.checkTermination(INPUT);
            System.err.println(term.name());
            // TODO we probably shouldn't rely on the return value like this?
            return (term != BlockCompressedInputStream.FileTermination.DEFECTIVE);
        }
        catch (IOException ioe) {
            throw new UserException("Exception reading terminator block of file: " + INPUT.getAbsolutePath());
        }
    }
}
