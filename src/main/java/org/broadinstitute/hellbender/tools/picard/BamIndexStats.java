package org.broadinstitute.hellbender.tools.picard;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;

import java.io.File;

/**
 * Command line program to print statistics from BAM index (.bai) file
 * Statistics include count of aligned and unaligned reads for each reference sequence
 * and a count of all records with no start coordinate.
 * Similar to the 'samtools idxstats' command.
 *
 * @author Martha Borkan
 */
@CommandLineProgramProperties(
        usage = "Generates BAM index statistics, including the number of aligned and unaligned SAMRecords for each reference sequence, " +
                "and the number of SAMRecords with no coordinate." +
                "Input BAM file must have a corresponding index file.\n",
        usageShort = "Generates index statistics from a BAM file",
        programGroup = ReadProgramGroup.class
)
public class BamIndexStats extends PicardCommandLineProgram {

    private static final Log log = Log.getInstance(BamIndexStats.class);

    @Argument(shortName= StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc="A BAM file to process.")
    public File INPUT;

    /**
     * Main method for the program.  Checks that input file is present and
     * readable, then iterates through the index printing meta data to stdout.
     */
    protected Object doWork() {

        if (INPUT.getName().endsWith(BAMIndex.BAMIndexSuffix))
            log.warn("INPUT should be BAM file not index file");
        IOUtil.assertFileIsReadable(INPUT);
        BAMIndexMetaData.printIndexStats(INPUT);

        return null;
    }
}
