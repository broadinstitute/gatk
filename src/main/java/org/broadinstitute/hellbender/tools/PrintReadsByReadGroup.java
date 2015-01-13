package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.Option;
import org.broadinstitute.hellbender.cmdline.StandardOptionDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

@CommandLineProgramProperties(
        usage = "Print reads by read group.",
        usageShort = "Print reads by read group",
        programGroup = ReadProgramGroup.class
)
public class PrintReadsByReadGroup extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "The SAM or BAM or CRAM file to print reads by read group.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "The directory to write SAM or BAM or CRAM files by read group.")
    public File OUTPUT_DIRECTORY = new File("");

    @Override
    protected Object doWork() {
        printReadsByReadGroup();
        return null;
    }

    /**
     * This is factored out of doWork only for unit testing.
     */
    protected void printReadsByReadGroup() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertDirectoryIsWritable(OUTPUT_DIRECTORY);
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        Map<String, SAMFileWriter> outs = createWriters(in);
        for (final SAMRecord rec : in) {
            outs.get(rec.getReadGroup().getReadGroupId()).addAlignment(rec);
        }
        CloserUtil.close(in);
        outs.values().forEach(CloserUtil::close);
    }

    private Map<String, SAMFileWriter> createWriters(final SamReader in) {
        final Map<String, SAMFileWriter> outs = new HashMap<>();

        final SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();

        final SAMFileHeader samFileHeaderIn = in.getFileHeader();
        final String extension = FilenameUtils.getExtension(INPUT.getName());

        samFileHeaderIn.getReadGroups().forEach(samReadGroupRecord -> {
            final SAMFileHeader samFileHeaderOut = samFileHeaderIn.clone();
            samFileHeaderOut.setReadGroups(Collections.singletonList(samReadGroupRecord));

            final String sample = samReadGroupRecord.getSample();
            final String readGroupId = samReadGroupRecord.getReadGroupId();
            final File outFile = new File(OUTPUT_DIRECTORY, sample + "." + readGroupId + "." + extension);

            outs.put(readGroupId, samFileWriterFactory.makeWriter(samFileHeaderOut, true, outFile, REFERENCE_SEQUENCE));
        });

        return outs;
    }
}
