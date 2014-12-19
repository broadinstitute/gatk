package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
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
    public Path OUTPUT_DIRECTORY = Paths.get("");

    public static void main(final String[] args) {
        new PrintReadsByReadGroup().instanceMain(args);
    }

    @Override
    protected int doWork() {
        printReadsByReadGroup();
        return 0;
    }

    /**
     * This is factored out of doWork only for unit testing.
     */
    protected void printReadsByReadGroup() {
        IOUtil.assertFileIsReadable(INPUT);
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
        samFileWriterFactory.setCreateIndex(true);
        samFileWriterFactory.setUseAsyncIo(true);

        final SAMFileHeader samFileHeaderIn = in.getFileHeader();
        final String[] tokens = INPUT.getName().split("\\.");
        final String extension = tokens[tokens.length - 1];

        samFileHeaderIn.getReadGroups().forEach(samReadGroupRecord -> {
            final SAMFileHeader samFileHeaderOut = samFileHeaderIn.clone();
            samFileHeaderOut.setReadGroups(Collections.singletonList(samReadGroupRecord));

            final String sample = samReadGroupRecord.getSample();
            final String readGroupId = samReadGroupRecord.getReadGroupId();
            final File outFile = OUTPUT_DIRECTORY.resolve(sample + "." + readGroupId + "." + extension).toFile();

            outs.put(readGroupId, samFileWriterFactory.makeWriter(samFileHeaderOut, true, outFile, REFERENCE_SEQUENCE));
        });

        return outs;
    }
}
