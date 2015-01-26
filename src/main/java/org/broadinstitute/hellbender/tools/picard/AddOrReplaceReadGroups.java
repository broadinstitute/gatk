package org.broadinstitute.hellbender.tools.picard;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;

import java.io.File;
import java.util.Arrays;

/**
 * Replaces read groups in a BAM file
 *
 * @author mdepristo
 */
@CommandLineProgramProperties(
        usage = "Replaces all read groups in the INPUT file with a single new read group and assigns " +
                "all reads to this read group in the OUTPUT BAM",
        usageShort = "Replaces read groups in a BAM or SAM file with a single new read group",
        programGroup = ReadProgramGroup.class
)
public class AddOrReplaceReadGroups extends PicardCommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input file (bam or sam).")
    public File INPUT = null;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (bam or sam).")
    public File OUTPUT = null;

    @Option(shortName = StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, optional = true,
            doc = "Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT.")
    public SAMFileHeader.SortOrder SORT_ORDER;

    @Option(shortName = "ID", doc = "Read Group ID")
    public String RGID = "1";

    @Option(shortName = "LB", doc = "Read Group Library")
    public String RGLB;

    @Option(shortName = "PL", doc = "Read Group platform (e.g. illumina, solid)")
    public String RGPL;

    @Option(shortName = "PU", doc = "Read Group platform unit (eg. run barcode)")
    public String RGPU;

    @Option(shortName = "SM", doc = "Read Group sample name")
    public String RGSM;

    @Option(shortName = "CN", doc = "Read Group sequencing center name", optional = true)
    public String RGCN;

    @Option(shortName = "DS", doc = "Read Group description", optional = true)
    public String RGDS;

    @Option(shortName = "DT", doc = "Read Group run date", optional = true)
    public Iso8601Date RGDT;

    @Option(shortName = "PI", doc = "Read Group predicted insert size", optional = true)
    public Integer RGPI;

    private final Log log = Log.getInstance(AddOrReplaceReadGroups.class);

    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        // create the read group we'll be using
        final SAMReadGroupRecord rg = new SAMReadGroupRecord(RGID);
        rg.setLibrary(RGLB);
        rg.setPlatform(RGPL);
        rg.setSample(RGSM);
        rg.setPlatformUnit(RGPU);
        if (RGCN != null) rg.setSequencingCenter(RGCN);
        if (RGDS != null) rg.setDescription(RGDS);
        if (RGDT != null) rg.setRunDate(RGDT);
        if (RGPI != null) rg.setPredictedMedianInsertSize(RGPI);

        log.info(String.format("Created read group ID=%s PL=%s LB=%s SM=%s%n", rg.getId(), rg.getPlatform(), rg.getLibrary(), rg.getSample()));

        // create the new header and output file
        final SAMFileHeader inHeader = in.getFileHeader();
        final SAMFileHeader outHeader = inHeader.clone();
        outHeader.setReadGroups(Arrays.asList(rg));
        if (SORT_ORDER != null) outHeader.setSortOrder(SORT_ORDER);

        final SAMFileWriter outWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(outHeader,
                outHeader.getSortOrder() == inHeader.getSortOrder(),
                OUTPUT);

        final ProgressLogger progress = new ProgressLogger(log);
        for (final SAMRecord read : in) {
            read.setAttribute(SAMTag.RG.name(), RGID);
            outWriter.addAlignment(read);
            progress.record(read);
        }

        // cleanup
        CloserUtil.close(in);
        outWriter.close();
        return null;
    }
}
