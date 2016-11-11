package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Iso8601Date;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;
import java.util.Arrays;

/**
 * Replaces read groups in a BAM file
 *
 * @author mdepristo
 */
@CommandLineProgramProperties(
        summary = "Replaces all read groups in the input file with a single new read group and assigns " +
                "all reads to this read group in the output SAM/BAM/CRAM",
        oneLineSummary = "Replaces read groups in a SAM/BAM/CRAM file with a single new read group",
        programGroup = ReadProgramGroup.class
)
public final class AddOrReplaceReadGroups extends PicardCommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
         doc = "Input SAM/BAM/CRAM file.")
    public File INPUT = null;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
        doc = "Output SAM/BAM/CRAM file.")
    public File OUTPUT = null;

    @Argument(shortName = StandardArgumentDefinitions.SORT_ORDER_SHORT_NAME, optional = true,
            doc = "Optional sort order to output in. If not supplied output is in the same order as input.")
    public SAMFileHeader.SortOrder SORT_ORDER;

    public static final String RGID_SHORT_NAME = "ID";
    @Argument(shortName = RGID_SHORT_NAME, doc = "Read Group ID")
    public String RGID = "1";

    public static final String RGLB_SHORT_NAME = "LB";
    @Argument(shortName = RGLB_SHORT_NAME, doc = "Read Group Library")
    public String RGLB;

    public static final String RGPL_SHORT_NAME = "PL";
    @Argument(shortName = RGPL_SHORT_NAME, doc = "Read Group platform (e.g. illumina, solid)")
    public String RGPL;

    public static final String RGPU_SHORT_NAME = "PU";
    @Argument(shortName = RGPU_SHORT_NAME, doc = "Read Group platform unit (eg. run barcode)")
    public String RGPU;

    public static final String RGSM_SHORT_NAME = "SM";
    @Argument(shortName = RGSM_SHORT_NAME, doc = "Read Group sample name")
    public String RGSM;

    public static final String RGCN_SHORT_NAME = "CN";
    @Argument(shortName = RGCN_SHORT_NAME, doc = "Read Group sequencing center name", optional = true)
    public String RGCN;

    public static final String RGDS_SHORT_NAME = "DS";
    @Argument(shortName = RGDS_SHORT_NAME, doc = "Read Group description", optional = true)
    public String RGDS;

    public static final String RGDT_SHORT_NAME = "DT";
    @Argument(shortName = RGDT_SHORT_NAME, doc = "Read Group run date", optional = true)
    public Iso8601Date RGDT;

    public static final String RGPI_SHORT_NAME = "PI";
    @Argument(shortName = RGPI_SHORT_NAME, doc = "Read Group predicted insert size", optional = true)
    public Integer RGPI;

    public static final String RGPG_SHORT_NAME = "PG";
    @Argument(shortName = RGPG_SHORT_NAME, doc = "Read Group program group", optional = true)
    public String RGPG;

    public static final String RGPM_SHORT_NAME = "PM";
    @Argument(shortName = RGPM_SHORT_NAME, doc = "Read Group platform model", optional = true)
    public String RGPM;

    @Override
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
        if (RGPG != null) rg.setProgramGroup(RGPG);
        if (RGPM != null) rg.setPlatformModel(RGPM);

        logger.info(String.format("Created read group ID=%s PL=%s LB=%s SM=%s%n", rg.getId(), rg.getPlatform(), rg.getLibrary(), rg.getSample()));

        // create the new header and output file
        final SAMFileHeader inHeader = in.getFileHeader();
        final SAMFileHeader outHeader = ReadUtils.cloneSAMFileHeader(inHeader);
        outHeader.setReadGroups(Arrays.asList(rg));
        if (SORT_ORDER != null) outHeader.setSortOrder(SORT_ORDER);

        try (final SAMFileWriter outWriter = createSAMWriter(
                OUTPUT,
                REFERENCE_SEQUENCE,
                outHeader,
                outHeader.getSortOrder() == inHeader.getSortOrder())) {

            final ProgressLogger progress = new ProgressLogger(logger);
            for (final SAMRecord read : in) {
                read.setAttribute(SAMTag.RG.name(), RGID);
                outWriter.addAlignment(read);
                progress.record(read);
            }
        } finally {
            // cleanup
            CloserUtil.close(in);
        }
        return null;
    }
}
