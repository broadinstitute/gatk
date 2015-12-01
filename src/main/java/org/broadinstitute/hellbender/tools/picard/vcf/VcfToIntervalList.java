package org.broadinstitute.hellbender.tools.picard.vcf;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;

import java.io.File;

/**
 * Creates an interval list from a VCF
 *
 * @author ggrant@broadinstitute.org
 */

@CommandLineProgramProperties(
        summary = "Converts a VCF file to a Picard Interval List",
        oneLineSummary = "Converts a VCF file to a Picard Interval List",
        programGroup = VariantProgramGroup.class
)
public final class VcfToIntervalList extends PicardCommandLineProgram {
    // The following attributes define the command-line arguments

    @Argument(doc = "The VCF input file. The file format is determined by file extension.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName= StandardArgumentDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "The output Picard Interval List")
    public File OUTPUT;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final IntervalList intervalList = VCFFileReader.fromVcf(INPUT);
        // Sort and write the output
        intervalList.uniqued().write(OUTPUT);
        return null;
    }
}
