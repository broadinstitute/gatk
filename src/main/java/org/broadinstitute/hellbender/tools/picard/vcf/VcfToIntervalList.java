package org.broadinstitute.hellbender.tools.picard.vcf;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
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
        usage = "Converts a VCF or BCF file to a Picard Interval List.",
        usageShort = "Converts a VCF or BCF file to a Picard Interval List.",
        programGroup = VariantProgramGroup.class
)
public final class VcfToIntervalList extends PicardCommandLineProgram {
    // The following attributes define the command-line arguments
    public static final Log LOG = Log.getInstance(VcfToIntervalList.class);

    @Argument(doc = "The BCF or VCF input file. The file format is determined by file extension.", shortName= StandardArgumentDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "The output Picard Interval List")
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