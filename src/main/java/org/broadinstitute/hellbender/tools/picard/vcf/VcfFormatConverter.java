package org.broadinstitute.hellbender.tools.picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;

/**
 * Converts an ASCII VCF file to a binary BCF or vice versa.
 *
 * @author jgentry@broadinstitute.org
 */
@CommandLineProgramProperties(
        usage = "Convert a VCF file to a BCF file, or BCF to VCF.\n" +
                "Input and output formats are determined by file extension.",
        usageShort = "Converts a VCF file to a BCF file, or BCF to VCF",
        programGroup = VariantProgramGroup.class
)
public final class VcfFormatConverter extends PicardCommandLineProgram {
    // The following attributes define the command-line arguments
    public static final Log LOG = Log.getInstance(VcfFormatConverter.class);

    @Argument(doc = "The BCF or VCF input file. The file format is determined by file extension.", shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "The BCF or VCF output file. The file format is determined by file extension.", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

	@Argument(doc = "Fail if an index is not available for the input VCF/BCF")
	public Boolean REQUIRE_INDEX = true;

	public VcfFormatConverter() {
		this.CREATE_INDEX = true;
	}

    @Override
    protected Object doWork() {
        final ProgressLogger progress = new ProgressLogger(LOG, 10000);
        
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

	    final VCFFileReader reader = new VCFFileReader(INPUT, REQUIRE_INDEX);
	    final VCFHeader header = new VCFHeader(reader.getFileHeader());
	    final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
	    if (CREATE_INDEX && sequenceDictionary == null) {
		    throw new UserException("A sequence dictionary must be available in the input file when creating indexed output.");
	    }

        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setOutputFile(OUTPUT)
                .setReferenceDictionary(sequenceDictionary);
        if (CREATE_INDEX)
            builder.setOption(Options.INDEX_ON_THE_FLY);
        else
            builder.unsetOption(Options.INDEX_ON_THE_FLY);
        final VariantContextWriter writer = builder.build();
        writer.writeHeader(header);
	    final CloseableIterator<VariantContext> iterator = reader.iterator();

	    while (iterator.hasNext()) {
		    final VariantContext context = iterator.next();
            writer.add(context);
            progress.record(context.getContig(), context.getStart());
        }

	    CloserUtil.close(iterator);
	    CloserUtil.close(reader);
        writer.close();

        return null;
    }
}
