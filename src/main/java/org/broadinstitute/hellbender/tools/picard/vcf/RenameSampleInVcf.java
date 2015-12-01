package org.broadinstitute.hellbender.tools.picard.vcf;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
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

import java.io.File;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Rename a sample within a VCF.",
        oneLineSummary = "Rename a sample within a VCF",
        programGroup = VariantProgramGroup.class
)
public final class RenameSampleInVcf extends PicardCommandLineProgram {
    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "Input single sample VCF.")
    public File INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output single sample VCF.")
    public File OUTPUT;

    @Argument(doc = "Existing name of sample in VCF; if provided, asserts that that is the name of the extant sample name", optional = true)
    public String OLD_SAMPLE_NAME = null;

    @Argument(doc = "New name to give sample in output VCF.")
    public String NEW_SAMPLE_NAME;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        try (final VCFFileReader in = new VCFFileReader(INPUT)) {
            final VCFHeader header = in.getFileHeader();

            if (header.getGenotypeSamples().size() > 1) {
                throw new IllegalArgumentException("Input VCF must be single-sample.");
            }

            if (OLD_SAMPLE_NAME != null && !OLD_SAMPLE_NAME.equals(header.getGenotypeSamples().get(0))) {
                throw new IllegalArgumentException("Input VCF did not contain expected sample. Contained: " + header.getGenotypeSamples().get(0));
            }

            final EnumSet<Options> options = EnumSet.copyOf(VariantContextWriterBuilder.DEFAULT_OPTIONS);
            if (CREATE_INDEX) options.add(Options.INDEX_ON_THE_FLY);
            else options.remove(Options.INDEX_ON_THE_FLY);

            final VCFHeader outHeader = new VCFHeader(header.getMetaDataInInputOrder(), CollectionUtil.makeList(NEW_SAMPLE_NAME));
            try (final VariantContextWriter out = new VariantContextWriterBuilder().setOutputFile(OUTPUT)
                    .setReferenceDictionary(outHeader.getSequenceDictionary()).setOptions(options).build()) {
                out.writeHeader(outHeader);

                for (final VariantContext ctx : in) {
                    out.add(ctx);
                }
            }
        }

        return null;
    }
}
