package org.broadinstitute.hellbender.tools.picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
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

/**
 * Takes a VCF file and a Sequence Dictionary (from a variety of file types) and updates the Sequence Dictionary in VCF.
 *
 * @author George Grant
 *
 */
@CommandLineProgramProperties(
        usage = "Takes a VCF and a second file that contains a sequence dictionary and updates the VCF with the new sequence dictionary.",
        usageShort = "Takes a VCF and a second file that contains a sequence dictionary and updates the VCF with the new sequence dictionary.",
        programGroup = VariantProgramGroup.class
)
public final class UpdateVcfSequenceDictionary extends PicardCommandLineProgram {
     @Argument(shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "Input VCF")
    public File INPUT;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output VCF to be written.")
    public File OUTPUT;

    @Argument(shortName = StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, doc = "A Sequence Dictionary (can be read from one of the " +
            "following file types (SAM, BAM, VCF, BCF, Interval List, Fasta, or Dict)")
    public File SEQUENCE_DICTIONARY;

    private final Log log = Log.getInstance(UpdateVcfSequenceDictionary.class);

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SAMSequenceDictionary samSequenceDictionary = SAMSequenceDictionaryExtractor.extractDictionary(SEQUENCE_DICTIONARY);

        final VCFFileReader fileReader = new VCFFileReader(INPUT, false);
        final VCFHeader fileHeader = fileReader.getFileHeader();

        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setReferenceDictionary(samSequenceDictionary)
                .clearOptions();
        if (CREATE_INDEX)
            builder.setOption(Options.INDEX_ON_THE_FLY);

        final VariantContextWriter vcfWriter = builder.setOutputFile(OUTPUT).build();
        fileHeader.setSequenceDictionary(samSequenceDictionary);
        vcfWriter.writeHeader(fileHeader);

        final ProgressLogger progress = new ProgressLogger(log, 10000);
        final CloseableIterator<VariantContext> iterator = fileReader.iterator();
        while (iterator.hasNext()) {
            final VariantContext context = iterator.next();
            vcfWriter.add(context);
            progress.record(context.getContig(), context.getStart());
        }

        CloserUtil.close(iterator);
        CloserUtil.close(fileReader);
        vcfWriter.close();

        return null;
    }
}
