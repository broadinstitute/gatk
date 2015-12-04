package org.broadinstitute.hellbender.tools.picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
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
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;

/**
 * Splits the input VCF file into two, one for indels and one for SNPs. The headers of the two output
 * files will be identical.
 * <p/>
 * An index file is created for the output file by default. Using an output file name with a ".gz"
 * extension will create gzip-compressed output.
 */
@CommandLineProgramProperties(
        summary = "Splits an input VCF file into two VCF files, one for indel records and one for SNPs. The" +
                "headers of the two output files will be identical. An index file is created and a" +
                "sequence dictionary is required by default.",
        oneLineSummary = "Splits an input VCF file into two VCF files",
        programGroup = VariantProgramGroup.class
)
public final class SplitVcfs extends PicardCommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, 
            doc = "The VCF input file")
    public File INPUT;

    @Argument(doc = "The VCF file to which SNP records should be written. The file format is determined by file extension.")
    public File SNP_OUTPUT;

    @Argument(doc = "The VCF file to which indel records should be written. The file format is determined by file extension.")
    public File INDEL_OUTPUT;

    @Argument(shortName = "D", doc = "The index sequence dictionary to use instead of the sequence dictionaries in the input files", optional = true)
    public File SEQUENCE_DICTIONARY;

    @Argument(doc = "If true an exception will be thrown if an event type other than SNP or indel is encountered")
    public Boolean STRICT = true;

    public SplitVcfs() {
        this.CREATE_INDEX = true;
    }

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        final ProgressLogger progress = new ProgressLogger(logger, 10000);

        final VCFFileReader fileReader = new VCFFileReader(INPUT);
        final VCFHeader fileHeader = fileReader.getFileHeader();

        final SAMSequenceDictionary sequenceDictionary =
                SEQUENCE_DICTIONARY != null
                        ? SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).getFileHeader(SEQUENCE_DICTIONARY).getSequenceDictionary()
                        : fileHeader.getSequenceDictionary();
        if (CREATE_INDEX && sequenceDictionary == null) {
            throw new UserException("A sequence dictionary must be available (either through the input file or by setting it explicitly) when creating indexed output.");
        }

        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setReferenceDictionary(sequenceDictionary)
                .clearOptions();
        if (CREATE_INDEX)
            builder.setOption(Options.INDEX_ON_THE_FLY);

        try (final VariantContextWriter snpWriter = builder.setOutputFile(SNP_OUTPUT).build();
             final VariantContextWriter indelWriter = builder.setOutputFile(INDEL_OUTPUT).build()) {
            snpWriter.writeHeader(fileHeader);
            indelWriter.writeHeader(fileHeader);

            int incorrectVariantCount = 0;

            final CloseableIterator<VariantContext> iterator = fileReader.iterator();
            while (iterator.hasNext()) {
                final VariantContext context = iterator.next();
                if (context.isIndel()) indelWriter.add(context);
                else if (context.isSNP()) snpWriter.add(context);
                else {
                    if (STRICT) throw new IllegalStateException("Found a record with type " + context.getType().name());
                    else incorrectVariantCount++;
                }

                progress.record(context.getContig(), context.getStart());
            }

            if (incorrectVariantCount > 0) {
                logger.debug("Found " + incorrectVariantCount + " records that didn't match SNP or INDEL");
            }

            CloserUtil.close(iterator);
            CloserUtil.close(fileReader);
        }

        return null;
    }
}
