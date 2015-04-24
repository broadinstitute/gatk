package org.broadinstitute.hellbender.tools.picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFUtils;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.util.*;

/**
 * Combines multiple VCF files into a single file. Input files must be sorted by their contigs
 * and, within contigs, by start position. Throws IllegalArgumentException if the contig lists
 * are not present in the input files, are not identical or if the sample lists are not the
 * same; this class uses the GATK to merge headers, which may throw exceptions if the headers
 * cannot be merged. See VCFUtils.smartMergeHeaders for details.
 * <p/>
 * An index file is created for the output file by default. Using an output file name with a
 * ".gz" extension will create gzip-compressed output.
 */
@CommandLineProgramProperties(
        usage = "Merges multiple VCF or BCF files into one VCF file. Input files must be sorted by their contigs " +
                "and, within contigs, by start position. The input files must have the same sample and " +
                "contig lists. An index file is created and a sequence dictionary is required by default.",
        usageShort = "Merges multiple VCF or BCF files into one VCF file or BCF",
        programGroup = VariantProgramGroup.class
)
public final class MergeVcfs extends PicardCommandLineProgram {

    @Argument(shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc="VCF or BCF input files File format is determined by file extension.", optional = false)
    public List<File> INPUT;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "The merged VCF or BCF file. File format is determined by file extension.")
    public File OUTPUT;

    @Argument(shortName = "D", doc = "The index sequence dictionary to use instead of the sequence dictionary in the input file", optional = true)
    public File SEQUENCE_DICTIONARY;

    private final Log log = Log.getInstance(MergeVcfs.class);

    public MergeVcfs() {
        this.CREATE_INDEX = true;
    }

    @Override
    protected Object doWork() {
        final ProgressLogger progress = new ProgressLogger(log, 10000);
        final List<String> sampleList = new ArrayList<String>();
        final Collection<CloseableIterator<VariantContext>> iteratorCollection = new ArrayList<CloseableIterator<VariantContext>>(INPUT.size());
        final Collection<VCFHeader> headers = new HashSet<VCFHeader>(INPUT.size());

        VariantContextComparator variantContextComparator = null;
        SAMSequenceDictionary sequenceDictionary = null;

        if (SEQUENCE_DICTIONARY != null) {
            sequenceDictionary = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(SEQUENCE_DICTIONARY).getFileHeader().getSequenceDictionary();
        }

        for (final File file : INPUT) {
            IOUtil.assertFileIsReadable(file);
            final VCFFileReader fileReader = new VCFFileReader(file, false);
            final VCFHeader fileHeader = fileReader.getFileHeader();

            if (variantContextComparator == null) {
                variantContextComparator = fileHeader.getVCFRecordComparator();
            } else {
                if (!variantContextComparator.isCompatible(fileHeader.getContigLines())) {
                    throw new IllegalArgumentException(
                            "The contig entries in input file " + file.getAbsolutePath() + " are not compatible with the others.");
                }
            }

            if (sequenceDictionary == null) sequenceDictionary = fileHeader.getSequenceDictionary();

            if (sampleList.isEmpty()) {
                sampleList.addAll(fileHeader.getSampleNamesInOrder());
            } else {
                if (!sampleList.equals(fileHeader.getSampleNamesInOrder())) {
                    throw new IllegalArgumentException("Input file " + file.getAbsolutePath() + " has sample entries that don't match the other files.");
                }
            }

            headers.add(fileHeader);
            iteratorCollection.add(fileReader.iterator());
        }

        if (CREATE_INDEX && sequenceDictionary == null) {
            throw new UserException("A sequence dictionary must be available (either through the input file or by setting it explicitly) when creating indexed output.");
        }

        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setOutputFile(OUTPUT)
                .setReferenceDictionary(sequenceDictionary)
                .clearOptions();
        if (CREATE_INDEX) {
            builder.setOption(Options.INDEX_ON_THE_FLY);
        }
        final VariantContextWriter writer = builder.build();

        writer.writeHeader(new VCFHeader(VCFUtils.smartMergeHeaders(headers, false), sampleList));

        final MergingIterator<VariantContext> mergingIterator = new MergingIterator<VariantContext>(variantContextComparator, iteratorCollection);
        while (mergingIterator.hasNext()) {
            final VariantContext context = mergingIterator.next();
            writer.add(context);
            progress.record(context.getContig(), context.getStart());
        }

        CloserUtil.close(mergingIterator);
        writer.close();
        return null;
    }
}
