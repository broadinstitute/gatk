package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.PloidyTable;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.util.Set;

public abstract class SVMergingWalker extends MultiVariantWalker {
    public static final String PLOIDY_TABLE_LONG_NAME = "ploidy-table";
    public static final String VARIANT_PREFIX_LONG_NAME = "variant-prefix";
    public static final String OMIT_MEMBERS_LONG_NAME = "omit-members";
    public static final String DEFAULT_NO_CALL_LONG_NAME = "default-no-call";
    public static final String MAX_RECORDS_IN_RAM_LONG_NAME = "max-records-in-ram";

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected GATKPath outputFile;

    /**
     * Expected format is tab-delimited and contains a header with the first column SAMPLE and remaining columns
     * contig names. Each row corresponds to a sample, with the sample ID in the first column and contig ploidy
     * integers in their respective columns.
     */
    @Argument(
            doc = "Sample ploidy table (.tsv)",
            fullName = PLOIDY_TABLE_LONG_NAME
    )
    protected GATKPath ploidyTablePath;

    @Argument(
            doc = "If supplied, generate variant IDs with this prefix",
            fullName = VARIANT_PREFIX_LONG_NAME,
            optional = true
    )
    protected String variantPrefix = null;

    @Argument(
            doc = "Omit cluster member ID annotations",
            fullName = OMIT_MEMBERS_LONG_NAME,
            optional = true
    )
    protected boolean omitMembers = false;

    /**
     * Default genotypes are assigned when they cannot be inferred from the inputs, such as when VCFs with different
     * variants and samples are provided.
     */
    @Argument(fullName = DEFAULT_NO_CALL_LONG_NAME,
            doc = "Default to no-call GT (e.g. ./.) instead of reference alleles (e.g. 0/0) when a genotype is not" +
                    " available",
            optional = true
    )
    protected boolean defaultNoCall = false;

    @Argument(fullName = MAX_RECORDS_IN_RAM_LONG_NAME,
            doc = "When writing VCF files that need to be sorted, this will specify the number of records stored in " +
                    "RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a " +
                    "VCF file, and increases the amount of RAM needed.",
            optional=true)
    public int maxRecordsInRam = 10000;

    protected SAMSequenceDictionary dictionary;
    protected ReferenceSequenceFile reference;
    protected PloidyTable ploidyTable;
    protected SortingCollection<VariantContext> sortingBuffer;
    protected VariantContextWriter writer;
    protected VCFHeader header;
    protected Set<String> samples;
    protected String currentContig;
    protected int numVariantsBuilt = 0;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        reference = ReferenceUtils.createReferenceReader(referenceArguments.getReferenceSpecifier());
        dictionary = reference.getSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        ploidyTable = new PloidyTable(ploidyTablePath.toPath());
        samples = getSamplesForVariants();
        writer = createVCFWriter(outputFile);
        header = createHeader();
        writer.writeHeader(header);
        sortingBuffer = SortingCollection.newInstance(
                VariantContext.class,
                new VCFRecordCodec(header, true),
                header.getVCFRecordComparator(),
                maxRecordsInRam,
                tmpDir.toPath());
    }

    @Override
    public Object onTraversalSuccess() {
        for (final VariantContext variant : sortingBuffer) {
            writer.add(variant);
        }
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if (sortingBuffer != null) {
            sortingBuffer.cleanup();
        }
        if (writer != null) {
            writer.close();
        }
    }


    protected VCFHeader createHeader() {
        final VCFHeader header = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder(), samples);
        header.setSequenceDictionary(dictionary);

        if (!omitMembers) {
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY,
                    VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Cluster variant ids"));
        }

        return header;
    }

    protected void write(final SVCallRecord call) {
        sortingBuffer.add(buildVariantContext(call));
    }

    protected VariantContext buildVariantContext(final SVCallRecord call) {
        // Add genotypes for missing samples
        final GenotypesContext filledGenotypes = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(
                call, samples, !defaultNoCall, ploidyTable, header);

        // Assign new variant ID
        final String newId = variantPrefix == null ? call.getId() : String.format("%s%08x", variantPrefix, numVariantsBuilt++);

        // Build new variant
        final SVCallRecord finalCall = new SVCallRecord(newId, call.getContigA(), call.getPositionA(), call.getStrandA(),
                call.getContigB(), call.getPositionB(), call.getStrandB(), call.getType(), call.getComplexSubtype(),
                call.getComplexEventIntervals(), call.getLength(), call.getEvidence(), call.getAlgorithms(), call.getAlleles(), filledGenotypes,
                call.getAttributes(), call.getFilters(), call.getLog10PError(), dictionary);
        final VariantContextBuilder builder = SVCallRecordUtils.getVariantBuilder(finalCall);
        if (omitMembers) {
            builder.rmAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY);
        }
        return builder.make();
    }
}
