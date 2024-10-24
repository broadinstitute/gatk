package org.broadinstitute.hellbender.tools.sv.cluster;

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
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.walkers.sv.JointGermlineCNVSegmentation;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.util.Set;

import static org.broadinstitute.hellbender.tools.walkers.sv.JointGermlineCNVSegmentation.BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME;
import static org.broadinstitute.hellbender.tools.walkers.sv.JointGermlineCNVSegmentation.FLAG_FIELD_LOGIC_LONG_NAME;

/***
 * Base class for tools that a simple interface for utilizing {@link SVClusterEngine}. It handles input/output easily,
 * including output sorting with spilling to disk to avoid excessive memory usage.
 */
public abstract class SVClusterWalker extends MultiVariantWalker {
    public static final String PLOIDY_TABLE_LONG_NAME = "ploidy-table";
    public static final String VARIANT_PREFIX_LONG_NAME = "variant-prefix";
    public static final String ENABLE_CNV_LONG_NAME = "enable-cnv";
    public static final String ALGORITHM_LONG_NAME = "algorithm";
    public static final String FAST_MODE_LONG_NAME = "fast-mode";
    public static final String OMIT_MEMBERS_LONG_NAME = "omit-members";
    public static final String DEFAULT_NO_CALL_LONG_NAME = "default-no-call";
    public static final String MAX_RECORDS_IN_RAM_LONG_NAME = "max-records-in-ram";

    /**
     * The enum Cluster algorithm.
     */
    public enum CLUSTER_ALGORITHM {
        /**
         * Defragment cnv cluster algorithm. Not supported with stratification.
         */
        DEFRAGMENT_CNV,
        /**
         * Single linkage cluster algorithm.
         */
        SINGLE_LINKAGE,
        /**
         * Max clique cluster algorithm.
         */
        MAX_CLIQUE
    }

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

    /**
     * When enabled, DEL and DUP variants will be clustered together. The resulting records with have an SVTYPE of CNV.
     */
    @Argument(
            doc = "Enable clustering DEL/DUP variants together as CNVs (does not apply to CNV defragmentation)",
            fullName = ENABLE_CNV_LONG_NAME,
            optional = true
    )
    protected boolean enableCnv = false;

    /**
     * Results in substantial space and time costs for large sample sets by clearing genotypes that are not needed for
     * clustering, but any associated annotation fields will be set to null in the output.
     */
    @Argument(
            doc = "Fast mode. Drops hom-ref and missing genotype fields and emits them as missing.",
            fullName = FAST_MODE_LONG_NAME,
            optional = true
    )
    protected boolean fastMode = false;

    @Argument(
            doc = "Omit cluster member ID annotations",
            fullName = OMIT_MEMBERS_LONG_NAME,
            optional = true
    )
    protected boolean omitMembers = false;

    @Argument(fullName = BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME,
            doc = "Strategy to use for choosing a representative value for a breakpoint cluster.",
            optional = true)
    protected CanonicalSVCollapser.BreakpointSummaryStrategy breakpointSummaryStrategy =
            CanonicalSVCollapser.BreakpointSummaryStrategy.REPRESENTATIVE;

    @Argument(fullName = JointGermlineCNVSegmentation.ALT_ALLELE_SUMMARY_STRATEGY_LONG_NAME,
            doc = "Strategy to use for choosing a representative alt allele for non-CNV biallelic sites with " +
                    "different subtypes.",
            optional = true)
    protected CanonicalSVCollapser.AltAlleleSummaryStrategy altAlleleSummaryStrategy =
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE;

    @Argument(fullName = FLAG_FIELD_LOGIC_LONG_NAME,
            doc = "Logic for collapsing Flag type INFO and FORMAT fields",
            optional = true)
    protected CanonicalSVCollapser.FlagFieldLogic flagFieldLogic = CanonicalSVCollapser.FlagFieldLogic.OR;

    @Argument(fullName = ALGORITHM_LONG_NAME,
            doc = "Clustering algorithm",
            optional = true
    )
    protected CLUSTER_ALGORITHM algorithm = CLUSTER_ALGORITHM.SINGLE_LINKAGE;

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
        currentContig = null;
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

    /**
     * Subclasses should override this method
     */
    public abstract void applyRecord(final SVCallRecord record);

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        SVCallRecord call = SVCallRecordUtils.create(variant, dictionary);
        if (fastMode && call.getType() != GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            // Strip out non-carrier genotypes to save memory and compute
            // Don't do for multi-allelic CNVs since carrier status can't be determined
            final GenotypesContext filteredGenotypes = GenotypesContext.copy(call.getCarrierGenotypeList());
            call = SVCallRecordUtils.copyCallWithNewGenotypes(call, filteredGenotypes);
        }
        // Update current contig
        if (!call.getContigA().equals(currentContig)) {
            currentContig = call.getContigA();
            logger.info("Processing contig " + currentContig + "...");
        }
        applyRecord(call);
    }

    protected VCFHeader createHeader() {
        final VCFHeader header = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder(), samples);
        header.setSequenceDictionary(dictionary);

        // Required info lines
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        header.addMetaDataLine(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVLEN));
        header.addMetaDataLine(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVTYPE));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.END2_ATTRIBUTE, 1,
                VCFHeaderLineType.Integer, "Second position"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, 1,
                VCFHeaderLineType.String, "Second contig"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.STRANDS_ATTRIBUTE, 1,
                VCFHeaderLineType.String, "First and second strands"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE,
                VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Source algorithms"));
        if (!omitMembers) {
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY,
                    VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Cluster variant ids"));
        }
        // Required format lines
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
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
