package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.walkers.sv.JointGermlineCNVSegmentation;
import org.broadinstitute.hellbender.tools.walkers.sv.SVMergingWalker;

import static org.broadinstitute.hellbender.tools.walkers.sv.JointGermlineCNVSegmentation.BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME;
import static org.broadinstitute.hellbender.tools.walkers.sv.JointGermlineCNVSegmentation.FLAG_FIELD_LOGIC_LONG_NAME;

/***
 * Base class for tools that a simple interface for utilizing {@link SVClusterEngine}. It handles input/output easily,
 * including output sorting with spilling to disk to avoid excessive memory usage.
 */
public abstract class SVClusterWalker extends SVMergingWalker {
    public static final String ENABLE_CNV_LONG_NAME = "enable-cnv";
    public static final String ALGORITHM_LONG_NAME = "algorithm";
    public static final String FAST_MODE_LONG_NAME = "fast-mode";

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

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        currentContig = null;
    }

    @Override
    protected VCFHeader createHeader() {
        final VCFHeader header = super.createHeader();

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

        // Required format lines
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
        return header;
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



}
