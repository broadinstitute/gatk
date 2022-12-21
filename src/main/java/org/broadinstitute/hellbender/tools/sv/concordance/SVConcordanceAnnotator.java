package org.broadinstitute.hellbender.tools.sv.concordance;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVAlleleCounter;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import picard.vcf.*;

import java.util.*;

/**
 * Generates SV records annotated with concordance metrics given a pair of "evaluation" and "truth" SVs.
 *
 * Multi-allelic CNVs are annotated only with copy state concordance.
 */
public class SVConcordanceAnnotator {

    protected final Logger logger = LogManager.getLogger(this.getClass());

    private final GenotypeConcordanceScheme scheme;
    final boolean useTruthAf;

    /**
     * @param useTruthAf annotate truth allele counts using original truth record attributes (AF/AC/AN)
     */
    public SVConcordanceAnnotator(final boolean useTruthAf) {
        this.useTruthAf = useTruthAf;
        this.scheme = new SVGenotypeConcordanceScheme();
    }

    /**
     * Annotates the given evaluation record with genotype concordance metrics.
     */
    public SVCallRecord annotate(final ClosestSVFinder.ClosestPair pair) {
        final SVCallRecord evalRecord = pair.getEvalItem();
        final GenotypesContext evalGenotypes = evalRecord.getGenotypes();
        final SVCallRecord truthRecord = pair.getClosest();

        final ArrayList<Genotype> newGenotypes = new ArrayList<>(evalGenotypes.size());
        final GenotypeConcordanceCounts counts = new GenotypeConcordanceCounts();
        final boolean isCnv = evalRecord.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV;
        Integer numCnvMatches = 0;
        for (final String sample : evalGenotypes.getSampleNames()) {
            GenotypeBuilder builder = new GenotypeBuilder(evalGenotypes.get(sample));
            if (isCnv) {
                final Boolean result = cnvGenotypesMatch(sample, evalRecord, truthRecord);
                builder = builder.attribute(GATKSVVCFConstants.TRUTH_CN_EQUAL_FORMAT, result == null ? null : result.booleanValue() ? 1 : 0);
                if (result == null) {
                    // Null this metric if we encounter any null results
                    numCnvMatches = null;
                } else if (numCnvMatches != null && result.booleanValue()) {
                    numCnvMatches++;
                }
            } else {
                final GenotypeConcordanceStates.TruthAndCallStates states = getStates(sample, evalRecord, truthRecord);
                counts.increment(states);
                builder = builder.attribute(GenotypeConcordance.CONTINGENCY_STATE_TAG,
                        scheme.getContingencyStateString(states.truthState, states.callState));
            }
            newGenotypes.add(builder.make());
        }
        final SVCallRecord recordWithGenotypes = SVCallRecordUtils.copyCallWithNewGenotypes(evalRecord, GenotypesContext.create(newGenotypes));
        final Map<String, Object> attributes = new HashMap<>(recordWithGenotypes.getAttributes());
        final ConcordanceState variantStatus = truthRecord == null ? ConcordanceState.FALSE_POSITIVE : ConcordanceState.TRUE_POSITIVE;
        final String closestVariantId = truthRecord == null ? null : truthRecord.getId();
        attributes.put(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, closestVariantId);
        attributes.put(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE, variantStatus.getAbbreviation());

        if (isCnv) {
            final Double cnvConcordance = numCnvMatches == null ? null : newGenotypes.isEmpty() ? Double.NaN : numCnvMatches / (double) newGenotypes.size();
            attributes.put(GATKSVVCFConstants.COPY_NUMBER_CONCORDANCE_INFO, cnvConcordance);
        } else if (truthRecord != null) {
            final GenotypeConcordanceSummaryMetrics metrics = new GenotypeConcordanceSummaryMetrics(VariantContext.Type.SYMBOLIC, counts, "truth", "eval", true);
            attributes.put(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, metrics.GENOTYPE_CONCORDANCE);
            attributes.put(GATKSVVCFConstants.NON_REF_GENOTYPE_CONCORDANCE_INFO, metrics.NON_REF_GENOTYPE_CONCORDANCE);
            attributes.put(GATKSVVCFConstants.HET_PPV_INFO, metrics.HET_PPV);
            attributes.put(GATKSVVCFConstants.HET_SENSITIVITY_INFO, metrics.HET_SENSITIVITY);
            attributes.put(GATKSVVCFConstants.HOMVAR_PPV_INFO, metrics.HOMVAR_PPV);
            attributes.put(GATKSVVCFConstants.HOMVAR_SENSITIVITY_INFO, metrics.HOMVAR_SENSITIVITY);
            attributes.put(GATKSVVCFConstants.VAR_PPV_INFO, metrics.VAR_PPV);
            attributes.put(GATKSVVCFConstants.VAR_SENSITIVITY_INFO, metrics.VAR_SENSITIVITY);
            attributes.put(GATKSVVCFConstants.VAR_SPECIFICITY_INFO, metrics.VAR_SPECIFICITY);
        }

        if (evalRecord.getType() != GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            // Compute allele frequency in eval
            final SVAlleleCounter counter = new SVAlleleCounter(evalRecord.getAltAlleles(), evalRecord.getGenotypes());
            attributes.put(VCFConstants.ALLELE_COUNT_KEY, counter.getCounts());
            attributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, counter.getFrequencies());
            attributes.put(VCFConstants.ALLELE_NUMBER_KEY, counter.getNumber());

            // Add in truth AF
            if (truthRecord == null) {
                attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO, null);
                attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO, null);
                attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO, null);
            } else if (useTruthAf) {
                // Use AF
                final Map<String, Object> truthAttr = truthRecord.getAttributes();
                attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO,
                        truthAttr.get(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO));
                attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO,
                        truthAttr.get(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO));
                attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO,
                        truthAttr.get(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO));
            } else {
                // Calculate truth AF
                final SVAlleleCounter truthCounter = new SVAlleleCounter(evalRecord.getAltAlleles(), truthRecord.getGenotypes());
                attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO, truthCounter.getCounts());
                attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO, truthCounter.getFrequencies());
                attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO, truthCounter.getNumber());
            }
        }

        return SVCallRecordUtils.copyCallWithNewAttributes(recordWithGenotypes, attributes);
    }

    /**
     * Get truth/call states for the genotypes of the given sample
     */
    private GenotypeConcordanceStates.TruthAndCallStates getStates(final String sample,
                                                                   final SVCallRecord eval,
                                                                   final SVCallRecord truth) {
        final List<Allele> altAlleles = eval.getAltAlleles();
        if (altAlleles.size() > 1) {
            throw new IllegalArgumentException("Record " + eval.getId() + " is multiallelic but this is not supported");
        }
        final Genotype evalGenotype = eval.getGenotypes().get(sample);
        final Genotype truthGenotype = truth == null ? null : truth.getGenotypes().get(sample);
        final GenotypeConcordanceStates.TruthState truthState = getTruthState(truthGenotype);
        final GenotypeConcordanceStates.CallState callState = getEvalState(evalGenotype);
        return new GenotypeConcordanceStates.TruthAndCallStates(truthState, callState);
    }

    /**
     * Returns whether the copy state of the given sample's genotype matches. Only use for multi-allelic CNVs.
     */
    public Boolean cnvGenotypesMatch(final String sample, final SVCallRecord eval, final SVCallRecord truth) {
        Utils.nonNull(sample);
        Utils.nonNull(eval);
        Utils.validateArg(eval.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, "Expected CNV evaluation variant");
        if (truth == null) {
            return null;
        }
        final Genotype evalGenotype = eval.getGenotypes().get(sample);
        final int evalCopyNumber = VariantContextGetters.getAttributeAsInt(evalGenotype, GATKSVVCFConstants.COPY_NUMBER_FORMAT, -1);
        final Genotype truthGenotype = truth.getGenotypes().get(sample);
        if (evalGenotype.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT) != truthGenotype.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) {
            throw new IllegalArgumentException("One genotype for sample " + sample + " has CN but the other does not");
        }
        final int copyNumber = VariantContextGetters.getAttributeAsInt(truthGenotype, GATKSVVCFConstants.COPY_NUMBER_FORMAT, -1);
        return copyNumber == evalCopyNumber;
    }

    @VisibleForTesting
    protected GenotypeConcordanceStates.TruthState getTruthState(final Genotype g) {
        // Treat non-existent genotype as hom-ref
        if (g == null || g.isHomRef()) {
            return GenotypeConcordanceStates.TruthState.HOM_REF;
        } else if (g.isHet()) {
            return GenotypeConcordanceStates.TruthState.HET_REF_VAR1;
        } else if (g.isHomVar()) {
            return GenotypeConcordanceStates.TruthState.HOM_VAR1;
        } else if (g.isNoCall() || g.getPloidy() == 0) {
            return GenotypeConcordanceStates.TruthState.NO_CALL;
        } else if (g.isMixed()) {
            return GenotypeConcordanceStates.TruthState.IS_MIXED;
        } else {
            throw new IllegalArgumentException("Could not determine truth state for genotype: " + g);
        }
    }

    @VisibleForTesting
    protected GenotypeConcordanceStates.CallState getEvalState(final Genotype g) {
        Utils.nonNull(g);
        if (g.isHomRef()) {
            return GenotypeConcordanceStates.CallState.HOM_REF;
        } else if (g.isHet()) {
            return GenotypeConcordanceStates.CallState.HET_REF_VAR1;
        } else if (g.isHomVar()) {
            return GenotypeConcordanceStates.CallState.HOM_VAR1;
        } else if (g.isNoCall() || g.getPloidy() == 0) {
            return GenotypeConcordanceStates.CallState.NO_CALL;
        } else if (g.isMixed()) {
            return GenotypeConcordanceStates.CallState.IS_MIXED;
        } else {
            throw new IllegalArgumentException("Could not determine eval state for genotype: " + g);
        }
    }

    /**
     * Based on {@link GA4GHSchemeWithMissingAsHomRef}. Unused rows have been removed for simplicity.
     */
    private class SVGenotypeConcordanceScheme extends GenotypeConcordanceScheme {

        @Override
        protected void initiateScheme() {
            /** ROW STATE                                               MISSING       HOM_REF       HET_REF_VAR1       HET_VAR1_VAR2        HOM_VAR1        NO_CALL        LOW_GQ     LOW_DP     VC_FILTERED GT_FILTERED IS_MIXED    **/
            addRow(GenotypeConcordanceStates.CallState.MISSING,         TN_ONLY,      TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.HOM_REF,         TN_ONLY,      TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.HET_REF_VAR1,    FP_TN,        FP_TN,        TP_TN,             TP_FN,               TP_FN,          EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.HET_REF_VAR2,    FP_TN,        FP_TN,        FP_TN_FN,          TP_FP_FN,            FP_FN,          EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.HET_VAR1_VAR2,   FP_ONLY,      FP_ONLY,      TP_FP,             TP_ONLY,             TP_FP_FN,       EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.HOM_VAR1,        FP_ONLY,      FP_ONLY,      TP_FP,             TP_FN,               TP_ONLY,        EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.HOM_VAR2,        FP_ONLY,      FP_ONLY,      FP_FN,             TP_FN,               FP_FN,          EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.NO_CALL,         EMPTY,        EMPTY,        EMPTY,             EMPTY,               EMPTY,          EMPTY,         NA,        NA,        NA,        NA,        NA);
            addRow(GenotypeConcordanceStates.CallState.IS_MIXED,        EMPTY,        EMPTY,        EMPTY,             EMPTY,               EMPTY,          EMPTY,         NA,        NA,        NA,        NA,        NA);
        }
    }
}
