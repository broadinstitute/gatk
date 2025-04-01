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
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import picard.vcf.*;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Generates SV records annotated with concordance metrics given a pair of "evaluation" and "truth" SVs.
 *
 * Multi-allelic CNVs are annotated only with copy state concordance.
 */
public class SVConcordanceAnnotator {

    protected final Logger logger = LogManager.getLogger(this.getClass());

    private final GenotypeConcordanceScheme scheme;
    private final Set<String> samples;

    /**
     * Default constructor where all eval record samples will be used for concordance.
     */
    public SVConcordanceAnnotator() {
        this(null);
    }

    /**
     * Option to enable annotating variants with overlap metrics.
     */
    public SVConcordanceAnnotator(final Set<String> samples) {
        this.samples = samples;
        this.scheme = new SVGenotypeConcordanceScheme();
    }

    /**
     * Annotates the given evaluation record with genotype concordance metrics.
     */
    public SVCallRecord annotate(final ClosestSVFinder.ClosestPair pair) {
        final SVCallRecord evalRecord = pair.getEvalItem();
        final GenotypesContext evalGenotypes = evalRecord.getGenotypes();
        final SVCallRecord truthRecord = pair.getClosest();
        final List<Allele> evalAltList = evalRecord.getAltAlleles();
        Utils.validate(evalAltList.size() == 1, "Eval record must have exactly 1 alt allele");
        final Allele evalAlt = evalAltList.get(0);
        final Allele truthAlt;
        if (truthRecord == null) {
            truthAlt = null;
        } else {
            final List<Allele> truthAltList = truthRecord.getAltAlleles();
            Utils.validate(truthAltList.size() == 1, "Truth record must have exactly 1 alt allele");
            truthAlt = truthAltList.get(0);
        }

        final ArrayList<Genotype> newGenotypes = new ArrayList<>(evalGenotypes.size());
        final GenotypeConcordanceCounts counts = new GenotypeConcordanceCounts();
        final boolean isMultiallelicCnv = evalRecord.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV;
        int numCnvMatches = 0;
        int numValidCnvComparisons = 0;
        for (final String sample : evalGenotypes.getSampleNames()) {
            GenotypeBuilder builder = new GenotypeBuilder(evalGenotypes.get(sample));
            if (samples == null || samples.contains(sample)) {
                if (isMultiallelicCnv) {
                    final Boolean result = copyNumbersMatch(sample, evalRecord, truthRecord);
                    builder = builder.attribute(GATKSVVCFConstants.TRUTH_CN_EQUAL_FORMAT, result == null ? null : result.booleanValue() ? 1 : 0);
                    if (result != null) {
                        numValidCnvComparisons++;
                        if (result.booleanValue()) {
                            numCnvMatches++;
                        }
                    }
                }
                Genotype evalGenotype = evalRecord.getGenotypes().get(sample);
                Genotype truthGenotype = truthRecord == null ? null : truthRecord.getGenotypes().get(sample);
                final GenotypeConcordanceStates.TruthAndCallStates states = evalRecord.isSimpleCNV() ? getCNVStates(evalGenotype, truthGenotype, evalAlt, truthAlt, evalRecord.getRefAllele()) : getStates(evalGenotype, truthGenotype);
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

        if (isMultiallelicCnv) {
            final Double cnvConcordance = numValidCnvComparisons == 0 ? null : numCnvMatches / (double) numValidCnvComparisons;
            attributes.put(GATKSVVCFConstants.COPY_NUMBER_CONCORDANCE_INFO, cnvConcordance);
        }
        if (truthRecord != null) {
            final GenotypeConcordanceSummaryMetrics metrics = new GenotypeConcordanceSummaryMetrics(VariantContext.Type.SYMBOLIC, counts, "truth", "eval", true);
            attributes.put(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, Double.isNaN(metrics.GENOTYPE_CONCORDANCE) ? null : metrics.GENOTYPE_CONCORDANCE);
            attributes.put(GATKSVVCFConstants.NON_REF_GENOTYPE_CONCORDANCE_INFO, Double.isNaN(metrics.NON_REF_GENOTYPE_CONCORDANCE) ? null : metrics.NON_REF_GENOTYPE_CONCORDANCE);
            attributes.put(GATKSVVCFConstants.HET_PPV_INFO, Double.isNaN(metrics.HET_PPV) ? null : metrics.HET_PPV);
            attributes.put(GATKSVVCFConstants.HET_SENSITIVITY_INFO, Double.isNaN(metrics.HET_SENSITIVITY) ? null : metrics.HET_SENSITIVITY);
            attributes.put(GATKSVVCFConstants.HOMVAR_PPV_INFO, Double.isNaN(metrics.HOMVAR_PPV) ? null : metrics.HOMVAR_PPV);
            attributes.put(GATKSVVCFConstants.HOMVAR_SENSITIVITY_INFO, Double.isNaN(metrics.HOMVAR_SENSITIVITY) ? null : metrics.HOMVAR_SENSITIVITY);
            attributes.put(GATKSVVCFConstants.VAR_PPV_INFO, Double.isNaN(metrics.VAR_PPV) ? null : metrics.VAR_PPV);
            attributes.put(GATKSVVCFConstants.VAR_SENSITIVITY_INFO, Double.isNaN(metrics.VAR_SENSITIVITY) ? null : metrics.VAR_SENSITIVITY);
            attributes.put(GATKSVVCFConstants.VAR_SPECIFICITY_INFO, Double.isNaN(metrics.VAR_SPECIFICITY) ? null : metrics.VAR_SPECIFICITY);
        }
        if (truthRecord != null) {
            final CanonicalSVLinkage.CanonicalLinkageResult linkageResult = pair.getLinkage();
            attributes.put(GATKSVVCFConstants.TRUTH_RECIPROCAL_OVERLAP_INFO, linkageResult.getReciprocalOverlap());
            attributes.put(GATKSVVCFConstants.TRUTH_SIZE_SIMILARITY_INFO, linkageResult.getSizeSimilarity());
            attributes.put(GATKSVVCFConstants.TRUTH_DISTANCE_START_INFO, linkageResult.getBreakpointDistance1());
            attributes.put(GATKSVVCFConstants.TRUTH_DISTANCE_END_INFO, linkageResult.getBreakpointDistance2());
        }

        if (evalRecord.getType() != GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            if (!evalRecord.getAllSamples().isEmpty() && !hasAlleleFrequencyAnnotations(evalRecord)) {
                // Compute allele frequency in eval
                final SVAlleleCounter counter = new SVAlleleCounter(evalRecord.getAltAlleles(), evalRecord.getGenotypes());
                attributes.put(VCFConstants.ALLELE_COUNT_KEY, counter.getCounts());
                attributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, counter.getFrequencies());
                attributes.put(VCFConstants.ALLELE_NUMBER_KEY, counter.getNumber());
            }

            // Add in truth AF
            if (truthRecord == null) {
                attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO, null);
                attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO, null);
                attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO, null);
            } else {
                if (hasAlleleFrequencyAnnotations(truthRecord)) {
                    // Use AF
                    final Map<String, Object> truthAttr = truthRecord.getAttributes();
                    attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO,
                            truthAttr.get(VCFConstants.ALLELE_COUNT_KEY));
                    attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO,
                            truthAttr.get(VCFConstants.ALLELE_FREQUENCY_KEY));
                    attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO,
                            truthAttr.get(VCFConstants.ALLELE_NUMBER_KEY));
                } else {
                    // Calculate truth AF
                    final SVAlleleCounter truthCounter = new SVAlleleCounter(evalRecord.getAltAlleles(), truthRecord.getGenotypes());
                    attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO, truthCounter.getCounts());
                    attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO, truthCounter.getFrequencies());
                    attributes.put(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO, truthCounter.getNumber());
                }
            }
        }

        return SVCallRecordUtils.copyCallWithNewAttributes(recordWithGenotypes, attributes);
    }

    /**
     * Returns whether the copy state of the given sample's genotype matches. Only use for multi-allelic CNVs.
     */
    protected Boolean copyNumbersMatch(final String sample, final SVCallRecord eval, final SVCallRecord truth) {
        Utils.nonNull(sample);
        if (eval == null || truth == null) {
            return null;
        }
        final Genotype evalGenotype = eval.getGenotypes().get(sample);
        final Genotype truthGenotype = truth.getGenotypes().get(sample);
        if (evalGenotype == null || truthGenotype == null) {
            return null;
        }
        final boolean evalHasCn = evalGenotype.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT);
        final boolean truthHasCn = truthGenotype.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT);
        if (!evalHasCn || !truthHasCn) {
            return null;
        }
        final int evalCopyNumber = VariantContextGetters.getAttributeAsInt(evalGenotype, GATKSVVCFConstants.COPY_NUMBER_FORMAT, -1);
        final int truthCopyNumber = VariantContextGetters.getAttributeAsInt(truthGenotype, GATKSVVCFConstants.COPY_NUMBER_FORMAT, -1);
        return truthCopyNumber == evalCopyNumber;
    }

    private boolean hasAlleleFrequencyAnnotations(final SVCallRecord record) {
        Utils.nonNull(record);
        final Map<String, Object> attr = record.getAttributes();
        return (attr.containsKey(VCFConstants.ALLELE_COUNT_KEY) && attr.get(VCFConstants.ALLELE_COUNT_KEY) != null)
                        && (attr.containsKey(VCFConstants.ALLELE_FREQUENCY_KEY) && attr.get(VCFConstants.ALLELE_FREQUENCY_KEY) != null)
                        && (attr.containsKey(VCFConstants.ALLELE_NUMBER_KEY) && attr.get(VCFConstants.ALLELE_NUMBER_KEY) != null);
    }

    /**
     * Get truth/call states for the genotypes of the given sample
     */
    private GenotypeConcordanceStates.TruthAndCallStates getStates(final Genotype evalGenotype, final Genotype truthGenotype) {
        final GenotypeConcordanceStates.TruthState truthState = getTruthState(truthGenotype);
        final GenotypeConcordanceStates.CallState callState = getEvalState(evalGenotype);
        return new GenotypeConcordanceStates.TruthAndCallStates(truthState, callState);
    }

    /**
     * Returns whether the copy state of the given sample's genotype matches. Only use for multi-allelic CNVs.
     */
    @VisibleForTesting
    protected GenotypeConcordanceStates.TruthAndCallStates getCNVStates(final Genotype evalGenotype, final Genotype truthGenotype,
                                                                        final Allele evalAlt, final Allele truthAlt, final Allele refAllele) {
        Utils.nonNull(evalGenotype);
        final List<Allele> evalAlleles = getCNVAlleles(evalGenotype, evalAlt, refAllele);
        List<Allele> truthAlleles = getCNVAlleles(truthGenotype, truthAlt, refAllele);
        final List<Allele> evalAlts = evalAlleles.stream().filter(a -> a == Allele.SV_SIMPLE_DEL || a == Allele.SV_SIMPLE_DUP).distinct().collect(Collectors.toUnmodifiableList());
        final List<Allele> truthAlts = truthAlleles.stream().filter(a -> a == Allele.SV_SIMPLE_DEL || a == Allele.SV_SIMPLE_DUP).distinct().collect(Collectors.toUnmodifiableList());
        Utils.validate(evalAlts.size() < 2, "Unexpected eval alts found");
        Utils.validate(truthAlts.size() < 2, "Unexpected truth alts found");
        if (!evalAlts.isEmpty() && !truthAlts.isEmpty() && !evalAlts.equals(truthAlts)) {
            // Mixed alt allele case - set truth allele to reference
            truthAlleles = replaceAltAllele(truthAlleles, truthAlts.get(0), refAllele);
        }
        final Genotype newEvalGenotype = new GenotypeBuilder(evalGenotype).alleles(evalAlleles).make();
        final Genotype newTruthGenotype = truthGenotype == null ? null : new GenotypeBuilder(truthGenotype).alleles(truthAlleles).make();
        return getStates(newEvalGenotype, newTruthGenotype);
    }

    private static List<Allele> replaceAltAllele(final List<Allele> list, final Allele search, final Allele replace) {
        Utils.nonNull(search);
        final List<Allele> newList = new ArrayList<>(list.size());
        for (final Allele a : list) {
            if (search.equals(a)) {
                newList.add(replace);
            } else {
                newList.add(a);
            }
        }
        return newList;
    }

    private static List<Allele> getCNVAlleles(final Genotype genotype, final Allele altAllele, final Allele refAllele) {
        if (genotype == null) {
            return Collections.emptyList();
        } else if (!altAllele.equals(Allele.SV_SIMPLE_CNV)) {
            return genotype.getAlleles();
        } else {
            return getCNVGenotypeAllelesFromCopyNumber(genotype, refAllele);
        }
    }


    /**
     * Generates genotype alleles, i.e. for the GT field, for CNVs (DEL and/or DUP). Multi-allelics result in
     * {@link Allele#NO_CALL} alleles.
     * @param genotype  genotype
     * @param refAllele   reference allele
     * @return  alleles for the sample at the site
     * @throws {@link IllegalArgumentException} if the alt allele(s) are not CNV(s)
     */
    private static List<Allele> getCNVGenotypeAllelesFromCopyNumber(final Genotype genotype, final Allele refAllele) {
        Utils.nonNull(refAllele);
        final int expectedCopyNumber = VariantContextGetters.getAttributeAsInt(genotype, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1);
        if (genotype.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT) == null) {
            // Had conflicting genotypes
            return Collections.nCopies(expectedCopyNumber, Allele.NO_CALL);
        }
        final int copyNumber = VariantContextGetters.getAttributeAsInt(genotype, GATKSVVCFConstants.COPY_NUMBER_FORMAT, -1);
        Utils.validateArg(copyNumber >= 0, "Invalid negative copy number: " + copyNumber);
        Utils.validateArg(expectedCopyNumber >= 0, "Invalid expected copy number: " + expectedCopyNumber);
        if (copyNumber == expectedCopyNumber) {
            // Most common in practice - use faster method
            return Collections.nCopies(expectedCopyNumber, refAllele);
        } else if (copyNumber < expectedCopyNumber) {
            // deletion case
            final int numAlt = expectedCopyNumber - copyNumber;
            return CanonicalSVCollapser.makeBiallelicList(Allele.SV_SIMPLE_DEL, refAllele, numAlt, expectedCopyNumber);
        } else {
            // duplication case
            final int numAlt = copyNumber - expectedCopyNumber;
            return CanonicalSVCollapser.makeBiallelicList(Allele.SV_SIMPLE_DUP, refAllele, numAlt, expectedCopyNumber);
        }
    }

    @VisibleForTesting
    protected GenotypeConcordanceStates.TruthState getTruthState(final Genotype g) {
        // Treat non-existent genotype as hom-ref
        if (g == null) {
            return GenotypeConcordanceStates.TruthState.NO_CALL;
        } else if (g.isHomRef()) {
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
