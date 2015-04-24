package org.broadinstitute.hellbender.tools.picard.vcf.concordance;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * Class that holds detail metrics about Genotype Concordance
 *
 * @author George Grant
 */
public final class GenotypeConcordanceDetailMetrics extends MetricBase {
    /** The type of the event (i.e. either SNP or INDEL) */
    public VariantContext.Type VARIANT_TYPE;

    /** The name of the 'truth' sample */
    public String TRUTH_SAMPLE;

    /** The name of the 'call' sample */
    public String CALL_SAMPLE;

    /** The state of the 'truth' sample (i.e. HOM_REF, HET_REF_VAR1, HET_VAR1_VAR2...) */
    public GenotypeConcordanceStates.TruthState TRUTH_STATE;

    /** The state of the 'call' sample (i.e. HOM_REF, HET_REF_VAR1...) */
    public GenotypeConcordanceStates.CallState CALL_STATE;

    /** The number of events of type TRUTH_STATE and CALL_STATE for the EVENT_TYPE and SAMPLEs */
    public long COUNT;

    /** The list of contingency table values (TP, TN, FP, FN) that are deduced from the truth/call state comparison, given the reference.
     *
     * In general, we are comparing two sets of alleles.  Therefore, we can have zero or more contingency table values represented in one comparison.  For example, if the truthset is
     * a heterozygous call with both alleles non-reference (HET_VAR1_VAR2), and the callset is a heterozygous call with both alleles non-reference with one of the alternate alleles
     * matching an alternate allele in the callset, we would have a true positive, false positive, and false negative.  The true positive is from the matching alternate alleles, the
     * false positive is the alternate allele found in the callset but not found in the truthset, and the false negative is the alternate in the truthset not found in the callset.
     *
     * We also include a true negative in cases where the reference allele is found in both the truthset and callset.
     */
    public String CONTINGENCY_VALUES;
}
