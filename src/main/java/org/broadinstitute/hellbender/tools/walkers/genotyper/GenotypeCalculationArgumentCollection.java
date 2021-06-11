package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.io.Serializable;

public final class GenotypeCalculationArgumentCollection implements Serializable, Cloneable {
    private static final long serialVersionUID = 1L;

    public static final String SUPPORTING_CALLSET_LONG_NAME = "population-callset";
    public static final String SUPPORTING_CALLSET_SHORT_NAME = "population";
    public static final String NUM_REF_SAMPLES_LONG_NAME = "num-reference-samples-if-no-call";
    public static final String CALL_CONFIDENCE_LONG_NAME = "standard-min-confidence-threshold-for-calling";
    public static final String CALL_CONFIDENCE_SHORT_NAME = "stand-call-conf";
    public static final String MAX_ALTERNATE_ALLELES_LONG_NAME = "max-alternate-alleles";
    public static final String MAX_GENOTYPE_COUNT_LONG_NAME = "max-genotype-count";
    public static final String SAMPLE_PLOIDY_SHORT_NAME = "ploidy";
    public static final String SAMPLE_PLOIDY_LONG_NAME = "sample-ploidy";

    public static final double DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING = 30.0;
    public static final int DEFAULT_MAX_ALTERNATE_ALLELES = 6;
    public static final int DEFAULT_MAX_GENOTYPE_COUNT = 1024;

    @Argument(fullName="use-posteriors-to-calculate-qual", shortName="gp-qual", optional = true, doc = "if available, use the genotype posterior probabilities to calculate the site QUAL")
    public boolean usePosteriorProbabilitiesToCalculateQual = false;

    /**
     * Creates a GenotypeCalculationArgumentCollection with default values.
     */
    public GenotypeCalculationArgumentCollection() {}

    /**
     * Creates a new GenotypeCalculationArgumentCollection with the values from other instance.
     * <p>
     *     Changes in direct field members of the returned object won't affect the values in the original argument
     *     collection.
     * </p>
     */
    @Override
    public GenotypeCalculationArgumentCollection clone() {
        try {
            return (GenotypeCalculationArgumentCollection) super.clone();
        } catch (final CloneNotSupportedException e) {
            throw new GATKException("this line of code should not be reached");
        }
    }

    @Advanced
    @Argument(fullName = "dont-use-dragstr-priors",
              doc      = "Forfeit the use of the DRAGstr model to calculate genotype priors. " +
                         "This argument does not have any effect in the absence of DRAGstr model parameters (--dragstr-model-params)", optional = true)
    public boolean dontUseDragstrPriors = false;

    /**
     * As of version 4.1.0.0, this argument is no longer needed because the new qual score is now on by default. See GATK 3.3 release notes for more details.
     */
    @Deprecated
    @Argument(fullName = "use-new-qual-calculator", shortName = "new-qual", doc = "Use the new AF model instead of the so-called exact model", optional = true)
    public boolean useNewAFCalculator = true;

    /**
     * Depending on the value of the --max_alternate_alleles argument, we may genotype only a fraction of the alleles being sent on for genotyping.
     * Using this argument instructs the genotyper to annotate (in the INFO field) the number of alternate alleles that were originally discovered at the site.
     */
    @Argument(fullName = "annotate-with-num-discovered-alleles", doc = "If provided, we will annotate records with the number of alternate alleles that were discovered (but not necessarily genotyped) at a given site", optional = true)
    public boolean ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED = false;

    /**
     * The expected heterozygosity value used to compute prior probability that a locus is non-reference.
     *
     * The default priors are for provided for humans:
     *
     * het = 1e-3
     *
     * which means that the probability of N samples being hom-ref at a site is:
     *
     * 1 - sum_i_2N (het / i)
     *
     * Note that heterozygosity as used here is the population genetics concept:
     *
     * http://en.wikipedia.org/wiki/Zygosity#Heterozygosity_in_population_genetics
     *
     * That is, a hets value of 0.01 implies that two randomly chosen chromosomes from the population of organisms
     * would differ from each other (one being A and the other B) at a rate of 1 in 100 bp.
     *
     * Note that this quantity has nothing to do with the likelihood of any given sample having a heterozygous genotype,
     * which in the GATK is purely determined by the probability of the observed data P(D | AB) under the model that there
     * may be a AB het genotype.  The posterior probability of this AB genotype would use the het prior, but the GATK
     * only uses this posterior probability in determining the prob. that a site is polymorphic.  So changing the
     * het parameters only increases the chance that a site will be called non-reference across all samples, but
     * doesn't actually change the output genotype likelihoods at all, as these aren't posterior probabilities at all.
     *
     * The quantity that changes whether the GATK considers the possibility of a het genotype at all is the ploidy,
     * which determines how many chromosomes each individual in the species carries.
     */
    @Argument(fullName = "heterozygosity", doc = "Heterozygosity value used to compute prior likelihoods for any locus.  See the GATKDocs for full details on the meaning of this population genetics concept", optional = true)
    public Double snpHeterozygosity = HomoSapiensConstants.SNP_HETEROZYGOSITY;

    /**
     * This argument informs the prior probability of having an indel at a site.
     */
    @Argument(fullName = "indel-heterozygosity", doc = "Heterozygosity for indel calling.  See the GATKDocs for heterozygosity for full details on the meaning of this population genetics concept", optional = true)
    public double indelHeterozygosity = HomoSapiensConstants.INDEL_HETEROZYGOSITY;

    /**
     * The standard deviation of the distribution of alt allele fractions.  The above heterozygosity parameters give the
     * *mean* of this distribution; this parameter gives its spread.
     */
    @Argument(fullName = "heterozygosity-stdev", doc = "Standard deviation of heterozygosity for SNP and indel calling.", optional = true)
    public double heterozygosityStandardDeviation = 0.01;

    /**
     * The minimum phred-scaled confidence threshold at which variants should be called. Only variant sites with QUAL equal
     * or greater than this threshold will be called. Note that since version 3.7, we no longer differentiate high confidence
     * from low confidence calls at the calling step. The default call confidence threshold is set low intentionally to achieve
     * high sensitivity, which will allow false positive calls as a side effect. Be sure to perform some kind of filtering after
     * calling to reduce the amount of false positives in your final callset. Note that when HaplotypeCaller is used in GVCF mode
     * (using either -ERC GVCF or -ERC BP_RESOLUTION) the call threshold is automatically set to zero. Call confidence thresholding
     * will then be performed in the subsequent GenotypeGVCFs command.
     *
     * Note that the default was changed from 10.0 to 30.0 in version 4.1.0.0 to accompany the switch to use the the new quality score by default.
     */
    @Argument(fullName = CALL_CONFIDENCE_LONG_NAME, shortName = CALL_CONFIDENCE_SHORT_NAME, doc = "The minimum phred-scaled confidence threshold at which variants should be called", optional = true)
    public double STANDARD_CONFIDENCE_FOR_CALLING = DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING;

    /**
     * If there are more than this number of alternate alleles presented to the genotyper (either through discovery or GENOTYPE_GIVEN ALLELES),
     * then only this many alleles will be used.  Note that genotyping sites with many alternate alleles is both CPU and memory intensive and it
     * scales exponentially based on the number of alternate alleles.  Unless there is a good reason to change the default value, we highly recommend
     * that you not play around with this parameter.
     *
     * See also {@link #MAX_GENOTYPE_COUNT}.
     */
    @Advanced
    @Argument(fullName = MAX_ALTERNATE_ALLELES_LONG_NAME, doc = "Maximum number of alternate alleles to genotype", optional = true)
    public int MAX_ALTERNATE_ALLELES = DEFAULT_MAX_ALTERNATE_ALLELES;

    /**
     * If there are more than this number of genotypes at a locus presented to the genotyper, then only this many genotypes will be used.
     * The possible genotypes are simply different ways of partitioning alleles given a specific ploidy asumption.
     * Therefore, we remove genotypes from consideration by removing alternate alleles that are the least well supported.
     * The estimate of allele support is based on the ranking of the candidate haplotypes coming out of the graph building step.
     * Note that the reference allele is always kept.
     *
     * Note that genotyping sites with large genotype counts is both CPU and memory intensive.
     * Unless there is a good reason to change the default value, we highly recommend that you not play around with this parameter.
     *
     * The maximum number of alternative alleles used in the genotyping step will be the lesser of the two:
     * 1. the largest number of alt alleles, given ploidy, that yields a genotype count no higher than {@link #MAX_GENOTYPE_COUNT}
     * 2. the value of {@link #MAX_ALTERNATE_ALLELES}
     *
     * See also {@link #MAX_ALTERNATE_ALLELES}.
     */
    @Advanced
    @Argument(fullName = MAX_GENOTYPE_COUNT_LONG_NAME, doc = "Maximum number of genotypes to consider at any site", optional = true)
    public int MAX_GENOTYPE_COUNT = DEFAULT_MAX_GENOTYPE_COUNT;

    /**
     *   Sample ploidy - equivalent to number of chromosomes per pool. In pooled experiments this should be = # of samples in pool * individual sample ploidy
     */
    @Argument(shortName = SAMPLE_PLOIDY_SHORT_NAME, fullName = SAMPLE_PLOIDY_LONG_NAME, doc="Ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy).", optional=true)
    public int samplePloidy = HomoSapiensConstants.DEFAULT_PLOIDY;

    /**
     * Supporting external panel. Allele counts from this panel (taken from AC,AN or MLEAC,AN or raw genotypes) will
     * be used to inform the frequency distribution underlying the genotype priors. These files must be VCF 4.2 spec or later.
     * Note that unlike CalculateGenotypePosteriors, HaplotypeCaller only allows one supporting callset.
     */
    @Argument(fullName=SUPPORTING_CALLSET_LONG_NAME, shortName=SUPPORTING_CALLSET_SHORT_NAME, doc="Callset to use in calculating genotype priors", optional=true)
    public FeatureInput<VariantContext> supportVariants = null;

    /**
     * When a variant is not seen in any panel, this argument controls whether to infer (and with what effective strength)
     * that only reference alleles were observed at that site. E.g. "If not seen in 1000Genomes, treat it as AC=0,
     * AN=2000".
     */
    @Argument(fullName= NUM_REF_SAMPLES_LONG_NAME,doc="Number of hom-ref genotypes to infer at sites not present in a panel",optional=true)
    public int numRefIfMissing = 0;

    @Argument(fullName= "genotype-assignment-method", shortName = "gam", doc = "How we assign genotypes", optional = true)
    public GenotypeAssignmentMethod genotypeAssignmentMethod = GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN;
}
