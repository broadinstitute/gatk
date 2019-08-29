package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Describes the results of the AFCalc
 *
 * Only the bare essentials are represented here, as all AFCalc models must return meaningful results for
 * all of these fields.
 *
 * Note that all of the values -- i.e. priors -- are checked now that they are meaningful, which means
 * that users of this code can rely on the values coming out of these functions.
 */
public final class AFCalculationResult {
    // In GVCF mode the STANDARD_CONFIDENCE_FOR_CALLING is 0 by default, and it's nice having this easily-interpretable
    // threshold that says "call anything with any evidence at all."  The problem is that *everything* has at least some evidence,
    // so this would end up putting every site, or at least too many sites, in the gvcf.  Thus this parameter is in place to say
    // that "0" really means "epsilon."
    private static final double EPSILON = 1.0e-10;

    private final double log10PosteriorOfNoVariant;

    private final Map<Allele, Double> log10pRefByAllele;

    /**
     * The AC values for all ALT alleles at the MLE
     */
    private final int[] alleleCountsOfMLE;

    /**
     * The list of alleles actually used in computing the AF
     */
    private final List<Allele> allelesUsedInGenotyping;

    /**
     * Create a results object capability of storing results for calls with up to maxAltAlleles
     */
    public AFCalculationResult(final int[] alleleCountsOfMLE,
                               final List<Allele> allelesUsedInGenotyping,
                               final double log10PosteriorOfNoVariant,
                               final Map<Allele, Double> log10pRefByAllele) {
        Utils.nonNull(alleleCountsOfMLE, "alleleCountsOfMLE cannot be null");
        Utils.nonNull(log10pRefByAllele, "log10pRefByAllele cannot be null");
        Utils.nonNull(allelesUsedInGenotyping, "allelesUsedInGenotyping cannot be null");
        Utils.validateArg(MathUtils.isValidLog10Probability(log10PosteriorOfNoVariant), "log10 posterior must be a valid log probability");

        if ( allelesUsedInGenotyping.isEmpty() ) {
            throw new IllegalArgumentException("allelesUsedInGenotyping must be non-null list of at least 1 value " + allelesUsedInGenotyping);
        }
        if ( alleleCountsOfMLE.length != allelesUsedInGenotyping.size() - 1) {
            throw new IllegalArgumentException("alleleCountsOfMLE.length " + alleleCountsOfMLE.length + " != allelesUsedInGenotyping.size() " + allelesUsedInGenotyping.size());
        }
        if ( log10pRefByAllele.size() != allelesUsedInGenotyping.size() - 1 ) {
            throw new IllegalArgumentException("log10pRefByAllele has the wrong number of elements: log10pRefByAllele " + log10pRefByAllele + " but allelesUsedInGenotyping " + allelesUsedInGenotyping);
        }
        if ( ! allelesUsedInGenotyping.containsAll(log10pRefByAllele.keySet()) ) {
            throw new IllegalArgumentException("log10pRefByAllele doesn't contain all of the alleles used in genotyping: log10pRefByAllele " + log10pRefByAllele + " but allelesUsedInGenotyping " + allelesUsedInGenotyping);
        }

        //make defensive copies of all arguments
        this.alleleCountsOfMLE = alleleCountsOfMLE.clone();
        this.allelesUsedInGenotyping = Collections.unmodifiableList(new ArrayList<>(allelesUsedInGenotyping));

        this.log10PosteriorOfNoVariant = log10PosteriorOfNoVariant;
        this.log10pRefByAllele = Collections.unmodifiableMap(new LinkedHashMap<>(log10pRefByAllele));
    }

    /**
     * Returns a vector with maxAltAlleles values containing AC values at the MLE
     *
     * The values of the ACs for this call are stored in the getAllelesUsedInGenotyping order,
     * starting from index 0 (i.e., the first alt allele is at 0).  The vector is always
     * maxAltAlleles in length, and so only the first getAllelesUsedInGenotyping.size() - 1 values
     * are meaningful.
     *
     * This method returns a copy of the internally-stored array.
     *
     * @return a vector with allele counts, not all of which may be meaningful
     */
    public int[] getAlleleCountsOfMLE() {
        return alleleCountsOfMLE.clone();
    }

    /**
     * Returns the AC of allele a la #getAlleleCountsOfMLE
     *
     * @param allele the allele whose AC we want to know.  Error if its not in allelesUsedInGenotyping
     * @throws IllegalStateException if allele isn't in allelesUsedInGenotyping
     * @return the AC of allele
     */
    public int getAlleleCountAtMLE(final Allele allele) {
        Utils.nonNull(allele);
        Utils.validate( allele.isNonReference(), () -> "Cannot get the alt allele index for reference allele " + allele);
        final int indexInAllAllelesIncludingRef = allelesUsedInGenotyping.indexOf(allele);
        Utils.validateArg(indexInAllAllelesIncludingRef != -1, () -> "could not find allele " + allele + " in " + allelesUsedInGenotyping);
        final int indexInAltAlleles = indexInAllAllelesIncludingRef - 1;
        return alleleCountsOfMLE[indexInAltAlleles];
    }

    /**
     * Get the list of alleles actually used in genotyping, which may be smaller than the actual list of alleles requested
     *
     * @return a non-empty list of alleles used during genotyping, the first of which is the reference allele
     */
    public List<Allele> getAllelesUsedInGenotyping() {
        return allelesUsedInGenotyping;
    }

    public double log10ProbOnlyRefAlleleExists() {
        return log10PosteriorOfNoVariant;
    }

    public double log10ProbVariantPresent() {
        return MathUtils.log10OneMinusPow10(log10PosteriorOfNoVariant);
    }

    @Override
    public String toString() {
        final List<String> byAllele = new LinkedList<>();
        for ( final Allele a : allelesUsedInGenotyping) {
            if (a.isNonReference()) {
                byAllele.add(String.format("%s => MLE %d / posterior %.2f", a, getAlleleCountAtMLE(a), getLog10PosteriorOfAlleleAbsent(a)));
            }
        }
        return String.format("AFCalc%n\t\tlog10PosteriorOfVariant=%.2f%n\t\t%s", log10ProbVariantPresent(), Utils.join("\n\t\t", byAllele));
    }

    /**
     * Are we confident that an allele is present
     */
    public boolean passesThreshold(final Allele allele, final double phredScaleQualThreshold) {
        Utils.nonNull(allele);
        return getLog10PosteriorOfAlleleAbsent(allele) + EPSILON < QualityUtils.qualToErrorProbLog10(phredScaleQualThreshold);
    }

    /**
     * Returns the log10 probability that allele is not segregating
     *
     * Note that this function is p not segregating so that we can store
     * internally the log10 value of AF == 0, which grows very quickly
     * negative and yet has sufficient resolution for high confidence tests.
     * For example, if log10pRef == -100, not an unreasonably high number,
     * if we tried to store log10pNonRef we'd be looking at 1 - 10^-100, which
     * quickly underflows to 1.  So the logic here is backward from what
     * you really want (the p of segregating) but we do that for numerical
     * reasons
     *
     * Unlike the sites-level annotation, this calculation is specific to allele, and can be
     * used to separately determine how much evidence there is that allele is independently
     * segregating as opposed to the site being polymorphic with any allele.  In the bi-allelic
     * case these are obviously the same but for multiple alt alleles there can be lots of
     * evidence for one allele but not so much for any other allele
     *
     * @param allele the allele we're interested in, must be in getAllelesUsedInGenotyping
     * @return the log10 probability that allele is not segregating at this site
     */
    public double getLog10PosteriorOfAlleleAbsent(final Allele allele) {
        Utils.nonNull(allele);
        final Double log10pNonRef = log10pRefByAllele.get(allele);
        Utils.nonNull(log10pNonRef, "Unknown allele " + allele);
        return log10pNonRef;
    }
}