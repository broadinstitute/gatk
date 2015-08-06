package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import org.broadinstitute.hellbender.utils.Utils;

import java.lang.reflect.Modifier;
import java.util.function.Supplier;

/**
 * Enumeration of usable AF calculation, their constraints (i.e. ploidy).
 *
 * Note: this is an enum so that it can be used by the CLI system as an argument.
 *
 * Note that the order these occur in the enum is the order of preference, so
 * the first value is taken over the second when multiple calculations satisfy
 * the needs of the request (i.e., considering ploidy).
 */
public enum AFCalculatorImplementation {

    /** default implementation */
    EXACT_INDEPENDENT(IndependentAllelesDiploidExactAFCalculator::new, 2),

    /** reference implementation of multi-allelic EXACT model.  Extremely slow for many alternate alleles */
    EXACT_REFERENCE(ReferenceDiploidExactAFCalculator::new, 2),

    /** original biallelic exact model, for testing only */
    EXACT_ORIGINAL(OriginalDiploidExactAFCalculator::new, 2, 2),

    /** implementation that supports any sample ploidy.  Currently not available for the HaplotypeCaller */
    EXACT_GENERAL_PLOIDY(GeneralPloidyExactAFCalculator::new);

    /**
     * Special max alt allele count indicating that this maximum is in fact unbound (can be anything).
     */
    public static final int UNBOUND_ALTERNATIVE_ALLELE_COUNT = -1;

    /**
     * Special ploidy constant that indicates that in fact the ploidy is unbound (can be anything).
     */
    public static final int UNBOUND_PLOIDY = -1;

    /**
     * Maximum number of supported alternative alleles.
     */
    private final int maxAltAlleles;

    /**
     * Supported ploidy.
     *
     * This is equal to {@link #UNBOUND_PLOIDY} if the class can handle any ploidy.
     */
    private final int requiredPloidy;

    /**
     * Reference to the default implementation.
     */
    public static final AFCalculatorImplementation DEFAULT = EXACT_INDEPENDENT;

    private final Supplier<AFCalculator> afCalculatorSupplier;

    /**
     * Constructs a new instance given all its properties
     * @param afCalculatorSupplier the calculator class that realizes this implementation.
     * @param requiredPloidy the required ploidy; zero or greater or {@link #UNBOUND_PLOIDY} to indicate that any ploidy is supported.
     * @param maxAltAlleles the maximum alternative alleles; zero or greater or {@link #UNBOUND_ALTERNATIVE_ALLELE_COUNT} to indicate that any maximum number of alternative alleles is supported.
     */
    AFCalculatorImplementation(final Supplier<AFCalculator> afCalculatorSupplier, final int requiredPloidy, final int maxAltAlleles) {
        Utils.nonNull(afCalculatorSupplier);
        this.afCalculatorSupplier = afCalculatorSupplier;
        this.requiredPloidy = requiredPloidy;
        this.maxAltAlleles = maxAltAlleles;
    }

    /**
     * Constructs a new instance leaving ploidy and max-allele count unbound.
     * @param afCalculatorSupplier the calculator class that realizes this implementation.
     */
    AFCalculatorImplementation(final Supplier<AFCalculator> afCalculatorSupplier) {
        this(afCalculatorSupplier,UNBOUND_PLOIDY, UNBOUND_ALTERNATIVE_ALLELE_COUNT);
    }

    /** Constructs a new instance leaving max-allele count unbound.
     * @param afCalculatorSupplier the calculator class that realizes this implementation.
     * @param requiredPloidy the required ploidy; zero or greater or {@link #UNBOUND_PLOIDY} to indicate that any ploidy is supported.
     */
    AFCalculatorImplementation(final Supplier<AFCalculator> afCalculatorSupplier, final int requiredPloidy) {
        this(afCalculatorSupplier,requiredPloidy,UNBOUND_PLOIDY);
    }

    /**
     * Checks whether a given ploidy and max alternative alleles combination is supported or not.
     * @param requestedPloidy the targeted ploidy.
     * @param requestedMaxAltAlleles the targeted max alternative alleles.
     * @return {@code true} iff this calculator implementation satisfies both requirements.
     */
    public boolean usableForParams(final int requestedPloidy, final int requestedMaxAltAlleles) {
        return (requiredPloidy == UNBOUND_PLOIDY || requiredPloidy == requestedPloidy)
                && (maxAltAlleles == UNBOUND_ALTERNATIVE_ALLELE_COUNT || maxAltAlleles >= requestedMaxAltAlleles);
    }

    public AFCalculator newInstance() { return this.afCalculatorSupplier.get(); }

    /**
     * Returns the best (fastest) model give the required ploidy and alternative allele count.
     *
     * @param requiredPloidy required ploidy
     * @param requiredAlternativeAlleleCount required alternative allele count.
     * @param preferred a preferred mode if any. A {@code null} indicate that we should be try to use the default instead.
     * @return never {@code null}
     */
    public static AFCalculatorImplementation bestValue(final int requiredPloidy, final int requiredAlternativeAlleleCount, final AFCalculatorImplementation preferred) {
        final AFCalculatorImplementation preferredValue = preferred == null ? DEFAULT : preferred;
        if (preferredValue.usableForParams(requiredPloidy,requiredAlternativeAlleleCount)) {
            return preferredValue;
        }
        if (EXACT_INDEPENDENT.usableForParams(requiredPloidy,requiredAlternativeAlleleCount)) {
            return EXACT_INDEPENDENT;
        }
        if (EXACT_REFERENCE.usableForParams(requiredPloidy,requiredAlternativeAlleleCount)) {  //TODO: this seems to be dead code. EXACT_REFERENCE will always lose to EXACT_INDEPENDENT.
            return EXACT_REFERENCE;
        }
        return EXACT_GENERAL_PLOIDY;
    }

    /**
     * Returns the value that corresponds to a given implementation calculator class.
     *
     * @param clazz the target class.
     *
     * @throws IllegalArgumentException if {@code clazz} is {@code null} or if it is abstract.
     *
     * @return never {@code null}.
     */
    public static AFCalculatorImplementation fromCalculatorClass(final Class<? extends AFCalculator> clazz) {
        Utils.nonNull(clazz, "input class cannot be null");
        Utils.validateArg(!Modifier.isAbstract(clazz.getModifiers()), "class " + clazz.getCanonicalName() + " should not be abstract");

        //Using iteration instead of a static map to avoid static state.
        for (final AFCalculatorImplementation impl : AFCalculatorImplementation.values()){
            if (clazz.equals(impl.newInstance().getClass())){
                return impl;
            }
        }
        throw new IllegalArgumentException("Attempt to retrieve AFCalculatorImplementation instance from a non-registered calculator class " + clazz.getName());
    }
}
