package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingEngine;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public final class AFCalculatorTestBuilder {
    final static Allele A = Allele.create("A", true);
    final static Allele C = Allele.create("C");
    final static Allele G = Allele.create("G");
    final static Allele T = Allele.create("T");
    final static Allele AA = Allele.create("AA");
    final static Allele AT = Allele.create("AT");
    final static Allele AG = Allele.create("AG");

    static int sampleNameCounter = 0;

    final int nSamples;
    final int numAltAlleles;
    final AFCalculatorImplementation modelType;
    final PriorType priorType;

    public AFCalculatorTestBuilder(final int nSamples, final int numAltAlleles,
                                   final AFCalculatorImplementation modelType, final PriorType priorType) {
        this.nSamples = nSamples;
        this.numAltAlleles = numAltAlleles;
        this.modelType = modelType;
        this.priorType = priorType;
    }

    @Override
    public String toString() {
        return String.format("AFCalcTestBuilder nSamples=%d nAlts=%d model=%s prior=%s", nSamples, numAltAlleles, modelType, priorType);
    }

    public enum PriorType {
        flat,
        human
    }

    public int getNumAltAlleles() {
        return numAltAlleles;
    }

    public int getnSamples() {
        return nSamples;
    }

    public AFCalculator makeModel() {
        return createCalculator(modelType, nSamples, getNumAltAlleles(), HomoSapiensConstants.DEFAULT_PLOIDY);
    }

    /**
     * Create a new AFCalc
     *
     * @param implementation the calculation to use
     * @param nSamples the number of samples we'll be using
     * @param maxAltAlleles the max. alt alleles to consider for SNPs
     * @param ploidy the sample ploidy.  Must be consistent with the implementation
     *
     * @return an initialized AFCalc
     */
    private static AFCalculator createCalculator(final AFCalculatorImplementation implementation, final int nSamples, final int maxAltAlleles, final int ploidy) {
        if ( implementation == null ) {
            throw new IllegalArgumentException("Calculation cannot be null");
        }
        if ( nSamples < 0 ) {
            throw new IllegalArgumentException("nSamples must be greater than zero " + nSamples);
        }
        if ( maxAltAlleles < 1 ) {
            throw new IllegalArgumentException("maxAltAlleles must be greater than zero " + maxAltAlleles);
        }
        if ( ploidy < 1 ) {
            throw new IllegalArgumentException("sample ploidy must be greater than zero " + ploidy);
        }

        if ( ! implementation.usableForParams(ploidy, maxAltAlleles) ) {
            throw new IllegalArgumentException("AFCalc " + implementation + " does not support requested ploidy " + ploidy);
        }

        return implementation.newInstance();
    }

    public double[] makePriors() {
        final int nPriorValues = 2*nSamples+1;

        switch ( priorType ) {
            case flat:
                return MathUtils.normalizeLog10(new double[nPriorValues]);  // flat priors

            //TODO break dependency with human... avoid special reference to this species.
            case human:
                final double[] humanPriors = new double[nPriorValues];
                GenotypingEngine.computeAlleleFrequencyPriors(nPriorValues - 1, humanPriors, 0.001, new ArrayList<>());
                return humanPriors;
            default:
                throw new RuntimeException("Unexpected type " + priorType);
        }
    }

    public VariantContext makeACTest(final List<Integer> ACs, final int nNonInformative, final int nonTypePL) {
        return makeACTest(ArrayUtils.toPrimitive(ACs.toArray(new Integer[]{})), nNonInformative, nonTypePL);
    }

    public VariantContext makeACTest(final int[] ACs, final int nNonInformative, final int nonTypePL) {
        final int nChrom = nSamples * 2;

        final int[] nhet = new int[numAltAlleles];
        final int[] nhomvar = new int[numAltAlleles];

        for ( int i = 0; i < ACs.length; i++ ) {
            final double p = ACs[i] / (1.0 * nChrom);
            nhomvar[i] = (int) Math.floor((nSamples - nNonInformative) * p * p);
            nhet[i] = ACs[i] - 2 * nhomvar[i];

            if ( nhet[i] < 0 )
                throw new IllegalStateException("Bug! nhet[i] < 0");
        }

        final long calcAC = MathUtils.sum(nhet) + 2 * MathUtils.sum(nhomvar);
        if ( calcAC != MathUtils.sum(ACs) )
            throw new IllegalStateException("calculated AC " + calcAC + " not equal to desired AC " + Utils.join(",", ACs));

        return makeACTest(nhet, nhomvar, nNonInformative, nonTypePL);
    }

    public VariantContext makeACTest(final int[] nhet, final int[] nhomvar, final int nNonInformative, final int nonTypePL) {
        List<Genotype> samples = new ArrayList<>(nSamples);

        for ( int altI = 0; altI < nhet.length; altI++ ) {
            for ( int i = 0; i < nhet[altI]; i++ )
                samples.add(makePL(GenotypeType.HET, nonTypePL, altI+1));
            for ( int i = 0; i < nhomvar[altI]; i++ )
                samples.add(makePL(GenotypeType.HOM_VAR, nonTypePL, altI+1));
        }

        final Genotype nonInformative = makeNonInformative();
        samples.addAll(Collections.nCopies(nNonInformative, nonInformative));

        final int nRef = Math.max((int) (nSamples - nNonInformative - MathUtils.sum(nhet) - MathUtils.sum(nhomvar)), 0);
        samples.addAll(Collections.nCopies(nRef, makePL(GenotypeType.HOM_REF, nonTypePL, 0)));

        samples = samples.subList(0, nSamples);

        if ( samples.size() > nSamples )
            throw new IllegalStateException("too many samples");

        VariantContextBuilder vcb = new VariantContextBuilder("x", "1", 1, 1, getAlleles());
        vcb.genotypes(samples);
        return vcb.make();
    }

    public List<Allele> getAlleles() {
        return Arrays.asList(A, C, G, T, AA, AT, AG).subList(0, numAltAlleles+1);
    }

    public List<Allele> getAlleles(final GenotypeType type, final int altI) {
        switch (type) {
            case HOM_REF: return Arrays.asList(getAlleles().get(0), getAlleles().get(0));
            case HET:     return Arrays.asList(getAlleles().get(0), getAlleles().get(altI));
            case HOM_VAR: return Arrays.asList(getAlleles().get(altI), getAlleles().get(altI));
            default: throw new IllegalArgumentException("Unexpected type " + type);
        }
    }

    public Genotype makePL(final List<Allele> expectedGT, int ... pls) {
        GenotypeBuilder gb = new GenotypeBuilder("sample" + sampleNameCounter++);
        gb.alleles(expectedGT);
        gb.PL(pls);
        return gb.make();
    }

    private int numPLs() {
        return GenotypeLikelihoods.numLikelihoods(numAltAlleles + 1, 2);
    }

    public Genotype makeNonInformative() {
        final int[] nonInformativePLs = new int[GenotypeLikelihoods.numLikelihoods(numAltAlleles, 2)];
        return makePL(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), nonInformativePLs);
    }

    public Genotype makePL(final GenotypeType type, final int nonTypePL, final int altI) {
        GenotypeBuilder gb = new GenotypeBuilder("sample" + sampleNameCounter++);
        gb.alleles(getAlleles(type, altI));

        final int[] pls = new int[numPLs()];
        Arrays.fill(pls, nonTypePL);

        int index = 0;
        switch ( type ) {
            case HOM_REF: index = GenotypeLikelihoods.calculatePLindex(0, 0); break;
            case HET:     index = GenotypeLikelihoods.calculatePLindex(0, altI); break;
            case HOM_VAR: index = GenotypeLikelihoods.calculatePLindex(altI, altI); break;
        }
        pls[index] = 0;
        gb.PL(pls);

        return gb.make();
    }
}