package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

@Analysis(description = "Counts different classes of variants in the sample")
public class CountVariants extends VariantEvaluator implements StandardEval {
    public CountVariants(VariantEvalEngine engine) {
        super(engine);
    }

    // the following fields are in output order:

    // basic counts on various rates found
    @DataPoint(description = "Number of processed loci", format = "%d")
    public long nProcessedLoci = 0;
    @DataPoint(description = "Number of called loci", format = "%d")
    public long nCalledLoci = 0;
    @DataPoint(description = "Number of reference loci", format = "%d")
    public long nRefLoci = 0;
    @DataPoint(description = "Number of variant loci", format = "%d")
    public long nVariantLoci = 0;

    // the following two calculations get set in the finalizeEvaluation
    @DataPoint(description = "Variants per loci rate", format = "%.8f")
    public double variantRate = 0;
    @DataPoint(description = "Number of variants per base", format = "%.8f")
    public double variantRatePerBp = 0;

    @DataPoint(description = "Number of snp loci", format = "%d")
    public long nSNPs = 0;
    @DataPoint(description = "Number of mnp loci", format = "%d")
    public long nMNPs = 0;
    @DataPoint(description = "Number of insertions", format = "%d")
    public long nInsertions = 0;
    @DataPoint(description = "Number of deletions", format = "%d")
    public long nDeletions = 0;
    @DataPoint(description = "Number of complex indels", format = "%d")
    public long nComplex = 0;
    @DataPoint(description = "Number of symbolic events", format = "%d")
    public long nSymbolic = 0;

    @DataPoint(description = "Number of mixed loci (loci that can't be classified as a SNP, Indel or MNP)", format = "%d")
    public long nMixed = 0;

    @DataPoint(description = "Number of no calls loci", format = "%d")
    public long nNoCalls = 0;
    @DataPoint(description = "Number of het loci", format = "%d")
    public long nHets = 0;
    @DataPoint(description = "Number of hom ref loci", format = "%d")
    public long nHomRef = 0;
    @DataPoint(description = "Number of hom var loci", format = "%d")
    public long nHomVar = 0;
    @DataPoint(description = "Number of singletons", format = "%d")
    public long nSingletons = 0;
    @DataPoint(description = "Number of derived homozygotes", format = "%d")
    public long nHomDerived = 0;

    // calculations that get set in the finalizeEvaluation method
    @DataPoint(description = "heterozygosity per locus rate", format = "%.2e")
    public double heterozygosity = 0;
    @DataPoint(description = "heterozygosity per base pair", format = "%.2f")
    public double heterozygosityPerBp = 0;
    @DataPoint(description = "heterozygosity to homozygosity ratio", format = "%.2f")
    public double hetHomRatio = 0;
    @DataPoint(description = "indel rate (insertion count + deletion count)", format = "%.2e")
    public double indelRate = 0;
    @DataPoint(description = "indel rate per base pair", format = "%.2f")
    public double indelRatePerBp = 0;
    @DataPoint(description = "insertion  to deletion ratio", format = "%.2f")
    public double insertionDeletionRatio = 0;
    
    private double perLocusRate(long n) {
        return rate(n, nProcessedLoci);
    }

    private long perLocusRInverseRate(long n) {
        return inverseRate(n, nProcessedLoci);
    }


    public int getComparisonOrder() {
        return 1;   // we only need to see each eval track
    }

    @Override
    public void update1(final VariantContext vc, final VariantEvalContext context) {
        nCalledLoci++;

        // Note from Eric:
        // This is really not correct.  What we really want here is a polymorphic vs. monomorphic count (i.e. on the Genotypes).
        // So in order to maintain consistency with the previous implementation (and the intention of the original author), I've
        // added in a proxy check for monomorphic status here.
        // Protect against the case when vc only has no-calls too - can happen if we stratify by sample and sample as a single no-call.
       if ( getEngine().getVariantEvalArgs().ignoreAC0Sites() && vc.isMonomorphicInSamples() ) {
            nRefLoci++;
        } else {
             switch (vc.getType()) {
                case NO_VARIATION:
                    // shouldn't get here
                    break;
                case SNP:
                    nVariantLoci++;
                    nSNPs++;
                    if (variantWasSingleton(vc)) nSingletons++;
                    break;
                case MNP:
                    nVariantLoci++;
                    nMNPs++;
                    if (variantWasSingleton(vc)) nSingletons++;
                    break;
                case INDEL:
                    nVariantLoci++;
                    if (vc.isSimpleInsertion())
                        nInsertions++;
                    else if (vc.isSimpleDeletion())
                        nDeletions++;
                    else
                        nComplex++;
                    break;
                case MIXED:
                    nVariantLoci++;
                    nMixed++;
                    break;
                case SYMBOLIC:
                    nSymbolic++;
                    break;
                default:
                    throw new GATKException("Unexpected VariantContext type " + vc.getType());
            }
        }

        // these operations are ordered to ensure that we don't get the base string of the ref unless we need it
        final String aaStr = vc.hasAttribute("ANCESTRALALLELE") ? vc.getAttributeAsString("ANCESTRALALLELE", null).toUpperCase() : null;
        final String refStr = aaStr != null ? vc.getReference().getBaseString().toUpperCase() : null;

        // ref  aa  alt  class
        // A    C   A    der homozygote
        // A    C   C    anc homozygote

        // A    A   A    ref homozygote
        // A    A   C
        // A    C   A
        // A    C   C

        for (final Genotype g : vc.getGenotypes()) {
            final String altStr = vc.getAlternateAlleles().size() > 0 ? vc.getAlternateAllele(0).getBaseString().toUpperCase() : null;

            switch (g.getType()) {
                case NO_CALL:
                    nNoCalls++;
                    break;
                case HOM_REF:
                    nHomRef++;

                    if ( aaStr != null && altStr != null && !refStr.equalsIgnoreCase(aaStr) ) {
                        nHomDerived++;
                    }

                    break;
                case HET:
                    nHets++;
                    break;
                case HOM_VAR:
                    nHomVar++;

                    if ( aaStr != null && altStr != null && !altStr.equalsIgnoreCase(aaStr) ) {
                        nHomDerived++;
                    }

                    break;
                case MIXED:
                    break;
                case UNAVAILABLE:
                    break;
                default:
                    throw new GATKException("BUG: Unexpected genotype type: " + g);
            }
        }
    }

    public void finalizeEvaluation() {
        nProcessedLoci = getEngine().getnProcessedLoci();
        variantRate = perLocusRate(nVariantLoci);
        variantRatePerBp = perLocusRInverseRate(nVariantLoci);
        heterozygosity = perLocusRate(nHets);
        heterozygosityPerBp = perLocusRInverseRate(nHets);
        hetHomRatio = ratio(nHets, nHomVar);
        indelRate = perLocusRate(nDeletions + nInsertions + nComplex);
        indelRatePerBp = perLocusRInverseRate(nDeletions + nInsertions + nComplex);
        insertionDeletionRatio = ratio(nInsertions, nDeletions);
    }

    @Override
    public boolean requiresTerritoryToBeSpecified() {
        return true;
    }
}