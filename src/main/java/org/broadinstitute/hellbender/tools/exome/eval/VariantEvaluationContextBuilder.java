package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.exome.CopyNumberTriStateAllele;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link VariantEvaluationContext} building helping class.
 */
public class VariantEvaluationContextBuilder extends VariantContextBuilder {

    private int[] truthAC = new int[CopyNumberTriStateAllele.ALL_ALLELES.size()];
    private int[] callsAC = new int[CopyNumberTriStateAllele.ALL_ALLELES.size()];
    private VariantContext truth;
    private List<VariantContext> calls;

    public VariantEvaluationContextBuilder() {
        super();
    }

    public VariantEvaluationContextBuilder(final VariantContext init) {
        super(init);
        alleles(init.getAlleles());
        genotypes(init.getGenotypes());
        if (init instanceof VariantEvaluationContext) {
            final VariantEvaluationContext casted = (VariantEvaluationContext) init;
            evidence(casted.getTruthVariantContext(), casted.getCallsVariantContexts());
        }
    }

    @Override
    public VariantEvaluationContextBuilder alleles(final Collection<Allele> alleles) {
        ParamUtils.noNulls(alleles, "the input cannot contain null alleles");
        validateAlleleStrings(alleles.stream().map(Allele::getDisplayString).collect(Collectors.toList()));
        if (!alleles.iterator().next().isReference()) {
            throw new IllegalArgumentException("the first allele must be a reference allele");
        }
        super.alleles(alleles);
        return this;
    }

    private void validateAlleleStrings(final List<String> allelesString) {
        ParamUtils.noNulls(allelesString, "the input cannot contain null strings");
        if (allelesString.isEmpty()) {
            throw new IllegalArgumentException("the allele list cannot be empty");
        } if (!allelesString.get(0).equals(CopyNumberTriStateAllele.REF.getDisplayString())) {
            throw new IllegalArgumentException("the first allele must be the reference allele: " + CopyNumberTriStateAllele.REF);
        }
        final Set<CopyNumberTriStateAllele> allelesFound = new HashSet<>();
        allelesFound.add(CopyNumberTriStateAllele.REF);
        for (int i = 1; i < allelesString.size(); i++) {
            final String alleleString = allelesString.get(i);
            final CopyNumberTriStateAllele allele;
            try {
                allele = CopyNumberTriStateAllele.valueOf(Allele.create(alleleString, false));
            } catch(final IllegalArgumentException ex) {
                throw new IllegalArgumentException(String.format("unknown allele with string: '%s'", alleleString));
            }
            if (!allelesFound.add(allele)) {
                throw new IllegalArgumentException("an allele cannot be listed twice.");
            }
        }
    }

    @Override
    public VariantEvaluationContextBuilder alleles(final List<String> alleleStrings) {
        validateAlleleStrings(alleleStrings);
        super.alleles(alleleStrings);
        return this;
    }

    @Override
    public VariantEvaluationContextBuilder alleles(final String... alleleStrings) {
        validateAlleleStrings(Arrays.asList(alleleStrings));
        super.alleles(alleleStrings);
        return this;
    }

    @Override
    public List<Allele> getAlleles() {
        final List<Allele> result = super.getAlleles();
        if (result == null) {
            throw new IllegalStateException("you must set the alleles before calling getAlleles");
        }
        return result;
    }

    @Override
    public VariantEvaluationContextBuilder genotypes(final GenotypesContext genotypes) {
        super.genotypes(genotypes);
        calculateAFs(genotypes);
        return this;
    }

    @Override
    public VariantEvaluationContextBuilder genotypesNoValidation(final GenotypesContext genotypes) {
        throw new UnsupportedOperationException("non-validated genotype assignation is not available");
    }

    @Override
    public VariantEvaluationContextBuilder genotypes(final Collection<Genotype> genotypes) {
        super.genotypes(genotypes);
        calculateAFs(genotypes);
        return this;
    }

    @Override
    public VariantEvaluationContextBuilder genotypes(final Genotype ... genotypes) {
        super.genotypes(genotypes);
        calculateAFs(Arrays.asList(genotypes));
        return this;
    }

    @Override
    public VariantEvaluationContextBuilder noGenotypes() {
        super.noGenotypes();
        Arrays.fill(truthAC, 0);
        Arrays.fill(callsAC, 0);
        return this;
    }

    private void calculateAFs(final Iterable<Genotype> genotypes) {
        final int[] truthAC = new int[CopyNumberTriStateAllele.ALL_ALLELES.size()];
        final int[] callsAC = new int[CopyNumberTriStateAllele.ALL_ALLELES.size()];
        for (final Genotype genotype : genotypes) {
            final List<Allele> alleles = genotype.getAlleles();
            if (alleles.size() > 1) {
                throw new GATKException("unexpected CNV genotype ploidy: " + alleles.size());
            } else if (!alleles.isEmpty() && alleles.get(0).isCalled()) {
                final int index = CopyNumberTriStateAllele.valueOf(alleles.get(0)).index();
                callsAC[index]++;
            }
            final String truthGT = String.valueOf(genotype.getExtendedAttribute(VariantEvaluationContext.TRUTH_GENOTYPE_KEY, VCFConstants.MISSING_VALUE_v4));
            final int truthAlleleIndex = truthGT.equals(VCFConstants.MISSING_VALUE_v4) ? -1 : Integer.parseInt(truthGT);
            if (truthAlleleIndex >= 0) {
                final List<Allele> contextAlleles = getAlleles();
                if (truthAlleleIndex >= contextAlleles.size()) {
                    throw new GATKException("unexpected CNV truth genotype makes reference to a non-existent allele: " + truthGT);
                }
                truthAC[CopyNumberTriStateAllele.valueOf(contextAlleles.get(truthAlleleIndex)).index()]++;
            }
            genotype.getAllele(0);
        }
        this.truthAC = truthAC;
        this.callsAC = callsAC;
    }

    @Override
    public VariantEvaluationContext make(final boolean leaveItModifiable) {
        final List<Allele> contextAlleles = getAlleles();
        final int[] alleleToCNVAlleleOrdinal = contextAlleles.stream()
                .mapToInt(a -> CopyNumberTriStateAllele.valueOf(a).index()).toArray();
        final int callsAN = (int) MathUtils.sum(this.callsAC);
        final int truthAN = (int) MathUtils.sum(this.truthAC);
        attribute(VariantEvaluationContext.CALLS_ALLELE_NUMBER_KEY, callsAN);
        attribute(VariantEvaluationContext.TRUTH_ALLELE_NUMBER_KEY, truthAN);
        if (callsAN > 0) {
            final double[] callsAF = IntStream.range(1, contextAlleles.size())
                    .mapToDouble(i -> callsAC[alleleToCNVAlleleOrdinal[i]] / (double) callsAN).toArray();
            attribute(VariantEvaluationContext.CALLS_ALLELE_FREQUENCY_KEY, callsAF);
        }
        if (truthAN > 0) {
            final double[] truthAF = IntStream.range(1, contextAlleles.size())
                    .mapToDouble(i -> truthAC[alleleToCNVAlleleOrdinal[i]] / (double) truthAN).toArray();
            attribute(VariantEvaluationContext.TRUTH_ALLELE_FREQUENCY_KEY, truthAF);
        }
        final VariantEvaluationContext context = new VariantEvaluationContext(super.make(leaveItModifiable));
        context.setEvidence(truth, calls);
        return context;
    }

    @Override
    public VariantEvaluationContext make() {
        return make(false);
    }

    public void evidence(final VariantContext truth, final List<VariantContext> calls) {
        this.truth = truth;
        this.calls = calls;
    }
}
