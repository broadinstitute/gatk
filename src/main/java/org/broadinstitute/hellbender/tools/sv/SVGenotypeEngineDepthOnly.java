package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.NaturalLogUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SVGenotypeEngineDepthOnly extends SVGenotypeEngine {

    public List<VCFHeaderLine> getHeaderLines() {
        final List<VCFHeaderLine> lines = new ArrayList<>();
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY, true));
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY, true));
        lines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_PL_KEY, true));
        lines.add(new VCFFormatHeaderLine(GATKSVVCFConstants.COPY_NUMBER_FIELD, 1, VCFHeaderLineType.Integer, "Copy number"));
        return lines;
    }

    public VariantContext genotypeVariant(final VariantContext variant) {
        final StructuralVariantType svType = variant.getStructuralVariantType();
        Utils.validate(SVGenotypeEngine.CNV_TYPES.contains(svType), "Variant " + variant.getID() + " is not a CNV");
        Utils.validate(SVGenotypeEngineFromModel.isDepthOnlyVariant(variant), "Variant " + variant.getID() + " is not depth-only");
        final List<Genotype> genotypes = variant.getGenotypes().stream()
                .map(g -> buildGenotype(g, svType))
                .collect(Collectors.toList());
        final VariantContextBuilder builder = new VariantContextBuilder(variant);
        builder.genotypes(genotypes);
        return builder.make();
    }

    private Genotype buildGenotype(final Genotype genotype, final StructuralVariantType svType) {
        final CopyNumberPosteriorDistribution dist = getCopyNumberStatePosterior(genotype);
        final int neutralCopyState = SVGenotypeEngine.getNeutralCopyNumber(genotype);
        // Return original genotype if information is missing
        if (dist == null) {
            return genotype;
        }

        final List<IntegerCopyNumberState> states = dist.getIntegerCopyNumberStateList();
        final double[] copyStatePosterior = IntStream.range(0, states.size())
                .mapToDouble(i -> dist.getCopyNumberPosterior(states.get(i))).toArray();
        NaturalLogUtils.normalizeFromLogToLinearSpace(copyStatePosterior);
        final double[] genotypeProbs = getGenotypePosterior(copyStatePosterior, neutralCopyState, svType);
        return genotypeFromGivenProbs(genotype, svType, genotypeProbs);
    }

    public static CopyNumberPosteriorDistribution getCopyNumberStatePosterior(final Genotype genotype) {
        if (!genotype.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY)) {
            return null;
        }
        final int[] copyStateQualsRaw = VariantContextGetters.getAttributeAsIntArray(genotype, GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY, null, 0);
        final double[] copyStateQuals = IntStream.of(copyStateQualsRaw).mapToDouble(x -> NaturalLogUtils.qualToLogErrorProb(x)).toArray();
        NaturalLogUtils.normalizeLog(copyStateQuals);

        final Map<IntegerCopyNumberState, Double> posteriorsMap = new HashMap<>(SVUtils.hashMapCapacity(copyStateQuals.length));
        for (int i = 0; i < copyStateQuals.length; i++) {
            posteriorsMap.put(new IntegerCopyNumberState(i), copyStateQuals[i]);
        }
        return new CopyNumberPosteriorDistribution(posteriorsMap);
    }

    public static double[] getGenotypePosterior(final double[] copyStatePosterior, final int neutralCopyState,
                                                final StructuralVariantType svType) {
        final double[] genotypePosterior;
        if (svType.equals(StructuralVariantType.DEL)) {
            genotypePosterior = getDeletionGenotypePosterior(neutralCopyState, copyStatePosterior);
        } else if (svType.equals(StructuralVariantType.DUP)) {
            genotypePosterior = getDuplicationGenotypePosterior(neutralCopyState, copyStatePosterior);
        } else {
            throw new IllegalArgumentException("SV type was neither " + StructuralVariantType.DEL.name() + " nor " + StructuralVariantType.DUP.name());
        }
        return MathUtils.normalizeSumToOne(genotypePosterior); // Fixes precision errors
    }

    private static double[] getDeletionGenotypePosterior(final int neutralCopyState, final double[] copyStatePosterior) {
        final int numGenotypes = neutralCopyState + 1;
        final double[] genotypePosteriors = new double[numGenotypes];
        double copyNeutralPosterior = 0;
        for (int i = neutralCopyState; i < copyStatePosterior.length; i++) {
            copyNeutralPosterior += copyStatePosterior[i];
        }
        genotypePosteriors[0] = copyNeutralPosterior;
        int idx = 1;
        for (int i = neutralCopyState - 1; i >= 0; i--) {
            genotypePosteriors[idx] = copyStatePosterior[i];
            idx++;
        }
        return genotypePosteriors;
    }

    private static double[] getDuplicationGenotypePosterior(final int neutralCopyState, final double[] copyStatePosterior) {
        final int numGenotypes = copyStatePosterior.length - neutralCopyState;
        final double[] genotypePosteriors = new double[numGenotypes];
        double copyNeutralPosterior = 0;
        for (int i = 0; i <= neutralCopyState; i++) {
            copyNeutralPosterior += copyStatePosterior[i];
        }
        genotypePosteriors[0] = copyNeutralPosterior;
        int idx = 1;
        for (int i = neutralCopyState + 1; i < copyStatePosterior.length; i++) {
            genotypePosteriors[idx] = copyStatePosterior[i];
            idx++;
        }
        return genotypePosteriors;
    }
}
