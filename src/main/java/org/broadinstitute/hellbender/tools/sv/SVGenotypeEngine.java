package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

public class SVGenotypeEngine {

    public static final Allele DEFAULT_REF_ALLELE = Allele.REF_N;
    public static final Allele BND_SYMBOLIC_ALLELE = Allele.create("<" + StructuralVariantType.BND.name() + ">", false);
    public static final List<StructuralVariantType> CNV_TYPES = Lists.newArrayList(StructuralVariantType.DEL, StructuralVariantType.DUP);

    public VariantContext genotypeVariant(final VariantContext variant) {
        final StructuralVariantType svType = variant.getStructuralVariantType();
        final List<Genotype> genotypes = variant.getGenotypes().stream()
                .map(g -> buildGenotypeUsingPL(g, svType))
                .collect(Collectors.toList());
        final VariantContextBuilder builder = new VariantContextBuilder(variant);
        builder.genotypes(genotypes);
        return builder.make();
    }

    protected Genotype buildGenotypeUsingPL(final Genotype genotype, final StructuralVariantType svType) {
        Utils.validateArg(genotype.hasPL(), "PL not assigned to this genotype");
        final double[] genotypeLikelihoods = IntStream.of(genotype.getPL()).mapToDouble(Double::valueOf).toArray();
        return buildGenotypeFromGivenLikelihoods(genotype, svType, genotypeLikelihoods);
    }

    protected Genotype genotypeFromGivenProbs(final Genotype genotype, final StructuralVariantType svType,
                                            final double[] genotypeProbs) {
        final int numStates = genotypeProbs.length;
        final double[] genotypeLikelihoods = new double[numStates];
        for (int i = 0; i < numStates; i++) {
            genotypeLikelihoods[i] = QualityUtils.phredScaleErrorRate(genotypeProbs[i]);
        }
        return buildGenotypeFromGivenLikelihoods(genotype, svType, genotypeLikelihoods);
    }

    protected Genotype buildGenotypeFromGivenLikelihoods(final Genotype genotype, final StructuralVariantType svType,
                                                         final double[] genotypeLikelihoods) {

        // Resize array to be consistent with neutral copy number, if necessary
        final int neutralCopyState = getNeutralCopyNumber(genotype);

        final int maxGenotypeStates = getMaxPossibleGenotypeStates(neutralCopyState, svType);
        final double[] resizedLikelihoods;
        if (genotypeLikelihoods.length > maxGenotypeStates) {
            resizedLikelihoods = resizeGenotypeLikelihoods(genotypeLikelihoods, maxGenotypeStates);
        } else {
            resizedLikelihoods = genotypeLikelihoods;
        }

        final int genotypeIndex = MathUtils.minElementIndex(resizedLikelihoods);
        final double genotypeQuality = neutralCopyState == 0 ? Integer.MAX_VALUE : MathUtils.secondSmallestMinusSmallest(resizedLikelihoods, 0);
        return buildGenotypeFromGivenQual(genotype, svType, resizedLikelihoods, genotypeIndex, genotypeQuality, neutralCopyState);
    }

    private static double[] resizeGenotypeLikelihoods(final double[] genotypeLikelihoods, final int newSize) {
        Utils.validateArg(newSize < genotypeLikelihoods.length, "New size must be smaller");
        final double[] resizedLikelihoods = new double[newSize];
        final int maxIndex = newSize - 1;
        for (int i = 0; i < maxIndex; i++) {
            resizedLikelihoods[i] = genotypeLikelihoods[i];
        }
        final double[] extraLikelihoods = new double[genotypeLikelihoods.length - maxIndex];
        for (int i = 0; i < extraLikelihoods.length; i++) {
            extraLikelihoods[i] = genotypeLikelihoods[maxIndex + i];
        }
        resizedLikelihoods[maxIndex] = QualityUtils.phredSum(extraLikelihoods);
        return resizedLikelihoods;
    }

    private static int getMaxPossibleGenotypeStates(final int neutralCopyState, final StructuralVariantType svType) {
        if (svType.equals(StructuralVariantType.DEL)
                || svType.equals(StructuralVariantType.INS)
                || svType.equals(StructuralVariantType.INV)
                || svType.equals(StructuralVariantType.BND)) {
            return neutralCopyState + 1;
        } else if (svType.equals(StructuralVariantType.DUP)) {
            return neutralCopyState == 0 ? 1 : Integer.MAX_VALUE;
        }
        throw new UserException.BadInput("Unsupported SVTYPE: " + svType);
    }

    private static Genotype buildGenotypeFromGivenQual(final Genotype genotype, final StructuralVariantType svType,
                                                       final double[] genotypeLikelihoods, final int genotypeIndex,
                                                       final double genotypeQuality, final int neutralCopyState) {
        final Allele altAllele = Allele.create("<" + svType.name() + ">", false);
        final Allele refAllele = DEFAULT_REF_ALLELE;
        // TODO: multi-allelic sites
        final List<Allele> alleles = getAlleles(svType, genotypeIndex, neutralCopyState, refAllele, altAllele);
        final Integer copyNumber = getGenotypeCopyNumber(svType, genotypeIndex, neutralCopyState);
        final int genotypeQualityInt = (int) Math.round(genotypeQuality);
        final int[] genotypeLikelihoodsInt = DoubleStream.of(genotypeLikelihoods).mapToInt(x -> (int) Math.round(x)).toArray();

        final GenotypeBuilder builder = new GenotypeBuilder(genotype);
        builder.alleles(alleles);
        builder.attribute(GATKSVVCFConstants.COPY_NUMBER_FIELD, copyNumber);
        builder.GQ(genotypeQualityInt);
        builder.PL(genotypeLikelihoodsInt);
        return builder.make();
    }

    private static Integer getGenotypeCopyNumber(final StructuralVariantType svType, final int genotypeIndex, final int neutralCopyState) {
        if (svType.equals(StructuralVariantType.INS)
                || svType.equals(StructuralVariantType.INV)
                || svType.equals(StructuralVariantType.BND)) {
            return null;
        } else if (svType.equals(StructuralVariantType.DEL)) {
            return Math.max(neutralCopyState - genotypeIndex, 0);
        } else if (svType.equals(StructuralVariantType.DUP)) {
            return neutralCopyState == 0 ? 0 : neutralCopyState + genotypeIndex;
        }
        throw new UserException.BadInput("Unsupported SVTYPE: " + svType);
    }

    private static List<Allele> getAlleles(final StructuralVariantType svType, final int genotypeIndex,
                                           final int neutralCopyState, final Allele refAllele, final Allele altAllele) {
        if (svType.equals(StructuralVariantType.DEL)
                || svType.equals(StructuralVariantType.INS)
                || svType.equals(StructuralVariantType.INV)
                || svType.equals(StructuralVariantType.BND)) {
            return getNonDupAlleles(genotypeIndex, neutralCopyState, refAllele, altAllele);
        } else if (svType.equals(StructuralVariantType.DUP)) {
            return getDupAlleles(neutralCopyState);
        }
        throw new UserException.BadInput("Unsupported SVTYPE: " + svType);
    }

    private static List<Allele> getNonDupAlleles(final int genotypeIndex, final int neutralCopyState,
                                                 final Allele refAllele, final Allele altAllele) {
        final int alleleCount = getAlleleCount(neutralCopyState);
        final int numAltAlleles = Math.min(genotypeIndex, neutralCopyState);
        final List<Allele> alleles = new ArrayList<>(neutralCopyState);
        for (int i = 0; i < alleleCount - numAltAlleles; i++) {
            alleles.add(refAllele);
        }
        for (int i = 0; i < numAltAlleles; i++) {
            alleles.add(altAllele);
        }
        return alleles;
    }

    private static List<Allele> getDupAlleles(final int neutralCopyState) {
        final int alleleCount = getAlleleCount(neutralCopyState);
        final List<Allele> alleles = new ArrayList<>(alleleCount);
        for (int i = 0; i < alleleCount; i++) {
            alleles.add(Allele.NO_CALL);
        }
        return alleles;
    }

    private static int getAlleleCount(final int neutralCopyState) {
        // Can't have mixed rows with some samples missing genotypes, so clamp to a minimum of 1
        return Math.max(1, neutralCopyState);
    }

    public static double calculateLog10PNoError(final Collection<Genotype> genotypes) {
        double log10ProbNoVariant = 0;
        for (final Genotype genotype : genotypes) {
            final int[] genotypeQuals = genotype.getPL();
            if (genotypeQuals.length > 0) {
                log10ProbNoVariant += QualityUtils.qualToErrorProbLog10(genotypeQuals[0]);
            }
        }
        return log10ProbNoVariant;
    }

    public static int getNeutralCopyNumber(final Genotype genotype) {
        Utils.validateArg(genotype.hasExtendedAttribute(GATKSVVCFConstants.NEUTRAL_COPY_NUMBER_KEY), "Missing required genotype field " + GATKSVVCFConstants.NEUTRAL_COPY_NUMBER_KEY);
        return VariantContextGetters.getAttributeAsInt(genotype, GATKSVVCFConstants.NEUTRAL_COPY_NUMBER_KEY, 0);
    }
}
