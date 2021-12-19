package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealMatrixChangingVisitor;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticLikelihoodsEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.SubsettedLikelihoodMatrix;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Total depth of coverage per sample and over all samples (DP)")
public final class AllelePseudoDepth implements GenotypeAnnotation {

    private static DecimalFormat DEPTH_FORMAT = new DecimalFormat("#.##");
    private static DecimalFormat FRACTION_FORMAT = new DecimalFormat("#.####");

    public static final VCFFormatHeaderLine DEPTH_HEADER_LINE =
            new VCFFormatHeaderLine(GATKVCFConstants.PSEUDO_DEPTH_KEY, VCFHeaderLineCount.R, VCFHeaderLineType.Float, "Allele depth based on Dirichlet posterior pseudo-counts");
    public static final VCFFormatHeaderLine FRACTION_HEADER_LINE =
            new VCFFormatHeaderLine(GATKVCFConstants.PSEUDO_FRACTION_KEY, VCFHeaderLineCount.R, VCFHeaderLineType.Float, "Allele Fraction based on Dirichlet posterior pseudo-counts");

    @Argument(fullName = "dirichlet-prior-pseudo-count",
              doc = "Pseudo-count used as prior for all alleles. The default is 1.0 resulting in a flat prior")
    @Advanced
    public double prior = 1.0;

    @Argument(fullName = "dirichlet-keep-prior-in-count",
              doc = "By default we don't keep the prior use in the output counts ase it makes it easier to interpret" +
                      "this quantity as the number of supporting reads specially in low depth sites. We this toggled the prior is included")
    public boolean keepPriorInCount = false;

    @Argument(fullName = "pseudo-count-weight-decay-rate",
              doc = "A what rate the weight of a read decreases base on its informativeness; e.g. 1.0 is " +
                      "linear decay (default), 2.0 is for quadratic decay",
              minValue = 0.0)
    public double weightDecay = 1.0;

    private static final List<String> KEYS = Arrays.asList(GATKVCFConstants.PSEUDO_DEPTH_KEY, GATKVCFConstants.PSEUDO_FRACTION_KEY);
    private static final List<VCFCompoundHeaderLine> HEADER_LINES = Arrays.asList(DEPTH_HEADER_LINE, FRACTION_HEADER_LINE);

    private Int2ObjectMap<double[]> priorPseudoCounts = new Int2ObjectOpenHashMap<>();
    @Override
    public List<VCFCompoundHeaderLine> getDescriptions() {
        return HEADER_LINES;
    }

    @Override
    public List<String> getKeyNames() {
        return KEYS;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void annotate(final ReferenceContext ref, final VariantContext vc,
                         final Genotype g, final GenotypeBuilder gb,
                         final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        if (likelihoods == null) {
            return;
        }
        final List<Allele> allelesToEmit = vc.getAlleles();
        if (allelesToEmit.size() <= 1) {
            return;
        }
        final double[] prior = composePriorPseudoCounts(allelesToEmit.size());
        final LikelihoodMatrix<GATKRead, Allele> sampleMatrix = likelihoods.sampleMatrix(likelihoods.indexOfSample(g.getSampleName()));
        final LikelihoodMatrix<GATKRead, Allele> sampleMatrixForAlleles =
                allelesToEmit.size() == likelihoods.numberOfAlleles()
                        ? sampleMatrix : new SubsettedLikelihoodMatrix<>(sampleMatrix, allelesToEmit);

        final double[] posteriors;
        if (sampleMatrix.evidence().isEmpty()) {
            posteriors = prior;
        } else {
            final RealMatrix lkMatrix = composeInputLikelihoodMatrix(likelihoods, sampleMatrixForAlleles);
            final double[] weights = calculateWeights(lkMatrix);
            posteriors = SomaticLikelihoodsEngine.alleleFractionsPosterior(lkMatrix, prior, weights);
        }

        final double[] frequencies = MathUtils.normalizeSumToOne(posteriors);
        if (!keepPriorInCount) {
            for (int i = 0; i < posteriors.length; i++) {
                posteriors[i] -= prior[i];
            }
        }
        gb.attribute(GATKVCFConstants.PSEUDO_DEPTH_KEY, Arrays.stream(posteriors)
                .mapToObj(DEPTH_FORMAT::format)
                .collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        gb.attribute(GATKVCFConstants.PSEUDO_FRACTION_KEY,Arrays.stream(frequencies)
                .mapToObj(FRACTION_FORMAT::format)
                .collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
    }

    private RealMatrix composeInputLikelihoodMatrix(AlleleLikelihoods<GATKRead, Allele> likelihoods, LikelihoodMatrix<GATKRead, Allele> sampleMatrixForAlleles) {
        final RealMatrix lkMatrix;
        sampleMatrixForAlleles.asRealMatrix().copy();
        if (!likelihoods.isNaturalLog()) {
            lkMatrix = sampleMatrixForAlleles.asRealMatrix().copy();
            final RealMatrixChangingVisitor log10ToLnTransformer = new DefaultRealMatrixChangingVisitor() {
                public double visit(int row, int column, double value) {
                    return Math.max(value, -.1 * sampleMatrixForAlleles.evidence().get(row).getMappingQuality()) * MathUtils.LOG_10;
                }
            };
            lkMatrix.walkInOptimizedOrder(log10ToLnTransformer);
        } else {
            lkMatrix = sampleMatrixForAlleles.asRealMatrix();
        }
        return lkMatrix;
    }

    private double[] calculateWeights(RealMatrix lkMatrix) {
        if (weightDecay == 0) {
            return null;
        } else if (weightDecay < 0.0) {
            throw new IllegalArgumentException("the weight decay must be 0 or greater");
        } else {
            final double[] result = new double[lkMatrix.getColumnDimension()];
            for (int i = 0; i < result.length; i++) {
                double best = lkMatrix.getEntry(0, i);
                double secondBest = Double.NEGATIVE_INFINITY;
                for (int j = 1; j < lkMatrix.getRowDimension(); j++) {
                    final double val = lkMatrix.getEntry(j, i);
                    if (val > best) {
                        secondBest = best;
                        best = val;
                    } else if (val > secondBest) {
                        secondBest = val;
                    }
                }
                final double secondAdjusted = secondBest - best;
                result[i] = 1 - Math.pow(10, secondAdjusted);
                if (weightDecay != 1.0) {
                    result[i] = Math.pow(result[i], weightDecay);
                }
            }
            return result;
        }
    }

    /**
     * Composes prior pseudo-count arrays on demand.
     *
     * @param numberOfAlleles
     * @return never {@code null.}
     */
    private double[] composePriorPseudoCounts(final int numberOfAlleles) {
        final double[] cached = priorPseudoCounts.get(numberOfAlleles);
        if (cached == null) {
            final double[] result = new double[numberOfAlleles];
            Arrays.fill(result, prior);
            priorPseudoCounts.put(numberOfAlleles, result);
            return result;
        } else {
            return cached;
        }
    }
}
