package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class OrientationBiasUtils {
    private static final Logger logger = LogManager.getLogger(OrientationBiasUtils.class);



    public static Double getGenotypeDouble(final Genotype g, final String fieldName, final double defaultValue) {
        final Object o2ExtendedAttribute = g.getExtendedAttribute(fieldName, VCFConstants.MISSING_VALUE_v4);
        final String o2PArtifactAsString = String.valueOf(o2ExtendedAttribute);
        if ((o2PArtifactAsString == null) || (o2PArtifactAsString.equals(VCFConstants.MISSING_VALUE_v4))) {
            return defaultValue;
        }
        return Double.parseDouble(o2PArtifactAsString);
    }
    public static Integer getGenotypeInteger(final Genotype g, final String fieldName, final int defaultValue) {
        final Object genotypeExtendedAttribute = g.getExtendedAttribute(fieldName, VCFConstants.MISSING_VALUE_v4);
        final String genotypePArtifactAsString = String.valueOf(genotypeExtendedAttribute);
        if ((genotypePArtifactAsString == null) || (genotypePArtifactAsString.equals(VCFConstants.MISSING_VALUE_v4))) {
            return defaultValue;
        }
        return Integer.parseInt(genotypePArtifactAsString);
    }

    /** Complement of the artifact mode is NOT considered.
     *
     *  Must be diploid genotype!
     *
     *  The genotype can contain indels, be a reference genotype, etc.  That is all considered.
     *
     * @param g Never {@code null}
     * @param transition Never {@code null}
     * @return whether the given genotype falls into the given artifact mode.  Complement of the artifact mode is NOT considered.
     */
    public static boolean isGenotypeInTransition(final Genotype g, final Transition transition) {
        Utils.nonNull(g, "Genotype cannot be null");
        Utils.nonNull(transition, "Artifact mode cannot be null");

        if (g.getAlleles().size() != 2) {
            throw new UserException.BadInput("Cannot run this method on non-diploid genotypes.");
        }
        return g.getAllele(0).getDisplayString().equals(transition.getLeft().toString()) &&
                g.getAllele(1).getDisplayString().equals(transition.getRight().toString());
    }

    /** See {@link #isGenotypeInTransition}, except that this will take into account complements.
     *
     * @param g See {@link #isGenotypeInTransition(Genotype, Transition)}
     * @param transition See {@link #isGenotypeInTransition(Genotype, Transition)}
     * @return See {@link #isGenotypeInTransition(Genotype, Transition)}
     */
    public static boolean isGenotypeInTransitionWithComplement(final Genotype g, final Transition transition) {
        Utils.nonNull(g, "Genotype cannot be null");
        Utils.nonNull(transition, "Transition cannot be null");

        if (g.getAlleles().size() != 2) {
            throw new UserException.BadInput("Cannot run this method on non-diploid genotypes.");
        }
        final boolean isInTransition = g.getAllele(0).getDisplayString().equals(transition.getLeft().toString()) &&
                g.getAllele(1).getDisplayString().equals(transition.getRight().toString());

        if (isInTransition) {
            return true;
        }

        final Transition transitionComplement = createReverseComplement(transition);
        return g.getAllele(0).getDisplayString().equals(transitionComplement.getLeft().toString()) &&
                g.getAllele(1).getDisplayString().equals(transitionComplement.getRight().toString());

    }

    /** Is this genotype in any of the specified artifact modes (or complements)
     *  See {@link #isGenotypeInTransitionWithComplement(Genotype, Transition)}
     *
     * @param g See {@link #isGenotypeInTransition(Genotype, Transition)}
     * @param transitions Collection to check.  Does not include complements.  Never {@code null}
     * @return whether the given genotype falls into any of the given artifact modes.  Complements of the artifact modes ARE considered.
     */
    public static boolean isGenotypeInTransitionsWithComplement(final Genotype g, final Collection<Transition> transitions) {
        Utils.nonNull(g, "Genotype cannot be null.");

        // transition null values are checked in the next call
        return transitions.stream()
                .anyMatch(am -> isGenotypeInTransitionWithComplement(g, am));
    }

    /**
     *  Create new Transitions that are complements of the given collection.
     *
     * @param transitions Never {@code null}
     * @return list of the artifact mode complements.
     */
    public static List<Transition> createReverseComplementTransitions(final Collection<Transition> transitions) {
        Utils.nonNull(transitions, "Transitions cannot be null.");
        return transitions.stream()
                .map(am -> createReverseComplement(am))
                .collect(Collectors.toList());
    }

    /** Create the complement of the given artifact mode.
     *
     * @param transition artifact mode to create a complement for.  Never {@code null}
     * @return new artifact mode that is the complement
     */
    public static Transition createReverseComplement(final Transition transition) {
        Utils.nonNull(transition, "Transition cannot be null.");
        return Transition.of((char) BaseUtils.simpleComplement((byte) transition.getLeft().charValue()),
                (char) BaseUtils.simpleComplement((byte) transition.getRight().charValue()));
    }

    /**
     *  Write output file that is a summary file with orientation bias filter counts.
     *  The preAdapterQ score for the complement of relevant artifact modes will be reported where applicable.
     *
     * @param sampleTransitionsWithoutComplement List of pairs (sampleName, artifact mode).  The complements of the artifact mode are not included.  These are the artifact modes that were considered for the sample.  Never {@code null}
     * @param variantContexts Variants to write to the file.  Must already be annotated with orientation bias results, Never {@code null}
     * @param preAdapterScoreMap Transition to preAdapterQ score.  Never {@code null}
     * @param outFile Destination location.  Never {@code null}
     */
    public static void writeOrientationBiasSummaryTable(final List<Pair<String, Transition>> sampleTransitionsWithoutComplement, final List<VariantContext> variantContexts,
                                                        final Map<Transition, Double> preAdapterScoreMap,
                                                        final File outFile) {
        Utils.nonNull(sampleTransitionsWithoutComplement);
        Utils.nonNull(variantContexts);
        Utils.nonNull(preAdapterScoreMap);
        Utils.nonNull(outFile);

        // This code is inefficient, since it loops through all variant contexts for each sampleTransition pair.
        try (final TableWriter<Pair<String, Transition> > writer =
                     TableUtils.writer(outFile, OrientationBiasFilterSummaryTableColumn.COLUMNS,
                             //lambda for creating DataLine with sampleName and segment fields
                             (sampleTransitionPair, dataLine) -> {
                                 dataLine.append(sampleTransitionPair.getLeft())
                                 .append(sampleTransitionPair.getRight().toString() + "(" +
                                         OrientationBiasUtils.createReverseComplement(sampleTransitionPair.getRight()).toString() + ")")
                                 .append(preAdapterScoreMap.getOrDefault(sampleTransitionPair.getRight(), OrientationBiasFilterer.PRE_ADAPTER_METRIC_NOT_ARTIFACT_SCORE))
                                 .append(calculateNumTransition(sampleTransitionPair.getLeft(), variantContexts, sampleTransitionPair.getRight()))
                                 .append(calculateNumTransitionFilteredByOrientationBias(sampleTransitionPair.getLeft(), variantContexts, sampleTransitionPair.getRight()))
                                 .append(calculateNumNotTransition(sampleTransitionPair.getLeft(), variantContexts, sampleTransitionPair.getRight()))
                                 .append(calculateUnfilteredNonRefGenotypeCount(variantContexts, sampleTransitionPair.getLeft()));
                             }
                     )) {
            for (final Pair<String, Transition> sampleTransition : sampleTransitionsWithoutComplement) {
                writer.writeRecord(Utils.nonNull(sampleTransition, "List of variantContexts contains a null."));
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e);
        }
    }

    /** Includes complements.  Excludes filtered variant contexts.  Excludes genotypes that were filtered by something other than the orientation bias filter. */
    @VisibleForTesting
    static long calculateNumTransition(final String sampleName, final List<VariantContext> variantContexts, final Transition transition) {
        final Transition complement = OrientationBiasUtils.createReverseComplement(transition);
        return getNumArtifactGenotypeStream(sampleName, variantContexts, transition, complement)
                .filter(g -> (g.getFilters() == null) || (g.getFilters().equals(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT)))
                .count();
    }

    /** Includes complements.  Excludes filtered variant contexts.  Excludes genotypes that were filtered by something other than the orientation bias filter. */
    @VisibleForTesting
    static long calculateNumNotTransition(final String sampleName, final List<VariantContext> variantContexts, final Transition transition) {
        final Transition complement = OrientationBiasUtils.createReverseComplement(transition);
        return getGenotypeStream(sampleName, variantContexts)
                .filter(g -> !OrientationBiasUtils.isGenotypeInTransition(g, complement) && !OrientationBiasUtils.isGenotypeInTransition(g, transition))
                .count();
    }

    private static Stream<Genotype> getNumArtifactGenotypeStream(String sampleName, List<VariantContext> variantContexts, Transition transition, Transition complement) {
        return getGenotypeStream(sampleName, variantContexts)
                .filter(g -> OrientationBiasUtils.isGenotypeInTransition(g, complement) || OrientationBiasUtils.isGenotypeInTransition(g, transition));
    }

    private static Stream<Genotype> getGenotypeStream(String sampleName, List<VariantContext> variantContexts) {
        return variantContexts.stream()
                .filter(vc -> !vc.isFiltered())
                .map(vc -> vc.getGenotype(sampleName));
    }

    /** Includes complements */
    @VisibleForTesting
    static long calculateNumTransitionFilteredByOrientationBias(final String sampleName, final List<VariantContext> variantContexts, final Transition transition) {
        final Transition complement = OrientationBiasUtils.createReverseComplement(transition);
        return getNumArtifactGenotypeStream(sampleName, variantContexts, transition, complement)
                .filter(g -> (g != null) && (g.getFilters() != null)
                        && g.getFilters().contains(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT))
                .count();
    }

    /**
     * @param variants Never {@code null}
     * @param sampleName Never {@code null}.  sampleName to base counts.
     * @return a count of the non-ref, unfiltered genotypes for the given sample name in the given variant contexts.
     */
    public static long calculateUnfilteredNonRefGenotypeCount(final List<VariantContext> variants, final String sampleName) {
        Utils.nonNull(variants);
        Utils.nonNull(sampleName);

        return getGenotypeStream(sampleName, variants)
                .filter(g -> !g.isFiltered())
                .filter(g -> !g.getAllele(0).basesMatch(g.getAllele(1)))
                .count();
    }

    /** Create an updated genotype string when trying to add a filter value.
     *
     * @param existingFilterValue filter string (presumably) from a genotype
     * @param newFilterToAdd new filter string to append.  Cannot be {@code null}.
     * @return Properly formatted genotype filter string
     */
    public static String addFilterToGenotype(final String existingFilterValue, final String newFilterToAdd) {

        Utils.nonNull(newFilterToAdd);

        if ((existingFilterValue == null) || (existingFilterValue.trim().length() == 0)
                || (existingFilterValue.equals(VCFConstants.UNFILTERED))
                || (existingFilterValue.equals(VCFConstants.PASSES_FILTERS_v4))) {
            return newFilterToAdd;
        } else if (existingFilterValue.length() > 0) {
            return existingFilterValue + VCFConstants.FILTER_CODE_SEPARATOR + newFilterToAdd;
        } else {
            final String appendedFilterString = existingFilterValue + VCFConstants.FILTER_CODE_SEPARATOR + newFilterToAdd;
            logger.warn("Existing genotype filter could be incorrect: " + existingFilterValue + " ... Proceeding with " + appendedFilterString + " ...");
            return appendedFilterString;
        }
    }
}
