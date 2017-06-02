package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.Transition;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Perform operations for the orientation bias filter.
 */
public class OrientationBiasFilterer {

    private static final Logger logger = LogManager.getLogger(OrientationBiasFilterer.class);

    public static final Double PRE_ADAPTER_METRIC_NOT_ARTIFACT_SCORE = 100.0;

    /**
     * mode of binomial distribution of pair orientation for the artifact alt reads.
     */
    public static final double BIAS_P = 0.96;

    /** Adds the annotations that can be created for individual variants.  The annotations placed on the genotypes here will be used by the filter.
     * @param vc the variant to prepare for the orientation bias filter. Never {@code null}
     * @param relevantTransitionsWithoutComplement the SNV artifact modes that we want to filter against.  Do not include the complement.  Never {@code null}
     * @param preAdapterQScoreMap mapping from Artifact mode to the preAdapterQ score.  Never {@code null}
     * @return updated VariantContext with new genotypes already populated.
     */
    public static VariantContext annotateVariantContextWithPreprocessingValues(final VariantContext vc, final SortedSet<Transition> relevantTransitionsWithoutComplement, final Map<Transition, Double> preAdapterQScoreMap) {

        Utils.nonNull(vc);
        Utils.nonNull(relevantTransitionsWithoutComplement);
        Utils.nonNull(preAdapterQScoreMap);

        final List<Transition> relevantTransitionsComplement = OrientationBiasUtils.createReverseComplementTransitions(relevantTransitionsWithoutComplement);

        final VariantContextBuilder vcb = new VariantContextBuilder(vc);
        final GenotypesContext genotypesContext = vc.getGenotypes();

        final List<Genotype> newGenotypes = new ArrayList<>();
        for (Genotype genotype: genotypesContext.iterateInSampleNameOrder()) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype);
            genotypeBuilder.attribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_ARTIFACT_MODE, String.valueOf(false));
            genotypeBuilder.attribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_RC_ARTIFACT_MODE, String.valueOf(false));

            final List<Allele> alleles = genotype.getAlleles();
            if (genotype.getPloidy() != 2) {
                logger.warn("No action required:  This tool will skip non-diploid sites.  Saw GT: " + genotype.getGenotypeString() + " at " + vc.toStringWithoutGenotypes());
            }

            // Get the reference allele as a String and make sure that there is only one ref allele and that it is length
            //  one, which would indicate that it could be a part of a SNP/SNV
            final List<String> refAlleles = alleles.stream().filter(a -> a.isReference()).map(a -> a.getBaseString()).collect(Collectors.toList());
            if (((refAlleles.size() == 1) && (refAlleles.get(0).length() == 1))) {
                final Character refAllele = (char) refAlleles.get(0).getBytes()[0];

                // Since we only look at the first alt allele on a site, we do not need a for loop over all non-ref alleles, e.g. for (int i = 1; i < alleles.size(); i++) {
                final Allele allele = genotype.getAllele(1);
                if (allele.isCalled() && allele.isNonReference() && !allele.equals(Allele.SPAN_DEL)
                        && allele.getBaseString().length() == 1) {

                    final Transition genotypeMode = Transition.transitionOf(refAllele, allele.getBaseString().charAt(0));
                    final boolean isRelevantArtifact = relevantTransitionsWithoutComplement.contains(genotypeMode);
                    final boolean isRelevantArtifactComplement = relevantTransitionsComplement.contains(genotypeMode);

                    genotypeBuilder.attribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_ARTIFACT_MODE, String.valueOf(isRelevantArtifact));
                    genotypeBuilder.attribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_RC_ARTIFACT_MODE, String.valueOf(isRelevantArtifactComplement));

                    genotypeBuilder.attribute(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_FIELD_NAME, preAdapterQScoreMap.getOrDefault(genotypeMode, PRE_ADAPTER_METRIC_NOT_ARTIFACT_SCORE));
                    genotypeBuilder.attribute(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_RC_FIELD_NAME, preAdapterQScoreMap.getOrDefault(genotypeMode.complement(), PRE_ADAPTER_METRIC_NOT_ARTIFACT_SCORE));

                    // FOB for the complement is ALT_F2R1/(ALT_F2R1 + ALT_F1R2) and for the actual artifact mode is ALT_F1R2/(ALT_F2R1 + ALT_F1R2)
                    // FOB is the fraction of alt reads indicating orientation bias error (taking into account artifact mode complement).
                    if (isRelevantArtifact || isRelevantArtifactComplement) {
                        final int totalAltAlleleCount = genotype.hasAD() ? genotype.getAD()[1] : 0;
                        final double fob = calculateFob(genotype, isRelevantArtifact);
                        genotypeBuilder.attribute(OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, ArtifactStatisticsScorer.calculateArtifactPValue(totalAltAlleleCount, (int) Math.round(fob * totalAltAlleleCount), BIAS_P));
                        genotypeBuilder.attribute(OrientationBiasFilterConstants.FOB, fob);
                    } else {
                        genotypeBuilder.attribute(OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, VCFConstants.EMPTY_ALLELE);
                        genotypeBuilder.attribute(OrientationBiasFilterConstants.FOB, VCFConstants.EMPTY_ALLELE);
                    }
                }
            }

            final Genotype newGenotype = genotypeBuilder.make();
            newGenotypes.add(newGenotype);
        }
        vcb.genotypes(newGenotypes);
        return vcb.make();
    }

    private static double calculateFob(final Genotype genotype, final boolean isRelevantArtifact) {
        final int altF2R1 = OrientationBiasUtils.getGenotypeInteger(genotype, GATKVCFConstants.OXOG_ALT_F2R1_KEY, 0);
        final int altF1R2 = OrientationBiasUtils.getGenotypeInteger(genotype, GATKVCFConstants.OXOG_ALT_F1R2_KEY, 0);
        return (isRelevantArtifact ? altF1R2 : altF2R1) / (double) (altF1R2 + altF2R1);
    }

    /** Filter genotypes with orientation bias filter while trying to keep the false discovery rate for the sample below the threshold specified.
     *
     * @param fdrThreshold Maximum FDR filter should allow.
     * @param relevantTransitionsWithoutComplements Transitions to filter.  Do not include complements. Never {@code null}.
     * @param preAdapterQAnnotatedVariants Variant contexts that have already been annotated by {@link OrientationBiasFilterer#annotateVariantContextWithPreprocessingValues(VariantContext, SortedSet, Map)} Never {@code null}.
     * @param preAdapterQScoreMap Mapping from Transition to the preAdapterQ score.  relevantTransitions should be included as keys, but not required. Never {@code null}.
     * @return The same variant contexts with the genotypes filter field populated with orientation bias filtering results.  If the variant contexts have no samples, then the variant contexts are returned without any annotation.
     */
    public static List<VariantContext> annotateVariantContextsWithFilterResults(final double fdrThreshold, final SortedSet<Transition> relevantTransitionsWithoutComplements, final List<VariantContext> preAdapterQAnnotatedVariants,
                                                                                final Map<Transition, Double> preAdapterQScoreMap) {

        // This holds the relevantArtifact modes and the complements, since we need to filter both.
        final SortedSet<Transition> relevantTransitions = new TreeSet<>();
        relevantTransitions.addAll(relevantTransitionsWithoutComplements);
        relevantTransitions.addAll(OrientationBiasUtils.createReverseComplementTransitions(relevantTransitionsWithoutComplements));

        if (preAdapterQAnnotatedVariants.size() == 0) {
            logger.info("No samples found in this file.  NO FILTERING BEING DONE.");
            return preAdapterQAnnotatedVariants;
        }

        final List<String> sampleNames = preAdapterQAnnotatedVariants.get(0).getSampleNamesOrderedByName();
        final Map<String, SortedMap<Genotype, VariantContext>> sampleNameToVariants = createSampleToGenotypeVariantContextSortedMap(sampleNames, preAdapterQAnnotatedVariants);

        // This map will hold all updated genotypes (across samples)
        final Map<VariantContext, List<Genotype>> newGenotypes = new HashMap<>();

        // For each sample, perform the actual filtering.
        for (final String sampleName: sampleNames) {

            // Count the total number of unfiltered genotypes for this sample, regardless of artifact mode or not.
            //  I.e. only remove filtered or ref/ref genotypes and count the rest.
            final long unfilteredGenotypeCount = OrientationBiasUtils.calculateUnfilteredNonRefGenotypeCount(preAdapterQAnnotatedVariants, sampleName);

            final SortedMap<Genotype, VariantContext> genotypesToConsiderForFiltering = sampleNameToVariants.get(sampleName);

            // Save some time, especially for the normal sample
            if (genotypesToConsiderForFiltering.keySet().size() == 0) {
                logger.info(sampleName + ": Nothing to filter.");
                continue;
            }

            // Populate counts of artifact modes.  These are genotypes that we can potentially cut.
            final Map<Transition, Long> transitionCount = createTransitionCountMap(relevantTransitions, genotypesToConsiderForFiltering);

            // Global number of artifacts to cut (not adjusted for preAdapterQ)
            final Map<Transition, Long> transitionNumToCut = createTransitionToNumCutPrePreAdapterQ(fdrThreshold, sampleName, unfilteredGenotypeCount, genotypesToConsiderForFiltering, transitionCount);

            // Adjust the artifact mode to cut based on the preAdapterQ score from picard
            for (final Transition transition : transitionNumToCut.keySet()) {
                final Transition modeOrReverseComplement = relevantTransitionsWithoutComplements.contains(transition) ? transition : transition.complement();
                final double suppression = ArtifactStatisticsScorer.calculateSuppressionFactorFromPreAdapterQ(preAdapterQScoreMap.get(modeOrReverseComplement));

                transitionNumToCut.put(transition, Math.round(transitionNumToCut.get(transition) * suppression));
                logger.info(sampleName + ": Cutting (" + transition + ") post-preAdapterQ: " + transitionNumToCut.get(transition));
            }

            // Add filtering results to genotypes and store a pair of genotype variant context, so that we can update variant contexts later.
            logger.info(sampleName + ": Adding orientation bias filter results to genotypes...");

            final Map<Transition, Long> transitionCutSoFar = new HashMap<>();
            relevantTransitions.stream().forEach(transition -> transitionCutSoFar.put(transition, 0L));
            for (final Genotype genotype : genotypesToConsiderForFiltering.keySet()) {
                final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype);
                // Since we only have the ALT_F2R1 and ALT_F1R2 counts for the first alt allele, we will not consider a site where transition is not artifact mode in the first alt allele.
                final Transition transition = Transition.transitionOf(genotype.getAllele(0).getBaseString().charAt(0), genotype.getAllele(1).getBaseString().charAt(0));

                if (!transitionNumToCut.keySet().contains(transition)) {
                    logger.warn("Have to skip genotype: " + genotype + " since it does not have the artifact mode in the first alt allele.  Total alleles: " + genotype.getAlleles().size());
                } else {

                    final Double pValue = OrientationBiasUtils.getGenotypeDouble(genotype, OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, 0.0);
                    final Double fractionOfReadsSupportingOrientationBias = OrientationBiasUtils.getGenotypeDouble(genotype, OrientationBiasFilterConstants.FOB, 0.0);
                    if (transitionCutSoFar.get(transition) < transitionNumToCut.get(transition)) {
                        final String updatedFilter = OrientationBiasUtils.addFilterToGenotype(genotype.getFilters(), OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT);
                        genotypeBuilder.filter(updatedFilter);
                        transitionCutSoFar.put(transition, transitionCutSoFar.get(transition) + 1);
                        logger.info("Cutting: " + genotype.getSampleName() + " " + genotype.getAllele(0) + " " + genotype.getAllele(1)
                                + " p=" + pValue + " Fob=" + fractionOfReadsSupportingOrientationBias);
                    } else {
                        // No need to do anything for the genotype filter, so just log it.
                        logger.info("Passing: " + genotype.getSampleName() + " " + genotype.getAllele(0) + " " + genotype.getAllele(1)
                                + " p=" + pValue + " Fob=" + fractionOfReadsSupportingOrientationBias);
                    }
                }

                newGenotypes.computeIfAbsent(genotypesToConsiderForFiltering.get(genotype), v -> new ArrayList<>()).add(genotypeBuilder.make());
            }
        }

        // Create the final variants with the modified genotypes
        logger.info("Updating genotypes and creating final list of variants...");
        final List<VariantContext> finalVariants = new ArrayList<>();
        final ProgressMeter resultProgressMeter = new ProgressMeter();
        resultProgressMeter.start();
        for (final VariantContext vc : preAdapterQAnnotatedVariants) {
            if (newGenotypes.containsKey(vc)) {
                final GenotypesContext gcc = GenotypesContext.copy(vc.getGenotypes());
                final List<Genotype> newGenotypesForThisVariantContext = newGenotypes.get(vc);
                newGenotypesForThisVariantContext.forEach(gcc::replace);
                final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(vc).genotypes(gcc);
                if (newGenotypesForThisVariantContext.stream().anyMatch(g -> (g != null) && (g.getFilters() != null) && (g.getFilters().contains(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT)))) {
                    variantContextBuilder.filter(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT);
                }
                final VariantContext updatedVariantContext = variantContextBuilder.make();
                finalVariants.add(updatedVariantContext);
            } else {
                finalVariants.add(vc);
            }
            resultProgressMeter.update(new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd()));
        }
        resultProgressMeter.stop();
        return finalVariants;
    }


    private static Map<Transition, Long> createTransitionToNumCutPrePreAdapterQ(double fdrThresh, String sampleName, long unfilteredGenotypeCount, final SortedMap<Genotype, VariantContext> genotypesToConsiderForFiltering, final Map<Transition, Long> transitionCount) {
        final long allTransitionCount = transitionCount.values().stream().mapToLong(Long::longValue).sum();
        final int totalNumToCut = calculateTotalNumToCut(fdrThresh, unfilteredGenotypeCount, genotypesToConsiderForFiltering);

        logger.info(sampleName + ": Cutting (total) pre-preAdapterQ: " + String.valueOf(totalNumToCut));

        // Adjust the number to cut based on artifact mode
        final Map<Transition, Long> transitionNumToCut = new HashMap<>();
        transitionCount.keySet().stream().forEach(transition -> transitionNumToCut.put(transition, 0L));
        for (final Transition transition : transitionNumToCut.keySet()) {
            transitionNumToCut.put(transition, Long.valueOf(Math.round(totalNumToCut * transitionCount.get(transition) / allTransitionCount)));
            logger.info(sampleName + ": Cutting (" + transition + ") pre-preAdapterQ: " + transitionNumToCut.get(transition));
        }
        return transitionNumToCut;
    }

    private static int calculateTotalNumToCut(final double fdrThresh, final long unfilteredGenotypeCount, final SortedMap<Genotype, VariantContext> genotypesToConsiderForFiltering) {
        final List<Double> pArtifactScores = genotypesToConsiderForFiltering.keySet().stream()
                .map(g -> OrientationBiasUtils.getGenotypeDouble(g, OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, 0.0))
                .collect(Collectors.toList());

        // When doing the Benjamini-Hochberg procedure, we need to include the non-artifact mode SNVs to guarantee the
        //  FDR is kept to the specified threshold.
        //   Pad the list of pArtifact scores so that the the non-artifacts (w/ a p-value of zero) are included.
        final int numToPadZeroes = (int) unfilteredGenotypeCount - pArtifactScores.size();
        List<Double> finalPArtifactScores = new ArrayList<>();
        finalPArtifactScores.addAll(pArtifactScores);
        finalPArtifactScores.addAll(Collections.nCopies(numToPadZeroes, 0.0));

        return calculateTotalNumToCut(fdrThresh, unfilteredGenotypeCount, finalPArtifactScores);
    }

    /**
     *
     * @param fdrThreshold desired maximum FDR threshold.  Must be greater than 0
     * @param unfilteredGenotypeCount total number of unfiltered variants
     * @param pArtifactScoresIncludingNonArtifact sorted list (descending) of the pArtifact scores.  Should include zeros for nonArtifact variants.
     * @return total number of artifact mode variants to cut in order to be below the specified FDR threshold.
     */
    @VisibleForTesting
    static int calculateTotalNumToCut(final double fdrThreshold, final long unfilteredGenotypeCount, final List<Double> pArtifactScoresIncludingNonArtifact) {
        ParamUtils.isPositive(fdrThreshold, "FDR threshold must be positive and greater than zero.");

        // Benjamini-Hochberg procedure https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini.E2.80.93Hochberg_procedure
        return IntStream.range(0, pArtifactScoresIncludingNonArtifact.size())
                .filter(i -> pArtifactScoresIncludingNonArtifact.get(i) < fdrThreshold * (i + 1) / unfilteredGenotypeCount)
                .findFirst().orElse(pArtifactScoresIncludingNonArtifact.size() - 1);
    }

    private static Map<Transition, Long> createTransitionCountMap(SortedSet<Transition> relevantTransitions, SortedMap<Genotype, VariantContext> genotypesToConsiderForFiltering) {
        final Map<Transition, Long> transitionCount = new HashMap<>();
        relevantTransitions.stream().forEach(transition -> transitionCount.put(transition, 0L));
        for (final Genotype g : genotypesToConsiderForFiltering.keySet()) {
            relevantTransitions.stream()
                    .filter(transition -> OrientationBiasUtils.isGenotypeInTransition(g, transition))
                    .forEach(transition -> transitionCount.put(transition, transitionCount.get(transition) + 1));
        }
        return transitionCount;
    }

    /** Creates a map that includes the artifact mode complements.
     *
     * @param sampleNames Names of samples to generate the map.
     * @param variants The associated VariantContexts.  The given sample names should be included.
     * @return a mapping from the sampleNames to the a sorted (by p_artifact score) map that associates genotypes to their enclosing variant context.
     */
    public static Map<String, SortedMap<Genotype, VariantContext>> createSampleToGenotypeVariantContextSortedMap(final List<String> sampleNames, final Collection<VariantContext> variants) {

        // Sorts in reverse order (highest p_artifact goes first and will not allow anything to be equal
        //  unless they share the same reference)
        // Note the negative sign is to sort in reverse error.
        final Comparator<Genotype> genotypePArtifactComparator = Comparator
                .comparingDouble((Genotype g) -> -OrientationBiasUtils.getGenotypeDouble(g, OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, 0.0))
                .thenComparingInt(g -> g.hashCode());

        final Map<String, SortedMap<Genotype, VariantContext>> sampleNameToVariants = new HashMap<>();

        final ProgressMeter customProgressMeter = new ProgressMeter(0.1);
        customProgressMeter.start();

        // Populate a mapping of genotypes that we might want to filter to their variant context.
        //  Make sure that the keys are sorted by cumulative probability of being an artifact.
        for (final String sampleName : sampleNames) {
            final SortedMap<Genotype, VariantContext> genotypesToConsiderForFiltering = new TreeMap<>(genotypePArtifactComparator);
            for (final VariantContext vc : variants) {
                vc.getGenotypes(sampleName).stream()
                        .filter(g -> isFilteringCandidate(g, vc))
                        .forEach(genotype -> genotypesToConsiderForFiltering.put(genotype, vc));
                customProgressMeter.update(new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd()));
            }
            sampleNameToVariants.put(sampleName, genotypesToConsiderForFiltering);
        }
        customProgressMeter.stop();
        return sampleNameToVariants;
    }

    /**
     * Determine whether this genotype can be filtered by the orientation bias filter.
     *
     * @param genotype
     * @param vc The variant context that contains the given genotype
     * @return whether this genotype should be considered for filtering due to orientation bias.
     */
    private static boolean isFilteringCandidate(final Genotype genotype, final VariantContext vc) {
        return !vc.isFiltered()
                && !genotype.isFiltered()
                && (genotype.getAnyAttribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_ARTIFACT_MODE).equals(String.valueOf(true))
                    || genotype.getAnyAttribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_RC_ARTIFACT_MODE).equals(String.valueOf(true)));
    }

    /** Ingest the current VCF header and update it with the information necessary for the Orientation Bias filter to run.
     *
     * @param inputVCFHeader original header.  Never {@code null}
     * @param commandLine The command line used to run this tool.
     * @param transitions  Never {@code null}
     * @return updated VCF Header
     */
    public static VCFHeader createVCFHeader(final VCFHeader inputVCFHeader, final String commandLine, final List<String> transitions) {
        Utils.nonNull(inputVCFHeader);
        Utils.nonNull(transitions);

        // Setup header for output file
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.add(new VCFFormatHeaderLine(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_FIELD_NAME, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Measure (across entire bam file) of orientation bias for a given REF/ALT error."));
        headerLines.add(new VCFFormatHeaderLine(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_RC_FIELD_NAME, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Measure (across entire bam file) of orientation bias for the complement of a given REF/ALT error."));
        headerLines.add(new VCFFormatHeaderLine(OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Orientation bias p value for the given REF/ALT artifact or its complement."));
        headerLines.add(new VCFFormatHeaderLine(OrientationBiasFilterConstants.FOB, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Fraction of alt reads indicating orientation bias error (taking into account artifact mode complement)."));
        headerLines.add(new VCFFormatHeaderLine(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_ARTIFACT_MODE, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Whether the variant can be one of the given REF/ALT artifact modes."));
        headerLines.add(new VCFFormatHeaderLine(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_RC_ARTIFACT_MODE, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Whether the variant can be one of the given REF/ALT artifact mode complements."));
        headerLines.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_FILTER_KEY, 1, VCFHeaderLineType.String, "Genotype-level filter"));
        headerLines.add(new VCFFilterHeaderLine(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT, "Orientation bias (in one of the specified artifact mode(s) or complement) seen in one or more samples."));
        headerLines.add(new VCFSimpleHeaderLine("orientation_bias_artifact_modes", String.join("|", transitions), "The artifact modes that were used for orientation bias artifact filtering for this VCF"));
        headerLines.add(new VCFHeaderLine("command", commandLine));
        final SampleList samples = new IndexedSampleList(inputVCFHeader.getGenotypeSamples());
        final Set<String> sampleNameSet = samples.asSetOfSamples();
        return new VCFHeader(headerLines, sampleNameSet);
    }
}
