package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParamUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.genotyper.GenotypePriorCalculator;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.function.BiPredicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * HaplotypeCaller's genotyping strategy implementation.
 */
public class HaplotypeCallerGenotypingEngine extends GenotypingEngine<StandardCallerArgumentCollection> {

    private static final Logger logger = LogManager.getLogger(HaplotypeCallerGenotypingEngine.class);
    private static final OneShotLogger DRAGENConaminationWarning = new OneShotLogger(logger);

    private final GenotypingModel genotypingModel;

    private final PloidyModel ploidyModel;
    private final ReferenceConfidenceMode referenceConfidenceMode;
    protected final double snpHeterozygosity;
    protected final double indelHeterozygosity;

    private final int maxGenotypeCountToEnumerate;
    private final Map<Integer, Integer> practicalAlleleCountForPloidy = new HashMap<>();

    private final boolean doPhysicalPhasing;

    private final DragstrParams dragstrParams;

    private final HaplotypeCallerArgumentCollection hcArgs;

    /**
     * {@inheritDoc}
     * @param configuration {@inheritDoc}
     * @param samples {@inheritDoc}
     * @param doPhysicalPhasing whether to try physical phasing.
     */
    public HaplotypeCallerGenotypingEngine(final HaplotypeCallerArgumentCollection configuration, final SampleList samples, final boolean doPhysicalPhasing, final boolean applyBQD) {
        super(configuration.standardArgs, samples, false);
        hcArgs = configuration;
        this.doPhysicalPhasing = doPhysicalPhasing;
        ploidyModel = new HomogeneousPloidyModel(samples,configuration.standardArgs.genotypeArgs.samplePloidy);
        dragstrParams = DragstrParamUtils.parse(configuration.likelihoodArgs.dragstrParams);
        genotypingModel = hcArgs.applyBQD || hcArgs.applyFRD ?
                new DRAGENGenotypesModel(applyBQD, hcArgs.applyFRD, hcArgs.informativeReadOverlapMargin, hcArgs.maxEffectiveDepthAdjustment, dragstrParams) :
                new IndependentSampleGenotypesModel();
        maxGenotypeCountToEnumerate = configuration.standardArgs.genotypeArgs.maxGenotypeCount;
        referenceConfidenceMode = configuration.emitReferenceConfidence;
        snpHeterozygosity = configuration.standardArgs.genotypeArgs.snpHeterozygosity;
        indelHeterozygosity = configuration.standardArgs.genotypeArgs.indelHeterozygosity;
    }

    @Override
    protected String callSourceString() {
        return "HC_call";
    }

    @Override
    protected boolean forceKeepAllele(final Allele allele) {
        return allele.equals(Allele.NON_REF_ALLELE,false) || referenceConfidenceMode != ReferenceConfidenceMode.NONE;
    }

    /**
     * Main entry point of class - given a particular set of haplotypes, samples and reference context, compute
     * genotype likelihoods and assemble into a list of variant contexts and genomic events ready for calling
     *
     * The list of samples we're working with is obtained from the readLikelihoods
     *
     * @param haplotypes                             Haplotypes to assign likelihoods to
     * @param readLikelihoods                        Map from reads->(haplotypes,likelihoods)
     * @param perSampleFilteredReadList              Map from sample to reads that were filtered after assembly and before calculating per-read likelihoods.
     * @param ref                                    Reference bytes at active region
     * @param refLoc                                 Corresponding active region genome location
     * @param activeRegionWindow                     Active window
     * @param givenAlleles                Alleles to genotype
     * @param emitReferenceConfidence whether we should add a &lt;NON_REF&gt; alternative allele to the result variation contexts.
     * @param maxMnpDistance Phased substitutions separated by this distance or less are merged into MNPs.  More than
     *                       two substitutions occurring in the same alignment block (ie the same M/X/EQ CIGAR element)
     *                       are merged until a substitution is separated from the previous one by a greater distance.
     *                       That is, if maxMnpDistance = 1, substitutions at 10,11,12,14,15,17 are partitioned into a MNP
     *                       at 10-12, a MNP at 14-15, and a SNP at 17.  May not be negative.
     * @param withBamOut whether to annotate reads in readLikelihoods for future writing to bamout
     * @param suspiciousLocations locations where possible alternative noisy error is affecting the result
     *
     * @return                                       A CalledHaplotypes object containing a list of VC's with genotyped events and called haplotypes
     *
     */
    public CalledHaplotypes assignGenotypeLikelihoods(final List<Haplotype> haplotypes,
                                                      final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                      final Map<String, List<GATKRead>> perSampleFilteredReadList,
                                                      final byte[] ref,
                                                      final SimpleInterval refLoc,
                                                      final SimpleInterval activeRegionWindow,
                                                      final FeatureContext tracker,
                                                      final List<Event> givenAlleles,
                                                      final boolean emitReferenceConfidence,
                                                      final int maxMnpDistance,
                                                      final SAMFileHeader header,
                                                      final boolean withBamOut,
                                                      final Set<Integer> suspiciousLocations,
                                                      final AlleleLikelihoods<GATKRead, Haplotype> preFilteringAlleleLikelihoods) {
        // sanity check input arguments
        Utils.nonEmpty(haplotypes, "haplotypes input should be non-empty and non-null");
        Utils.validateArg(readLikelihoods != null && readLikelihoods.numberOfSamples() > 0, "readLikelihoods input should be non-empty and non-null");
        Utils.validateArg(ref != null && ref.length > 0, "ref bytes input should be non-empty and non-null");
        Utils.nonNull(refLoc, "refLoc must be non-null");
        Utils.validateArg(refLoc.size() == ref.length, " refLoc length must match ref bytes");
        Utils.nonNull(activeRegionWindow, "activeRegionWindow must be non-null");
        Utils.nonNull(givenAlleles, "givenAlleles must be non-null");
        Utils.validateArg(refLoc.contains(activeRegionWindow), "refLoc must contain activeRegionWindow");
        ParamUtils.isPositiveOrZero(maxMnpDistance, "maxMnpDistance may not be negative.");

        // update the haplotypes so we're ready to call, getting the ordered list of positions on the reference
        // that carry events among the haplotypes
        EventMap.buildEventMapsForHaplotypes(haplotypes, ref, refLoc, hcArgs.assemblerArgs.debugAssembly, maxMnpDistance);
        final SortedSet<Integer> eventStarts = EventMap.getEventStartPositions(haplotypes);

        // Walk along each position in the key set and create each event to be outputted
        final Set<Haplotype> calledHaplotypes = new HashSet<>();
        final List<VariantContext> returnCalls = new ArrayList<>();
        final int ploidy = configuration.genotypeArgs.samplePloidy;
        final List<Allele> noCallAlleles = GATKVariantContextUtils.noCallAlleles(ploidy);

        if (withBamOut) {
            //add annotations to reads for alignment regions and calling regions
            AssemblyBasedCallerUtils.annotateReadLikelihoodsWithRegions(readLikelihoods, activeRegionWindow);
        }

        // null if there is no potential uses of DRAGstr in this region.
        final DragstrReferenceAnalyzer dragstrs = constructDragstrReferenceSTRAnalyzerIfNecessary(haplotypes, ref, refLoc, eventStarts);

        final BiPredicate<GATKRead, SimpleInterval> readQualifiesForGenotypingPredicate = composeReadQualifiesForGenotypingPredicate(hcArgs);

        // haplo-genotype posteriors by joint detection interval for exact joint detection genotyping within determined intervals
        // note: this will be empty is haplotypes are not partially determined
        final OverlapDetector<Pair<GenotypingLikelihoods<Haplotype>, GenotypingLikelihoods<Haplotype>>> haploGTLikelihoodAndPosteriorOverlapDetector =
                computeHaploGenotypeLikelihoodsAndPosteriors(readLikelihoods, ref, refLoc, activeRegionWindow, header, ploidy, dragstrs, readQualifiesForGenotypingPredicate);

        for( final int loc : eventStarts ) {
            if( loc < activeRegionWindow.getStart() || loc > activeRegionWindow.getEnd() ) {
                continue;
            }

            final List<VariantContext> eventsAtThisLoc = AssemblyBasedCallerUtils.getVariantsFromActiveHaplotypes(loc,
                    haplotypes, !hcArgs.disableSpanningEventGenotyping);

            final List<VariantContext> eventsAtThisLocWithSpanDelsReplaced = replaceSpanDels(eventsAtThisLoc,
                    Allele.create(ref[loc - refLoc.getStart()], true), loc);

            VariantContext mergedVC = AssemblyBasedCallerUtils.makeMergedVariantContext(eventsAtThisLocWithSpanDelsReplaced);

            if( mergedVC == null ) {
                continue;
            }

            final Map<Allele, List<Haplotype>> alleleMapper = AssemblyBasedCallerUtils.createAlleleMapper(mergedVC, loc, haplotypes, !hcArgs.disableSpanningEventGenotyping);

            if( hcArgs.assemblerArgs.debugAssembly && logger != null ) {
                logger.info("Genotyping event at " + loc + " with alleles = " + mergedVC.getAlleles());
            }

            mergedVC = removeAltAllelesIfTooManyGenotypes(ploidy, alleleMapper, mergedVC);

            int mergedAllelesListSizeBeforePossibleTrimming = mergedVC.getAlleles().size();

            // Note: read-allele likelihoods are not used in joint detection genotyping, but they are always used for annotations
            AlleleLikelihoods<GATKRead, Allele> readAlleleLikelihoods = readLikelihoods.marginalize(alleleMapper);
            final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
            final SimpleInterval variantCallingRelevantOverlap = new SimpleInterval(mergedVC).expandWithinContig(hcArgs.informativeReadOverlapMargin, sequenceDictionary);
            readAlleleLikelihoods.retainEvidence(read -> readQualifiesForGenotypingPredicate.test(read, variantCallingRelevantOverlap));
            readAlleleLikelihoods.setVariantCallingSubsetUsed(variantCallingRelevantOverlap);

            if (configuration.isSampleContaminationPresent()) {
                // This warrants future evaluation as to the best way to handle disqualified reads
                if (hcArgs.applyBQD || hcArgs.applyFRD) {
                    DRAGENConaminationWarning.warn("\\n=============================================================================" +
                            "Sample contamination specified with FRD/BQD enabled. Contamination calling with either BQD or FRD genotyping models enabled is currently unsupported and may produce unexpected results. Use at your own risk." +
                            "\n=============================================================================");
                }

                readAlleleLikelihoods.contaminationDownsampling(configuration.getSampleContamination());
            }

            // note: this debugging does not pertain to joint detection, where read-allele likelihoods are never computed
            if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
                HaplotypeCallerGenotypingDebugger.println("\n=============================================================================");
                HaplotypeCallerGenotypingDebugger.println("Event at: " + mergedVC + " with " + readAlleleLikelihoods.evidenceCount() + " reads and "+readAlleleLikelihoods.filteredSampleEvidence(0).size()+" disqualified");
                HaplotypeCallerGenotypingDebugger.println("=============================================================================");
                HaplotypeCallerGenotypingDebugger.println("haplotype alleles key:");
                for (Map.Entry<Allele, List<Haplotype>> allele : alleleMapper.entrySet()) {
                    HaplotypeCallerGenotypingDebugger.println("Allele: "+allele.getKey()+" Haps: "+allele.getValue().stream().map(readLikelihoods::indexOfAllele).map(i -> Integer.toString(i)).collect(Collectors.joining(", ")));
                }
                HaplotypeCallerGenotypingDebugger.println("Read-allele matrix:");
                String allele_string = readAlleleLikelihoods.alleles().stream().map(al -> al.toString()).collect(Collectors.joining(" "));
                HaplotypeCallerGenotypingDebugger.println(allele_string);
                for (int sn = 0 ; sn < readAlleleLikelihoods.numberOfSamples(); sn++){
                    for (int evn = 0 ; evn < readAlleleLikelihoods.sampleEvidence(sn).size(); evn++) {
                        String outputStr = "read: " + readLikelihoods.sampleEvidence(sn).indexOf(readAlleleLikelihoods.sampleEvidence(sn).get(evn)) + " " + readAlleleLikelihoods.sampleEvidence(sn).get(evn).getName();

                        for (Allele curAllele : readAlleleLikelihoods.alleles()) {
                            int idx = readAlleleLikelihoods.indexOfAllele(curAllele);
                            outputStr = outputStr + " " + readAlleleLikelihoods.sampleMatrix(sn).get(idx, evn);
                        }
                        HaplotypeCallerGenotypingDebugger.println(outputStr);
                    }
                }

                HaplotypeCallerGenotypingDebugger.println("Normalized Read-Allele matrix:");
                for (int sn = 0 ; sn < readAlleleLikelihoods.numberOfSamples(); sn++){
                    for (int evn = 0 ; evn < readAlleleLikelihoods.sampleEvidence(sn).size(); evn++) {
                        String outputStr = "read: " + readLikelihoods.sampleEvidence(sn).indexOf(readAlleleLikelihoods.sampleEvidence(sn).get(evn)) + " " + readAlleleLikelihoods.sampleEvidence(sn).get(evn).getName();

                        double max = Double.NEGATIVE_INFINITY;
                        for (Allele curAllele : readAlleleLikelihoods.alleles()) {
                            int idx = readAlleleLikelihoods.indexOfAllele(curAllele);
                            max = Math.max(readAlleleLikelihoods.sampleMatrix(sn).get(idx, evn), max);
                        }

                        for (Allele curAllele : readAlleleLikelihoods.alleles()) {
                            int idx = readAlleleLikelihoods.indexOfAllele(curAllele);
                            outputStr = outputStr + " " + (readAlleleLikelihoods.sampleMatrix(sn).get(idx, evn) - max);
                        }
                        HaplotypeCallerGenotypingDebugger.println(outputStr);
                    }
                }
            }   // end of read-allele likelihoods debugging messages

            // TODO: does this work if we have the haplo-genotype override below?
            if (emitReferenceConfidence) {
                mergedVC = ReferenceConfidenceUtils.addNonRefSymbolicAllele(mergedVC);
                readAlleleLikelihoods.addNonReferenceAllele(Allele.NON_REF_ALLELE);
                mergedAllelesListSizeBeforePossibleTrimming++;
            }

            // use joint detection haplo-genotype posteriors if available
            final List<Pair<GenotypingLikelihoods<Haplotype>, GenotypingLikelihoods<Haplotype>>> haploGTLikelihoodsAndPosteriors = new ArrayList<>(haploGTLikelihoodAndPosteriorOverlapDetector.getOverlaps(mergedVC));
            Utils.validate(haploGTLikelihoodsAndPosteriors.size() < 2, "At most one set of haplotype genotype posteriors should overlap this variant.");
            final Optional<GenotypingLikelihoods<Allele>> genotypeLikelihoodOverride = (haploGTLikelihoodsAndPosteriors.size() == 1 && ploidy == 2) ?
                    Optional.of(computeAlleleGenotypesFromHaploGenotypes(ploidy, noCallAlleles, loc, mergedVC, haploGTLikelihoodsAndPosteriors.get(0).getLeft())) : Optional.empty();
            final Optional<GenotypingLikelihoods<Allele>> genotypePosteriorOverride = (haploGTLikelihoodsAndPosteriors.size() == 1 && ploidy == 2) ?
                    Optional.of(computeAlleleGenotypesFromHaploGenotypes(ploidy, noCallAlleles, loc, mergedVC, haploGTLikelihoodsAndPosteriors.get(0).getRight())) : Optional.empty();

            // these genotypes have the PLs calculated and filled out but are otherwise empty (no-call alleles, no annotations)
            final GenotypesContext genotypes = calculateGLsForThisEvent(readAlleleLikelihoods, mergedVC, noCallAlleles, ref,
                    loc - refLoc.getStart(), dragstrs, genotypeLikelihoodOverride, genotypePosteriorOverride);

            // In joint detection we have already accounted for *haplotype* priors, no need to double-count priors here!
            final GenotypePriorCalculator gpc = genotypeLikelihoodOverride.isPresent() ? GenotypePriorCalculator.flatPrior() :
                    resolveGenotypePriorCalculator(dragstrs, loc - refLoc.getStart() + 1, snpHeterozygosity, indelHeterozygosity);
            final VariantContext call = calculateGenotypes(new VariantContextBuilder(mergedVC).genotypes(genotypes).make(), gpc, givenAlleles);

            // TODO: looks like even in joint detection we will need the read-allele likelihoods for annotations
            if( call != null ) {
                readAlleleLikelihoods = prepareReadAlleleLikelihoodsForAnnotation(readLikelihoods, perSampleFilteredReadList,
                        emitReferenceConfidence, alleleMapper, readAlleleLikelihoods, call, variantCallingRelevantOverlap);

                VariantContext annotatedCall = makeAnnotatedCall(ref, refLoc, tracker, header, mergedVC,
                        mergedAllelesListSizeBeforePossibleTrimming, readAlleleLikelihoods, call, annotationEngine, preFilteringAlleleLikelihoods);

                if (suspiciousLocations.contains(loc)){
                    annotatedCall.getCommonInfo().putAttribute(GATKVCFConstants.POSSIBLE_FP_ADJACENT_TP_KEY, true);
                }


                if (dragstrs != null && GATKVariantContextUtils.containsInlineIndel(annotatedCall)) {
                    final int strOffset = loc - refLoc.getStart() + 1;
                    annotatedCall = DragstrVariantContextAnnotations.annotateVariantContextWithDragstrParametersUsed(annotatedCall, dragstrParams, dragstrs.period(strOffset), dragstrs.repeatLength(strOffset));
                }
                returnCalls.add( annotatedCall );

                if (withBamOut) {
                    AssemblyBasedCallerUtils.annotateReadLikelihoodsWithSupportedAlleles(call, readAlleleLikelihoods);
                }
                // maintain the set of all called haplotypes
                call.getAlleles().stream().map(alleleMapper::get).filter(Objects::nonNull).forEach(calledHaplotypes::addAll);
            }
        }

        final List<VariantContext> phasedCalls = doPhysicalPhasing ? AssemblyBasedCallerUtils.phaseCalls(returnCalls, calledHaplotypes) : returnCalls;
        return new CalledHaplotypes(phasedCalls, calledHaplotypes);
    }

    /**
     *
     * @param ploidy
     * @param noCallAlleles
     * @param loc
     * @param mergedVC
     * @param haplotypePosteriors
     * @return
     */
    private GenotypingLikelihoods<Allele> computeAlleleGenotypesFromHaploGenotypes(int ploidy, List<Allele> noCallAlleles, int loc, VariantContext mergedVC, GenotypingLikelihoods<Haplotype> haplotypePosteriors) {
        final List<Haplotype> determinedHaplotypes = haplotypePosteriors.asListOfAlleles(); // this only contains haplotypes whose determined span contains this locus
        final Ordering<Haplotype> haplotypeOrdering = Ordering.explicit(determinedHaplotypes);
        final Map<Allele, List<Haplotype>> alleleToHaplotypes = AssemblyBasedCallerUtils.createAlleleMapper(mergedVC, loc, determinedHaplotypes, !hcArgs.disableSpanningEventGenotyping);

        // map from haplotype pair to corresponding haplo-genotype index (in the canonical GL order)
        final Map<List<Haplotype>, Integer> haploPairToGTIndex = Utils.stream(GenotypeAlleleCounts.iterable(ploidy, haplotypePosteriors.numberOfAlleles()))
                .collect(Collectors.toMap(gac -> gac.asAlleleList(determinedHaplotypes), gac -> gac.index()));

        // map from allele-genotype index to list of compatible haplo-genotype indices
        final Map<Integer, int[]> alleleGTIndexToHaploGTIndices = new HashMap<>();
        for (final GenotypeAlleleCounts gac : GenotypeAlleleCounts.iterable(ploidy, mergedVC.getNAlleles())) {
            final List<Allele> alleleList = gac.asAlleleList(mergedVC.getAlleles());
            final int[] haploGTIndices = Lists.cartesianProduct(alleleList.stream().map(alleleToHaplotypes::get).toList()).stream()
                            .map(haplotypeOrdering::sortedCopy).mapToInt(haploPairToGTIndex::get).toArray();
            alleleGTIndexToHaploGTIndices.put(gac.index(), haploGTIndices);
        }

        // for each sample, traverse the allele-based GACs and add up the log-space posteriors of all corresponding haplotype-based posterior GLs
        final List<GenotypeLikelihoods> genotypeLikelihoodsInSampleOrder = new ArrayList<>();
        for (int sampleIndex = 0; sampleIndex < haplotypePosteriors.numberOfSamples(); sampleIndex++) {
            final double[] log10HaploGenotypePosteriors = haplotypePosteriors.sampleLikelihoods(sampleIndex).getAsVector();

            final double[] jointDetectionPosteriorGLs = IntStream.range(0, GenotypeIndexCalculator.genotypeCount(ploidy, mergedVC.getNAlleles()))
                    .mapToObj(alleleGTIndex -> alleleGTIndexToHaploGTIndices.get(alleleGTIndex))
                    .map(haploGTIndices -> MathUtils.applyToArray(haploGTIndices, idx -> log10HaploGenotypePosteriors[idx]))
                    .mapToDouble(MathUtils::log10SumLog10).toArray();

            genotypeLikelihoodsInSampleOrder.add(GenotypeLikelihoods.fromLog10Likelihoods(jointDetectionPosteriorGLs));
        }

        return new GenotypingLikelihoods<>(new IndexedAlleleList<>(mergedVC.getAlleles()), ploidyModel, genotypeLikelihoodsInSampleOrder);
    }

    // returns pairs of GenotypingLikelihoods of haplo-genotypes, the first being actual likelihoods and the second being
    // posteriors -- the same object with priors applied
    private OverlapDetector<Pair<GenotypingLikelihoods<Haplotype>, GenotypingLikelihoods<Haplotype>>> computeHaploGenotypeLikelihoodsAndPosteriors(AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods, byte[] ref, SimpleInterval refLoc, SimpleInterval activeRegionWindow, SAMFileHeader header, int ploidy, DragstrReferenceAnalyzer dragstrs, BiPredicate<GATKRead, SimpleInterval> readQualifiesForGenotypingPredicate) {
        final OverlapDetector<Pair<GenotypingLikelihoods<Haplotype>, GenotypingLikelihoods<Haplotype>>> haploGenotypeLikelihoodAndPosteriorOverlapDetector = new OverlapDetector<>(0,0);
        // if haplotypes are partially determined, do DRAGEN exact genotyping within each determined span
        final List<Haplotype> haplotypes = readLikelihoods.alleles();
        if (!haplotypes.isEmpty() && haplotypes.get(0).isPartiallyDetermined()) {

            final Map<SimpleInterval, List<Haplotype>> jdGroups = haplotypes.stream().collect(Collectors.groupingBy(h -> ((PartiallyDeterminedHaplotype) h).getDeterminedSpan()));
            final List<SimpleInterval> jdIntervals = jdGroups.keySet().stream().sorted(Comparator.comparingInt(SimpleInterval::getStart)).toList();

            for (final SimpleInterval jdInterval : jdIntervals) {
                if( !jdInterval.overlaps(activeRegionWindow)) {
                    continue;
                }

                final List<Haplotype> jdHaplotypes = jdGroups.get(jdInterval);
                final AlleleLikelihoods<GATKRead, Haplotype> jdReadLikelihoods = readLikelihoods.removeAllelesToSubset(jdHaplotypes);

                // TODO: this mimics the read subsetting around to a padded interval around the merged VC in the allele genotyping code below
                // TODO: do we also wish to subset to overlapping reads here?
                final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
                final SimpleInterval variantCallingRelevantOverlap = jdInterval.expandWithinContig(hcArgs.informativeReadOverlapMargin, sequenceDictionary);
                jdReadLikelihoods.retainEvidence(read -> readQualifiesForGenotypingPredicate.test(read, variantCallingRelevantOverlap));

                // NOTE: we do not want the DRAGEN genotypes model here because it includes FRD and BQD, which only apply once
                // we have genotypes based on alleles, not haplotypes
                final GenotypingLikelihoods<Haplotype> log10HaploGenotypeLikelihoods = new IndependentSampleGenotypesModel().calculateLikelihoods(jdReadLikelihoods,
                        new GenotypingData<>(ploidyModel, jdReadLikelihoods), ref, jdInterval.getStart() - refLoc.getStart(), dragstrs);

                // calculate heterozygosity for each event
                final List<Event> allEvents = jdHaplotypes.stream().flatMap(haplotype -> haplotype.getEventMap().getEvents().stream()).distinct().toList();
                final Map<Event, Double> log10EventHetPriors = new HashMap<>();
                for (final Event event : allEvents) {
                    final double log10SNPHetPrior = Math.log10(snpHeterozygosity);
                    if (event.isSNP()) {
                        log10EventHetPriors.put(event, log10SNPHetPrior);
                    } else {
                        // TODO: is this right?
                        final int pos = event.getStart() - refLoc.getStart() + 1;
                        final boolean noDragstr = hcArgs.likelihoodArgs.dragstrParams == null || hcArgs.standardArgs.genotypeArgs.dontUseDragstrPriors || dragstrs == null;
                        final double log10IndelHetPrior = noDragstr ? Math.log10(indelHeterozygosity) :
                                -.1 * dragstrParams.api(dragstrs.period(pos), dragstrs.repeatLength(pos));
                        log10EventHetPriors.put(event, event.isIndel() ? log10IndelHetPrior : Math.max(log10IndelHetPrior, log10SNPHetPrior));
                    }
                }

                // TODO: we compute this even for all ref haplotype -- that seems sketchy!!
                // calculate haplotype heterozygosities by adding up (in log space) their event heterozygosities
                final Map<Haplotype, Double> log10HaplotypeHetPriors = jdHaplotypes.stream().collect(Collectors.toMap(h->h,
                        h -> h.getEventMap().getEvents().stream().mapToDouble(log10EventHetPriors::get).sum()));

                // TODO: likewise: we are calculating this for the ref haplotype. . .
                // note that DRAGEN's prior for homozygosity is NOT the same as Hardy-Weinberg equilibrium!
                // TODO: should the DRAGEN het-hom ratio for haplotypes be the same as for alleles?
                final Map<Haplotype, Double> log10HaplotypeHomPriors = jdHaplotypes.stream().collect(Collectors.toMap(h->h,
                        h -> log10HaplotypeHetPriors.get(h) - Math.log10(hcArgs.likelihoodArgs.dragstrHetHomRatio)));

                // calculate diploid haplotype-based genotype priors
                final double[] log10HaploGenotypePriors = new double[GenotypeIndexCalculator.genotypeCount(ploidy, log10HaploGenotypeLikelihoods.numberOfAlleles())];

                for (final GenotypeAlleleCounts gac : GenotypeAlleleCounts.iterable(ploidy, log10HaploGenotypeLikelihoods.numberOfAlleles())) {
                    // note: the zeroth element here is not necessarily the ref-ref haplo-genotype
                    log10HaploGenotypePriors[gac.index()] = gac.sumOverAlleleIndicesAndCounts((allele, count) -> count == ploidy ? log10HaplotypeHomPriors.get(log10HaploGenotypeLikelihoods.getAllele(allele))
                            : count * log10HaplotypeHetPriors.get(log10HaploGenotypeLikelihoods.getAllele(allele)));
                }

                // add haplo-genotype priors to haplo-genotype likelihoods in log space and normalize to get posteriors
                final List<GenotypeLikelihoods> log10HaploGenotypePosteriorsInSampleOrder = IntStream.range(0, log10HaploGenotypeLikelihoods.numberOfSamples()).boxed()
                        .map(s -> MathUtils.normalizeLog10(MathUtils.ebeAdd(log10HaploGenotypePriors, log10HaploGenotypeLikelihoods.sampleLikelihoods(s).getAsVector())))
                        .map(GenotypeLikelihoods::fromLog10Likelihoods)
                        .toList();

                GenotypingLikelihoods<Haplotype> log10HaploGenotypePosteriors = new GenotypingLikelihoods<>(log10HaploGenotypeLikelihoods, ploidyModel, log10HaploGenotypePosteriorsInSampleOrder);
                haploGenotypeLikelihoodAndPosteriorOverlapDetector.addLhs(Pair.of(log10HaploGenotypeLikelihoods, log10HaploGenotypePosteriors), jdInterval);
            }
        }
        return haploGenotypeLikelihoodAndPosteriorOverlapDetector;
    }

    /**
     * If there is potential to use DRAGstr in the region based on the event map, then this method composes the
     *   STR finder.
     * @param haplotypes reconstructed haplotypes.
     * @param ref the reference haplotype sequence.
     * @param refLoc the interval of the sequence provided.
     * @param startPosKeySet the location of events within the haplotypes on the reference sequence.
     * @return {@code null} iff there is no chance that we would be using DRAGstr in this region based
     *    on the reconstructed haplotypes.
     */
    private  DragstrReferenceAnalyzer constructDragstrReferenceSTRAnalyzerIfNecessary(final List<Haplotype> haplotypes,
                                                                                      final byte[] ref,
                                                                                      final SimpleInterval refLoc,
                                                                                      final SortedSet<Integer> startPosKeySet) {
        if (isDragstrSTRAnalyzerNecessary(startPosKeySet, haplotypes)) {
            final int offset = startPosKeySet.first() - refLoc.getStart();
            final int to = startPosKeySet.last() - refLoc.getStart() + 2;
            // +2 = +1+1
            // where one +1 is because starPosKeySet indexes are 1-based and offset/to are 0-based.
            //   and the other +1 is because we need to analyze/include one base after each event position including the last.
            return  DragstrReferenceAnalyzer.of(ref, offset, to, dragstrParams.maximumPeriod());
        } else {
            return null;
        }
    }

    /**
     * Confirms whether there is the need to analyze the region's reference sequence for the presence of STRs.
     * <p>
     *     This is only the case when DRAGstr is activate, we are going to use their priors and there is some indel
     *     amongst the haplotypes.
     * </p>
     */
    private boolean isDragstrSTRAnalyzerNecessary(SortedSet<Integer> startPosKeySet, List<Haplotype> haplotypes) {
        return !startPosKeySet.isEmpty() && dragstrParams != null
                && !hcArgs.standardArgs.genotypeArgs.dontUseDragstrPriors &&
                haplotypes.stream()
                        .anyMatch(h -> h.getEventMap().getEvents().stream()
                                .anyMatch(e -> GATKVariantContextUtils.containsInlineIndel(e.refAllele(), Collections.singletonList(e.altAllele()))));
    }


    /**
     * Composes the appropriate test to determine if a read is to be retained for evidence/likelihood calculation for variants
     * located in a target region.
     * @param hcArgs configuration that may affect the criteria use to retain or filter-out reads.
     * @return never {@code null}.
     */
    private BiPredicate<GATKRead, SimpleInterval> composeReadQualifiesForGenotypingPredicate(final HaplotypeCallerArgumentCollection hcArgs) {
        if (hcArgs.applyBQD || hcArgs.applyFRD) {
            return (read, target) -> softUnclippedReadOverlapsInterval(read, target);
        } else {
            // NOTE: we must make this comparison in target -> read order because occasionally realignment/assembly produces
            // reads that consume no reference bases and this can cause them to overlap adjacent
            return (read, target) -> target.overlaps(read);
        }
    }

    /**
     * Checks whether a read's extended mapping region (unclipping soft-clips) overlaps a given target interval
     * even if it is only by one base. Adjacency is not good enough.
     * @param read the read to test.
     * @param target the interval to test.
     * @return {@code true} iff there is an overlap.
     */
    private boolean softUnclippedReadOverlapsInterval(final GATKRead read, final Locatable target) {
        return read.getContig().equalsIgnoreCase(target.getContig())
                && read.getSoftStart() <= target.getEnd()
                && read.getSoftEnd() >= target.getStart()
                && read.getSoftStart() <= read.getSoftEnd(); // is this possible, ever? this test was performed before extracting this method so we keep it just in case.
    }

    private GenotypePriorCalculator resolveGenotypePriorCalculator(final DragstrReferenceAnalyzer strs, final int pos,
                                                                   final double snpHeterozygosity, final double indelHeterozygosity) {
        if (hcArgs.likelihoodArgs.dragstrParams == null || hcArgs.standardArgs.genotypeArgs.dontUseDragstrPriors) {
            return GenotypePriorCalculator.assumingHW(Math.log10(snpHeterozygosity), Math.log10(indelHeterozygosity));
        } else if (strs == null) {
            return GenotypePriorCalculator.givenHetToHomRatio(Math.log10(snpHeterozygosity), Math.log10(indelHeterozygosity), Math.log10(Math.min(snpHeterozygosity, indelHeterozygosity)), hcArgs.likelihoodArgs.dragstrHetHomRatio);
        } else {
            final int period = strs.period(pos);
            final int repeats = strs.repeatLength(pos);
            return GenotypePriorCalculator.givenDragstrParams(dragstrParams, period, repeats, Math.log10(snpHeterozygosity), hcArgs.likelihoodArgs.dragstrHetHomRatio);
        }
    }

    @VisibleForTesting
    static public List<VariantContext> replaceSpanDels(final List<VariantContext> eventsAtThisLoc, final Allele refAllele, final int loc) {
        return eventsAtThisLoc.stream().map(vc -> replaceWithSpanDelVC(vc, refAllele, loc)).collect(Collectors.toList());
    }

    @VisibleForTesting
    static VariantContext replaceWithSpanDelVC(final VariantContext variantContext, final Allele refAllele, final int loc) {
        if (variantContext.getStart() == loc) {
            return variantContext;
        } else {
            VariantContextBuilder builder = new VariantContextBuilder(variantContext)
                    .start(loc)
                    .stop(loc)
                    .alleles(Arrays.asList(refAllele, Allele.SPAN_DEL))
                    .genotypes(GenotypesContext.NO_GENOTYPES);
            return builder.make();
        }

    }

    /**
     * If the number of alleles is so high that enumerating all possible genotypes is impractical, as determined by
     * {@link #maxGenotypeCountToEnumerate}, remove alt alleles from the input {@code alleleMapper} that are
     * not well supported by good-scored haplotypes.
     * Otherwise do nothing.
     *
     * Alleles kept are guaranteed to have higher precedence than those removed, where precedence is determined by
     * {@link AlleleScoredByHaplotypeScores}.
     *
     * After the remove operation, entries in map are guaranteed to have the same relative order as they were in the input map,
     * that is, entries will be only be removed but not not shifted relative to each other.
     *  @param ploidy        ploidy of the sample
     * @param alleleMapper  original allele to haplotype map
     */
    private VariantContext removeAltAllelesIfTooManyGenotypes(final int ploidy, final Map<Allele, List<Haplotype>> alleleMapper, final VariantContext mergedVC) {

        final int originalAlleleCount = alleleMapper.size();
        practicalAlleleCountForPloidy.putIfAbsent(ploidy, GenotypeIndexCalculator.computeMaxAcceptableAlleleCount(ploidy, maxGenotypeCountToEnumerate));
        final int practicalAlleleCount = practicalAlleleCountForPloidy.get(ploidy);

        if (originalAlleleCount > practicalAlleleCount) {
            final List<Allele> allelesToKeep = whichAllelesToKeepBasedonHapScores(alleleMapper, practicalAlleleCount);
            alleleMapper.keySet().retainAll(allelesToKeep);
            logger.warn(String.format("At position %s removed alt alleles where ploidy is %d and original allele count is %d, whereas after trimming the allele count becomes %d. Alleles kept are:%s",
                    new SimpleInterval(mergedVC).toString(), ploidy, originalAlleleCount, practicalAlleleCount, allelesToKeep));
            return removeExcessAltAllelesFromVC(mergedVC, allelesToKeep);
        } else {
            return mergedVC;
        }
    }

    /**
     * Returns a list of alleles that is a subset of the key set of input map {@code alleleMapper}.
     * The size of the returned list is min({@code desiredNumOfAlleles}, alleleMapper.size()).
     *
     * Alleles kept are guaranteed to have higher precedence than those removed, where precedence is determined by
     * {@link AlleleScoredByHaplotypeScores}.
     *
     * Entries in the returned list are guaranteed to have the same relative order as they were in the input map.
     *
     * @param alleleMapper          original allele to haplotype map
     * @param desiredNumOfAlleles   desired allele count, including ref allele
     */
    @VisibleForTesting
    static List<Allele> whichAllelesToKeepBasedonHapScores(final Map<Allele, List<Haplotype>> alleleMapper,
                                                           final int desiredNumOfAlleles) {

        if(alleleMapper.size() <= desiredNumOfAlleles){
            return alleleMapper.keySet().stream().collect(Collectors.toList());
        }

        final PriorityQueue<AlleleScoredByHaplotypeScores> alleleMaxPriorityQ = new PriorityQueue<>();
        for(final Allele allele : alleleMapper.keySet()){
            final List<Double> hapScores = alleleMapper.get(allele).stream().map(Haplotype::getScore).sorted().collect(Collectors.toList());
            final Double highestScore = hapScores.size() > 0 ? hapScores.get(hapScores.size()-1) : Double.NEGATIVE_INFINITY;
            final Double secondHighestScore = hapScores.size()>1 ? hapScores.get(hapScores.size()-2) : Double.NEGATIVE_INFINITY;

            alleleMaxPriorityQ.add(new AlleleScoredByHaplotypeScores(allele, highestScore, secondHighestScore));
        }

        final Set<Allele> allelesToRetain = new LinkedHashSet<>();
        while(allelesToRetain.size()<desiredNumOfAlleles){
            allelesToRetain.add(alleleMaxPriorityQ.poll().getAllele());
        }
        return alleleMapper.keySet().stream().filter(allelesToRetain::contains).collect(Collectors.toList());
    }

    /**
     * A utility class that provides ordering information, given best and second best haplotype scores.
     * If there's a tie between the two alleles when comparing their best haplotype score, the second best haplotype score
     * is used for breaking the tie. In the case that one allele doesn't have a second best allele, i.e. it has only one
     * supportive haplotype, its second best score is set as {@link Double#NEGATIVE_INFINITY}.
     * In the extremely unlikely cases that two alleles, having the same best haplotype score, neither have a second
     * best haplotype score, or the same second best haplotype score, the order is exactly the same as determined by
     * the Allele class's comparison
     */
    private static final class AlleleScoredByHaplotypeScores implements Comparable<AlleleScoredByHaplotypeScores>{
        private final Allele allele;
        private final Double bestHaplotypeScore;
        private final Double secondBestHaplotypeScore;

        public AlleleScoredByHaplotypeScores(final Allele allele, final Double bestHaplotypeScore, final Double secondBestHaplotypeScore){
            this.allele = allele;
            this.bestHaplotypeScore = bestHaplotypeScore;
            this.secondBestHaplotypeScore = secondBestHaplotypeScore;
        }

        @Override
        public int compareTo(final AlleleScoredByHaplotypeScores other) {

            if(allele.isReference() && other.allele.isNonReference()){
                return -1;
            } else if(allele.isNonReference() && other.allele.isReference()){
                return 1;
            } else if(bestHaplotypeScore > other.bestHaplotypeScore) {
                return -1;
            } else if (bestHaplotypeScore < other.bestHaplotypeScore) {
                return 1;
            } else if (!secondBestHaplotypeScore.equals(other.secondBestHaplotypeScore)) {
                return secondBestHaplotypeScore > other.secondBestHaplotypeScore ? -1 : 1;
            } else {
                return allele.compareTo(other.allele);
            }
        }

        public Allele getAllele(){
            return allele;
        }
    }

    /**
     * Returns an VC that is similar to {@code inputVC} in every aspect except that alleles not in {@code allelesToKeep}
     * are removed in the returned VC.
     * @throws IllegalArgumentException if 1) {@code allelesToKeep} is null or contains null elements; or
     *                                     2) {@code allelesToKeep} doesn't contain a reference allele; or
     *                                     3) {@code allelesToKeep} is not a subset of {@code inputVC.getAlleles()}
     */
    @VisibleForTesting
    static VariantContext removeExcessAltAllelesFromVC(final VariantContext inputVC, final Collection<Allele> allelesToKeep){
        Utils.validateArg(allelesToKeep!=null, "alleles to keep is null");
        Utils.validateArg(!allelesToKeep.contains(null), "alleles to keep contains null elements");
        Utils.validateArg(allelesToKeep.stream().anyMatch(Allele::isReference), "alleles to keep doesn't contain reference allele!");
        Utils.validateArg(inputVC.getAlleles().containsAll(allelesToKeep), "alleles to keep is not a subset of input VC alleles");
        if(inputVC.getAlleles().size() == allelesToKeep.size()) return inputVC;

        final VariantContextBuilder vcb = new VariantContextBuilder(inputVC);
        final List<Allele> originalList = inputVC.getAlleles();
        originalList.retainAll(allelesToKeep);
        vcb.alleles(originalList);
        return vcb.make();
    }

    @VisibleForTesting
    static protected VariantContext makeAnnotatedCall(byte[] ref, SimpleInterval refLoc, FeatureContext tracker, SAMFileHeader header, VariantContext mergedVC, int mergedAllelesListSizeBeforePossibleTrimming, AlleleLikelihoods<GATKRead, Allele> readAlleleLikelihoods, VariantContext call, VariantAnnotatorEngine annotationEngine, final AlleleLikelihoods<GATKRead, Haplotype> preFilteringAlleleLikelihoods) {
        final SimpleInterval locus = new SimpleInterval(mergedVC);
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        final SimpleInterval refLocInterval= new SimpleInterval(refLoc);
        final ReferenceDataSource refData = new ReferenceMemorySource(new ReferenceBases(ref, refLocInterval), sequenceDictionary);
        final ReferenceContext referenceContext = new ReferenceContext(refData, locus, refLocInterval);

        final VariantContext untrimmedResult =  annotationEngine.annotateContext(call, tracker, referenceContext, readAlleleLikelihoods, Optional.empty(), Optional.empty(), Optional.ofNullable(preFilteringAlleleLikelihoods), a -> true);

        // propagate the tag indicating that the VC was collapsed
        if ( mergedVC.getAttribute(AssemblyBasedCallerUtils.EXT_COLLAPSED_TAG) != null ) {
            untrimmedResult.getCommonInfo().putAttribute(AssemblyBasedCallerUtils.EXT_COLLAPSED_TAG,
                    mergedVC.getAttribute(AssemblyBasedCallerUtils.EXT_COLLAPSED_TAG));
        }
        // NOTE: We choose to reverseTrimAlleles() here as opposed to when we actually do the trimming because otherwise we would have to resolve
        //       the mismatching readAlleleLikelihoods object which is keyed to the old, possibly incorrectly trimmed alleles.
        return untrimmedResult.getAlleles().size() == mergedAllelesListSizeBeforePossibleTrimming ? untrimmedResult
                : GATKVariantContextUtils.reverseTrimAlleles(untrimmedResult);
    }

    /**
     * For a particular event described in inputVC, form PL vector for each sample by looking into allele read map and filling likelihood matrix for each allele
     * @param readLikelihoods          Allele map describing mapping from reads to alleles and corresponding likelihoods
     * @param mergedVC               Input VC with event to genotype
     * @parma glsOverride           GenotypingLikelihoods calculated elsewhere, such as from joint detection haplo-genotyping, in which case we omit
     *                              the GL calculation and only do additional modifications such as DRAGEN BQD and FRD
     * @return                       GenotypesContext object wrapping genotype objects with PLs
     */
    protected GenotypesContext calculateGLsForThisEvent(final AlleleLikelihoods<GATKRead, Allele> readLikelihoods,
                                                        final VariantContext mergedVC, final List<Allele> noCallAlleles,
                                                        final byte[] paddedReference, final int offsetForRefIntoEvent,
                                                        final DragstrReferenceAnalyzer dragstrs, Optional<GenotypingLikelihoods<Allele>> genotypeLikelihoodsOverride,
                                                        Optional<GenotypingLikelihoods<Allele>> genotypePosteriorsOverride) {
        Utils.nonNull(readLikelihoods, "readLikelihoods");
        Utils.nonNull(mergedVC, "mergedVC");
        final List<Allele> vcAlleles = mergedVC.getAlleles();
        final AlleleList<Allele> alleleList = readLikelihoods.numberOfAlleles() == vcAlleles.size() ? readLikelihoods : new IndexedAlleleList<>(vcAlleles);
        final GenotypingLikelihoods<Allele> likelihoods = genotypingModel.calculateLikelihoods(alleleList,new GenotypingData<>(ploidyModel,readLikelihoods),
                paddedReference, offsetForRefIntoEvent, dragstrs, genotypeLikelihoodsOverride, genotypePosteriorsOverride);
        final int sampleCount = samples.numberOfSamples();
        final GenotypesContext result = GenotypesContext.create(sampleCount);
        for (int s = 0; s < sampleCount; s++) {
            result.add(new GenotypeBuilder(samples.getSample(s)).alleles(noCallAlleles).PL(likelihoods.sampleLikelihoods(s).getAsPLs()).make());
        }
        return result;
    }

    /**
     * Returns the ploidy-model used by this genotyping engine.
     *
     * @return never {@code null}.
     */
    public PloidyModel getPloidyModel() {
        return ploidyModel;
    }

    // Builds the read-likelihoods collection to use for annotation considering user arguments and the collection
    // used for genotyping.
    private AlleleLikelihoods<GATKRead, Allele> prepareReadAlleleLikelihoodsForAnnotation(
            final AlleleLikelihoods<GATKRead, Haplotype> readHaplotypeLikelihoods,
            final Map<String, List<GATKRead>> perSampleFilteredReadList,
            final boolean emitReferenceConfidence,
            final Map<Allele, List<Haplotype>> alleleMapper,
            final AlleleLikelihoods<GATKRead, Allele> readAlleleLikelihoodsForGenotyping,
            final VariantContext call,
            final SimpleInterval relevantReadsOverlap) {

        final AlleleLikelihoods<GATKRead, Allele> readAlleleLikelihoodsForAnnotations;

        // We can reuse for annotation the likelihood for genotyping as long as there is no contamination filtering
        // or the user want to use the contamination filtered set for annotations.
        // Otherwise (else part) we need to do it again.
        if (hcArgs.useFilteredReadMapForAnnotations || !configuration.isSampleContaminationPresent()) {
            readAlleleLikelihoodsForAnnotations = readAlleleLikelihoodsForGenotyping;
            // the input likelihoods are supposed to have been filtered to only overlapping reads so no need to
            // do it again.
        } else {
            readAlleleLikelihoodsForAnnotations = readHaplotypeLikelihoods.marginalize(alleleMapper);
            readAlleleLikelihoodsForAnnotations.retainEvidence(relevantReadsOverlap::overlaps);
            if (emitReferenceConfidence) {
                readAlleleLikelihoodsForAnnotations.addNonReferenceAllele(Allele.NON_REF_ALLELE);
            }
        }

        if (call.getAlleles().size() != readAlleleLikelihoodsForAnnotations.numberOfAlleles()) {
            readAlleleLikelihoodsForAnnotations.updateNonRefAlleleLikelihoods(new IndexedAlleleList<>(call.getAlleles()));
        }

        // Skim the filtered map based on the location so that we do not add filtered read that are going to be removed
        // right after a few lines of code below.
        final Map<String, List<GATKRead>> overlappingFilteredReads = overlappingFilteredReads(perSampleFilteredReadList, relevantReadsOverlap);

        readAlleleLikelihoodsForAnnotations.addEvidence(overlappingFilteredReads,0);

        return readAlleleLikelihoodsForAnnotations;
    }

    private static Map<String, List<GATKRead>> overlappingFilteredReads(final Map<String, List<GATKRead>> perSampleFilteredReadList, final SimpleInterval loc) {
        final Map<String,List<GATKRead>> overlappingFilteredReads = new HashMap<>(perSampleFilteredReadList.size());

        for (final Map.Entry<String,List<GATKRead>> sampleEntry : perSampleFilteredReadList.entrySet()) {
            final List<GATKRead> originalList = sampleEntry.getValue();
            final String sample = sampleEntry.getKey();
            if (originalList == null || originalList.isEmpty()) {
                continue;
            }
            final List<GATKRead> newList = originalList.stream()
                    .filter(read -> read.overlaps(loc))
                    .collect(Collectors.toCollection(() -> new ArrayList<>(originalList.size())));

            if (!newList.isEmpty()) {
                overlappingFilteredReads.put(sample,newList);
            }
        }
        return overlappingFilteredReads;
    }
}