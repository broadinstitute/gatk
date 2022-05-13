package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandOddsRatio;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.InverseAllele;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.io.ComponentNameProvider;
import org.jgrapht.io.DOTExporter;
import org.jgrapht.io.IntegerComponentNameProvider;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Filtering haplotypes that contribute weak alleles to the genotyping.
 *
 * @author Ilya Soifer &lt;ilya.soifer@ultimagen.com&gt;
 * @author Yossi Farjoun &lt;farjoun@broadinstitute.org&gt;
 *
 */

public abstract class AlleleFiltering {

    final protected static Logger logger = LogManager.getLogger(AlleleFiltering.class);
    final protected AssemblyBasedCallerArgumentCollection assemblyArgs;
    final private OutputStreamWriter assemblyDebugOutStream;
    AlleleFiltering(final AssemblyBasedCallerArgumentCollection assemblyArgs, final OutputStreamWriter assemblyDebugOutStream){
        this.assemblyArgs = assemblyArgs;
        this.assemblyDebugOutStream = assemblyDebugOutStream;
    }

    /**
     * Finds alleles that are likely not contributing much to explaining the data and remove the haplotypes
     * that contribute them.
     *
     * The alleles from the active region are divided into clusters of alleles that likely "compete" with each
     * other, where compete means that they are the same allele up to a sequencing error although they might
     * be assigned to a different genomic location. In each cluster we iteratively calculate the quality of
     * each allele relative to other alleles in the cluster and remove the allele with the lowest quality.
     * We then also select in each cluster alleles with high SOR and remove them.
     *
     * Every haplotype that contributes a filtered allele is filtered out.
     *
     * @param readLikelihoods unfiltered read x haplotype likelihood matrix
     * @param activeWindowStart location of the active windows (assemblyResult.getPaddedReferenceLoc().getStart()
     * @param suspiciousLocations set of suspicious locations for further marking in genotyping
     * @return Subsetted read x haplotype where only the haplotypes that do not contribute filtered alleles show. Also
     * locations of filtered alleles on the genome added to `suspiciousLocations` list
     */

    public AlleleLikelihoods<GATKRead, Haplotype> filterAlleles(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                final int activeWindowStart, Set<Integer> suspiciousLocations){

        logger.debug("SHA:: filter alleles - start");
        final AlleleLikelihoods<GATKRead, Haplotype> subsettedReadLikelihoodsFinal = subsetHaplotypesByAlleles(readLikelihoods, assemblyArgs, activeWindowStart, suspiciousLocations);
        logger.debug("SHA:: filter alleles - end");

        readLikelihoods.setFilteredHaplotypeCount(readLikelihoods.numberOfAlleles() - subsettedReadLikelihoodsFinal.numberOfAlleles());

        if (assemblyDebugOutStream != null) {
            try {
                assemblyDebugOutStream.write("\nThere were " + subsettedReadLikelihoodsFinal.alleles().size() + " haplotypes found after subsetting by alleles. Here they are:\n");
                subsettedReadLikelihoodsFinal.alleles().forEach(h -> {
                    try {
                        assemblyDebugOutStream.write(h.toString());
                        assemblyDebugOutStream.append("\n");
                    } catch (IOException e) {
                        throw new UserException("Error writing to debug output stream", e);
                    }
                });
            } catch (IOException e) {
                throw new UserException("Error writing to debug output stream", e);
            }
        }

        return subsettedReadLikelihoodsFinal;
    }

    /**
     * Returns all alleles from haplotype
     * @param haplotype Input
     * @return set of AlleleAndContext
     */
    static private Set<AlleleAndContext> getAlleles(final Haplotype haplotype){
        final Collection<VariantContext> vcs = haplotype.getEventMap().getVariantContexts();

        return vcs.stream().flatMap(
                vc -> vc.getAlleles().stream().map(
                        al ->  new AlleleAndContext(vc.getContig(), vc.getStart(), al, vc.getReference()))
        ).collect(Collectors.toSet());
    }

    /**
     * Main function that filters haplotypes that contribute weak alleles
     * @param readLikelihoods read x haplotype matrix
     * @param assemblyArgs HaplotypeCaller/Mutect2 parameters
     * @param activeWindowStart Genomic location of the start
     * @param suspiciousLocations set of positions where the alleles are being filtered (modified)
     * @return read x haplotype matrix where the filtered haplotypes are removed
     * @throws IOException if output file can't be written
     */
    private AlleleLikelihoods<GATKRead, Haplotype> subsetHaplotypesByAlleles(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                         final AssemblyBasedCallerArgumentCollection assemblyArgs,
                                                                         final int activeWindowStart, Set<Integer> suspiciousLocations) {
        // 1. Collect all alleles in the active region
        final Set<Haplotype> disabledHaplotypes = new HashSet<>();
        final Map<Haplotype, Collection<AlleleAndContext>> haplotypeAlleleMap  = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);
        readLikelihoods.alleles().forEach(h -> getAlleles(h).stream().filter(al -> !al.isReference()).forEach(jh -> haplotypeAlleleMap.get(h).add(jh)));

        // 2. Split them into sets to genotype together. The goal is to cluster true allele with all its variants from
        // biased seq. error.
        // The alleles split into clusters of alleles that potentially interact (compete with each other for reads)
        // First we generate a graph with edge for each pair of alleles that do not occur in the same haplotype
        // Then we only keep the edges where the alleles are close or up to hmer indel from each other
        // the connected components of the graph are genotyped together
        final OccurrenceMatrix<Haplotype, AlleleAndContext> occm = new OccurrenceMatrix<>(haplotypeAlleleMap);
        List<Pair<AlleleAndContext, AlleleAndContext>> nonCoOcurringAlleles = occm.nonCoOcurringColumns();
        final List<Pair<AlleleAndContext, AlleleAndContext>> closeNonCoOccurringAlleles = filterByDistance(nonCoOcurringAlleles, 0, 3);
        nonCoOcurringAlleles = filterSameUpToHmerPairs(filterByDistance(nonCoOcurringAlleles,0,20),
                findReferenceHaplotype(readLikelihoods.alleles()), activeWindowStart);
        nonCoOcurringAlleles.addAll(closeNonCoOccurringAlleles);
        final List<Set<AlleleAndContext>> independentAlleles = occm.getIndependentSets(nonCoOcurringAlleles);

        // 3. For each cluster - remove weak alleles
        for (final Set<AlleleAndContext> alleleSet : independentAlleles) {

            // debugging - write the interaction map of the location (we will keep this function from the unused approach
            // where we only attempted to filter alleles that strongly affect an another allele's quality. This approach
            // failed to deliver a significant improvement and thus is not used.
            // interaction map is the graph of how much quality of each allele is improved when another allele is removed
            if (assemblyArgs.writeFilteringGraphs) {
                if (alleleSet.size() > 1 ) {
                    final List<AlleleAndContext> alleleSetAsList = new ArrayList<>(alleleSet);
                    final Map<AlleleAndContext, Integer> initialRPLsMap = new HashMap<>();
                    final DefaultDirectedWeightedGraph<AlleleAndContext, DefaultWeightedEdge> intm =
                            interactionMatrixToGraph(getInteractionMatrix(alleleSetAsList, haplotypeAlleleMap,
                                    readLikelihoods, initialRPLsMap), initialRPLsMap);
                    printInteractionGraph(intm, initialRPLsMap, alleleSet);
                }
            }

            boolean removedAlleles = true;
            final Set<Haplotype> activeHaplotypes = new HashSet<>(readLikelihoods.alleles());

            while (removedAlleles) {
                removedAlleles = false;
                // b. Marginalize: calculate quality of each allele relative to all other alleles
                logger.debug("GAL::start of iteration");
                final List<AlleleAndContext> activeAlleles = activeHaplotypes.stream()
                        .flatMap(h -> getAlleles(h).stream().filter(alleleSet::contains))
                        .distinct()
                        .collect(Collectors.toList());;

                final Map<AlleleAndContext, List<Haplotype>> alleleHaplotypeMap = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);
                readLikelihoods.alleles().stream().filter(activeHaplotypes::contains)
                        .forEach(h ->
                                    getAlleles(h).stream()
                                            .filter(alleleSet::contains)
                                            .filter(al -> !al.isReference())
                                            .forEach(jh -> alleleHaplotypeMap.get(jh).add(h))

                        );

                logger.debug("AHM::printout start");
                for (final AlleleAndContext al : alleleHaplotypeMap.keySet()) {
                    logger.debug("AHM::allele block ---> ");
                    for (final Allele h : alleleHaplotypeMap.get(al)) {
                        logger.debug(() -> String.format("AHM:: (%d) %s/%s: %s", al.getLoc(), al.getAllele().getBaseString(), al.getRefAllele().getBaseString(), h.getBaseString()));
                    }
                    logger.debug("AHM::allele block ---< ");

                }
                logger.debug("AHM::printout end");



                final List<AlleleLikelihoods<GATKRead, Allele>> alleleLikelihoods =
                        activeAlleles.stream().map(al -> getAlleleLikelihoodMatrix(readLikelihoods, al,
                                haplotypeAlleleMap, activeHaplotypes)).collect(Collectors.toList());
                //    c. Calculate SOR and RPL
                // Note that the QUAL is calculated as a PL, that is -10*log likelihood. This means that high PL is low quality allele
                final List<Integer> collectedRPLs = IntStream.range(0, activeAlleles.size()).mapToObj(i -> getAlleleLikelihoodVsInverse(alleleLikelihoods.get(i), activeAlleles.get(i))).collect(Collectors.toList());
                final List<Double> collectedSORs = IntStream.range(0, activeAlleles.size()).mapToObj(i -> getAlleleSOR(alleleLikelihoods.get(i), activeAlleles.get(i))).collect(Collectors.toList());

                //    d. Generate variants that are below SOR threshold and below RPL threshold
                final List<AlleleAndContext> filteringCandidates = identifyBadAlleles(collectedRPLs,
                                                                                collectedSORs,
                                                                                activeAlleles,
                                                                                assemblyArgs.prefilterQualThreshold,
                                                                                assemblyArgs.prefilterSorThreshold);


                //very weak candidates are filtered out in any case, even if they are alone (they will be filtered anyway even in the GVCF mode)
                // the very weak quality is hardcoded
                final List<AlleleAndContext> filteringCandidatesStringent = identifyBadAlleles(collectedRPLs,
                        collectedSORs,
                        activeAlleles,
                        1,
                        Integer.MAX_VALUE);


                //for now we just mark all locations with close alleles, one of which is weak.
                //We write them in suspiciousLocations and they will be then annotated as SUSP_NOISY... in the VCF
                if ((filteringCandidates.size() > 0 ) && (alleleSet.size()>0)) {
                    activeAlleles.forEach(laa -> suspiciousLocations.add(laa.getLoc()));
                }

                //    e. For every variant - calculate what is the effect of its deletion and if higher than threshold - delete and continue

                // (This is a currently disabled code from the approach that would disable only the candidates that strongly
                // affect other alleles
                //AlleleAndContext candidateToDisable = identifyStrongInteractingAllele(filteringCandidates,
                //        hcargs.prefilterQualThreshold, activeAlleles, collectedRPLs, readLikelihoods, haplotypeAlleleMap, alleleHaplotypeMap); )

                // if weak candidate had been identified - add its haplotypes into blacklist, remove the allele from the
                // current cluster and genotype again.
                if ((filteringCandidates.size()>0 && activeAlleles.size()>1) ||
                (activeAlleles.size()==1 && filteringCandidatesStringent.size()>0) ||
                        (filteringCandidates.size()>0 && this.assemblyArgs.filterLoneAlleles)) {

                    if ((filteringCandidatesStringent.size()>0) && (filteringCandidates.size() == 0 )) {
                        throw new GATKException.ShouldNeverReachHereException("The thresholds for stringent allele " +
                                "filtering should always be higher than for the relaxed one");
                    }

                    final AlleleAndContext candidateToDisable = filteringCandidates.get(0);
                    logger.debug(() -> String.format("GAL:: Remove %s", candidateToDisable.toString()));
                    removedAlleles = true;
                    final List<Haplotype> haplotypesToRemove = alleleHaplotypeMap.get(candidateToDisable);
                    disabledHaplotypes.addAll(haplotypesToRemove);
                    activeHaplotypes.removeAll(haplotypesToRemove);
                }
                logger.debug("GAL::end of iteration");

            }
        }

        // finalizing: remove all disabled genotypes
        logger.debug("----- SHA list of removed haplotypes start ----");
        for (Haplotype hap : disabledHaplotypes) {
            logger.debug(() -> String.format("SHA :: Removed haplotype : %s ", hap.toString()));
        }
        logger.debug("----- SHA list of removed haplotypes end ----");

        final Set<Haplotype> eventualAlleles = new HashSet<>();
        readLikelihoods.alleles().stream().filter(al -> !disabledHaplotypes.contains(al)).forEach(eventualAlleles::add);
        logger.debug("----- SHA list of remaining haplotypes start ----");
        for (Haplotype hap : eventualAlleles) {
            logger.debug(() -> String.format("SHA :: Remaining haplotype : %s ", hap.toString()));
        }
        logger.debug("----- SHA list of remaining haplotypes end ----");



        final AlleleLikelihoods<GATKRead, Haplotype> currentReadLikelihoods = readLikelihoods.removeAllelesToSubset(eventualAlleles);
        logger.debug("----- SHA list of remaining alleles start ----");
        final Set<AlleleAndContext> locAllele = new HashSet<>();
        currentReadLikelihoods.alleles().forEach(h -> getAlleles(h).stream().filter(al -> !al.isReference()).forEach(locAllele::add));
        for (final AlleleAndContext al: locAllele) {
            logger.debug(() -> String.format("---- SHA :: %s ", al.toString()));
        }
        logger.debug("----- SHA list of remaining alleles end ----");

        return currentReadLikelihoods;
    }


    /**
     * Finds a list of alleles that are candidate for removal in the order of precedence (first - the best candidate to be removed)
     *
     * @param collectedRPLs list of each allele qualities (collected by {@link AlleleFiltering#getAlleleLikelihoodVsInverse}
     * @param collectedSORs list of each allele SORs (collected by {@link AlleleFiltering#getAlleleSOR(AlleleLikelihoods, Allele)}
     * @param alleles list of alleles in the same order as in collectedRPLs/collectedSORs
     * @param qualThreshold only variants with quality below qualThreshold will be considered
     * @param sorThreshold only variants with SOR above threshold will be considered
     * @return list of alleles that can be removed
     */
    private List<AlleleAndContext> identifyBadAlleles(final List<Integer> collectedRPLs, final List<Double> collectedSORs,
                                                      final List<AlleleAndContext> alleles,
                                                      final double qualThreshold,
                                                      final double sorThreshold) {

        //collected RPLs are the -10*QUAL of the alleles. high RPL means low quality.
        // SORs are regular: high SOR - strongly biased
        final int[] rplsIndices = getSortedIndexList(collectedRPLs);
        final int[] sorIndices = rplsIndices; // the variants that have high sor are ordered according to their quality


        //this list will contain all alleles that should be filtered in the order of priority
        final List<AlleleAndContext> result = new ArrayList<>();
        final double THRESHOLD = -1 * qualThreshold; // quality threshold is like in GATK (GL) and we collected PL, so QUAL 30 will appear as -30.
        final double SOR_THRESHOLD = sorThreshold;

        //note that high value is a poor quality allele, so the worst allele is the highest collectedRPL
        //we first collect all allleles with low quality: from the lowest
        for (int i = rplsIndices.length-1 ; (i >= 0) && (collectedRPLs.get(rplsIndices[i])>THRESHOLD) ; i--) {
            result.add(alleles.get(rplsIndices[i]));
        }
        int rplCount = result.size();
        //we then add alleles with high SOR. Note that amongh all allleles with the SOR higher than the SOR_THRESHOLD
        //we will first filter the one with the lowest QUAL.
        logger.debug(() -> String.format("SHA:: Have %d candidates with low QUAL", rplCount));
        for (int i = sorIndices.length-1 ; (i >= 0) && (collectedSORs.get(sorIndices[i])>SOR_THRESHOLD) ; i--) {
            if (!result.contains(alleles.get(sorIndices[i]))) {
                result.add(alleles.get(sorIndices[i]));
            }
        }
        logger.debug(() -> String.format("SHA:: Have %d candidates with high SOR", result.size() - rplCount));
        return result;
    }


    /**
     * Generates from read x haplotype matrix a read x allele matrix with two alleles: Allele and ~Allele where Allele
     * is supported by all haplotypes that contain this allele and ~Allele is supported by all other haplotypes.
     * Similar to {@link AlleleLikelihoods#marginalize(Map)}.
     *
     * @param readLikelihoods read x haplotype matrix
     * @param allele allele to consider
     * @param haplotypeAlleleMap map between alleles and haplotypes
     * @param enabledHaplotypes list of haplotypes that are currently inactive
     * @return read x allele matrix
     */
    private AlleleLikelihoods<GATKRead, Allele> getAlleleLikelihoodMatrix(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                                     final AlleleAndContext allele,
                                                                                     final Map<Haplotype, Collection<AlleleAndContext>> haplotypeAlleleMap,
                                                                                     final Set<Haplotype> enabledHaplotypes
                                                                                     ){
        final Map<Allele,List<Haplotype>> alleleHaplotypeMap = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);

        final Allele notAllele= InverseAllele.of(allele.getAllele(), true);
        readLikelihoods.alleles().stream().filter(enabledHaplotypes::contains)
                .filter(h->haplotypeAlleleMap.get(h).contains(allele))
                .forEach(alleleHaplotypeMap.get(allele)::add);
        readLikelihoods.alleles().stream().filter(enabledHaplotypes::contains)
                .filter(h -> !haplotypeAlleleMap.get(h).contains(allele))
                .forEach(alleleHaplotypeMap.get(notAllele)::add);

        final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods = readLikelihoods.marginalize(alleleHaplotypeMap);
        logger.debug(() -> String.format("GALM: %s %d %d", allele.toString(), alleleHaplotypeMap.get(allele).size(), alleleHaplotypeMap.get(notAllele).size()));
        return alleleLikelihoods;
    }

    //functions to get allele likelihoods and SOR. Differ between the mutect and the HC implementations
    abstract int getAlleleLikelihoodVsInverse(final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods, final Allele allele);

    private double getAlleleSOR(final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods, final Allele allele) {
        final Allele notAllele = InverseAllele.of(allele, true);
        final int [][] contingency_table = StrandOddsRatio.getContingencyTableWrtAll(alleleLikelihoods, notAllele, Collections.singletonList(allele), 1);
        final double sor = StrandOddsRatio.calculateSOR(contingency_table);
        logger.debug(() -> String.format("GAS:: %s: %f (%d %d %d %d)", allele.toString(), sor, contingency_table[0][0], contingency_table[0][1], contingency_table[1][0], contingency_table[1][1]));
        return sor;

    }

    //filters pairs of alleles by distance
    private List<Pair<AlleleAndContext, AlleleAndContext>> filterByDistance(
            final List<Pair<AlleleAndContext, AlleleAndContext>> allelePairs,
            final int minDist, final int maxDist) {
        logger.debug(() -> String.format("FBD: input %d pairs ", allelePairs.size()));
        final List<Pair<AlleleAndContext, AlleleAndContext>> result = new ArrayList<>(allelePairs);
        result.removeIf(v -> Math.abs(v.getLeft().getLoc() - v.getRight().getLoc())>maxDist);
        result.removeIf(v -> Math.abs(v.getLeft().getLoc() - v.getRight().getLoc())<minDist);
        logger.debug(() -> String.format("FBD: output %d pairs ", allelePairs.size()));

        return result;
    }

    //filters pairs of alleles that are not same up to hmer indel
    private List<Pair<AlleleAndContext, AlleleAndContext>> filterSameUpToHmerPairs(final List<Pair<AlleleAndContext,
            AlleleAndContext>> allelePairs, final Haplotype refHaplotype, final int activeWindowStart) {

        final List<Pair<AlleleAndContext, AlleleAndContext>> result = new ArrayList<>();
        for (final Pair<AlleleAndContext, AlleleAndContext> allelePair: allelePairs) {
            final Pair<Haplotype, Haplotype> modifiedHaplotypes = new ImmutablePair<>(
                    refHaplotype.insertAllele(
                            allelePair.getLeft().getRefAllele(),
                            allelePair.getLeft().getAllele(),
                            allelePair.getLeft().getLoc()),
                    refHaplotype.insertAllele(
                            allelePair.getRight().getRefAllele(),
                            allelePair.getRight().getAllele(),
                            allelePair.getRight().getLoc()));

            if ( BaseUtils.equalUpToHmerChange(modifiedHaplotypes.getLeft().getBases(), modifiedHaplotypes.getRight().getBases()) ) {
                result.add(allelePair);
            }

        }

        return result;
    }


    // find (the) reference haplotype within a list of haplotypes
    static Haplotype findReferenceHaplotype(final List<Haplotype> haplotypeList) {
        for (final Haplotype h: haplotypeList ) {
            if (h.isReference()) {
                return h;
            }
        }
        return null;
    }

    // if alleles are different in length, return their the length of their (potentially) common prefix, otherwise return 0
    private int getCommonPrefixLength(final Allele al1, final Allele al2){
        if (al1.length()!=al2.length()){
            return Math.min(al1.length(), al2.length());
        } else {
            return 0;
        }
    }

    // sort an integer list
    private int[] getSortedIndexList(final List<Integer> values) {
        return IntStream.range(0, values.size()).
                mapToObj(i -> new ImmutablePair<>(i, values.get(i)))
                .sorted(Comparator.comparingInt( v -> (int)v.getRight()))
                .mapToInt(v-> v.getLeft()).toArray();

    }

    // the functions below are currently unused but kept for potential future uses.
    // The goal of these functions is to look at how one allele affects the other and make decisions
    // only for the alleles that really affect others. The approach did not currently work that well
    @SuppressWarnings("unused")
    private AlleleAndContext identifyStrongInteractingAllele(final List<AlleleAndContext> candidateList,
                                                             final float prefilterThreshold,
                                                             final List<AlleleAndContext> allAlleles,
                                                             final List<Integer> rpls,
                                                             final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                             final Map<Haplotype, Collection<AlleleAndContext>> haplotypeAlleleMap,
                                                             final Map<AlleleAndContext, List<Haplotype>> alleleHaplotypeMap
                                                              ){


        logger.debug("ISIA :: start");
        final Map<AlleleAndContext, Integer> initialRPLsMap = new HashMap<>();
        IntStream.range(0, allAlleles.size()).forEach(i -> initialRPLsMap.put(allAlleles.get(i), rpls.get(i)));

        for (final AlleleAndContext cand: candidateList){
            logger.debug(() -> String.format("ISIA :: test %s", cand.toString()));
            if ( initialRPLsMap.get(cand) > (-1)*prefilterThreshold){
                logger.debug( String.format("ISIA:: selected %s due to low QUAL", cand));
                return cand;
            }

            if (allAlleles.size() <=1) {
                return null;
            }

            final Map<AlleleAndContext, Integer> interactionVector = getInteractionVector(cand,
                    haplotypeAlleleMap, alleleHaplotypeMap, readLikelihoods, initialRPLsMap);
            for (final AlleleAndContext allele: interactionVector.keySet()){
                logger.debug(() -> String.format(" --- %s: %d", allele.toString(), initialRPLsMap.get(allele) - interactionVector.get(allele)));
                if (initialRPLsMap.get(allele) - interactionVector.get(allele) > prefilterThreshold ){
                    logger.debug(String.format("ISIA:: selected %s", cand));
                    return cand;
                }
            }
        }
        logger.debug("ISIA :: end");

        return null;

    }


    // function to calculate interactions matrix between the alleles
    private Map<AlleleAndContext, Map<AlleleAndContext, Integer>> getInteractionMatrix(
            final List<AlleleAndContext> alleles,
            final Map<Haplotype, Collection<AlleleAndContext>> haplotypeAlleleMap,
            final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
            final Map<AlleleAndContext, Integer> initialRPLsMap) {

        final Map<AlleleAndContext, List<Haplotype>> alleleHaplotypeMap = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);
        final Set<Haplotype> haplotypes = new HashSet<>(readLikelihoods.alleles());
        readLikelihoods.alleles().stream().forEach(h -> getAlleles(h).stream().filter(al -> alleles.contains(al)).filter(al -> !al.isReference()).forEach(
                        jh -> alleleHaplotypeMap.get(jh).add(h))
                );

        final List<AlleleAndContext> allAlleles = new ArrayList<>(alleleHaplotypeMap.keySet());

        final List<AlleleLikelihoods<GATKRead, Allele>> initialAlleleLikelihoods =
                allAlleles.stream().map(c -> getAlleleLikelihoodMatrix(readLikelihoods, c, haplotypeAlleleMap, haplotypes)).collect(Collectors.toList());

        final List<Integer> initialRPLs = IntStream.range(0, allAlleles.size()).mapToObj(i -> getAlleleLikelihoodVsInverse(initialAlleleLikelihoods.get(i),
                allAlleles.get(i))).collect(Collectors.toList());

        for (int i = 0 ; i < allAlleles.size(); i++) {
            initialRPLsMap.put(allAlleles.get(i), initialRPLs.get(i));
        }

        final Map<AlleleAndContext, Map<AlleleAndContext, Integer>> result = new HashMap<>();
        for ( final AlleleAndContext alleleToDisable : allAlleles) {
            Map<AlleleAndContext, Integer> rplsWithoutAlleleMap = getInteractionVector(alleleToDisable, haplotypeAlleleMap, alleleHaplotypeMap, readLikelihoods, initialRPLsMap);
            result.put(alleleToDisable, rplsWithoutAlleleMap);
        }

        return result;
    }

    // function to create interaction of a single allele with other alleles
    private Map<AlleleAndContext, Integer> getInteractionVector(
            final AlleleAndContext alleleToDisable,
            final Map<Haplotype, Collection<AlleleAndContext>> haplotypeAlleleMap,
            final Map<AlleleAndContext, List<Haplotype>> alleleHaplotypeMap,
            final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
            final Map<AlleleAndContext, Integer> initialRPLsMap) {


        final Set<AlleleAndContext> allAlleles = initialRPLsMap.keySet();
        final List<AlleleAndContext> allelesWithoutDisabledAllele = allAlleles.stream().filter(al -> al!=alleleToDisable).collect(Collectors.toList());
        final Set<Haplotype> haplotypes = haplotypeAlleleMap.keySet();
        final Set<Haplotype> haplotypesWithoutDisabledAllele = haplotypes.stream().filter( h -> !alleleHaplotypeMap.get(alleleToDisable).contains(h)).collect(Collectors.toSet());

        final List<AlleleLikelihoods<GATKRead, Allele>> disabledAlleleLikelihood =
                allelesWithoutDisabledAllele.stream().map(c -> getAlleleLikelihoodMatrix(readLikelihoods, c, haplotypeAlleleMap, haplotypesWithoutDisabledAllele)).collect(Collectors.toList());

        final List<Integer> rplsWithoutAllele = IntStream.range(0, allelesWithoutDisabledAllele.size()).mapToObj(i -> getAlleleLikelihoodVsInverse(disabledAlleleLikelihood.get(i),
                allelesWithoutDisabledAllele.get(i))).collect(Collectors.toList());

        final Map<AlleleAndContext, Integer> rplsWithoutAlleleMap = new HashMap<>();
        IntStream.range(0, allelesWithoutDisabledAllele.size()).forEach( i -> rplsWithoutAlleleMap.put(allelesWithoutDisabledAllele.get(i), rplsWithoutAllele.get(i)));

        return rplsWithoutAlleleMap;
    }


    private DefaultDirectedWeightedGraph<AlleleAndContext, DefaultWeightedEdge> interactionMatrixToGraph(final Map<AlleleAndContext, Map<AlleleAndContext, Integer>> interactionMatrix,
                                                                                                 final Map<AlleleAndContext, Integer> initialRPL ){
        final DefaultDirectedWeightedGraph<AlleleAndContext, DefaultWeightedEdge> result = new DefaultDirectedWeightedGraph<>(DefaultWeightedEdge.class);
        initialRPL.keySet().stream().forEach(x -> result.addVertex(x));


        for ( final AlleleAndContext loc1 : interactionMatrix.keySet() ) {
            for ( final AlleleAndContext loc2 : interactionMatrix.get(loc1).keySet()){
                final int diff = interactionMatrix.get(loc1).get(loc2) - initialRPL.get(loc2);
                if (diff < 0){
                    final DefaultWeightedEdge edge = result.addEdge(loc1, loc2);
                    result.setEdgeWeight(edge, diff);
                }
            }
        }
        return result;
    }

    //debug function - prints dot file with edges between the alleles that affect each other
    void printInteractionGraph(final DefaultDirectedWeightedGraph<AlleleAndContext, DefaultWeightedEdge>  intm,
                               final Map<AlleleAndContext, Integer> rpls,
                               final Set<AlleleAndContext> alleleSet ){
        final IntegerComponentNameProvider<AlleleAndContext> p1 = new IntegerComponentNameProvider<>();
        final ComponentNameProvider<AlleleAndContext> p2 = (v -> v.toString() + " = " + rpls.get(v)) ;
        final ComponentNameProvider<DefaultWeightedEdge> p4 = (e -> String.valueOf(intm.getEdgeWeight(e)));

        final DOTExporter<AlleleAndContext, DefaultWeightedEdge> dotExporter = new DOTExporter<>(p1, p2, p4,
                null, null);
        final String contig = alleleSet.iterator().next().getContig();
        final int rangeStart = alleleSet.stream().mapToInt(al -> al.getLoc()).min().getAsInt();
        final int rangeEnd = alleleSet.stream().mapToInt(al -> al.getLoc()).max().getAsInt();
        try {
            final Writer outfile = new FileWriter(String.format("allele.interaction.%s.%d-%d.dot", contig, rangeStart, rangeEnd));
            dotExporter.exportGraph(intm, outfile);
        }
        catch (IOException e) {
            throw new RuntimeException("Unable to write a DOT file" + String.format("allele.interaction.%s.%d-%d.dot", contig, rangeStart, rangeEnd));
        }
    }
}

