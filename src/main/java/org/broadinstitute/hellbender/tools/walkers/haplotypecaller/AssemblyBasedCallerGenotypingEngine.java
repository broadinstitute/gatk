package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * Created by davidben on 12/21/16.
 */
public abstract class AssemblyBasedCallerGenotypingEngine extends GenotypingEngine<AssemblyBasedCallerArgumentCollection> {

    protected static final int ALLELE_EXTENSION = 2;
    protected final boolean doPhysicalPhasing;
    private static final String phase01 = "0|1";
    private static final String phase10 = "1|0";

    /**
     * {@inheritDoc}
     * @param configuration {@inheritDoc}
     * @param samples {@inheritDoc}
     * @param doPhysicalPhasing whether to try physical phasing.
     */
    public AssemblyBasedCallerGenotypingEngine(final AssemblyBasedCallerArgumentCollection configuration, final SampleList samples,
                                           final AFCalculatorProvider afCalculatorProvider, final boolean doPhysicalPhasing) {
        super(configuration, samples, afCalculatorProvider, false);
        this.doPhysicalPhasing= doPhysicalPhasing;
    }

    @Override
    protected boolean forceSiteEmission() {
        return configuration.outputMode == OutputMode.EMIT_ALL_SITES || configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES;
    }

    /**
     * Carries the result of a call to #assignGenotypeLikelihoods
     */
    public static class CalledHaplotypes {
        private final List<VariantContext> calls;
        private final Set<Haplotype> calledHaplotypes;

        public CalledHaplotypes(final List<VariantContext> calls, final Set<Haplotype> calledHaplotypes) {
            this.calls = Utils.nonNull(calls, "calls cannot be null");
            this.calledHaplotypes = Utils.nonNull(calledHaplotypes, "calledHaplotypes cannot be null");
            Utils.validateArg(calls.isEmpty() == calledHaplotypes.isEmpty(), "Calls and calledHaplotypes should both be empty or both not but got calls=" + calls + " calledHaplotypes=" + calledHaplotypes);
        }

        /**
         * Get the list of calls made at this location
         * @return a non-null (but potentially empty) list of calls
         */
        public List<VariantContext> getCalls() {
            return calls;
        }

        /**
         * Get the set of haplotypes that we actually called (i.e., underlying one of the VCs in getCalls().
         * @return a non-null set of haplotypes
         */
        public Set<Haplotype> getCalledHaplotypes() {
            return calledHaplotypes;
        }
    }

    /**
     * Go through the haplotypes we assembled, and decompose them into their constituent variant contexts
     *
     * @param haplotypes the list of haplotypes we're working with
     * @param ref the reference bases (over the same interval as the haplotypes)
     * @param refLoc the span of the reference bases
     * @param activeAllelesToGenotype alleles we want to ensure are scheduled for genotyping (GGA mode)
     * @param maxMnpDistance Phased substitutions separated by this distance or less are merged into MNPs.  More than
     *                       two substitutions occurring in the same alignment block (ie the same M/X/EQ CIGAR element)
     *                       are merged until a substitution is separated from the previous one by a greater distance.
     *                       That is, if maxMnpDistance = 1, substitutions at positions 10,11,12,14,15,17 are partitioned into a MNP
     *                       at 10-12, a MNP at 14-15, and a SNP at 17.
     * @return never {@code null} but perhaps an empty list if there is no variants to report.
     */
    protected TreeSet<Integer> decomposeHaplotypesIntoVariantContexts(final List<Haplotype> haplotypes,
                                                                      final byte[] ref,
                                                                      final SimpleInterval refLoc,
                                                                      final List<VariantContext> activeAllelesToGenotype,
                                                                      final int maxMnpDistance) {
        final boolean inGGAMode = ! activeAllelesToGenotype.isEmpty();

        // Using the cigar from each called haplotype figure out what events need to be written out in a VCF file
        // IMPORTANT NOTE: This needs to be done even in GGA mode, as this method call has the side effect of setting the
        // event maps in the Haplotype objects!
        final TreeSet<Integer> startPosKeySet = EventMap.buildEventMapsForHaplotypes(haplotypes, ref, refLoc, configuration.debug, maxMnpDistance);

        if ( inGGAMode ) {
            startPosKeySet.clear();
            for( final VariantContext compVC : activeAllelesToGenotype ) {
                startPosKeySet.add(compVC.getStart());
            }
        }

        return startPosKeySet;
    }


    protected static class LocationAndAlleles {
        private final int loc;
        private final List<Allele> alleles;

        public LocationAndAlleles(final int loc, final List<Allele> alleles) {
            this.loc = loc;
            this.alleles = alleles;
        }

        public int getLoc() {
            return loc;
        }

        public List<Allele> getAlleles() {
            return alleles;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final LocationAndAlleles that = (LocationAndAlleles) o;

            if (loc != that.loc) return false;
            return alleles != null ? alleles.equals(that.alleles) : that.alleles == null;
        }

        @Override
        public int hashCode() {
            int result = loc;
            result = 31 * result + (alleles != null ? alleles.hashCode() : 0);
            return result;
        }

    }

    protected static List<VariantContext> getVCsAtThisLocation(final List<Haplotype> haplotypes,
                                                               final int loc,
                                                               final List<VariantContext> activeAllelesToGenotype) {
        return getVCsAtThisLocation(haplotypes, loc, activeAllelesToGenotype, false);
    }

    protected static List<VariantContext> getVCsAtThisLocation(final List<Haplotype> haplotypes,
                                                               final int loc,
                                                               final List<VariantContext> activeAllelesToGenotype,
                                                               final boolean includeSpanningEvents) {
        // the overlapping events to merge into a common reference view

                if (loc == 10008952) {
            int foo = loc;
        }

        final List<VariantContext> results = new ArrayList<>();
        final Set<LocationAndAlleles> uniqueLocationsAndAlleles = new HashSet<>();

        if( activeAllelesToGenotype.isEmpty() ) {
            haplotypes.stream()
                    .flatMap(h -> Utils.stream(h.getEventMap().getSpanningEvents(loc)))
                    .filter(Objects::nonNull)
                    .filter(v -> (includeSpanningEvents || v.getStart() == loc))
                    .forEach(v -> {
                        final LocationAndAlleles locationAndAlleles = new LocationAndAlleles(v.getStart(), v.getAlleles());
                        if (! uniqueLocationsAndAlleles.contains(locationAndAlleles)) {
                            uniqueLocationsAndAlleles.add(locationAndAlleles);
                            results.add(v);
                        }
                    });
        } else { // we are in GGA mode!
            int compCount = 0;
            for( final VariantContext compVC : activeAllelesToGenotype ) {
                if( compVC.getStart() <= loc && compVC.getEnd() >= loc) {
                    if (! (includeSpanningEvents || compVC.getStart() == loc)) {
                        continue;
                    }
                    int alleleCount = 0;
                    for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                        final List<Allele> alleleSet = Arrays.asList(compVC.getReference(), compAltAllele);

                        //TODO: this source name seems arbitrary and probably just has to be unique
                        //TODO: how about replace it by vcSourceName = String.parseInt(nameCounter++)?
                        final String vcSourceName = "Comp" + compCount + "Allele" + alleleCount;
                        // check if this event is already in the list of events due to a repeat in the input alleles track
                        final VariantContext candidateEventToAdd = new VariantContextBuilder(compVC).alleles(alleleSet).source(vcSourceName).make();

                        final LocationAndAlleles locationAndAlleles = new LocationAndAlleles(candidateEventToAdd.getStart(), candidateEventToAdd.getAlleles());
                        if (! uniqueLocationsAndAlleles.contains(locationAndAlleles)) {
                            uniqueLocationsAndAlleles.add(locationAndAlleles);
                            results.add(candidateEventToAdd);
                        }

                        alleleCount++;
                    }
                }
                compCount++;
            }
        }

        return results;
    }

    protected static Map<Allele, List<Haplotype>> createAlleleMapper(final List<VariantContext> eventsAtThisLoc,
                                                                     final VariantContext mergedVC,
                                                                     final int loc,
                                                                     final List<Haplotype> haplotypes) {
        return createAlleleMapper(eventsAtThisLoc, mergedVC, loc, haplotypes, null);
    }

    protected static Map<Allele, List<Haplotype>> createAlleleMapper(final List<VariantContext> eventsAtThisLoc,
                                                                     final VariantContext mergedVC,
                                                                     final int loc,
                                                                     final List<Haplotype> haplotypes,
                                                                     final List<VariantContext> activeAllelesToGenotype) {
        Utils.validateArg(haplotypes.size() > eventsAtThisLoc.size(), "expected haplotypes.size() >= eventsAtThisLoc.size() + 1");

        final Map<Allele, List<Haplotype>> result = new LinkedHashMap<>();

        final Allele ref = mergedVC.getReference();
        result.put(ref, new ArrayList<>());

        //Note: we can't use the alleles implied by eventsAtThisLoc because they are not yet merged to a common reference
        //For example, a homopolymer AAAAA reference with a single and double deletion would yield (i) AA* A and (ii) AAA* A
        //in eventsAtThisLoc, when in mergedVC it would yield AAA* AA A
        mergedVC.getAlternateAlleles().stream().filter(a -> !a.isSymbolic()).forEach(a -> result.put(a, new ArrayList<>()));

        if (loc == 10037037) {
            int foo = loc;
        }

        for (final Haplotype h : haplotypes) {

            final Iterator<VariantContext> spanningEvents = h.getEventMap().getSpanningEvents(loc);

            final boolean hasSpanningEvents = spanningEvents.hasNext();

            while (spanningEvents.hasNext()) {
                final VariantContext spanningEvent = spanningEvents.next();
                if (spanningEvent.getStart() == loc) {
                    //TODO: this only works if eventsAtThisLoc's and mergedVC's alleles are ordered identically.
                    //TODO: this is in fact the case (see AssemblyBasedCallerUtils::makeMergedVariantContext), but
                    //TODO: that's not good enough

                    //TODO: due to the ordering mentioned above the SPAN_DEL containing variants should come first in eventsAtThisLoc
                    int alleleIndex = 0;
                    for (int n = 0; n < eventsAtThisLoc.size(); n++) {
                        if (eventsAtThisLoc.get(n).getAlternateAllele(0).equals(Allele.SPAN_DEL)) {
                            alleleIndex = 1;
                            continue;
                        }
                        if (spanningEvent.hasSameAllelesAs(eventsAtThisLoc.get(n))) {
                            result.get(mergedVC.getAlternateAllele(alleleIndex)).add(h);
                            break;
                        }
                        alleleIndex++;
                    }

                } else {
                    if (activeAllelesToGenotype != null && activeAllelesToGenotype.size() > 0) {
                        // in HC GGA mode we need to check to make sure that spanning deletion
                        // events actually match one of the alleles we were given to genotype
                        for (VariantContext givenVC : activeAllelesToGenotype) {
                            if (givenVC.getStart() == spanningEvent.getStart()) {
                                for (Allele a : spanningEvent.getAlternateAlleles()) {
                                    if (givenVC.hasAlternateAllele(a)) {
                                        if (! result.containsKey(Allele.SPAN_DEL)) {
                                            result.put(Allele.SPAN_DEL, new ArrayList<>());
                                        }
                                        result.get(Allele.SPAN_DEL).add(h);
                                        break;
                                    }
                                }
                            }

                        }
                    } else {
                        if (! result.containsKey(Allele.SPAN_DEL)) {
                            result.put(Allele.SPAN_DEL, new ArrayList<>());
                        }
                        result.get(Allele.SPAN_DEL).add(h);
                    }
                    break;
                }
            }

            if (! hasSpanningEvents) {    //no events --> this haplotype supports the reference at this locus
                result.get(ref).add(h);
            }
        }
        return result;
    }

    // Builds the read-likelihoods collection to use for annotation considering user arguments and the collection
    // used for genotyping.
    protected ReadLikelihoods<Allele> prepareReadAlleleLikelihoodsForAnnotation(
            final ReadLikelihoods<Haplotype> readHaplotypeLikelihoods,
            final Map<String, List<GATKRead>> perSampleFilteredReadList,
            final boolean emitReferenceConfidence,
            final Map<Allele, List<Haplotype>> alleleMapper,
            final ReadLikelihoods<Allele> readAlleleLikelihoodsForGenotyping,
            final VariantContext call) {

        final ReadLikelihoods<Allele> readAlleleLikelihoodsForAnnotations;
        final SimpleInterval loc = new SimpleInterval(call);

        // We can reuse for annotation the likelihood for genotyping as long as there is no contamination filtering
        // or the user want to use the contamination filtered set for annotations.
        // Otherwise (else part) we need to do it again.
        if (configuration.useFilteredReadMapForAnnotations || !configuration.isSampleContaminationPresent()) {
            readAlleleLikelihoodsForAnnotations = readAlleleLikelihoodsForGenotyping;
            readAlleleLikelihoodsForAnnotations.filterToOnlyOverlappingReads(loc);
        } else {
            readAlleleLikelihoodsForAnnotations = readHaplotypeLikelihoods.marginalize(alleleMapper, loc);
            if (emitReferenceConfidence) {
                readAlleleLikelihoodsForAnnotations.addNonReferenceAllele(Allele.NON_REF_ALLELE);
            }
        }

        if (call.getAlleles().size() != readAlleleLikelihoodsForAnnotations.numberOfAlleles()) {
            readAlleleLikelihoodsForAnnotations.updateNonRefAlleleLikelihoods(new IndexedAlleleList<>(call.getAlleles()));
        }

        // Skim the filtered map based on the location so that we do not add filtered read that are going to be removed
        // right after a few lines of code below.
        final Map<String, List<GATKRead>> overlappingFilteredReads = overlappingFilteredReads(perSampleFilteredReadList, loc);

        readAlleleLikelihoodsForAnnotations.addReads(overlappingFilteredReads,0);

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

    /**
     * Tries to phase the individual alleles based on pairwise comparisons to the other alleles based on all called haplotypes
     *
     * @param calls             the list of called alleles
     * @param calledHaplotypes  the set of haplotypes used for calling
     * @return a non-null list which represents the possibly phased version of the calls
     */
    protected static List<VariantContext> phaseCalls(final List<VariantContext> calls, final Set<Haplotype> calledHaplotypes) {

        // construct a mapping from alternate allele to the set of haplotypes that contain that allele
        final Map<VariantContext, Set<Haplotype>> haplotypeMap = constructHaplotypeMapping(calls, calledHaplotypes);

        // construct a mapping from call to phase set ID
        final Map<VariantContext, Pair<Integer, String>> phaseSetMapping = new HashMap<>();
        final int uniqueCounterEndValue = constructPhaseSetMapping(calls, haplotypeMap, calledHaplotypes.size() - 1, phaseSetMapping);

        // we want to establish (potential) *groups* of phased variants, so we need to be smart when looking at pairwise phasing partners
        return constructPhaseGroups(calls, phaseSetMapping, uniqueCounterEndValue);
    }

    /**
     * Construct the mapping from alternate allele to the set of haplotypes that contain that allele
     *
     * @param originalCalls    the original unphased calls
     * @param calledHaplotypes  the set of haplotypes used for calling
     * @return non-null Map
     */
    protected static Map<VariantContext, Set<Haplotype>> constructHaplotypeMapping(final List<VariantContext> originalCalls,
                                                                                   final Set<Haplotype> calledHaplotypes) {
        final Map<VariantContext, Set<Haplotype>> haplotypeMap = new HashMap<>(originalCalls.size());
        for ( final VariantContext call : originalCalls ) {
            // don't try to phase if there is not exactly 1 alternate allele
            if ( ! isBiallelic(call) ) {
                haplotypeMap.put(call, Collections.<Haplotype>emptySet());
                continue;
            }

            // keep track of the haplotypes that contain this particular alternate allele
            final Allele alt = call.getAlternateAllele(0);
            final Predicate<VariantContext> hasThisAlt = vc -> (vc.getStart() == call.getStart() && vc.getAlternateAlleles().contains(alt)) ||
                    (Allele.SPAN_DEL.equals(alt) && vc.getStart() < call.getStart() && vc.getEnd() >= call.getStart());
            final Set<Haplotype> hapsWithAllele = calledHaplotypes.stream()
                    .filter(h -> h.getEventMap().getVariantContexts().stream().anyMatch(hasThisAlt))
                    .collect(Collectors.toCollection(HashSet<Haplotype>::new));

            haplotypeMap.put(call, hapsWithAllele);
        }

        return haplotypeMap;
    }


    /**
     * Construct the mapping from call (variant context) to phase set ID
     *
     * @param originalCalls    the original unphased calls
     * @param haplotypeMap     mapping from alternate allele to the set of haplotypes that contain that allele
     * @param totalAvailableHaplotypes the total number of possible haplotypes used in calling
     * @param phaseSetMapping  the map to populate in this method;
     *                         note that it is okay for this method NOT to populate the phaseSetMapping at all (e.g. in an impossible-to-phase situation)
     * @return the next incremental unique index
     */
    protected static int constructPhaseSetMapping(final List<VariantContext> originalCalls,
                                                  final Map<VariantContext, Set<Haplotype>> haplotypeMap,
                                                  final int totalAvailableHaplotypes,
                                                  final Map<VariantContext, Pair<Integer, String>> phaseSetMapping) {

        final int numCalls = originalCalls.size();
        int uniqueCounter = 0;

        // use the haplotype mapping to connect variants that are always/never present on the same haplotypes
        for ( int i = 0; i < numCalls - 1; i++ ) {
            final VariantContext call = originalCalls.get(i);
            final Set<Haplotype> haplotypesWithCall = haplotypeMap.get(call);
            if ( haplotypesWithCall.isEmpty() ) {
                continue;
            }

            final boolean callIsOnAllHaps = haplotypesWithCall.size() == totalAvailableHaplotypes;

            for ( int j = i+1; j < numCalls; j++ ) {
                final VariantContext comp = originalCalls.get(j);
                final Set<Haplotype> haplotypesWithComp = haplotypeMap.get(comp);
                if ( haplotypesWithComp.isEmpty() ) {
                    continue;
                }

                // if the variants are together on all haplotypes, record that fact.
                // another possibility is that one of the variants is on all possible haplotypes (i.e. it is homozygous).
                final boolean compIsOnAllHaps = haplotypesWithComp.size() == totalAvailableHaplotypes;
                if ( (haplotypesWithCall.size() == haplotypesWithComp.size() && haplotypesWithCall.containsAll(haplotypesWithComp)) || callIsOnAllHaps || compIsOnAllHaps ) {

                    // create a new group if these are the first entries
                    if ( ! phaseSetMapping.containsKey(call) ) {
                        // note that if the comp is already in the map then that is very bad because it means that there is
                        // another variant that is in phase with the comp but not with the call.  Since that's an un-phasable
                        // situation, we should abort if we encounter it.
                        if ( phaseSetMapping.containsKey(comp) ) {
                            phaseSetMapping.clear();
                            return 0;
                        }

                        // An important note: even for homozygous variants we are setting the phase as "0|1" here.
                        // We do this because we cannot possibly know for sure at this time that the genotype for this
                        // sample will actually be homozygous downstream: there are steps in the pipeline that are liable
                        // to change the genotypes.  Because we can't make those assumptions here, we have decided to output
                        // the phase as if the call is heterozygous and then "fix" it downstream as needed.
                        phaseSetMapping.put(call, Pair.of(uniqueCounter, phase01));
                        phaseSetMapping.put(comp, Pair.of(uniqueCounter, phase01));
                        uniqueCounter++;
                    }
                    // otherwise it's part of an existing group so use that group's unique ID
                    else if ( ! phaseSetMapping.containsKey(comp) ) {
                        final Pair<Integer, String> callPhase = phaseSetMapping.get(call);
                        phaseSetMapping.put(comp, Pair.of(callPhase.getLeft(), callPhase.getRight()));
                    }
                }
                // if the variants are apart on *all* haplotypes, record that fact
                else if ( haplotypesWithCall.size() + haplotypesWithComp.size() == totalAvailableHaplotypes ) {

                    final Set<Haplotype> intersection = new HashSet<>();
                    intersection.addAll(haplotypesWithCall);
                    intersection.retainAll(haplotypesWithComp);
                    if ( intersection.isEmpty() ) {
                        // create a new group if these are the first entries
                        if ( ! phaseSetMapping.containsKey(call) ) {
                            // note that if the comp is already in the map then that is very bad because it means that there is
                            // another variant that is in phase with the comp but not with the call.  Since that's an un-phasable
                            // situation, we should abort if we encounter it.
                            if ( phaseSetMapping.containsKey(comp) ) {
                                phaseSetMapping.clear();
                                return 0;
                            }

                            phaseSetMapping.put(call, Pair.of(uniqueCounter, phase01));
                            phaseSetMapping.put(comp, Pair.of(uniqueCounter, phase10));
                            uniqueCounter++;
                        }
                        // otherwise it's part of an existing group so use that group's unique ID
                        else if ( ! phaseSetMapping.containsKey(comp) ){
                            final Pair<Integer, String> callPhase = phaseSetMapping.get(call);
                            phaseSetMapping.put(comp, Pair.of(callPhase.getLeft(), callPhase.getRight().equals(phase01) ? phase10 : phase01));
                        }
                    }
                }
            }
        }

        return uniqueCounter;
    }

    /**
     * Assemble the phase groups together and update the original calls accordingly
     *
     * @param originalCalls    the original unphased calls
     * @param phaseSetMapping  mapping from call (variant context) to phase group ID
     * @param indexTo          last index (exclusive) of phase group IDs
     * @return a non-null list which represents the possibly phased version of the calls
     */
    protected static List<VariantContext> constructPhaseGroups(final List<VariantContext> originalCalls,
                                                               final Map<VariantContext, Pair<Integer, String>> phaseSetMapping,
                                                               final int indexTo) {
        final List<VariantContext> phasedCalls = new ArrayList<>(originalCalls);

        // if we managed to find any phased groups, update the VariantContexts
        for ( int count = 0; count < indexTo; count++ ) {
            // get all of the (indexes of the) calls that belong in this group (keeping them in the original order)
            final List<Integer> indexes = new ArrayList<>();
            for ( int index = 0; index < originalCalls.size(); index++ ) {
                final VariantContext call = originalCalls.get(index);
                if ( phaseSetMapping.containsKey(call) && phaseSetMapping.get(call).getLeft() == count ) {
                    indexes.add(index);
                }
            }
            if ( indexes.size() < 2 ) {
                throw new IllegalStateException("Somehow we have a group of phased variants that has fewer than 2 members");
            }

            // create a unique ID based on the leftmost one
            final String uniqueID = createUniqueID(originalCalls.get(indexes.get(0)));

            // update the VCs
            for ( final int index : indexes ) {
                final VariantContext originalCall = originalCalls.get(index);
                final VariantContext phasedCall = phaseVC(originalCall, uniqueID, phaseSetMapping.get(originalCall).getRight());
                phasedCalls.set(index, phasedCall);
            }
        }

        return phasedCalls;
    }

    /**
     * Is this variant bi-allelic?  This implementation is very much specific to this class so shouldn't be pulled out into a generalized place.
     *
     * @param vc the variant context
     * @return true if this variant context is bi-allelic, ignoring the NON-REF symbolic allele, false otherwise
     */
    private static boolean isBiallelic(final VariantContext vc) {
        return vc.isBiallelic() || (vc.getNAlleles() == 3 && vc.getAlternateAlleles().contains(Allele.NON_REF_ALLELE));
    }

    /**
     * Create a unique identifier given the variant context
     *
     * @param vc   the variant context
     * @return non-null String
     */
    private static String createUniqueID(final VariantContext vc) {
        return String.format("%d_%s_%s", vc.getStart(), vc.getReference().getDisplayString(), vc.getAlternateAllele(0).getDisplayString());
    }

    /**
     * Add physical phase information to the provided variant context
     *
     * @param vc   the variant context
     * @param ID   the ID to use
     * @param phaseGT the phase GT string to use
     * @return phased non-null variant context
     */
    private static VariantContext phaseVC(final VariantContext vc, final String ID, final String phaseGT) {
        final List<Genotype> phasedGenotypes = new ArrayList<>();
        for ( final Genotype g : vc.getGenotypes() ) {
            phasedGenotypes.add(new GenotypeBuilder(g).attribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, ID).attribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, phaseGT).make());
        }
        return new VariantContextBuilder(vc).genotypes(phasedGenotypes).make();
    }

}
