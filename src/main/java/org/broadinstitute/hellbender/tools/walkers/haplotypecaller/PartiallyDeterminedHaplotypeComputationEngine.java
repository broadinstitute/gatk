package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.*;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SmallBitSet;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.jgrapht.Graph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

import java.nio.ByteBuffer;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Class that manages the complicated steps involved in generating artificial haplotypes for the PDHMM:
 *
 * The primary method to this class is {@link #generatePDHaplotypes} which will be called with an existing AssemblyResultSet object
 * as well as the list of alleles that were found from the pileupcaller code. The method attempts to replace the existing haplotypes
 * from the provided result set object with new haplotypes that were artificially generated from merging the assembled alleles
 * with those found by the pileup caller. Crucially, the constructed haplotypes will be PartiallyDeterminedHaplotype objects
 * which have markers indicating what bases are "undetermined" and should not be penalized in the HMM when computing likelihoods
 * The method might hit one of a few heuristic barriers and chose to fallback on only the assembled haplotypes if the
 * processing would become too complicated.
 */
public class PartiallyDeterminedHaplotypeComputationEngine {
    final static int MAX_PD_HAPS_TO_GENERATE = 256*2; //(2048 is Illumina's #) (without optimizing the hmm to some degree this is probably unattainable)
    final static int MAX_BRANCH_PD_HAPS = 128; //(128 is Illumina's #)
    final static int MAX_VAR_IN_EVENT_GROUP = 17; // (20 is Illumina's #)

    //To make this somewhat cleaner of a port from Illumina, we have two base spaces. R and U space. R space is vcf coordinate space,
    //U is a 0->N (N = region size) based space where Insertions are +0.5 and deletions are + 1 from their original position

    // We use this comparator for haplotype construction to make it slightly easier to build compound haplotypes (i.e. snp+insertion/deletion at the same anchor base)
    public static final Comparator<Event> HAPLOTYPE_SNP_FIRST_COMPARATOR = Comparator.comparingInt(Event::getStart)
            // Decide arbitrarily so as not to accidentally throw away overlapping variants
            .thenComparingInt(e -> e.refAllele().length())
            .thenComparingInt(e -> e.altAllele().length())
            .thenComparing(Event::altAllele);


    /**
     * The workhorse method for the PDHMM branching code. This method is responsible for the lion's share of the work in taking a list of
     * assembly haplotypes/alleles and the alleles found in the ColumnwiseDetection code and converting them into the artificial haplotypes
     * that are necessary for variant calling. This code might fail due to having too complex a site, if that happens it will return
     * the input sourceSet object without anything being changed.
     *
     * Overall pattern is as follows:
     *  1. Generate the combined list of variants from assembly and pileup-caller that pass the heuristic filers
     *  2. Generate lists of variants that overlap in groups
     *  3. Use SmithWaterman to realign pairs of 2/3 events that contain an indel, recording pairs that realign in a better
     *     set of sub-events or result in the reference.
     *  4. Merge event groups that have equivalent events.
     *  5. Iterate over every start position and attempt to build determined variants:
     *      5a. Select an allele at a position (ref or one of the alts)
     *      5b. Ask each event group for branches of mutually exclusive events (overlapping or from the equivalent event code) for the determined variants.
     *      5c. Merge all event groups into distinct branches (combinatorially)
     *      5d. For each branch, build a PDHaplotype with the determined allele and the allowed variants for that branch.
     *
     * @param sourceSet AssemblyResultSet to be modified with the new haplotypes
     * @param badPileupEvents Pileup alleles that should be filtered if they are part of the assembly
     * @param goodPileupEvents Pileup alleles that pass the heuristics to be included in genotyping
     * @param aligner SmithWatermanAligner to use for filtering out equivalent event sets
     * @return unchanged assembly result set if failed, updated haplotypes otherwise
     */
    public static AssemblyResultSet generatePDHaplotypes(final AssemblyResultSet sourceSet,
                                                         final Set<Event> badPileupEvents,
                                                         final Collection<Event> goodPileupEvents,
                                                         final SmithWatermanAligner aligner,
                                                         final AssemblyBasedCallerArgumentCollection args) {
        final Haplotype referenceHaplotype = sourceSet.getReferenceHaplotype();
        final Locatable callingSpan = sourceSet.getRegionForGenotyping().getSpan();

        final PileupDetectionArgumentCollection pileupArgs = args.pileupDetectionArgs;
        final boolean debug = pileupArgs.debugPileupStdout;

        final List<Event> eventsInOrder = makeFinalListOfEventsInOrder(sourceSet, badPileupEvents, goodPileupEvents, referenceHaplotype, pileupArgs, debug);

        // TODO this is where we filter out if indels > 32 (a heuristic known from DRAGEN that is not implemented here)
        SortedMap<Integer, List<Event>> eventsByStartPos = eventsInOrder.stream()
                .collect(Collectors.groupingBy(Event::getStart, TreeMap::new, Collectors.toList()));

        List<List<Event>> disallowedCombinations = smithWatermanRealignPairsOfVariantsForEquivalentEvents(referenceHaplotype, aligner, args.getHaplotypeToReferenceSWParameters(), debug, eventsInOrder);
        dragenDisallowedGroupsMessage(referenceHaplotype.getStart(), debug, disallowedCombinations);

        final List<EventGroup> eventGroups = getEventGroupClusters(eventsInOrder, disallowedCombinations);
        // if any of our merged event groups is too large, abort.
        if (eventGroups == null) {
            Utils.printIf(debug, () -> "Found event group with too many variants! Aborting haplotype building");
            return sourceSet;
        }
        Utils.printIf(debug,() -> "Event groups after merging:\n"+eventGroups.stream().map(eg -> eg.toDisplayString(referenceHaplotype.getStart())).collect(Collectors.joining("\n")));

        Set<Haplotype> outputHaplotypes = new LinkedHashSet<>();
        if (pileupArgs.determinePDHaps) {
            // only add the reference haplotype if we are producing regular haplotype objects (not PD haplotypes for the haplotype alleles)
            outputHaplotypes.add(referenceHaplotype);
        }

        /*
          There are several nested loops here:
          Layer 1: iterate over all determined event positions
          Layer 2: iterate over all alleles at that position, including the reference allele unless we are making determined haplotypes, to set as
                    the determined allele.  Other positions are treated as undetermined,
          Layer 3: make partially determined haplotypes based on each set of branch exclusions
         */
        for (final int determinedLocus : eventsByStartPos.keySet()) {   // it's a SortedMap -- iterating over its keyset is okay!
            final List<Event> allEventsHere = eventsByStartPos.get(determinedLocus);
            Utils.printIf(debug, () -> "working with variants: " + allEventsHere + " at position " + determinedLocus);

            if (!Range.closed(callingSpan.getStart(), callingSpan.getEnd()).contains(determinedLocus)) {
                Utils.printIf(debug, () -> "Skipping determined hap construction! Outside of span: " + callingSpan);
                continue;
            }

            for (int determinedAlleleIndex = (pileupArgs.determinePDHaps?0:-1); determinedAlleleIndex < allEventsHere.size(); determinedAlleleIndex++) { //note -1 for I here corresponds to the reference allele at this site
                final boolean determinedAlleleIsRef = determinedAlleleIndex == -1;
                final Set<Event> determinedEvents = determinedAlleleIsRef ? Set.of() : Set.of(allEventsHere.get(determinedAlleleIndex));
                final Event determinedEventToTest = allEventsHere.get(determinedAlleleIsRef ? 0 : determinedAlleleIndex);
                Utils.printIf(debug, () -> "Working with allele at site: "+(determinedAlleleIsRef? "[ref:"+(determinedLocus-referenceHaplotype.getStart())+"]" : PartiallyDeterminedHaplotype.getDRAGENDebugEventString(referenceHaplotype.getStart()).apply(determinedEventToTest)));

                List<Set<Event>> branchExcludeAlleles = new ArrayList<>();
                branchExcludeAlleles.add(new HashSet<>()); // Add the null branch (assuming no exclusions)

                /* Note for future posterity:
                 * An assembly region could potentially have any number of (within some limitations) of event groups. When we are constructing
                 * haplotypes out of the assembled variants we want to take the dot product of the branches for each set of event groups that
                 * we find. I.E. if an event group with mutex variants (B,C) requires two branches for Variant A and Variant A also leads to two branches in
                 * another event group with mutex variants (D,E). Then we would want to ultimately generate the branches A,B,D -> A,B,E -> A,C,D -> A,C,E.
                 * This is why we iterate over branchExcludeAlleles internally here.
                 */
                for(EventGroup group : eventGroups ) {
                    if (group.causesBranching()) {
                        List<Set<Event>> branchingSets = group.setsForBranching(allEventsHere, determinedEvents, true);
                        // Combinatorially expand the branches as necessary
                        List<Set<Event>> newBranchesToAdd = new ArrayList<>();
                        for (Set<Event> excludedVars : branchExcludeAlleles) {
                            //For every exclude group, fork it by each subset we have:
                            for (int i = 1; i < branchingSets.size(); i++) { //NOTE: iterate starting at 1 here because we special case that branch at the end
                                newBranchesToAdd.add(Sets.union(excludedVars, branchingSets.get(i)).immutableCopy());
                            }
                            // Be careful since this event group might have returned nothing
                            if (!branchingSets.isEmpty()) {
                                excludedVars.addAll(branchingSets.get(0));
                            }
                        }
                        branchExcludeAlleles.addAll(newBranchesToAdd);

                        if (branchExcludeAlleles.size() > MAX_BRANCH_PD_HAPS) {
                            Utils.printIf(debug, () -> "Found too many branches for variants at: " + determinedEventToTest.getStart() + " aborting and falling back to Assembly Variants!");
                            return sourceSet;
                        }
                    }
                }

                branchExcludeAllelesMessage(referenceHaplotype, debug, eventsInOrder, branchExcludeAlleles);


                // from each set of branch exclusions make a single PD haplotype or a combinatorial number of determined haplotypes
                for (Set<Event> excludeEvents : branchExcludeAlleles) {

                    List<Haplotype> branchHaps = new ArrayList<>();
                    List<Event> newBranch = new ArrayList<>();

                    // Handle the simple case of making PD haplotypes
                    if (!pileupArgs.determinePDHaps) {
                        for (final int locus : eventsByStartPos.keySet()) {
                            if (locus == determinedLocus) {
                                newBranch.add(determinedEventToTest);
                            } else {
                                // We know here that nothing illegally overlaps because there are no groups.
                                // Also exclude any events that overlap the determined allele since we cant construct them (also this stops compound alleles from being formed)
                                // NOTE: it is important that we allow reference alleles to overlap undetermined variants as it leads to mismatches against DRAGEN otherwise.
                                eventsByStartPos.get(locus).stream().filter(vc -> !excludeEvents.contains(vc)).forEach(newBranch::add);
                            }
                        }
                        newBranch.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
                        PartiallyDeterminedHaplotype newPDHaplotypeFromEvents = createNewPDHaplotypeFromEvents(referenceHaplotype, determinedEvents, determinedLocus, newBranch, allEventsHere);
                        branchHaps.add(newPDHaplotypeFromEvents);
                    } else {
                        // TODO currently this approach doesn't properly handle a bunch of duplicate events...
                        // If we are producing determined bases, then we want to enforce that every new event at least has THIS event as a variant.
                        List<List<Event>> growingEventGroups = new ArrayList<>();
                        growingEventGroups.add(new ArrayList<>());
                        for (final int locus : eventsByStartPos.keySet()) {
                            // Iterate through the growing combinatorial expansion of haps, split it into either having or not having the variant.
                            if (growingEventGroups.size() > MAX_BRANCH_PD_HAPS) {
                                Utils.printIf(debug, () -> "Too many branch haplotypes ["+growingEventGroups.size()+"] generated from site, falling back on assembly variants!");
                                return sourceSet;
                            } else if (locus == determinedLocus) {
                                growingEventGroups.forEach(group -> group.add(determinedEventToTest));
                            } else if (determinedLocus < locus) {
                                eventsByStartPos.get(locus).stream()
                                        .filter(event -> !excludeEvents.contains(event))
                                        .flatMap(event -> growingEventGroups.stream().map(group -> growEventGroup(group, event)))
                                        .forEach(growingEventGroups::add);
                            }
                        }

                        for (List<Event> subset : growingEventGroups) {
                            subset.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
                            Utils.printIf(debug, () -> "Constructing Hap From Events:"+ formatEventsLikeDragenLogs(subset,  referenceHaplotype.getStart()));
                            branchHaps.add(constructHaplotypeFromEvents(referenceHaplotype, subset, true));
                        }
                    }
                    branchHaplotypesDebugMessage(referenceHaplotype, debug, excludeEvents, branchHaps);

                    outputHaplotypes.addAll(branchHaps);
                    if (outputHaplotypes.size() > MAX_PD_HAPS_TO_GENERATE) {
                        Utils.printIf(debug,() -> "Too many Haps ["+outputHaplotypes.size()+"] generated at this site! Aborting!");
                        return sourceSet;
                    }
                }
            }
        }

        if (outputHaplotypes.size() > MAX_PD_HAPS_TO_GENERATE) {
            Utils.printIf(debug,() -> "Too many branch haplotypes found, aborting ["+outputHaplotypes.size()+"]");
            return sourceSet;
        }
        sourceSet.storeAssemblyHaplotypes();

        // TODO this is an entirely unnecessary step that can be done away with but i leave in because it makes debugging against previous versions much easier.
        final List<Haplotype> result = outputHaplotypes.stream().sorted(Comparator.comparing(Haplotype::getBaseString)).toList();
        sourceSet.replaceAllHaplotypes(result);
        Utils.printIf(debug, () -> "Constructed Haps for Branch"+sourceSet.getHaplotypeList().stream().map(Haplotype::toString).collect(Collectors.joining("\n")));
        sourceSet.setPartiallyDeterminedMode(!pileupArgs.determinePDHaps);
        Utils.printIf(debug, () -> "Returning "+outputHaplotypes.size()+" to the HMM");
        return sourceSet;
    }

    /**
     * Starting from all events from an AssemblyResultSet, remove "bad" events from the pileups that should be filtered
     * according to heuristic in PileupBasedAlleles, then add "good" events from the pileups, excluding SNPs that are too close
     * to assembled indels.
     *
     * @param sourceSet             AssemblyResultSet from which the original set of assembled events comes
     * @param badPileupEvents       Events found in pileups that fail certain heuristic conditions and need to be removed
     * @param goodPileupEvents      Events found in pileups that pass ceetain other heuristic and are to be added
     * @param referenceHaplotype    Only used for a debug message to match DRAGEN
     * @param pileupArgs            Arguments containing, in particular, the distance that "good" pileup SNPs must be
     *                              From any assembled indels in order to be added
     * @param debug                 Controls whether to output certain debugging messages
     * @return                      A final sorted list of events
     */
    private static List<Event> makeFinalListOfEventsInOrder(final AssemblyResultSet sourceSet, final Set<Event> badPileupEvents,
                                                            final Collection<Event> goodPileupEvents, final Haplotype referenceHaplotype,
                                                            final PileupDetectionArgumentCollection pileupArgs, final boolean debug) {
        // start with all assembled events, using maxMnpDistance = 0 because we currently don't support MNPs here,
        // then remove bad pileup events and add good pileup events other than SNPs too close to indels
        removeBadPileupEventsMessage(debug, sourceSet, badPileupEvents);
        final Set<Event> passingEvents = sourceSet.getVariationEvents(0).stream()
                .filter(event -> !badPileupEvents.contains(event))
                .collect(Collectors.toSet());
        final List<Event> indels = passingEvents.stream().filter(Event::isIndel).toList();
        goodPileupEvents.stream().filter(event -> !passingEvents.contains(event)).filter(event -> event.isIndel() ||
                        indels.stream().noneMatch(indel -> event.withinDistanceOf(indel, pileupArgs.snpAdjacentToAssemblyIndel)))
                .forEach(passingEvents::add);
        final List<Event> eventsInOrder = passingEvents.stream().sorted(HAPLOTYPE_SNP_FIRST_COMPARATOR).toList();
        finalEventsListMessage(referenceHaplotype.getStart(), debug, eventsInOrder);
        return eventsInOrder;
    }

    /**
     * Helper method that handles one of the Heuristics baked by DRAGEN into this artificial haplotype generation code.
     *
     * To help mitigate the risk of generating combinatorial haplotypes with SNPs/Indels that that might or might not add
     * up to equivalent events, DRAGEN enforces that events are NOT allowed to be present in the same haplotype if they
     * (when run through smith waterman) add up to other events that were found at this assembly region.
     *
     * To cut down on the complexity of the task; we (and DRAGEN) follow this procedure:
     * 1. look at all sets of 2 variants where at least one is an indel and none overlap.
     *    a) for each set construct an artificial haplotype with only those two variants on it
     *    b) Smith Waterman align it against the reference to generate the cheapest cigar string representation
     *    c) Construct the event map for the new artificial haplotype, if any events in the new event map are in our list of variants
     *       but are NOT the constituent events that were used to construct the haplotype then disallow the pair
     * 2. Look at all sets of 3 variants that do not contain disallowed pairs found in step 1.
     *    a-b-c) repeat steps 1a,1b,and 1c on the 3 event sets
     *
     * @return A list of lists of variant contexts that correspond to disallowed groups. This list may be empty if none are found.
     */
    private static List<List<Event>> smithWatermanRealignPairsOfVariantsForEquivalentEvents(Haplotype referenceHaplotype, SmithWatermanAligner aligner, SWParameters swParameters, boolean debug, List<Event> eventsInOrder) {
        List<List<Event>> disallowedPairs = new ArrayList<>();

        //Iterate over all 2 element permutations in which one element is an indel and test for alignments
        for (int i = 0; i < eventsInOrder.size(); i++) {
            final Event firstEvent = eventsInOrder.get(i);
            if (firstEvent.isIndel()) {
                // For every indel, make every 2-3 element subset (without overlapping) of variants to test for equivalency
                for (int j = 0; j < eventsInOrder.size(); j++) {
                    final Event secondEvent = eventsInOrder.get(j);
                    // Don't compare the event to itself, to overlapping events, or to indels I've already examined me (to prevent double counting)
                    if (j != i && !eventsOverlapForPDHapsCode(firstEvent, secondEvent) && ((!secondEvent.isIndel()) || j > i)) {
                        final List<Event> events = new ArrayList<>(Arrays.asList(firstEvent, secondEvent));
                        events.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
                        Utils.printIf(debug, () -> "Testing events: "+ formatEventsLikeDragenLogs(events,  referenceHaplotype.getStart()));
                        if (constructArtificialHaplotypeAndTestEquivalentEvents(referenceHaplotype, aligner, swParameters, eventsInOrder, events, debug)) {
                            disallowedPairs.add(events);
                        }
                    }
                }
            }
        }

        //TODO NOTE: there are some discrepancies with the iteration over 3x variants in some complicated cases involving
        //TODO       lots of transitively disallowed pairs. Hopefully this is a minor effect.
        //Now iterate over all 3 element pairs and make sure none of the
        for (int i = 0; i < eventsInOrder.size(); i++) {
            final Event firstEvent = eventsInOrder.get(i);
            if (firstEvent.isIndel()) {
                // For every indel, make every 2-3 element subset (without overlapping) of variants to test for equivalency
                for (int j = 0; j < eventsInOrder.size(); j++) {
                    final Event secondEvent = eventsInOrder.get(j);
                    // Don't compare the event to itself, to overlapping events, or to indels I've already examined me (to prevent double counting)
                    if (j != i && !eventsOverlapForPDHapsCode(firstEvent, secondEvent) && ((!secondEvent.isIndel()) || j > i)) {
                        // if i and j area already disallowed keep going
                        if (disallowedPairs.stream().anyMatch(p -> p.contains(firstEvent) && p.contains(secondEvent))) {
                            continue;
                        }
                        final List<Event> events = new ArrayList<>(Arrays.asList(firstEvent, secondEvent));
                        // If our 2 element arrays weren't inequivalent, test subsets of 3 including this:
                        for (int k = j+1; k < eventsInOrder.size(); k++) {
                            final Event thirdEvent = eventsInOrder.get(k);
                            if (k != i && !eventsOverlapForPDHapsCode(thirdEvent, firstEvent) && !eventsOverlapForPDHapsCode(thirdEvent, secondEvent)) {
                                // if k and j or k and i are disallowed, keep looking
                                if (disallowedPairs.stream().anyMatch(p -> (p.contains(firstEvent) && p.contains(thirdEvent)) || (p.contains(secondEvent) && p.contains(thirdEvent)))) {
                                    continue;
                                }
                                List<Event> subList = new ArrayList<>(events);
                                subList.add(thirdEvent);
                                subList.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
                                Utils.printIf(debug,() ->"Testing events: " + formatEventsLikeDragenLogs(subList,  referenceHaplotype.getStart()));
                                if (constructArtificialHaplotypeAndTestEquivalentEvents(referenceHaplotype, aligner, swParameters, eventsInOrder, subList, debug)) {
                                    disallowedPairs.add(subList);
                                }
                            }
                        }
                    }
                }
            }
        }

        return disallowedPairs;
    }

    /**
     * Partition events into the largest possible clusters such that events in distinct clusters are mutually compatible
     * i.e. can exist on the same haplotype.
     *
     * Incompatibility comes in two types, the first and more obvious being overlap.  The second is defined by the
     * heuristic in which we inject two or three events into the reference haplotype and realign the resulting two- or three-event
     * haplotype against the reference haplotype via Smith-Waterman.  If the resulting reduced representation contains other events
     * that have already been found the two or three events are mutually forbidden.
     *
     * To find the desired partition we calculate the connected components of an undirected graph whose edges correspond to
     * either type of incompatibility.  That is, overlapping events and Smith-Waterman-forbidden pairs get an edge; Smith-Waterman-forbidden
     * trios get three edges, one for every two events in the trio.
     *
     * @param eventsInOrder all events in a calling region
     * @param swForbiddenPairsAndTrios list of lists of two or three non-overlapping events that are mutually excluded according
     *                                to the Smith-Waterman heuristic in which injecting them into the reference haplotype and
     *                                simplifying via Smith Waterman yields other events
     */
    private static List<EventGroup> getEventGroupClusters(List<Event> eventsInOrder, List<List<Event>> swForbiddenPairsAndTrios) {
        final List<List<Event>> allMutexes = new ArrayList<>(swForbiddenPairsAndTrios);

        // edges due to overlapping position
        for (int e1 = 0; e1 < eventsInOrder.size(); e1++) {
            final Event event1 = eventsInOrder.get(e1);
            for (int e2 = e1 + 1; e2 < eventsInOrder.size() && eventsInOrder.get(e2).getStart() <= event1.getEnd() + 1; e2++) {
                final Event event2 = eventsInOrder.get(e2);
                if (eventsOverlapForPDHapsCode(event1, event2)) {
                    allMutexes.add(List.of(event1, event2));
                }
            }
        }

        // for each mutex, add edges 0-1, 1-2. . . to put all mutex events in one connected component
        final Graph<Event, DefaultEdge> graph = new SimpleGraph<>(DefaultEdge.class);
        eventsInOrder.forEach(graph::addVertex);
        allMutexes.forEach(mutex -> new IndexRange(0, mutex.size() - 1).forEach(n -> graph.addEdge(mutex.get(n), mutex.get(n+1))));

        final List<Set<Event>> components = new ConnectivityInspector<>(graph).connectedSets();
        return components.stream().anyMatch(comp -> comp.size() > MAX_VAR_IN_EVENT_GROUP) ? null :
                components.stream().map(component -> new EventGroup(component, allMutexes)).toList();
    }

    /**
     * Overlaps method to handle indels and snps correctly. Specifically for this branching codes purposes,
     * indels don't overlap on their anchor bases and insertions don't overlap anything except deletions spanning them or other insertions
     * at the same base.
     *
     */
    @VisibleForTesting
    static boolean eventsOverlapForPDHapsCode(Event e1, Event e2){
        if (!e1.getContig().equals(e2.getContig())) {
            return false;
        }

        double end1 = e1.getEnd() + (e1.isSimpleInsertion() ? 0.5 : 0);
        double end2 = e2.getEnd() + (e2.isSimpleInsertion() ? 0.5 : 0);
        return !(dragenStart(e1) > end2 || dragenStart(e2) > end1);
    }


    /**
     * This method is the helper that manages the actual SW alignment and testing of a group of variants vs the reference haplotype.
     *
     * The method is as follows, construct and artificial haplotype of the provided events, then realign it vs the reference and test
     * if any of the resulting variants are present in the inputs (but doesn't match)
     *
     * NOTE: as per DRAGEN implementation, a set is considered invalid if we re-SmithWaterman align and we get:
     * 1) A different result hap from the ref + alleles we added
     * 2) At least one of the resulting events is NOT in the initial set of alleles we added to the ref
     * 3) Was at least one of the resulting alleles a variant we discovered in the pileups/assembly
     *
     * #3 means that if we add alleles to the ref and realign to get a different representation that was NOT found in the assemblies
     *    we don't consider the events to be a problem, we only care in the even that they add up to a variant that is going to be
     *    in our PDhaplotypes anyway!
     *
     * @return true if we SHOULD NOT allow the eventsToTest alleles to appear as alleles together in determined haplotypes
     */
    @VisibleForTesting
    private static boolean constructArtificialHaplotypeAndTestEquivalentEvents(Haplotype referenceHaplotype, SmithWatermanAligner aligner, SWParameters swParameters, Collection<Event> events, List<Event> eventsToTest, boolean debug) {
        final Haplotype realignHap = constructHaplotypeFromEvents(referenceHaplotype, eventsToTest, false);
        //Special case to capture events that equal the reference (and thus have empty event maps).
        if (Arrays.equals(realignHap.getBases(), referenceHaplotype.getBases())) {
            Utils.printIf(debug,()->"Events add up to the reference! disallowing pair");
            return true;
        }
        //ALIGN!
        realignHap.setCigar(CigarUtils.calculateCigar(referenceHaplotype.getBases(), realignHap.getBases(), aligner, swParameters, SWOverhangStrategy.INDEL));
        EventMap.buildEventMapsForHaplotypes(Collections.singletonList(realignHap), referenceHaplotype.getBases(), referenceHaplotype.getGenomeLocation(), false,0);
        //TODO this differs from DRAGEN, specifically in DRAGEN they do the realignment and compare the SmithWatermanScore with the SWscore of the
        //TODO non-realigned haplotypes, then they only call an event set equivalent if the score for the new found alignment is less than the score
        //TODO for the existing ones. Since we are simply realigning and checking the event map outputs its possible that we consider events to
        //TODO be disallowed that have an equal SmithWaterman score to the original but a different (but equivalent) variant representation.
        //TODO This is likely a minor effect on the overall correctness.
        final boolean wasEquivalentEvent = realignHap.getEventMap().getEvents().stream()
                // Are there any variants NOT in our initial list
                .filter(event -> eventsToTest.stream().noneMatch(event::equals))
                // Do any of variants (that were not in our set of 2-3 targets) appear in our overall list of alleles
                .anyMatch(event -> events.stream().anyMatch(event::equals));
        Utils.printIf(debug, () -> formatEventsLikeDragenLogs(realignHap.getEventMap().getEvents(),  referenceHaplotype.getStart(), "\n"));
        Utils.printIf(debug && wasEquivalentEvent,()->"Events mismatched!");

        return wasEquivalentEvent;
    }
    /**
     * NOTE: this accepts multiple alleles stacked up at the same base (assuming the order is SNP -> INDEL)
     * NOTE: However this class does NOT accept multiple SNPS overlapping or SNPs overlapping deletions
     */
    @VisibleForTesting
    public static Haplotype constructHaplotypeFromEvents(final Haplotype refHap, final List<Event> events, final boolean setEventMap) {
        Utils.validate(refHap.isReference() && refHap.getCigar().numCigarElements() == 1, "This is not a valid base haplotype for construction");

        //ASSERT that everything is fully overlapping the reference.
        events.stream().forEach(v -> Utils.validate(refHap.contains(v), () -> "Provided Variant Context"+v+"doesn't overlap haplotype "+refHap));

        final int refStart = refHap.getStart();
        int lastPositionAdded = refStart;   // genomic coordinate

        byte[] refbases = refHap.getBases();
        CigarBuilder runningCigar = new CigarBuilder();
        final int resultHaplotypeLength = refHap.length() + events.stream().mapToInt(e -> e.altAllele().length() - e.refAllele().length()).sum();
        final ByteBuffer newHapBases = ByteBuffer.allocate(resultHaplotypeLength);

        Utils.validate(IntStream.range(0, events.size()-1).allMatch(n ->
                events.get(n).getEnd() < actualStartExcludingInitialIndelBase(events.get(n+1))), () -> "PD event list: " + events + " is out of order.");

        //ASSUME sorted for now
        for (Event event : events) {
            Allele refAllele = event.refAllele();
            Allele altAllele = event.altAllele();

            final int actualStart = actualStartExcludingInitialIndelBase(event);   // the +1 for indels accounts for the initial alt=ref dummy base

            int basesBeforeNextEvent = actualStart - lastPositionAdded;
            runningCigar.add(new CigarElement(basesBeforeNextEvent, CigarOperator.M));

            final int altRefLengthDiff = altAllele.length() - refAllele.length();
            final CigarElement eventElement = altRefLengthDiff == 0 ? new CigarElement(refAllele.length(), CigarOperator.X) :
                    new CigarElement(Math.abs(altRefLengthDiff), altRefLengthDiff < 0 ? CigarOperator.D : CigarOperator.I);
            runningCigar.add(eventElement);

            // bases before the event, including the dummy initial indel base, followed by event bases, excluding the dummy indel base
            newHapBases.put(ArrayUtils.subarray(refbases, lastPositionAdded - refStart, actualStart - refStart));
            newHapBases.put(altRefLengthDiff == 0 ? altAllele.getBases() : basesAfterFirst(altAllele));

            lastPositionAdded = event.getEnd() + 1; //TODO this is probably not set for future reference
        }

        // bases after the last event
        int refStartIndex = lastPositionAdded - refStart;
        newHapBases.put(ArrayUtils.subarray(refbases, refStartIndex, refbases.length));
        runningCigar.add(new CigarElement(refbases.length - refStartIndex, CigarOperator.M));

        final Haplotype result = new Haplotype(newHapBases.array(), false, refHap, runningCigar.make());
        if (setEventMap) {
            EventMap.buildEventMapsForHaplotypes(List.of(result), refHap.getBases(), refHap, false,0);
            // NOTE: we set this AFTER generating the event maps because the event map code above is being generated from the ref hap so this offset will cause out of bounds errors
            result.setAlignmentStartHapwrtRef(refHap.getAlignmentStartHapwrtRef());
        }
        return result;
    }

    private static int actualStartExcludingInitialIndelBase(final Event event) {
        return event.getStart() + (event.isIndel() ? 1 : 0);
    }

    // skip the dummy intial base of an indel
    private static byte[] basesAfterFirst(final Allele altAllele) {
        return Arrays.copyOfRange(altAllele.getBases(), 1, altAllele.length());
    }

    /**
     * Construct a PD haplotype from scratch
     *
     * Generally we are constructing a new haplotype with all the reference bases for SNP events and with the longest possible allele for INDEL events.
     * For deletions, we extend the haplotype by the ref length
     *
     * NOTE: we assume each provided VC is in start position order, and that if there are overlapping SNPs and indels that the SNPs come first
     */
    @VisibleForTesting
    //TODO When we implement JointDetection we will need to allow multiple eventWithVariants to be prsent...
    static PartiallyDeterminedHaplotype createNewPDHaplotypeFromEvents(final Haplotype refHap, final Set<Event> determinedEvents, final int determinedLocus,
                                                                       final List<Event> constituentEvents, final List<Event> allEventsAtDeterminedLocus) {
        Utils.validate(refHap.isReference() && refHap.getCigar().numCigarElements() == 1, "This is not a valid base haplotype for construction");
        final boolean refIsDetermined = determinedEvents.isEmpty();

        //TODO add a more stringent check that the format of constituentEvents works
        int refStart = refHap.getStart();
        int lastPositionAdded = refStart;

        byte[] refBasesToAddTo = refHap.getBases();
        CigarBuilder runningCigar = new CigarBuilder(false); // NOTE: in some incredibly rare edge cases involving the legacy assembly region trimmer a deletion can hang past the edge of an active window.
        final int upperBoundOnResultLength = refHap.length() + constituentEvents.stream().mapToInt(e -> Math.max(e.altAllele().length(), e.refAllele().length()) - 1).sum();
        final ByteBuffer newHapBases = ByteBuffer.allocate(upperBoundOnResultLength);
        final ByteBuffer pdBytes = ByteBuffer.allocate(upperBoundOnResultLength);

        //ASSUME sorted for now
        // use the reverse list to save myself figuring out cigars for right now
        boolean lastEventWasSnp = false;
        for (Event event : constituentEvents) {

            final int actualStart = event.getStart() + (event.isIndel() ? 1 : 0);   // the +1 for indels accounts for the initial alt=ref dummy base
            int basesBeforeNextEvent = actualStart - lastPositionAdded;

            // Special case for two SNPs at the same position -- merge them into a single undetermined SNP via bitwise OR
            if (basesBeforeNextEvent == -1 && event.isSNP() && lastEventWasSnp) {
                final byte byteForThisSnp = PartiallyDeterminedHaplotype.getPDBytesForHaplotypes(event.refAllele(), event.altAllele())[0];
                final int bufferPositionOfLastSnp = pdBytes.position() - 1;
                final byte byteForLastSnp = pdBytes.get(bufferPositionOfLastSnp);
                pdBytes.put(bufferPositionOfLastSnp,(byte) (byteForLastSnp | byteForThisSnp) );
                continue;
            } else if (event.getStart() == determinedLocus && refIsDetermined) {
                // Ref alleles (even if they overlap undetermined events) should be skipped
                continue;
            }

            Utils.validate(basesBeforeNextEvent >= 0, () -> "Event " + event + " is out of order in PD event list: " + constituentEvents + ".");

            Allele refAllele = event.refAllele();
            Allele altAllele = event.altAllele();
            final int altRefLengthDiff = altAllele.length() - refAllele.length();

            boolean isInsertion = altRefLengthDiff > 0; // If it's an insertion we flip to "ADD" the bases to the ref.
            final boolean isEvent = determinedEvents.stream().anyMatch(event::equals);
            runningCigar.add(new CigarElement(basesBeforeNextEvent, CigarOperator.M));

            // Figure out the cigar element to add:
            // - If we are in the ref, simply add the cigar corresponding to the allele we are using
            if (event.isSNP()) {    // SNPs are straightforward, whether or not it's the special event
                runningCigar.add(new CigarElement(refAllele.length(), refIsDetermined || !isEvent ? CigarOperator.M : CigarOperator.X));
            } else if (isEvent) {   // special indel
                // The event is considered a deletion of the longer allele, regardless of which is ref
                // we subtract 1 for the dummy initial indel base.
                final int elementLength = isInsertion && refIsDetermined ? 0 : Math.max(refAllele.length(), altAllele.length()) - 1;
                final CigarOperator operator = isInsertion ? CigarOperator.I : (refIsDetermined ? CigarOperator.M : CigarOperator.D);
                runningCigar.add(new CigarElement(elementLength, operator));
            } else {    // non-special indel. Insertions are treated as such; deletions become matches
                runningCigar.add(new CigarElement(Math.abs(altRefLengthDiff), altRefLengthDiff > 0 ? CigarOperator.I : CigarOperator.M));
            }

            // bases before the event, including the dummy initial indel base, followed by event bases, excluding the dummy indel base
            newHapBases.put(ArrayUtils.subarray(refBasesToAddTo, lastPositionAdded - refStart, actualStart - refStart)); // bases before the variant
            pdBytes.put(new byte[actualStart - lastPositionAdded]); // bases before the variant -- all zeroes

            // Now add event bases.  If this is the blessed variant, add the ref or alt as appropriate
            // Otherwise make sure we are adding the longest allele (for indels) or the ref allele for snps.
            final boolean alleleToUseIsRef = (isEvent && refIsDetermined) || (!isEvent && altRefLengthDiff <= 0);
            final Allele alleleToUse = alleleToUseIsRef ? refAllele : altAllele;
            final Allele otherAllele = alleleToUseIsRef ? altAllele : refAllele;
            final byte[] basesToAdd = event.isIndel() ? basesAfterFirst(alleleToUse) : alleleToUse.getBases();
            newHapBases.put(basesToAdd); // refbases added
            pdBytes.put(isEvent ? new byte[basesToAdd.length] :
                    PartiallyDeterminedHaplotype.getPDBytesForHaplotypes(alleleToUse, otherAllele)); // refbases added

            lastPositionAdded = event.getEnd() + 1; //TODO this is probably not set for future reference
            lastEventWasSnp = event.isSNP();
        }

        // Finish off the haplotype with the final bases
        int refStartIndex = lastPositionAdded - refStart;
        newHapBases.put(ArrayUtils.subarray(refBasesToAddTo, refStartIndex, refBasesToAddTo.length));
        pdBytes.put(new byte[refBasesToAddTo.length - refStartIndex]);
        runningCigar.add(new CigarElement(refBasesToAddTo.length - refStartIndex, CigarOperator.M));

        return new PartiallyDeterminedHaplotype(
                new Haplotype(ArrayUtils.subarray(newHapBases.array(), 0, newHapBases.position()), false, refHap.getGenomeLocation(), runningCigar.make()),
                ArrayUtils.subarray(pdBytes.array(), 0, pdBytes.position()),
                constituentEvents,
                determinedEvents,
                runningCigar.make(),
                determinedLocus,
                allEventsAtDeterminedLocus,
                refHap.getAlignmentStartHapwrtRef());

    }

    // A helper class for managing mutually exclusive event clusters and the logic around forming valid events vs each other.
    private static class EventGroup {
        private final ImmutableList<Event> eventsInOrder;
        private final ImmutableMap<Event, Integer> eventIndices;

        // There are BitSets and SmallBitSets going around.  The SmallBitSets represent sets of events -- mutexes and other subsets
        // of the event group.  This BitSet is a set of sets of events!  Each element/bit represents a subset of events.  If the bit is true
        // then the corresponding subset of events is allowed i.e. doesn't contain a mutex.  For example, if the 11th bit is true we know
        // that the subset containing events 0, 1, and 3 (11 = 2^3 + 2^1 + 2^0) is allowed because neither {0,1}, {0,3}, {1,3}, or {0,1,3}
        // are mutually excluded.
        private final BitSet allowedSubsets;

        // Optimization to save ourselves recomputing the subsets at every point its necessary to do so.
        List<Set<Event>> cachedEventSets = null;

        /**
         * @param events all events in this event group
         * @param mutexPairsAndTrios all mutexes of two or three events, defined by overlap and by the Smith-Waterman heuristic,
         *                           that determined the event groups in this calling region
         */
        public EventGroup(final Collection<Event> events, List<List<Event>> mutexPairsAndTrios) {
            Utils.validate(events.size() <= MAX_VAR_IN_EVENT_GROUP, () -> "Too many events (" + events.size() + ") for populating bitset.");
            eventsInOrder = events.stream().sorted(HAPLOTYPE_SNP_FIRST_COMPARATOR).collect(ImmutableList.toImmutableList());
            eventIndices = IntStream.range(0, events.size()).boxed().collect(ImmutableMap.toImmutableMap(eventsInOrder::get, n -> n));
            allowedSubsets = new BitSet(1 << eventsInOrder.size());

            final List<List<Event>> overlappingMutexes = mutexPairsAndTrios.stream()
                            .filter(mutext -> mutext.stream().anyMatch(eventIndices::containsKey)).toList();
            for (final List<Event> mutex : overlappingMutexes) {
                Utils.validate(mutex.stream().allMatch(eventIndices::containsKey), () -> "Mutex group " + mutex + " only partially overlaps event group " + this);
            }
            computeAllowedSubsets(overlappingMutexes);
        }

        /**
         * Compute the BitSet of allowed subsets of events i.e. subsets of events that do not contain a mutex group.
         *
         * NOTE: we can use SmallBitSets because we limit ourselves to at most 22 variants per group.
         *
         * @param mutexPairsAndTrios The groups of mutually forbidden events that generated all the event groups in this calling region.
         */


         private void computeAllowedSubsets(List<List<Event>> mutexPairsAndTrios) {
             if (eventsInOrder.size() < 2) {
                 return;
             }

             // initialize all events as being allowed and then disallow them in turn .
             allowedSubsets.set(1, 1 << eventsInOrder.size());

             // make SmallBitSet of the event indices of each mutex
             final List<SmallBitSet> mutexes = mutexPairsAndTrios.stream().map(mutex -> new SmallBitSet(mutex.stream().map(eventIndices::get).toList())).toList();

             // Now forbid all subsets that contain forbidden combinations
             //TODO This method is potentially very inefficient! We don't technically have to iterate over every i,
             //TODO we know there is an optimization involving minimizing the number of checks necessary here by iterating
             //TODO using the bitmask values themselves for the loop
             if (!mutexes.isEmpty()) {
                 for (final SmallBitSet subset = new SmallBitSet().increment(); !subset.hasElementGreaterThan(eventsInOrder.size()); subset.increment()) {
                     if (mutexes.stream().anyMatch(subset::contains)) {
                         allowedSubsets.set(subset.index(), false);
                     }
                 }
             }
         }

        /**
         * This method handles the logic involved in getting all of the allowed subsets of alleles for this event group.
         *
         * @param disallowSubsets
         * @return
         */
        public List<Set<Event>> setsForBranching(final List<Event> locusEvents, final Set<Event> determinedEvents, final boolean disallowSubsets) {
            final SmallBitSet locusOverlapSet = overlapSet(locusEvents);
            final SmallBitSet determinedOverlapSet = overlapSet(determinedEvents);

            // Special case (if we are determining bases outside of this mutex cluster we can reuse the work from previous iterations)
            if (locusOverlapSet.isEmpty() && cachedEventSets != null) {
                return cachedEventSets;
            }

            final List<SmallBitSet> allowedAndDetermined = new ArrayList<>();
            // Iterate from the full set (containing every event) to the empty set (no events), which lets us output the largest possible subsets
            // NOTE: we skip over 0 here since that corresponds to ref-only events, handle those externally to this code
            for (final SmallBitSet subset = SmallBitSet.fullSet(eventsInOrder.size()); !subset.isEmpty(); subset.decrement()) {
                if (allowedSubsets.get(subset.index()) && subset.intersection(locusOverlapSet).equals(determinedOverlapSet)) {
                    // Only check for subsets if we need to
                    if (!disallowSubsets || allowedAndDetermined.stream().noneMatch(group -> group.contains(subset))) {
                        allowedAndDetermined.add(subset.copy());    // copy subset since the decrement() mutates it in-place
                    }
                }
            }

            // Now that we have all the mutex groups, unpack them into lists of variants
            List<Set<Event>> output = new ArrayList<>();
            for (SmallBitSet grp : allowedAndDetermined) {
                Set<Event> newGrp = new HashSet<>();
                for (int i = 0; i < eventsInOrder.size(); i++) {
                    if (!grp.get(i)) {
                        newGrp.add(eventsInOrder.get(i));
                    }
                }
                output.add(newGrp);
            }
            // Cache the result
            if(locusOverlapSet.isEmpty()) {
                cachedEventSets = Collections.unmodifiableList(output);
            }
            return output;
        }

        // create the SmallBitSet of those elements from some collection of events that overlap this EventGroup
        private SmallBitSet overlapSet(final Collection<Event> events) {
            return new SmallBitSet(events.stream().map(e -> eventIndices.getOrDefault(e, -1)).filter(n -> n != -1).toList());
        }

        public boolean causesBranching() {
            return eventsInOrder.size() > 1;
        }

        //Print The event group in Illumina indexed ordering:
        public String toDisplayString(int startPos) {
            return "EventGroup: " + formatEventsLikeDragenLogs(eventsInOrder, startPos);
        }

        public int size() { return eventsInOrder.size(); }
    }

    private static List<Event> growEventGroup(final List<Event> group, final Event event) {
        final List<Event> result = new ArrayList<>(group);
        result.add(event);
        return result;
    }

    // To match DRAGEN we must define a modified start position for indels, which is used to determine overlaps when creating event groups
    private static double dragenStart(final Event event) {
        return event.getStart() + (event.isIndel() ? (event.isSimpleDeletion() ? 1 : 0.5) : 0);
    }

    /**
     * The remaining methods in this class are for debugging only, trying to produce output logs as similar as possible
     * to DRAGEN's for the sake of comparison and functional equivalence.
     */

    private static String formatEventsLikeDragenLogs(final Collection<Event> events, final int refStart) {
        return formatEventsLikeDragenLogs(events, refStart,"->");
    }

    private static String formatEventsLikeDragenLogs(final Collection<Event> events, final int refStart, final CharSequence delimiter) {
        return events.stream()
                .map(PartiallyDeterminedHaplotype.getDRAGENDebugEventString(refStart))
                .collect(Collectors.joining(delimiter));
    }

    private static void removeBadPileupEventsMessage(final boolean debug, final AssemblyResultSet assemblyResultSet, final Set<Event> badPileupEvents) {
        if (debug) {
            final Set<Event> intersection = Sets.intersection(assemblyResultSet.getVariationEvents(0), badPileupEvents);
            Utils.printIf(!intersection.isEmpty(), () ->"Removing assembly variant due to columnwise heuristics: " + intersection);
        }
    }

    private static void finalEventsListMessage(final int refStart, final boolean debug, final Collection<Event> eventsInOrder) {
        Utils.printIf(debug, () -> "Variants to PDHapDetermination:\n" + formatEventsLikeDragenLogs(eventsInOrder, refStart));
    }

    private static void dragenDisallowedGroupsMessage(int refStart, boolean debug, List<List<Event>> disallowedPairs) {
        Utils.printIf(debug, () -> "disallowed groups:" + disallowedPairs.stream().map(group -> formatEventsLikeDragenLogs(group, refStart)).collect(Collectors.joining("\n")));
    }

    private static void branchExcludeAllelesMessage(Haplotype referenceHaplotype, boolean debug, Collection<Event> eventsInOrder, List<Set<Event>> branchExcludeAlleles) {
        if (debug) {
            System.out.println("Branches:");
            for (int i = 0; i < branchExcludeAlleles.size(); i++) {
                final int ifinal = i;
                System.out.println("Branch "+i+" VCs:");
                System.out.println("exclude:" + formatEventsLikeDragenLogs(branchExcludeAlleles.get(i), referenceHaplotype.getStart()));
                //to match dragen debug output for personal sanity
                final Collection<Event> included = eventsInOrder.stream().filter(variantContext -> !branchExcludeAlleles.get(ifinal).contains(variantContext)).collect(Collectors.toList());
                System.out.println("include:" + formatEventsLikeDragenLogs(included, referenceHaplotype.getStart()));
            }
        }
    }

    private static void branchHaplotypesDebugMessage(final Haplotype referenceHaplotype, final boolean debug, final Set<Event> excludeEvents, final List<Haplotype> branchHaps) {
        Utils.printIf(debug, () -> "Constructed Haps for Branch" + formatEventsLikeDragenLogs(excludeEvents,  referenceHaplotype.getStart(), ",") + ":");
        Utils.printIf(debug, () -> branchHaps.stream().map(h -> h.getCigar() + " " + h).collect(Collectors.joining("\n")));
    }

}
