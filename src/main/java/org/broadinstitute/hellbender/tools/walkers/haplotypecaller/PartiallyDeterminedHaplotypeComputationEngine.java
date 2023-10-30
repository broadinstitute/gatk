package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.*;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.Lazy;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang.math.DoubleRange;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

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

    // TODO: make these command line arguments?
    // TODO: is 3 repeats the right criterion for STR, regardless of unit length?
    final static int MAX_STR_UNIT_LENGTH = 6;
    final static int MIN_NUM_STR_REPEATS = 3;

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

        final List<SimpleInterval> strIntervals = findSTRs(referenceHaplotype);

        List<List<Event>> swMutexes = smithWatermanRealignPairsOfVariantsForEquivalentEvents(referenceHaplotype, aligner, args.getHaplotypeToReferenceSWParameters(), debug, eventsInOrder);
        dragenDisallowedGroupsMessage(referenceHaplotype.getStart(), debug, swMutexes);

        // TODO: add in command line argument to activate joint detection and only use the STR overlap detector if joint detection is turned on
        final List<EventGroup> eventGroups = getEventGroupClusters(eventsInOrder, swMutexes, strIntervals);
        // if any of our merged event groups is too large, abort.
        if (eventGroups == null) {
            Utils.printIf(debug, () -> "Found event group with too many variants! Aborting haplotype building");
            return sourceSet;
        }
        Utils.printIf(debug,() -> "Event groups after merging:\n"+eventGroups.stream().map(eg -> eg.toDisplayString(referenceHaplotype.getStart())).collect(Collectors.joining("\n")));

        /*
          Outer loop: iterate over which EventGroup to treat as determined
          Inner loop: iterate over subsets of the determined EventGroup, building PD haplotypes from these determined subsets
                      and branches of undetermined events drawn from other EventGroups

          To be clear, each EventGroup may have multiple determined subsets, and
         */
        final Set<Haplotype> outputHaplotypes = pileupArgs.useDeterminedHaplotypesDespitePdhmmMode ? Sets.newLinkedHashSet(List.of(referenceHaplotype)) : Sets.newLinkedHashSet();  // NOTE: output changes if this is not a LinkedHashSet!

        for (int determinedEventGroupIndex = 0; determinedEventGroupIndex < eventGroups.size(); determinedEventGroupIndex++) {
            final EventGroup determinedEventGroup = eventGroups.get(determinedEventGroupIndex);
            Utils.printIf(debug, () -> "working with determined EventGroup: " + determinedEventGroup.toDisplayString(referenceHaplotype.getStart()));

            // TODO: do we need this check?
            // if (!Range.closed(callingSpan.getStart(), callingSpan.getEnd()).contains(determinedLocus)) {
            //    Utils.printIf(debug, () -> "Skipping determined hap construction! Outside of span: " + callingSpan);
            //    continue;
            //}

            // note that the undetermined event branching from other EventGroups is independent of the determined subset
            final List<Set<Event>> undeterminedBranches = computeUndeterminedBranches(eventGroups, determinedEventGroupIndex);

            if (undeterminedBranches == null) {
                Utils.printIf(debug, () -> "Found too many branches, aborting and falling back to Assembly Variants!");
                return sourceSet;
            }

            // by default, the whole event group is the determined span, but we can optionally revert to pre-joint detection behavior
            // where only one event at a time is determined and the rest of the event group is undetermined
            final List<Pair<Set<Event>, SimpleInterval>> eventSetsAndDeterminedSpans = pileupArgs.disableJointDetection ? determinedEventGroup.determinedUndeterminedHybridSets() :
                    determinedEventGroup.determinedEventSets().stream().map(s -> Pair.of(s, determinedEventGroup.getSpan())).toList();
                    ;

            for (final Pair<Set<Event>, SimpleInterval> eventSetAndDeterminedSpan : eventSetsAndDeterminedSpans) {
                final Set<Event> eventsFromDeterminedGroup = eventSetAndDeterminedSpan.getLeft();
                final SimpleInterval determinedSpan = eventSetAndDeterminedSpan.getRight();
                Utils.printIf(debug, () -> "Constructing PD haplotypes with determined events (if no determined events, PD haplotypes match the reference over the determined span): " +
                        eventsFromDeterminedGroup.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugEventString(referenceHaplotype.getStart())).collect(Collectors.joining(", ")));

                // add the determined events to the undetermined branches
                final List<ImmutableSet<Event>> branches = undeterminedBranches.stream()
                        .map(branch -> Sets.union(branch, eventsFromDeterminedGroup).immutableCopy()).toList();
                branchExcludeAllelesMessage(referenceHaplotype, debug, eventsInOrder, branches);

                for (Set<Event> branch : branches) {
                    if (!pileupArgs.useDeterminedHaplotypesDespitePdhmmMode) {
                        final List<Event> eventsInPDHap = branch.stream().sorted(HAPLOTYPE_SNP_FIRST_COMPARATOR).toList();
                        // TODO: in the disable-joint-detection mode, don't use the span of the determined event group
                        PartiallyDeterminedHaplotype newPDHaplotype = createNewPDHaplotypeFromEvents(referenceHaplotype, eventsFromDeterminedGroup, determinedSpan, eventsInPDHap);
                        branchHaplotypesDebugMessage(referenceHaplotype, debug, branch, List.of(newPDHaplotype));
                        outputHaplotypes.add(newPDHaplotype);
                    } else {
                        // TODO currently this approach doesn't properly handle a bunch of duplicate events...
                        // Start with a single "seed" haplotype containing the determined event(s). At each downstream locus,
                        // every haplotype can grow by every event, or remain unchanged.  For example, if we have haplotypes
                        // H1 and H2 so far and reach a locus with events A and B our list grows to H1/ref, H2/ref, H1/A, H1/B, H2/A, H2/B.
                        List<List<Event>> fullyDeterminedHaplotypes = new ArrayList<>();
                        fullyDeterminedHaplotypes.add(new ArrayList<>(eventsFromDeterminedGroup)); // the seed haplotype

                        for (final int locus : eventsByStartPos.keySet()) {
                            // TODO: don't know if the logic on the next line is right, but this code path could probably be deleted anyway
                            if (determinedEventGroup.getSpan().getEnd() < locus) {
                                final List<List<Event>> children = eventsByStartPos.get(locus).stream().filter(branch::contains)
                                        .flatMap(event -> fullyDeterminedHaplotypes.stream().map(group -> growEventList(group, event))).toList();
                                children.forEach(child -> child.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR));
                                fullyDeterminedHaplotypes.addAll(children);

                                if (fullyDeterminedHaplotypes.size() > MAX_BRANCH_PD_HAPS) {
                                    Utils.printIf(debug, () -> "Too many branch haplotypes [" + fullyDeterminedHaplotypes.size() + "] generated from site, falling back on assembly variants!");
                                    return sourceSet;
                                }
                            }
                        }
                        fullyDeterminedHaplotypes.forEach(events -> Utils.printIf(debug, () -> "Constructing Haplotype From Events:" + formatEventsLikeDragenLogs(events, referenceHaplotype.getStart())));
                        final List<Haplotype> branchHaps = fullyDeterminedHaplotypes.stream()
                                .map(events -> constructHaplotypeFromEvents(referenceHaplotype, events, true)).toList();
                        branchHaplotypesDebugMessage(referenceHaplotype, debug, branch, branchHaps);
                        outputHaplotypes.addAll(branchHaps);
                    }

                    if (outputHaplotypes.size() > MAX_PD_HAPS_TO_GENERATE) {
                        Utils.printIf(debug,() -> "Too many branch haplotypes found, aborting ["+outputHaplotypes.size()+"]");
                        return sourceSet;
                    }
                }
            }
        }

        sourceSet.storeAssemblyHaplotypes();

        // TODO: Sorting haplotypes is unnecessary but makes debugging against previous versions much easier.
        final List<Haplotype> result = outputHaplotypes.stream().sorted(Comparator.comparing(Haplotype::getBaseString)).toList();
        sourceSet.replaceAllHaplotypes(result);
        Utils.printIf(debug, () -> "Constructed Haps for Branch"+sourceSet.getHaplotypeList().stream().map(Haplotype::toString).collect(Collectors.joining("\n")));

        // Set PDHMM flag to enable later PDHMM genotyping.  We only set it now because earlier code might fail and revert to non-PDHMM mode.
        sourceSet.setPartiallyDeterminedMode(!pileupArgs.useDeterminedHaplotypesDespitePdhmmMode);
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

    // TODO: should there be any padding at the edges of STRs?
    @VisibleForTesting
    static List<SimpleInterval> findSTRs(final Haplotype referenceHaplotype) {
        final List<SimpleInterval> strs = new ArrayList<>();
        final DragstrReferenceAnalyzer analyzer = DragstrReferenceAnalyzer.of(referenceHaplotype.getBases(), 0, referenceHaplotype.length(), MAX_STR_UNIT_LENGTH);

        final int refStart = referenceHaplotype.getStart();
        final String refContig = referenceHaplotype.getContig();
        int n = 0;
        while (n < referenceHaplotype.length()) {
            final int numRepeats = analyzer.repeatLength(n);
            if (numRepeats >= MIN_NUM_STR_REPEATS) {  // found an STR starting at n
                final int strLength = analyzer.repeatUnit(n).length * numRepeats;
                strs.add(new SimpleInterval(refContig, refStart + n, refStart + n + strLength));
                n += strLength; // jump to end of this STR
            } else {
                n++;
            }
        }

        return strs;
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
     * Partition events into the smallest possible contiguous clusters such that events in different clusters are compatible i.e. can
     * exist on the same haplotype.
     *
     * More precisely, given a collection of Events, a list of tuples of events that are mutually excluded, and (optionally)
     * a set of STR intervals, the resulting partition into EventGroups satisfies the following:
     *  1) Each Event is assigned to exactly one EventGroup.
     *  2) Events that overlap are assigned to the same EventGroup.
     *  3) Events that overlap the same STR interval are assigned to the same EventGroup.
     *  4) Events that belong to the same mutex are assigned to the same EventGroup.
     *  5) Each EventGroup contains every Event in its span.  That is, between the first Event start and the last Event end
     *     of an EventGroup there are no Events assigned to another EventGroup.
     *  6) No EventGroup can be split into proper subsets while satisfying 2-5.
     *
     * In the common case that all events are compatible and don't overlap the same STRs -- they are non-overlapping SNPs, for example --
     * the result is a bunch of singleton event groups.
     *
     * The input mutexes are defined by the heuristic in which we inject two or three events into the reference haplotype and
     * realign the resulting two- or three-event haplotype against the reference haplotype via Smith-Waterman.  If the resulting
     * reduced representation contains other events that have already been found the two or three events are mutually forbidden.
     *
     * To find the desired partition we imagine lighting up genomic intervals:
     *  1) light up the span of every STR
     *  2) light up the span from the leftmost start to the rightmost end of every mutex
     *  3) light up the span from the start to the end of each Event
     * Each resulting contiguous lit-up interval defines an EventGroup.
     *
     * Note: The resulting partition is the same as DRAGEN's except for the extremely rare (perhaps impossible) edge case in which
     * DRAGEN would yield a non-contiguous EventGroup.  This could *only* occur if a Smith-Waterman mutex somehow happened
     * without an STR (which itself is extremely dubious) and the Events of that mutex spanned an intervening non-overlapping
     * Event that was compatible with all the Events of the mutex.
     *
     * @param eventsInOrder all events in a calling region
     * @param swMutexes list of lists of two or three non-overlapping events that are mutually excluded according
     *                                to the Smith-Waterman heuristic in which injecting them into the reference haplotype and
     *                                simplifying via Smith Waterman yields other events
     */
    @VisibleForTesting
    static List<EventGroup> getEventGroupClusters(List<Event> eventsInOrder, List<List<Event>> swMutexes, List<SimpleInterval> strIntervals) {
        // In addition to the Smith-Waterman mutexes, we find the mutexes due to overlap
        // Note that overlapping SNPs are allowed on the undetermined parts of PD haplotypes.  Here we include overlapping SNPS
        // as mutexes and let the EventGroup class handle them.
        final List<List<Event>> overlapMutexes = new ArrayList<>();
        for (int e1 = 0; e1 < eventsInOrder.size(); e1++) {
            final Event event1 = eventsInOrder.get(e1);
            for (int e2 = e1 + 1; e2 < eventsInOrder.size() && eventsInOrder.get(e2).getStart() <= event1.getEnd() + 1; e2++) {
                final Event event2 = eventsInOrder.get(e2);

                // TODO: we need to differentiate overlap mutexes from SW mutexes in EventGroup constructor
                if (event1.overlaps(event2)) {  // the EventGroup constructor will sort out when the stricter condition of DRAGEN overlap is relevant
                    overlapMutexes.add(List.of(event1, event2));
                }
            }
        }

        final List <DoubleRange> allIntervals = new ArrayList<>();
        strIntervals.forEach(str -> allIntervals.add(new DoubleRange(str.getStart(), str.getEnd())));   // light up STRs
        eventsInOrder.forEach(event -> allIntervals.add(new DoubleRange(event.getStart(), dragenEnd(event))));    // light up each Event from VCF start to DRAGEN end
        swMutexes.forEach(tuple -> allIntervals.add(unionOfDragenAndVCFSpans(tuple)));  // light up each mutex
        allIntervals.sort(Comparator.comparingDouble(DoubleRange::getMinimumDouble));   // no need to break ties here -- the end result is unique regardless
        allIntervals.add(new DoubleRange(Double.MAX_VALUE));    // add a dummy interval at the end with no possible overlap to simplify the loop below

        final List<EventGroup> eventGroups = new ArrayList<>();
        double start = allIntervals.get(0).getMinimumDouble();
        double end = allIntervals.get(0).getMaximumDouble();
        int eventIndex = 0;
        final List<Event> eventsInDragenOrder = eventsInOrder.stream().sorted(Comparator.comparingDouble(e -> dragenStart(e))).toList();

        for (final DoubleRange interval : allIntervals) {
            if (interval.getMinimumDouble() <= end) {    // contiguous with previous span.  Note that DoubleRanges are closed (inclusive) on both ends
                end = Math.max(end, interval.getMaximumDouble());
            } else {    // add the previous span and begin a new one.  Note that the dummy last interval at triggers processing of the last EventGroup
                final DoubleRange span = new DoubleRange(start, end);
                final Set<Event> eventsInSpan = new HashSet<>();
                while (eventIndex < eventsInDragenOrder.size() && span.overlapsRange(new DoubleRange(eventsInDragenOrder.get(eventIndex).getStart(), dragenEnd(eventsInDragenOrder.get(eventIndex))))) {
                    eventsInSpan.add(eventsInDragenOrder.get(eventIndex));
                    eventIndex++;
                }
                if (!eventsInSpan.isEmpty()) {  // the interval could be an STR with no Events within
                    if (eventsInSpan.size() > MAX_VAR_IN_EVENT_GROUP) {
                        return null;
                    }
                    try {
                        eventGroups.add(new EventGroup(eventsInSpan, swMutexes, overlapMutexes));
                    } catch (Exception e) {
                        System.out.println("The exception was: " + e);
                        System.out.println("Events in order: " + eventsInOrder);
                        System.out.println("Events in DRAGEN order: " + eventsInDragenOrder);
                        System.out.println("Event groups found so far: " + eventGroups);
                        System.out.println("All intervals in order: " + allIntervals);
                        System.out.println("current span: " + span);
                        System.out.println("Events in the current span: " + eventsInSpan);
                        System.out.println("All SW mutexes here: " + swMutexes);
                        System.out.println("All STR intervals here: " + strIntervals);
                        throw e;
                    }
                }

                start = interval.getMinimumDouble();
                end = interval.getMaximumDouble();
            }
        }

        return eventGroups;
    }


    // same as above without the STR overlap detector
    @VisibleForTesting
    static List<EventGroup> getEventGroupClusters(List<Event> eventsInOrder, List<List<Event>> swForbiddenPairsAndTrios) {
        return getEventGroupClusters(eventsInOrder, swForbiddenPairsAndTrios, List.of());
    }

    /**
     * Compute the branches (sets of mutually allowed Events) that define partially determined haplotypes over their undetermined parts.
     *
     * Each undetermined event group has one or more 1) maximal subsets that 2) contain no mutually-exclusive pairs or trios of events.
     * When constructing PD haplotypes we take all possible unions of these subsets over all undetermined event groups.
     *
     * Returns null if the combinatorial explosion of branches becomes too large.
     */
    @VisibleForTesting
    static List<Set<Event>> computeUndeterminedBranches(final List<EventGroup> eventGroups, final int determinedEventGroupIndex) {
        List<Set<Event>> branches = new ArrayList<>();
        branches.add(new HashSet<>());   // start with a single empty branch

        for (int eventGroupIndex = 0; eventGroupIndex < eventGroups.size(); eventGroupIndex++) {
            if (eventGroupIndex == determinedEventGroupIndex) {
                continue;
            }
            final EventGroup group = eventGroups.get(eventGroupIndex);

            final List<Set<Event>> setsToAdd = group.undeterminedEventSets();

            // Take every possible union of existing branches and this event group's sets to add. As an optimization
            // we add the 0th set's elements in-place, and append unions with the 1st, 2nd etc sets to add.
            final List<HashSet<Event>> extraBranches = setsToAdd.size() < 2 ? List.of() : branches.stream()
                    .flatMap(branch -> setsToAdd.stream().skip(1).map(setToAdd -> Sets.newHashSet(Sets.union(branch, setToAdd))))
                    .toList();

            // add the 0th exclusion set in-place
            if (!setsToAdd.isEmpty()) {
                branches.forEach(branch -> branch.addAll(setsToAdd.get(0)));
            }
            branches.addAll(extraBranches);

            if (branches.size() > MAX_BRANCH_PD_HAPS) {
                return null;
            }
        }
        return branches;
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
    static PartiallyDeterminedHaplotype createNewPDHaplotypeFromEvents(final Haplotype refHap, final Set<Event> determinedEvents, final SimpleInterval determinedSpan,
                                                                       final List<Event> constituentEvents) {
        Utils.validate(refHap.isReference() && refHap.getCigar().numCigarElements() == 1, "This is not a valid base haplotype for construction");

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
            }

            Utils.validate(basesBeforeNextEvent >= 0, () -> "Event " + event + " is out of order in PD event list: " + constituentEvents + ".");

            Allele refAllele = event.refAllele();
            Allele altAllele = event.altAllele();
            final int altRefLengthDiff = altAllele.length() - refAllele.length();

            boolean isInsertion = altRefLengthDiff > 0; // If it's an insertion we flip to "ADD" the bases to the ref.
            final boolean isDeterminedEvent = determinedEvents.stream().anyMatch(event::equals);    // TODO: could this be replaced by simply a call to contains?
            runningCigar.add(new CigarElement(basesBeforeNextEvent, CigarOperator.M));

            // Figure out the cigar element to add:
            // - If we are in the ref, simply add the cigar corresponding to the allele we are using
            if (event.isSNP()) {    // SNPs are straightforward, whether or not it's the special event
                runningCigar.add(new CigarElement(refAllele.length(), !isDeterminedEvent ? CigarOperator.M : CigarOperator.X));
            } else if (isDeterminedEvent) {
                // TODO: this comment seems wrong -- it seems like for determined events there is no hocus pocus where we treat insertions as deletions
                // The event is considered a deletion of the longer allele, regardless of which is ref
                // we subtract 1 for the dummy initial indel base.
                final int elementLength = Math.max(refAllele.length(), altAllele.length()) - 1;
                final CigarOperator operator = isInsertion ? CigarOperator.I : CigarOperator.D;
                runningCigar.add(new CigarElement(elementLength, operator));
            } else {    // non-special indel. Insertions are treated as such; deletions become matches
                // in the undetermined part of PD haplotypes, indels are always treated as deletions relative to the longer allele,
                // regardless of which is ref.  Hence when we have an insertion the alt allele is treated as a match and the ref allele is
                // treated as a deletion
                // TODO: double-check this line
                runningCigar.add(new CigarElement(Math.abs(altRefLengthDiff), altRefLengthDiff > 0 ? CigarOperator.I : CigarOperator.M));
            }

            // bases before the event, including the dummy initial indel base, followed by event bases, excluding the dummy indel base
            newHapBases.put(ArrayUtils.subarray(refBasesToAddTo, lastPositionAdded - refStart, actualStart - refStart)); // bases before the variant
            pdBytes.put(new byte[actualStart - lastPositionAdded]); // bases before the variant -- all zeroes

            // Now add event bases.  If this is the blessed variant, add the ref or alt as appropriate
            // Otherwise make sure we are adding the longest allele (for indels) or the ref allele for snps.
            final boolean alleleToUseIsRef = (!isDeterminedEvent && altRefLengthDiff <= 0);
            final Allele alleleToUse = alleleToUseIsRef ? refAllele : altAllele;
            final Allele otherAllele = alleleToUseIsRef ? altAllele : refAllele;
            final byte[] basesToAdd = event.isIndel() ? basesAfterFirst(alleleToUse) : alleleToUse.getBases();
            newHapBases.put(basesToAdd); // refbases added
            pdBytes.put(isDeterminedEvent ? new byte[basesToAdd.length] :
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
                determinedSpan,
                refHap.getAlignmentStartHapwrtRef());

    }

    // A helper class for managing mutually exclusive event clusters and the logic around forming valid events vs each other.
    @VisibleForTesting
    static class EventGroup {
        private final ImmutableList<Event> eventsInOrder;
        private final ImmutableMap<Event, Integer> eventIndices;

        // There are BitSets and SmallBitSets going around.  The SmallBitSets represent sets of events -- mutexes and other subsets
        // of the event group.  This BitSet is a set of sets of events!  Each element/bit represents a subset of events.  If the bit is true
        // then the corresponding subset of events is allowed i.e. doesn't contain a mutex.  For example, if the 11th bit is true we know
        // that the subset containing events 0, 1, and 3 (11 = 2^3 + 2^1 + 2^0) is allowed because neither {0,1}, {0,3}, {1,3}, or {0,1,3}
        // are mutually excluded.
        private final BitSet allowedUndeterminedSubsets;

        // same as above, except that overlapping SNPs are not allowed to coexist in determined subsets of event
        private final BitSet allowedDeterminedSubsets;

        private final Lazy<List<Set<Event>>> lazyUndeterminedEventSets;

        private final SimpleInterval span;

        /**
         * @param events all events in this event group
         * @param swMutexes all mutexes of two or three events, defined by the Smith-Waterman heuristic,
         *                           that determined the event groups in this calling region
         */
        public EventGroup(final Collection<Event> events, final List<List<Event>> swMutexes, final List<List<Event>> overlapMutexes) {
            Utils.validate(events.size() <= MAX_VAR_IN_EVENT_GROUP, () -> "Too many events (" + events.size() + ") for populating bitset.");
            eventsInOrder = events.stream().sorted(HAPLOTYPE_SNP_FIRST_COMPARATOR).collect(ImmutableList.toImmutableList());
            eventIndices = IntStream.range(0, events.size()).boxed().collect(ImmutableMap.toImmutableMap(eventsInOrder::get, n -> n));

            final List<List<Event>> relevantSWMutexes = swMutexes.stream().filter(mutex -> eventIndices.containsKey(mutex.get(0))).toList();

            // the Smith-Waterman mutexes are non-negotiable and apply in both the determined and undetermined parts of PD haplotypes
            final List<List<Event>> mutexesForDeterminedEvents = new ArrayList<>(relevantSWMutexes);
            final List<List<Event>> mutexesForUndeterminedEvents = new ArrayList<>(relevantSWMutexes);

            // overlapping SNPs can't both be determined but they may coexist in the undetermined span by being merged into an undetermined SNP pseudo-allele
            // two indels or one indel and a SNP that overlap in the DRAGEN sense (excluding the dummy leading base etc) may never coexist
            for (final List<Event> mutex : overlapMutexes) {
                if (!eventIndices.containsKey(mutex.get(0))) {
                    continue;
                }
                if (mutex.get(0).isSNP() && mutex.get(1).isSNP()) { // if they are both SNPs and are in overlapMutexes, they must be at the same locus, no check needed
                    mutexesForDeterminedEvents.add(mutex);
                } else if (eventsOverlapForPDHapsCode(mutex.get(0), mutex.get(1))) {    // the strict DRAGEN sense of overlap is a strict mutex, except for undetermined SNPs as above
                    mutexesForDeterminedEvents.add(mutex);
                    mutexesForUndeterminedEvents.add(mutex);
                } else {    // VCF overlap but not DRAGEN overlap -- they may coexist behind the scenes in the undetermined part of PD haplotypes
                    mutexesForDeterminedEvents.add(mutex);
                }
            }

            for (final List<Event> mutex : mutexesForDeterminedEvents) {
                Utils.validate(mutex.stream().allMatch(eventIndices::containsKey), () -> "Mutex group " + mutex + " only partially overlaps event group " + this);
            }

            allowedUndeterminedSubsets = computeAllowedSubsets(mutexesForUndeterminedEvents);
            allowedDeterminedSubsets = computeAllowedSubsets(mutexesForDeterminedEvents);
            lazyUndeterminedEventSets = new Lazy<>(() -> computeUndeterminedEventSetsForPDHaplotypes());

            final String contig = eventsInOrder.get(0).getContig();
            final int minStart = eventsInOrder.stream().mapToInt(Event::getStart).min().getAsInt();
            final int maxEnd = eventsInOrder.stream().mapToInt(Event::getEnd).max().getAsInt();
            span = new SimpleInterval(contig, minStart, maxEnd);
        }

        /**
         * Compute the BitSet of allowed subsets of events i.e. subsets of events that do not contain a mutex group.
         *
         * NOTE: we can use SmallBitSets because we limit ourselves to at most 22 variants per group.
         *
         * @param mutexPairsAndTrios The groups of mutually forbidden events that generated all the event groups in this calling region.
         */
         private BitSet computeAllowedSubsets(List<List<Event>> mutexPairsAndTrios) {
             final BitSet result = new BitSet(1 << eventsInOrder.size());
             // initialize all subsets as being allowed and then disallow any subset containing a mutex below
             result.set(0, 1 << eventsInOrder.size());

             if (eventsInOrder.size() < 2) {
                 return result;
             }

             // make SmallBitSet of the event indices of each mutex
             final List<SmallBitSet> mutexes = mutexPairsAndTrios.stream().map(mutex -> new SmallBitSet(mutex.stream().map(eventIndices::get).toList())).toList();

             // Now forbid all subsets that contain forbidden combinations
             //TODO This method is potentially very inefficient! We don't technically have to iterate over every i,
             //TODO we know there is an optimization involving minimizing the number of checks necessary here by iterating
             //TODO using the bitmask values themselves for the loop
             if (!mutexes.isEmpty()) {
                 for (final SmallBitSet subset = new SmallBitSet().increment(); !subset.hasElementGreaterThan(eventsInOrder.size()); subset.increment()) {
                     if (mutexes.stream().anyMatch(subset::contains)) {
                         result.set(subset.index(), false);
                     }
                 }
             }
             return result;
         }

        /**
         * Find a minimal and complete collection of allowed subsets of this EventGroup for the purposes of the
         * undetermined parts of PD haplotypes.
         *
         * That is, we find every subset of this EventGroup that satisfies
         *
         * 1) contains no two or three mutually exclusive elements i.e. is an allowed subset as computed in this object's constructor
         * 2) is not a subset of any larger subset that satisfies 1) and 2)
         *
         * The resulting subsets need not be disjoint.
         *
         * Because these subsets contain no mutexes, one can sensibly make partially-determined haplotypes out of them,
         * with each subset entailing a unique set of undetermined alleles.  Property 2) above ensures that we make as
         * few partially-determined haplotypes as possible.
         *
         * @return
         */
        private List<Set<Event>> computeUndeterminedEventSetsForPDHaplotypes() {
            if (eventsInOrder.size() == 1) {
                return Collections.singletonList(Set.of(eventsInOrder.get(0))); // if only one event, there are no mutexes
            }

            final List<SmallBitSet> maximalAllowedSubsets = new ArrayList<>();

            // TODO: I think this comment may be wrong
            // Iterate from the full set (containing every event) to the empty set (no events), which lets us output the largest possible subsets
            // NOTE: we skip over 0 here since that corresponds to ref-only events, handle those externally to this code
            for (final SmallBitSet subset = SmallBitSet.fullSet(eventsInOrder.size()); !subset.isEmpty(); subset.decrement()) {
                if (allowedUndeterminedSubsets.get(subset.index()) && maximalAllowedSubsets.stream().noneMatch(group -> group.contains(subset))) {
                    maximalAllowedSubsets.add(subset.copy());    // copy subset since the decrement() mutates it in-place
                }
            }

            // Now that we have all the maximal sets, map bitset indices to Events
            return maximalAllowedSubsets.stream().map(bitset -> bitset.asSet(eventsInOrder)).toList();
        }

        public List<Set<Event>> undeterminedEventSets() {
            return lazyUndeterminedEventSets.get();
        }

        /**
         * Find every allowed subset of this EventGroup for the purposes of the determined part of PD haplotypes.  This simply
         * means every subset (including the empty and full subsets) that does not contain a mutex.
         *
         * @return
         */
        public List<Set<Event>> determinedEventSets() {
            final List<Set<Event>> result = new ArrayList<>();
            for (final SmallBitSet subset = new SmallBitSet(); !subset.hasElementGreaterThan(eventsInOrder.size()); subset.increment()) {
                if (allowedDeterminedSubsets.get(subset.index())) {
                    result.add(subset.asSet(eventsInOrder));
                }
            }

            return result;
        }

        /**
         * Pre-joint detection behavior where individual events are set as determined and the rest of the event group is
         * undetermined.
         * @return a list of branches, each branch a mix of determined and undetermined events from this EventGroup paired with
         * a SimpleInterval denoting which span is considered determined
         */
        public List<Pair<Set<Event>, SimpleInterval>> determinedUndeterminedHybridSets() {
            final List<Pair<Set<Event>, SimpleInterval>> result = new ArrayList<>();

            // case where ref is determined -- every event at a given locus shares the same undetermined sets
            for (final int locus : eventsInOrder.stream().map(Event::getStart).collect(Collectors.toSet())) {
                // TODO: should this instead be all events that overlap the locus?
                final List<Event> eventsAtDeterminedLocus = eventsInOrder.stream().filter(event -> event.getStart() == locus).toList();
                final SimpleInterval determinedSpan = IntervalUtils.getSpanningInterval(eventsAtDeterminedLocus);
                final SmallBitSet locusEvents = new SmallBitSet(eventsAtDeterminedLocus.stream().map(eventIndices::get).toList());

                final List<SmallBitSet> maximalAllowedSubsets = new ArrayList<>();

                // Iterate from the full set (containing every event) to the empty set (no events), which lets us output the largest possible subsets
                // conditions 1) allowed 2) doesn't overlap the determined locus 3) maximal
                for (final SmallBitSet subset = SmallBitSet.fullSet(eventsInOrder.size()); !subset.isEmpty(); subset.decrement()) {
                    if (allowedUndeterminedSubsets.get(subset.index()) && subset.intersection(locusEvents).isEmpty() &&
                            maximalAllowedSubsets.stream().noneMatch(group -> group.contains(subset))) {
                        maximalAllowedSubsets.add(subset.copy());    // copy subset since the decrement() mutates it in-place
                    }
                }

                maximalAllowedSubsets.forEach(mas -> result.add(Pair.of(mas.asSet(eventsInOrder), determinedSpan)));
            }

            // case where one event at a time is determined
            for (final Event determinedEvent : eventsInOrder) {
                final int determinedEventIndex = eventIndices.get(determinedEvent);
                // here locusEvents exclude the determinedEvent
                final SmallBitSet locusEvents = new SmallBitSet(eventsInOrder.stream().filter(determinedEvent::overlaps).map(eventIndices::get)
                        .filter(i -> i != determinedEventIndex).toList());

                final List<SmallBitSet> maximalAllowedSubsets = new ArrayList<>();

                // Iterate from the full set (containing every event) to the empty set (no events), which lets us output the largest possible subsets
                // conditions 1) allowed 2) contains the determined event but nothing else at that locus 3) maximal
                for (final SmallBitSet subset = SmallBitSet.fullSet(eventsInOrder.size()); !subset.isEmpty(); subset.decrement()) {
                    if (subset.get(determinedEventIndex) && allowedUndeterminedSubsets.get(subset.index()) && subset.intersection(locusEvents).isEmpty() &&
                            maximalAllowedSubsets.stream().noneMatch(group -> group.contains(subset))) {
                        maximalAllowedSubsets.add(subset.copy());    // copy subset since the decrement() mutates it in-place
                    }
                }
                maximalAllowedSubsets.forEach(mas -> result.add(Pair.of(mas.asSet(eventsInOrder), new SimpleInterval(determinedEvent))));
            }

            // Now that we have all the maximal sets, map bitset indices to Events
            return result;
        }

        public SimpleInterval getSpan() { return span; }

        //Print The event group in Illumina indexed ordering:
        public String toDisplayString(int startPos) {
            return "EventGroup: " + formatEventsLikeDragenLogs(eventsInOrder, startPos);
        }

        public int size() { return eventsInOrder.size(); }

        @VisibleForTesting
        List<Event> eventsInOrderForTesting() { return eventsInOrder; }
    }

    private static List<Event> growEventList(final List<Event> eventList, final Event event) {
        final List<Event> result = new ArrayList<>(eventList);
        result.add(event);
        return result;
    }

    // To match DRAGEN we must define a modified start position for indels, which is used to determine overlaps when creating event groups
    private static double dragenStart(final Event event) {
        return event.getStart() + (event.isIndel() ? (event.isSimpleDeletion() ? 1 : 0.5) : 0);
    }

    private static double dragenEnd(final Event event) {
        return event.getEnd() + (event.isSimpleInsertion() ? 0.5 : 0);
    }

    // the union of VCF span and DRAGEN span uses the VCF start (which can be before the DRAGEN start) and the DRAGEN end (which can be after the VCF end)
    private static final DoubleRange unionOfDragenAndVCFSpans(final List<Event> tuple) {
        final int size = tuple.size();
        Utils.validate(size == 2 || size == 3, "this is supposed to be a pair or a trio");
        final double start01 = Math.min(tuple.get(0).getStart(), tuple.get(1).getStart());
        final double end01 = Math.max(dragenEnd(tuple.get(0)), dragenEnd(tuple.get(1)));
        return size == 2 ? new DoubleRange(start01, end01) :
                new DoubleRange(Math.min(start01, tuple.get(2).getStart()), Math.max(end01, dragenEnd(tuple.get(2))));
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

    private static void branchExcludeAllelesMessage(Haplotype referenceHaplotype, boolean debug, Collection<Event> eventsInOrder, List<? extends Set<Event>> branches) {
        if (debug) {
            System.out.println("Branches:");
            for (int i = 0; i < branches.size(); i++) {
                System.out.println("Branch " + i + " VCs:");
                final Set<Event> included = branches.get(i);
                final Set<Event> excluded = eventsInOrder.stream().filter(event -> !included.contains(event)).collect(Collectors.toSet());
                System.out.println("exclude:" + formatEventsLikeDragenLogs(excluded, referenceHaplotype.getStart()));
                //to match dragen debug output for personal sanity
                System.out.println("include:" + formatEventsLikeDragenLogs(included, referenceHaplotype.getStart()));
            }
        }
    }

    private static void branchHaplotypesDebugMessage(final Haplotype referenceHaplotype, final boolean debug, final Set<Event> excludeEvents, final List<Haplotype> branchHaps) {
        Utils.printIf(debug, () -> "Constructed Haps for Branch" + formatEventsLikeDragenLogs(excludeEvents,  referenceHaplotype.getStart(), ",") + ":");
        Utils.printIf(debug, () -> branchHaps.stream().map(h -> h.getCigar() + " " + h).collect(Collectors.joining("\n")));
    }

}
