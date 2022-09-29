package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.Tuple;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * TODO
 */
public class PartiallyDeterminedHaplotypeComputationEngine {
    final static int MAX_PD_HAPS_TO_GENERATE = 256*2; //(2048 is illuminas #) (without optimizing the hmm to some degree this is probably unattainable)
    final static int MAX_BRANCH_PD_HAPS = 128; //(128 is illuminas #)
    final static int MAX_VAR_IN_EVENT_GROUP = 15; // (20 is illuminas #)

    //To make this somewhat cleaner of a port from Illumina, we have two base spaces. R and U space. R space is vcf coordinate space,
    //U is a 0->N (N = region size) based space where Insertions are +0.5 and deletions are + 1 from their original position

    public static final Comparator<VariantContext> HAPLOTYPE_SNP_FIRST_COMPARATOR = Comparator.comparingInt(VariantContext::getStart)
            // Decide arbitrarily so as not to accidentally throw away overlapping variants
            .thenComparingInt(vc -> vc.getReference().length())
            .thenComparingInt(vc -> vc.getAlternateAllele(0).length())
            .thenComparing(vc -> vc.getAlternateAllele(0));


    /**
     * TODO
     */
    public static AssemblyResultSet generatePDHaplotypes(final AssemblyResultSet sourceSet,
                                                         final Haplotype referenceHaplotype,
                                                         final SortedSet<VariantContext> assemblyVariants,
                                                         final List<VariantContext> pileupAllelesFoundShouldFilter,
                                                         final List<VariantContext> pileupAllelesPassingFilters,
                                                         final int snpAdjacentToIndelLimit,
                                                         final SmithWatermanAligner aligner,
                                                         final SWParameters swParameters,
                                                         final boolean makeDeterminedHapsInstead,
                                                         final boolean debugSite) {
        List<Haplotype> output = new ArrayList<>();

        final TreeSet<VariantContext> variantsInOrder = new TreeSet<>(
                HAPLOTYPE_SNP_FIRST_COMPARATOR);

        // First we filter the assembly variants based on badness from the graph
        for (VariantContext delVariant : pileupAllelesFoundShouldFilter) {

            List<VariantContext> variantsToRemove = assemblyVariants.stream().filter(
                    v -> v.getStart() == delVariant.getStart() &&
                            delVariant.getReference().equals(v.getReference()) &&
                            delVariant.getAlternateAllele(0).equals(v.getAlternateAllele(0))).collect(Collectors.toList());

            if (!variantsToRemove.isEmpty()) {
                if (debugSite) System.out.println("Removing assembly variants due to columnwise heurisits: " + variantsToRemove);
                variantsToRemove.forEach(assemblyVariants::remove); //TODO don't blow up the original assemblys?
            }
        }

        // Ignore any snps from pileups that were close to indels
        final List<VariantContext> givenAllelesFiltered = pileupAllelesPassingFilters.stream()
                .filter(vc -> vc.isIndel() ||
                        assemblyVariants.stream().filter(VariantContext::isIndel).noneMatch(indel -> vc.withinDistanceOf(indel, snpAdjacentToIndelLimit))).collect(Collectors.toList());

        // TODO make sure this TREE-SET is properly comparing the VCs
        variantsInOrder.addAll(assemblyVariants);
        variantsInOrder.addAll(givenAllelesFiltered);

        if (debugSite) System.out.println("Variants to PDHapDetermination:\n"+
                variantsInOrder.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("\n")));

        // TODO this is where we filter out if indels > 32
        // TODO this is where we do eventTestSet work
        List<VariantContext> vcsAsList = new ArrayList<>(variantsInOrder);
        List<List<VariantContext>> dissalowedPairs = new LinkedList<>();

        // NOTE: we iterate over this several times and expect it to be sorted.
        Map<Double, List<VariantContext>> eventsByDRAGENCoordinates = new LinkedHashMap<>();
        SortedMap<Integer, List<VariantContext>> variantsByStartPos = new TreeMap<>();
        List<EventGroup> eventGroups = new ArrayList<>();
        int lastEventEnd = -1;
        //TODO figure out the end for indels correctly
        for (VariantContext vc : variantsInOrder) {
            // Break everything into independent groups (don't worry about transitivitiy right now)
            Double eventKey = vc.getStart() + (vc.isSimpleInsertion()? 0.5:0) + (vc.isSimpleDeletion()? 1 : 0) - referenceHaplotype.getStartPosition();
            eventsByDRAGENCoordinates.putIfAbsent(eventKey, new ArrayList<>());
            eventsByDRAGENCoordinates.get(eventKey).add(vc);
            variantsByStartPos.putIfAbsent(vc.getStart(), new ArrayList<>());
            variantsByStartPos.get(vc.getStart()).add(vc);
            if (debugSite) System.out.println("testing:"+vc);
            if (debugSite) System.out.println("EventKey:"+eventKey);
            if (eventKey <= lastEventEnd + 0.5) {
                eventGroups.get(eventGroups.size()-1).addEvent(vc);
            } else {
                eventGroups.add(new EventGroup(vc));
            }
            int newEnd = (int) (vc.getEnd() - referenceHaplotype.getStartPosition());
            if (debugSite) System.out.println("LastEventEnd:"+lastEventEnd+"    newEventEnd:"+newEnd);
            lastEventEnd = Math.max(newEnd, lastEventEnd);
        }
        //Print the event groups
        if (debugSite) eventsByDRAGENCoordinates.entrySet().stream().map(e -> {
            return String.format("%.1f", e.getKey()) + " -> " + e.getValue().stream()
                    .map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition()))
                    .collect(Collectors.joining(","));
        }).forEach(System.out::println);

        // Iterate over all events starting with all indels
        for (int i = 0; i < vcsAsList.size(); i++) {
            final VariantContext firstEvent = vcsAsList.get(i);
            if (firstEvent.isIndel()) {
                // For every indel, make every 2-3 element subset (without overlapping) of variants to test for equivalency
                for (int j = 0; j < vcsAsList.size(); j++) {
                    final VariantContext secondEvent = vcsAsList.get(j);
                    // Don't compare myslef, any overlappers to myself, or indels i've already examined me (to prevent double counting)
                    if (j != i && !eventsOverlapForPDHapsCode(firstEvent, secondEvent, true) && ((!secondEvent.isIndel()) || j > i)) {
                        final List<VariantContext> events = new ArrayList<>(Arrays.asList(firstEvent, secondEvent));
                        events.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
                        if (debugSite) System.out.println("Testing events: "+ events.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                        if (constructArtificialHaplotypeAndTestEquivalentEvents(referenceHaplotype, aligner, swParameters, variantsInOrder, events, debugSite)) {
                            dissalowedPairs.add(events);
                        } else {

                            // If our 2 element arrays weren't inequivalent, test subsets of 3 including this:
                            for (int k = j+1; k < vcsAsList.size(); k++) {
                                final VariantContext thirdEvent = vcsAsList.get(k);
                                if (k != i && !eventsOverlapForPDHapsCode(thirdEvent,firstEvent,true) && !eventsOverlapForPDHapsCode(thirdEvent,secondEvent,true) ) {
                                    List<VariantContext> subList = new ArrayList<>(events);
                                    subList.add(thirdEvent);
                                    subList.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
                                    if (debugSite) System.out.println("Testing events: "+ subList.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                                    if (constructArtificialHaplotypeAndTestEquivalentEvents(referenceHaplotype, aligner, swParameters, variantsInOrder, subList, debugSite)) {
                                        dissalowedPairs.add(subList);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (debugSite) {
            System.out.println("Dissalowed Variant pairs:");
            dissalowedPairs.stream().map(l -> l.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->"))).forEach(System.out::println);
        }

        if (debugSite) {
            System.out.println("Event groups before merging:\n"+eventGroups.stream().map(eg -> eg.toDisplayString((int)referenceHaplotype.getStartPosition())).collect(Collectors.joining("\n")));
        }
        //Now that we have the dissalowed groups, lets merge any of them from seperate groups:
        //TODO this is not an efficient way of doing this
        for (List<VariantContext> pair : dissalowedPairs) {
            EventGroup eventGrpLeft = null;
            for (VariantContext event : pair) {
                //TODO add exception if its not present for some reason
                EventGroup grpForEvent = eventGroups.stream().filter(grp -> grp.contains(event)).findFirst().get();
                // If the event isn't in the same event group as its predicessor, merge this group with that one and
                if (eventGrpLeft != grpForEvent) {
                    if (eventGrpLeft == null) {
                        eventGrpLeft = grpForEvent;
                    } else {
                        eventGrpLeft.mergeEvent(grpForEvent);
                        eventGroups.remove(grpForEvent);
                    }
                }
            }
        }
        if (debugSite) {
            System.out.println("Event groups after merging:\n"+eventGroups.stream().map(eg -> eg.toDisplayString((int)referenceHaplotype.getStartPosition())).collect(Collectors.joining("\n")));
        }

        //Now we have finished with the work of merging event groups transitively by position and mutually exclusiveness. Now every group should be entirely independant of one antoher:
        if (eventGroups.stream().map(eg -> eg.populateBitset(dissalowedPairs)).anyMatch(b->!b)) {
            // if any of our event groups is too large, abort.
            if (debugSite ) System.out.println("Found event group with too many variants! Aborting haplotype building");
            return sourceSet;
        };

        Set<Haplotype> outputHaplotypes = new HashSet<>();
        outputHaplotypes.add(referenceHaplotype);

        //Iterate over very VCF start position in R space
        List<Map.Entry<Integer, List<VariantContext>>> entriesRInOrder = new ArrayList<>(variantsByStartPos.entrySet());
        /**
         * Overall Loop:
         * Iterate over every cluster of variants with the same start position.
         */
        for (int indexOfDeterminedInR = 0; indexOfDeterminedInR < entriesRInOrder.size(); indexOfDeterminedInR++) {
            Map.Entry<Integer, List<VariantContext>> variantSiteGroup = entriesRInOrder.get(indexOfDeterminedInR);
            if (debugSite) System.out.println("working with variants of the group: " + variantSiteGroup);

            final List<VariantContext> determinedVariants = variantSiteGroup.getValue();

            /**
             * Determined Event Loop:
             * We iterate over every ref position and select single alleles (including ref) from that reference position (in R space) to be "determined"
             *
             * NOTE: we skip the reference allele in the event that we are making determined haplotypes instead of undetermined haplotypes
             */
            for (int IndexOfAllele = (makeDeterminedHapsInstead?0:-1); IndexOfAllele < determinedVariants.size(); IndexOfAllele++) { //note -1 for I here corresponds to the reference allele at this site
                if (debugSite) System.out.println("Working with allele at site: "+(IndexOfAllele ==-1? "[ref:"+(variantSiteGroup.getKey()-referenceHaplotype.getStartPosition())+"]" : PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int)referenceHaplotype.getStartPosition()).apply(determinedVariants.get(IndexOfAllele))));
                // This corresponds to the DRAGEN code for
                // 0 0
                // 0 1
                // 1 0
                final boolean isRef = IndexOfAllele == -1;
                final VariantContext determinedEventToTest = determinedVariants.get(isRef ? 0 : IndexOfAllele);

                /*
                 * Here we handle any of the necessary work to deal with the event groups and maybe forming compund branches out of the groups
                 */
                //Set Determined pairs:
                List<Tuple<VariantContext, Boolean>> determinedPairs = new ArrayList<>();
                for(int j = 0; j < determinedVariants.size(); j++) {
                    determinedPairs.add(new Tuple<>(determinedVariants.get(j), IndexOfAllele == j));
                }

                // Loop over eventGroups, have each of them return a list of VariantContexts
                List<Set<VariantContext>> branchExcludeAlleles = new ArrayList<>();
                branchExcludeAlleles.add(new HashSet<>()); // Add the null branch.

                for(EventGroup group : eventGroups ) {
                    if (group.causesBranching()) {
                        List<List<Tuple<VariantContext, Boolean>>> groupVCs = group.getVariantGroupsForEvent(determinedPairs, true);
                        // Combinatorially expand the branches as necessary
                        List<Set<VariantContext>> newBranchesToAdd = new ArrayList<>();
                        for (Set<VariantContext> excludedVars : branchExcludeAlleles) {
                            //For every exclude group, fork it by each subset we have:
                            for (int i = 1; i < groupVCs.size(); i++) {
                                Set<VariantContext> newSet = new HashSet<>(excludedVars);
                                groupVCs.get(i).stream().filter(t -> !t.b).forEach(t -> newSet.add(t.a));
                                newBranchesToAdd.add(newSet);
                            }
                            // Be careful since this event group might have returned nothing
                            if (!groupVCs.isEmpty()) {
                                groupVCs.get(0).stream().filter(t -> !t.b).forEach(t -> excludedVars.add(t.a));
                            }
                        }
                        branchExcludeAlleles.addAll(newBranchesToAdd);

                        if (branchExcludeAlleles.size() > MAX_BRANCH_PD_HAPS) {
                            if (debugSite ) System.out.println("Found too many branches for variants at: "+determinedEventToTest.getStart()+" aborting and falling back to Assembly Varinats!");
                            return sourceSet;
                        }
                    }
                }

                if (debugSite) {
                    System.out.println("Branches:");
                    for (int i = 0; i < branchExcludeAlleles.size(); i++) {
                        System.out.println("Branch "+i+" excluded VCs:");
                        System.out.println(branchExcludeAlleles.get(i).stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int)referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                    }
                }


                /**
                 * Now handle each branch independently of the others. (the logic is the same in every case except we must ensure that none of the excluded alleles get included when constructing haps.
                 */
                for (Set<VariantContext> excludeEvents : branchExcludeAlleles) {

                    List<Haplotype> branchHaps = new ArrayList<>();
                    // iterate over all of the
                    List<VariantContext> newBranch = new ArrayList<>();

                    // Handle the simple case of making PD haplotypes
                    if (!makeDeterminedHapsInstead) {
                        for (int secondRIndex = 0; secondRIndex < entriesRInOrder.size(); secondRIndex++) {
                            if (secondRIndex != indexOfDeterminedInR) {
                                // We know here that nothing illegally overlaps because there are no groups.
                                // Also exclude anything overlapping myself to protect
                                entriesRInOrder.get(secondRIndex).getValue().stream().filter(vc -> !excludeEvents.contains(vc)).filter(vc -> isRef && !vc.overlaps(determinedEventToTest)).forEach(newBranch::add);
                            } else {
                                newBranch.add(determinedEventToTest);
                            }
                        }
                        newBranch.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
                        PartiallyDeterminedHaplotype newPDHaplotypeFromEvents = createNewPDHaplotypeFromEvents(referenceHaplotype, determinedEventToTest, isRef, newBranch);
                        branchHaps.add(newPDHaplotypeFromEvents);

                    } else {
                        // TODO currently this approach doesn't properly handle a bunch of duplicate events...
                        // If we are producing determined bases, then we want to enforce that every new event at least has THIS event as a variant.
                        List<List<VariantContext>> variantGroupsCombinatorialExpansion = new ArrayList<>();
                        variantGroupsCombinatorialExpansion.add(new ArrayList<>());
                        // TODO, i'm reasonably sure its valid to say that indexOfDeterminedInR is the FIRST event in the combinatorial expansion, everything before that has been processed as true
                        for (int secondRIndex = indexOfDeterminedInR; secondRIndex < entriesRInOrder.size(); secondRIndex++) {
                            if (variantGroupsCombinatorialExpansion.size() > MAX_BRANCH_PD_HAPS) {
                                if(debugSite ) System.out.println("Too many branch haplotypes ["+variantGroupsCombinatorialExpansion.size()+"] generated from site, falling back on assebmly variants!");
                                return sourceSet;
                            }
                            // Iterate through the growing combinatorial expansion of haps, split it into either having or not having the variant.
                            if (secondRIndex == indexOfDeterminedInR) {
                                for (List<VariantContext> hclist : variantGroupsCombinatorialExpansion) {
                                    hclist.add(determinedEventToTest);
                                }
                            // Othewise make sure to include the combinatorial expansion of events at the other site
                            } else {
                                List<List<VariantContext>> hapsPerVCsAtRSite = new ArrayList<>();
                                for (VariantContext vc : entriesRInOrder.get(secondRIndex).getValue()) {
                                    for (List<VariantContext> hclist : variantGroupsCombinatorialExpansion) {
                                        if (!excludeEvents.contains(vc)) {
                                            List<VariantContext> newList = new ArrayList<>(hclist);
                                            newList.add(vc);
                                            hapsPerVCsAtRSite.add(newList);
                                        }
                                    }
                                }
                                //Add them after to prevent accidentally adding duplicates of events at a site
                                variantGroupsCombinatorialExpansion.addAll(hapsPerVCsAtRSite);
                            }
                        }

                        for (List<VariantContext> subset : variantGroupsCombinatorialExpansion) {
                            subset.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
                            if (debugSite) System.out.println("Construcing Hap From Events:"+ subset.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                            branchHaps.add(constructHaplotypeFromVariants(referenceHaplotype, subset, true));
                        }
                    }
                    // Add the branch haps to the results:
                    if (debugSite) {
                        System.out.println("Constructed Haps for Branch"+excludeEvents.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining(",")) + ":");
                        System.out.println(branchHaps.stream().map(h -> h.getCigar() + " " + h.toString()).collect(Collectors.joining("\n")));
                    }

                    outputHaplotypes.addAll(branchHaps);
                    if (outputHaplotypes.size() > MAX_PD_HAPS_TO_GENERATE) {
                        if (debugSite  ) System.out.println("Too many Haps ["+outputHaplotypes.size()+"] generated at this site! Aborting!");
                        return sourceSet;
                    }
                }
            }
        }

        if (outputHaplotypes.size() > MAX_PD_HAPS_TO_GENERATE) {
            if (debugSite ) System.out.println("Too many branch haplotypes found, aborting ["+outputHaplotypes.size()+"]");
            return sourceSet;
        }
        // TODO be careful here: The old haps still exist... WE must be careful to use the PD haps correctly
        sourceSet.storeAssemblyHaplotypes();
        sourceSet.replaceAllHaplotypes(outputHaplotypes);
        if (!makeDeterminedHapsInstead) {
            // Setting a boolean on the source-set to indicate to downstream code that we have PD haplotypes
            sourceSet.setPartiallyDeterminedMode();
        }
        if (debugSite ) System.out.println("Returning "+outputHaplotypes.size()+" to the HMM");
        return sourceSet;
    }

    /**
     * Overlaps method to handle indels and snps correctly
     *
     * @param snpsOverlap if false, don't ever evaluate snps as overlapping other snps (we do this because sometimes we need to construct artifical haps where we don't allow overlapping)
     */
    static boolean eventsOverlapForPDHapsCode(VariantContext vc1, VariantContext vc2, boolean snpsOverlap){
        if (!snpsOverlap && vc2.isSNP() && vc1.isSNP()) {
            return false;
        }
        if (!vc1.getContig().equals(vc2.getContig())) {
            return false;
        }
        double vc1start = vc1.isIndel() ? (vc1.isSimpleDeletion() ? vc1.getStart() + 1 : vc1.getStart() + 0.5) : vc1.getStart();
        double vc1end = vc1.isSimpleInsertion() ? vc1.getEnd() + 0.5 : vc1.getEnd();
        double vc2start = vc2.isIndel() ? (vc2.isSimpleDeletion() ? vc2.getStart() + 1 : vc2.getStart() + 0.5) : vc2.getStart();
        double vc2end = vc2.isSimpleInsertion() ? vc2.getEnd() + 0.5 : vc2.getEnd();

        //Pulled directly from CoordMath.java (not using here because of doubles)
        return (vc2start >= vc1start && vc2start <= vc1end) || (vc2end >=vc1start && vc2end <= vc1end) || vc1start >= vc2start && vc1end <= vc2end;
    }


    /**
     * This method is the helper that manages the actual SW alignment and testing of a group of variants vs the reference haplotype.
     *
     * The method is as follows, construct and artificial haplotype of the provided events, then realign it vs the reference and test
     * if any of the resulting variants are present in the inputs (but doesn't match)
     */
    @VisibleForTesting
    private static boolean constructArtificialHaplotypeAndTestEquivalentEvents(Haplotype referenceHaplotype, SmithWatermanAligner aligner, SWParameters swParameters, TreeSet<VariantContext> vcs, List<VariantContext> eventsToTest, boolean debugSite) {
        final Haplotype realignHap = constructHaplotypeFromVariants(referenceHaplotype, eventsToTest, false);
        //ALIGN!
        realignHap.setCigar(CigarUtils.calculateCigar(referenceHaplotype.getBases(), realignHap.getBases(), aligner, swParameters, SWOverhangStrategy.INDEL));
        EventMap.buildEventMapsForHaplotypes(Collections.singletonList(realignHap), referenceHaplotype.getBases(), referenceHaplotype.getGenomeLocation(), false,0);
        //TODO this isn't exactly how it works in DRAGEN, or at least we've lost the score here...
        final boolean wasEquivalentEvent = realignHap.getEventMap().getVariantContexts().stream().filter(eMapVC ->
                // Are there any variants NOT in our initial list
                eventsToTest.stream().noneMatch(v -> {
                    return doVariantsMatch(eMapVC, v);
                }))
                // Do any of those variants appear in our overall list of alleles
                .anyMatch(eMapVc -> vcs.stream().anyMatch(v -> {
                    return doVariantsMatch(eMapVc, v);
                }));
        if (debugSite) System.out.println(
                realignHap.getEventMap().getVariantContexts().stream()
                .map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition()))
                .collect(Collectors.joining("\n")));
        if (wasEquivalentEvent) {
            if (debugSite) System.out.println("Events mismatched!");
        }

        return wasEquivalentEvent;
    }

    // A helper method to assert that variant contexts in the event map match those outside of it.
    private static boolean doVariantsMatch(VariantContext eMapVC, VariantContext v) {
        return eMapVC.getStart() == v.getStart() &&
                eMapVC.getReference().equals(v.getReference()) &&
                eMapVC.getAlternateAllele(0).equals(v.getAlternateAllele(0)) &&
                eMapVC.getAlternateAlleles().size() == v.getAlternateAlleles().size();
    }

    //Helper method that does the work of computing what braches need to exist for later PDHap construction
    //TODO this will EVENTUALLY take as input the results fro smith waterman
    //TODO this approach will not eventually work whne we need to split groups
    private static List<List<VariantContext>> getBranchedHaplotypeGroups(final Map<Double, List<VariantContext>> eventGroups) {
        return recursiveHapGroupHelper(new ArrayList<>(eventGroups.entrySet()), -1, new ArrayList<>());
    }


    //TODO this is going to be totally and completely redone...
    private static List<List<VariantContext>> recursiveHapGroupHelper(final List<Map.Entry<Double, List<VariantContext>>> listOfGroups,
                                                                      final int positionOfLastAddedElement,
                                                                      final List<VariantContext> continueGroup) {
        List<List<VariantContext>> output = new ArrayList<>();
        for (int i = positionOfLastAddedElement + 1; i < listOfGroups.size(); i++) {
//            int endBase = continueGroup.isEmpty()? continueGroup.get(continueGroup.size()-1).getEnd() : -1;

            List<VariantContext> currentMutexGroup = listOfGroups.get(i).getValue();
            if (currentMutexGroup.size() > 1) {
                for (int j = 1; j < currentMutexGroup.size(); j++) {
                    List<VariantContext> subList = new ArrayList<>(continueGroup);
                    subList.add(currentMutexGroup.get(j));
                    output.addAll(recursiveHapGroupHelper(listOfGroups, i, subList));
                }
            }
            continueGroup.add(currentMutexGroup.get(0));
        }
        // If we reach here then we have constructed
        output.add(continueGroup);
        return output;
    }

    /**
     *
     *
     * NOTE: this accepts multiple alleles stacked up at the same base (assuming the order is SNP -> INDEL)
     * @param refHap
     * @param variantContexts
     * @return
     */
    @VisibleForTesting
    public static Haplotype constructHaplotypeFromVariants(final Haplotype refHap, final List<VariantContext> variantContexts, final boolean setEventMap) {
        //ASSERT that the base is ref and cool
        if (!refHap.isReference() || refHap.getCigar().numCigarElements() > 1) {
            throw new RuntimeException("This is not a valid base haplotype for construction");
        }
        final long genomicStartPosition = refHap.getStartPosition();
        long refOffsetOfNextBaseToAdd = genomicStartPosition;

        byte[] refbases = refHap.getBases();
        CigarBuilder runningCigar = new CigarBuilder();
        byte[] newRefBases = {};

        //ASSUME sorted for now
        // use the reverse list to save myself figuring out cigars for right now
        for (VariantContext vc : variantContexts) {
            if (vc.getAlternateAlleles().size() > 1) {
                throw new RuntimeException("too may alt alleles");
            }
            Allele refAllele = vc.getReference();
            Allele altAllele = vc.getAlternateAllele(0);

            int intermediateRefStartPosition = (int) (refOffsetOfNextBaseToAdd - genomicStartPosition);
            int intermediateRefEndPosition = Math.toIntExact(vc.getStart() - genomicStartPosition);

            if ((vc.isIndel() && intermediateRefEndPosition - intermediateRefStartPosition < -1) || (!vc.isIndel() && intermediateRefEndPosition - intermediateRefStartPosition < 0)) {//todo clean this up
                throw new RuntimeException("Variant "+vc+" is out of order in the PD event list: "+variantContexts);
            }
            if (intermediateRefEndPosition - intermediateRefStartPosition > 0) { // Append the cigar element for the anchor base if necessary.
                runningCigar.add(new CigarElement(intermediateRefEndPosition - intermediateRefStartPosition, CigarOperator.M));
            }
            // Include the ref base for indel if the base immediately proceeding this event is not already tracked
            boolean includeRefBaseForIndel = vc.isIndel() && (intermediateRefStartPosition <= intermediateRefEndPosition);

            CigarElement newCigarElement;
            if (refAllele.length() == altAllele.length()) {
                newCigarElement = new CigarElement(refAllele.length(), CigarOperator.X);
            } else {
                // If the last base was filled by another event, don't attempt to fill in the indel ref base.
                if (includeRefBaseForIndel) {
                    runningCigar.add(new CigarElement(1, CigarOperator.M)); //When we add an indel we end up inserting a matching base
                }
                newCigarElement = new CigarElement(Math.abs(altAllele.length() - refAllele.length()),
                        refAllele.length() > altAllele.length() ?
                                CigarOperator.D : CigarOperator.I);
            }
            runningCigar.add(newCigarElement);

            if (intermediateRefEndPosition - intermediateRefStartPosition > 0) {
                newRefBases = ArrayUtils.addAll(newRefBases, ArrayUtils.subarray(refbases, intermediateRefStartPosition, (int) (vc.getStart() - genomicStartPosition))); // bases before the variant
            }
            // Handle the ref base for indels that exlcude their ref bases
            if (refAllele.length() != altAllele.length() && !includeRefBaseForIndel) {
                newRefBases = ArrayUtils.addAll(newRefBases, Arrays.copyOfRange(altAllele.getBases(),1, altAllele.length()));
            // else add the snp
            } else {
                newRefBases = ArrayUtils.addAll(newRefBases, altAllele.getBases()); // refbases added
            }
            refOffsetOfNextBaseToAdd = vc.getEnd() + 1; //TODO this is probably not set for future reference
        }

        // Finish off the haplotype with the final bases
        int refStartIndex = (int) (refOffsetOfNextBaseToAdd - genomicStartPosition);
        newRefBases = ArrayUtils.addAll(newRefBases, ArrayUtils.subarray(refbases, refStartIndex, refbases.length));
        runningCigar.add(new CigarElement(refbases.length - refStartIndex, CigarOperator.M));

        final Haplotype outHaplotype = new Haplotype(newRefBases, false, refHap.getGenomeLocation(), runningCigar.make());
        if (setEventMap) {
            EventMap.buildEventMapsForHaplotypes(Collections.singletonList(outHaplotype), refHap.getBases(), refHap.getGenomeLocation(), false,0);
            // NOTE: we set this AFTER generating the event maps because hte event map code above is being generated from the ref hap so this offset will cause out of bounds errors
            outHaplotype.setAlignmentStartHapwrtRef(refHap.getAlignmentStartHapwrtRef()); //TODO better logic here
        }
        return outHaplotype;
    }

    /**
     * Construct a PD haplotype from scratch
     *
     *
     * Generally we are constructing a new haplotype with all of the reference bases for SNP events and with the longest possible allele for INDEL events.
     * For deleteions, we extend the haplotype by the ref lenght
     *
     * NOTE: we assume each provided VC
     *
     * //TODO these both need guardrails against exceptions
     *
     * @param base
     * @param eventWithVariant
     * @param useRef
     * @param constituentEvents
     */
    @VisibleForTesting
    static PartiallyDeterminedHaplotype createNewPDHaplotypeFromEvents(Haplotype base, VariantContext eventWithVariant, boolean useRef, List<VariantContext> constituentEvents) {

        //ASSERT that the base is ref and cool
        if (!base.isReference() || base.getCigar().numCigarElements() > 1) {
            throw new RuntimeException("This is not a valid base haplotype for construction");
        }
        long genomicStartPosition = base.getStartPosition();
        long refOffsetOfNextBaseToAdd = genomicStartPosition;

        byte[] refBasesToAddTo = base.getBases();
        CigarBuilder runningCigar = new CigarBuilder();
        byte[] newHaplotypeBasees = {};
        byte[] pdBytes = {};

        //ASSUME sorted for now
        // use the reverse list to save myself figuring out cigars for right now
        for (VariantContext vc : constituentEvents) {
            int intermediateRefStartPosition = (int) (refOffsetOfNextBaseToAdd - genomicStartPosition);
            int intermediateRefEndPosition = Math.toIntExact(vc.getStart() - genomicStartPosition);

            // An extra special case if we are a SNP following a SNP
            if (vc.isSNP() && intermediateRefEndPosition - intermediateRefStartPosition == -1 && ((pdBytes[pdBytes.length-1] & PartiallyDeterminedHaplotype.SNP) != 0)  ) {
                byte[] array = PartiallyDeterminedHaplotype.getPDBytesForHaplotypes(vc.getReference(), vc.getAlternateAllele(0));
                pdBytes[pdBytes.length-1] = (byte) (pdBytes[pdBytes.length-1] | array[0]); // adding any partial bases if necessary
                continue;
            }

                //Check that we are allowed to add this event (and if not we are)
            if ((vc.isIndel() && intermediateRefEndPosition - intermediateRefStartPosition < -1) || (!vc.isIndel() && intermediateRefEndPosition - intermediateRefStartPosition < 0)) {//todo clean this up
                throw new RuntimeException("Variant "+vc+" is out of order in the PD event list: "+constituentEvents);
            }

            // Add the cigar for bases we skip over
            if (intermediateRefEndPosition - intermediateRefStartPosition > 0) {
                runningCigar.add(new CigarElement(intermediateRefEndPosition - intermediateRefStartPosition, CigarOperator.M));
            }

            // Include the ref base for indel if the base immediately proceeding this event is not already tracked
            boolean includeRefBaseForIndel = vc.isIndel() && (intermediateRefStartPosition <= intermediateRefEndPosition);

            // Determine the alleles to add
            Allele refAllele = vc.getReference();
            Allele altAllele = vc.getAlternateAllele(0);
            boolean shouldFlip = altAllele.length() > refAllele.length();
            boolean isEvent = false;
            byte[] basesToAdd;
            // If this is the blessed variant, add
            if (vc.getStart()==eventWithVariant.getStart()) {
                if (!useRef && eventWithVariant.getAlternateAlleles().size() > 1) {
                    throw new RuntimeException("the Blessed variant must be monoallelic");
                }
                isEvent = true;
                basesToAdd = useRef? refAllele.getBases() : altAllele.getBases();
            // Otherwise make sure we are adding the longest allele (for indels) or the ref allele for snps.
            } else {
                basesToAdd = refAllele.length() >= altAllele.length() ? refAllele.getBases() : altAllele.getBases();
            }

            // We want to always add the longest allele
            if (vc.isIndel() && !includeRefBaseForIndel) {
                basesToAdd = Arrays.copyOfRange(basesToAdd, 1, basesToAdd.length);
            }


            // Figure out the cigar to add:
            // - If we are in the ref, simply add the cigar corresponding to the allele we are using
            // -
            CigarElement newCigarElement;
            // if this is the event
            if (isEvent && useRef) {
                if (vc.isIndel() && !includeRefBaseForIndel) {
                    newCigarElement = new CigarElement(refAllele.length() - 1 , CigarOperator.M);
                } else {
                    newCigarElement = new CigarElement(refAllele.length(), CigarOperator.M);
                }
           // If we aren't in the blessed variant, add a match and make sure the array is set accordingly
            } else {
                if (refAllele.length()==altAllele.length()) {
                    newCigarElement = new CigarElement(refAllele.length() , CigarOperator.M);
                } else {
                    if (includeRefBaseForIndel) {
                        runningCigar.add(new CigarElement(1, CigarOperator.M));
                    }
                    if (shouldFlip) {
                        //When we add an indel we end up inserting a matching base
                        newCigarElement = new CigarElement(altAllele.length() - refAllele.length(), CigarOperator.I);
                    } else {
                        newCigarElement = new CigarElement(basesToAdd.length - 1, CigarOperator.M);
                    }
                }
            }
            runningCigar.add(newCigarElement);

            // Add ref basses up to this if necessary
            if (intermediateRefEndPosition - intermediateRefStartPosition > 0) {
                newHaplotypeBasees = ArrayUtils.addAll(newHaplotypeBasees, ArrayUtils.subarray(refBasesToAddTo, intermediateRefStartPosition, (int) (vc.getStart() - genomicStartPosition))); // bases before the variant
                pdBytes = ArrayUtils.addAll(pdBytes, new byte[vc.getStart() - (int)refOffsetOfNextBaseToAdd]); // bases before the variant
            }
            newHaplotypeBasees = ArrayUtils.addAll(newHaplotypeBasees, basesToAdd); // refbases added
            if (includeRefBaseForIndel) {
                pdBytes = ArrayUtils.add(pdBytes, (byte)0);
            }
            pdBytes = ArrayUtils.addAll(pdBytes, isEvent?
                    new byte[basesToAdd.length - (includeRefBaseForIndel?1:0)] :
                    PartiallyDeterminedHaplotype.getPDBytesForHaplotypes(shouldFlip?
                            altAllele :
                            refAllele,
                            shouldFlip? refAllele :
                                    altAllele)); // refbases added
            refOffsetOfNextBaseToAdd = vc.getEnd() + 1; //TODO this is probably not set for future reference
        }

        // Finish off the haplotype with the final bases
        int refStartIndex = (int) (refOffsetOfNextBaseToAdd - genomicStartPosition);
        newHaplotypeBasees = ArrayUtils.addAll(newHaplotypeBasees, ArrayUtils.subarray(refBasesToAddTo, refStartIndex, refBasesToAddTo.length));
        pdBytes = ArrayUtils.addAll(pdBytes, new byte[refBasesToAddTo.length - refStartIndex]);
        runningCigar.add(new CigarElement(refBasesToAddTo.length - refStartIndex, CigarOperator.M));

        return new PartiallyDeterminedHaplotype(
                new Haplotype(newHaplotypeBasees, false, base.getGenomeLocation(), runningCigar.make()),//TODO the GenomeLoc is now probably incorrect...
                false,
                pdBytes,
                constituentEvents,
                eventWithVariant,
                runningCigar.make());
    }

    // A helper class for managing mutually exclusive event clusters and the logic arround forming valid events vs eachother.
    private static class EventGroup {
        List<VariantContext> variantsInBitmapOrder;
        HashSet<VariantContext> variantContextSet; //TODO eventually map?
        //From Illumina (there is a LOT of math that will eventually go into these)/
        BitSet allowedEvents = null;

        // Optimizaiton to save ourselves recomputing the subsets at every point its necessary to do so.
        List<List<Tuple<VariantContext,Boolean>>> cachedEventLitsts = null;

        public EventGroup(final VariantContext variantContext) {
            variantsInBitmapOrder = new ArrayList<>();
            variantContextSet = new HashSet<>();
            variantsInBitmapOrder.add(variantContext);
            variantContextSet.add(variantContext);
        }
        public EventGroup() {
            variantsInBitmapOrder = new ArrayList<>();
            variantContextSet = new HashSet<>();
        }

        /**
         * This is the primary method for handling mutually exclusive events in this subgroup. This code amd methods comes directly from DRAGEN:
         *
         * Create a #Variants bitset to store valid pairings:
         *      - The index of each element corresponds to an enumerated subset of alleles in this group
         *      - Each bit in the index corresponds to the presence or absence of a variant from the vcf list.
         *          - For example with variants [A,B,C] the number 5 corresponds to subset [A,C]
         *      - A false in the bitset corresponds to a dissalowed pair.
         *      - NOTE: we can use 32bit ints for these bitshift operations by virtue of the fact that we limit ourselves to at most 22 variants per group.
         * Iterate through pairs of Variants that overlap and mark off any pairings including this.
         * Iterate through the mutex variants and ensure pairs containing all mutex variant groups are marked as true
         *
         * @param dissalowedEvents Pairs of events dissalowed
         * @return false if the event group is too large to process
         */
        public boolean populateBitset(List<List<VariantContext>> dissalowedEvents) {
            if (variantsInBitmapOrder.size() > MAX_VAR_IN_EVENT_GROUP) {
                return false;
            }
            if (variantsInBitmapOrder.size() < 2) {
                return true;
            }

            allowedEvents = new BitSet(variantsInBitmapOrder.size());
            allowedEvents.flip(1, 1 << variantsInBitmapOrder.size());
            // initialize all events as being allowed and then dissalow them in turn .

            // Ensure the list is in positional order before commencing.
            variantsInBitmapOrder.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
            List<Integer> bitmasks = new ArrayList<>();
            // Mark as disallowed all events that overlap eachother
            for (int i = 0; i < variantsInBitmapOrder.size(); i++) {
                for (int j = i+1; j < variantsInBitmapOrder.size(); j++) {
                    if (eventsOverlapForPDHapsCode(variantsInBitmapOrder.get(i), variantsInBitmapOrder.get(j), false)) {
                        bitmasks.add(1 << i | 1 << j);
                    }
                }
            }
            // mark as dissalowed any sets of variants from the bitmask.
            for (List<VariantContext> disallowed : dissalowedEvents) {
                //
                if (disallowed.stream().anyMatch(v -> variantContextSet.contains(v))){
                    int bitmask = 0;
                    for (VariantContext v : disallowed) {
                        int indexOfV = variantsInBitmapOrder.indexOf(v);
                        if (indexOfV < 0) {
                            throw new RuntimeException("Something went wrong in event group merging, variant "+v+" is missing from the event group despite being in a mutex pair: "+disallowed+"\n"+this);
                        }
                        bitmask += 1 << variantsInBitmapOrder.indexOf(v);
                    }
                    bitmasks.add(bitmask);
                }
            }

            // Now iterate through the list and dissalow all events with every bitmask
            if (!bitmasks.isEmpty()) {
                events:
                for (int i = 1; i < allowedEvents.length(); i++) {
                    for (final int mask : bitmasks) {
                        if ((i & mask) == mask) { // are the bits form the mask true?
                            allowedEvents.set(i, false);
                            continue events;
                            // Once i is set we don't need to keep checking bitmasks
                        }
                    }
                }
            }

            return true;
        }

        /**
         * This method handles the logic involved in getting all of the allowed subsets of alleles for this event group.
         *
         * @param dissalowSubsets
         * @return
         */
        public List<List<Tuple<VariantContext,Boolean>>> getVariantGroupsForEvent(final List<Tuple<VariantContext, Boolean>> eventsForMask, final boolean dissalowSubsets) {
            // If we are dealing with an external to this list event
            int eventMask = 0;
            int maskValues = 0;
            for(Tuple<VariantContext, Boolean> event : eventsForMask) {
                if (variantContextSet.contains(event.a)) {
                    int index = variantsInBitmapOrder.indexOf(event.a);
                    eventMask = eventMask | (1 << index);
                    maskValues = maskValues | ((event.b ? 1 : 0) << index);
                }
            }
            // Special case (if we are determining bases outside of this mutex cluster we can reuse the work from previous iterations)
            if (eventMask == 0 && cachedEventLitsts != null) {
                return cachedEventLitsts;
            }

            List<Integer> ints = new ArrayList<>();
            // Iterate from the BACK of the list (i.e. ~supersets -> subsets)
            // NOTE: we skip over 0 here since that corresponds to ref-only events, handle those externally to this code
            outerLoop:
            for (int i = allowedEvents.length(); i > 0; i--) {
                // If the event is allowed AND if we are looking for a particular event to be present or absent.
                if (allowedEvents.get(i) && (eventMask == 0 || ((i & eventMask) == maskValues))) {
                    // Only check for subsets if we need to
                    if (dissalowSubsets) {
                        for (Integer group : ints) {
                            // if the current i is a subset of an existing group
                            if ((i | group) == group) {
                                continue outerLoop;
                            }
                        }
                    }
                    ints.add(i);
                }
            }

            // Now that we have all the mutex groups, unpack them into lists of variants
            List<List<Tuple<VariantContext,Boolean>>> output = new ArrayList<>();
            for (Integer grp : ints) {
                List<Tuple<VariantContext,Boolean>> newGrp = new ArrayList<>();
                for (int i = 0; i < variantsInBitmapOrder.size(); i++) {
                    // if the corresponding bit is 1, set it as such, otherwise set it as 0.
                    newGrp.add(new Tuple<>(variantsInBitmapOrder.get(i), ((1<<i) & grp) != 0));
                }
                output.add(newGrp);
            }
            // Cache the result
            if(eventMask==0) {
                cachedEventLitsts = Collections.unmodifiableList(output);
            }
            return output;
        }

        public boolean causesBranching() {
            return variantsInBitmapOrder.size() > 1;
        }

        //Print The event group in Illumina indexed ordering:
        public String toDisplayString(int startPos) {
            return "EventGroup: " + variantsInBitmapOrder.stream().map(vc -> PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString(startPos).apply(vc)).collect(Collectors.joining("->"));
        }

        public boolean contains(final VariantContext event) {
            return variantContextSet.contains(event);
        }

        public void addEvent(final VariantContext event) {
            variantsInBitmapOrder.add(event);
            variantContextSet.add(event);
            allowedEvents = null;
        }

        public EventGroup mergeEvent(final EventGroup other) {
            variantsInBitmapOrder.addAll(other.variantsInBitmapOrder);
            variantContextSet.addAll(other.variantsInBitmapOrder);
            allowedEvents = null;
            return this;
        }
    }
}
