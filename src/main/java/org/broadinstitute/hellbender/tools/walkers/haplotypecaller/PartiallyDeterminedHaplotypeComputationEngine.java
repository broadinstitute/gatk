package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
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
import java.util.stream.Collectors;

/**
 * TODO
 */
public class PartiallyDeterminedHaplotypeComputationEngine {
    final static int MAX_PD_HAPS_TO_GENERATE = 256; //2048 TODO probably not...
    final static int MAX_BRANCH_PD_HAPS = 256; //128
    final static int MAX_VAR_IN_BRANCH = 8; //128

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

        final TreeSet<VariantContext> vcs = new TreeSet<>(
                AssemblyResultSet.HAPLOTYPE_VARIANT_CONTEXT_COMPARATOR);

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
        vcs.addAll(assemblyVariants);
        vcs.addAll(givenAllelesFiltered);

        if (debugSite) System.out.println("Variants to PDHapDetermination:\n"+
                vcs.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("\n")));

        // TODO this is where we filter out if indels > 32
        // TODO this is where we do eventTestSet work
        List<VariantContext> vcsAsList = new ArrayList<>(vcs);
        List<List<VariantContext>> dissalowedPairs = new LinkedList<>();

        // NOTE: we iterate over this several times and expect it to be sorted.
//        Map<Double, List<VariantContext>> eventGroups = new LinkedHashMap<>();
//        for (VariantContext vc : vcs) {
//            // Break everything into independent groups (don't worry about transitivitiy right now)
//            Double eventKey = vc.getStart() + (vc.isSimpleInsertion()? 0.5:0) + (vc.isSimpleDeletion()? 1 : 0) - referenceHaplotype.getStartPosition();
//            eventGroups.putIfAbsent(eventKey, new ArrayList<>());
//            eventGroups.get(eventKey).add(vc);
//        }
        // TODO lets break everything into event sets
        Map<Double, List<VariantContext>> eventGroups = new LinkedHashMap<>();
        int lastEventEnd = -1;
        Double lastEventStart = null;

        for (VariantContext vc : vcs) {
            // Break everything into independent groups (don't worry about transitivitiy right now)
            if (vc.getStart() > lastEventEnd) {
                Double key = new Double(vc.getStart());
                eventGroups.put(key, new ArrayList<VariantContext>(Collections.singletonList(vc)));
                lastEventStart = key;
                lastEventEnd = vc.getEnd();//TODO check this

            } else {
                eventGroups.get(lastEventStart).add(vc);
                lastEventEnd = Math.max(lastEventEnd, vc.getEnd());
            }
        }
        //Print the event groups
        if (debugSite) eventGroups.entrySet().stream().map(e -> {
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
                    if (j != i && !secondEvent.overlaps(firstEvent) && ((!secondEvent.isIndel()) || j > i)) {
                        final List<VariantContext> events = new ArrayList<>(Arrays.asList(firstEvent, secondEvent));
                        events.sort(AssemblyResultSet.HAPLOTYPE_VARIANT_CONTEXT_COMPARATOR);
                        if (debugSite) System.out.println("Testing events: "+ events.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                        if (constructArtificialHaplotypeAndTestEquivalentEvents(referenceHaplotype, aligner, swParameters, vcs, events, debugSite)) {
                            dissalowedPairs.add(events);
                        } else {

                            // If our 2 element arrays weren't inequivalent, test subsets of 3 including this:
                            for (int k = j+1; k < vcsAsList.size(); k++) {
                                final VariantContext thirdEvent = vcsAsList.get(k);
                                if (k != i && !thirdEvent.overlaps(firstEvent) && !thirdEvent.overlaps(secondEvent) ) {
                                    List<VariantContext> subList = new ArrayList<>(events);
                                    subList.add(thirdEvent);
                                    subList.sort(AssemblyResultSet.HAPLOTYPE_VARIANT_CONTEXT_COMPARATOR);
                                    if (debugSite) System.out.println("Testing events: "+ subList.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                                    if (constructArtificialHaplotypeAndTestEquivalentEvents(referenceHaplotype, aligner, swParameters, vcs, subList, debugSite)) {
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
        // now iterate over the linked hash map...
        //TODO get branches of compatable variants to treat in tandem


        //TODO clean this up to be more more corrects
        final List<List<VariantContext>> branchedHaplotypeGroups = getBranchedHaplotypeGroups(eventGroups);

        if (debugSite) System.out.println("Branches of PD Haps to construct:\n"+
                branchedHaplotypeGroups.stream().map(variantContexts ->
                variantContexts.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")))
                .collect(Collectors.joining("\n")));

        Set<Haplotype> outputHaplotypes = new HashSet<>();
        outputHaplotypes.add(referenceHaplotype);

        if (!makeDeterminedHapsInstead) {
            //Generate the PDHaplotypes based on each of the branches
            for (List<VariantContext> branch : branchedHaplotypeGroups) {
                if (debugSite) System.out.println("Handling Branch \n" + branch.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                List<PartiallyDeterminedHaplotype> branchHaps = new ArrayList<>();
                for (VariantContext variantContext : branch) {
                    // Add both the variant and the event
                    branchHaps.add(createNewPDHaplotypeFromEvents(referenceHaplotype, variantContext, true, branch));
                    branchHaps.add(createNewPDHaplotypeFromEvents(referenceHaplotype, variantContext, false, branch));
                }
                if (debugSite) System.out.println("Constructed PD Haps:" + branchHaps.stream().map(PartiallyDeterminedHaplotype::toString).collect(Collectors.joining("\n")));
                //CHECK for overwhelming branch haps
                if (branchHaps.size()>= MAX_BRANCH_PD_HAPS) {
                    if (debugSite) System.out.println("Too many branch haplotypes ["+branchHaps.size()+"] generated from branch: "+branch.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                    return sourceSet;
                }
                outputHaplotypes.addAll(branchHaps);
            }
        } else {

            //BIG TODO unexpand the haps...
            //GENERATE THE FULL EXPANSION OF HAPS!
            for (List<VariantContext> branch : branchedHaplotypeGroups) {
                if (debugSite) System.out.println("Handling Branch \n" + branch.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                List<List<VariantContext>> combinatorialExpansionHaps = new ArrayList<>();
                combinatorialExpansionHaps.add(new ArrayList<>());
                if (branch.size() >= MAX_VAR_IN_BRANCH) {
                    if (debugSite) System.out.println("Too many variants ["+branch.size()+"] generated from branch: "+branch.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                    return sourceSet;
                }
                for (VariantContext vc : branch) {
                    List<List<VariantContext>> hapsPerVC = new ArrayList<>();
                    for (List<VariantContext> hclist : combinatorialExpansionHaps) {
                        List<VariantContext> newList = new ArrayList<>(hclist);
                        newList.add(vc);
                        hapsPerVC.add(newList);
                    }
                    combinatorialExpansionHaps.addAll(hapsPerVC);
                }
                if (combinatorialExpansionHaps.size() >= MAX_BRANCH_PD_HAPS) {
                    if (debugSite) System.out.println("Too many branch haplotypes ["+combinatorialExpansionHaps.size()+"] generated from branch: "+branch.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                    return sourceSet;
                }
                List<Haplotype> branchHaps = new ArrayList<>();
                for (List<VariantContext> subset : combinatorialExpansionHaps) {
                    branchHaps.add(constructHaplotypeFromVariants(referenceHaplotype, subset, true));
                }
                if (debugSite) System.out.println("Constructed PD Haps:" + branchHaps.stream().map(Haplotype::toString).collect(Collectors.joining("\n")));
                //CHECK for overwhelming branch haps

                outputHaplotypes.addAll(branchHaps);
            }
        }

        if (outputHaplotypes.size() > MAX_PD_HAPS_TO_GENERATE) {
            if (debugSite) System.out.println("Too many branch haplotypes found, aborting ["+outputHaplotypes.size()+"]");
            return sourceSet;
        }
        sourceSet.replaceAllHaplotypes(outputHaplotypes);

        return sourceSet;
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

    // A helper method
    private static boolean doVariantsMatch(VariantContext eMapVC, VariantContext v) {
        return eMapVC.getStart() == v.getStart() &&
                eMapVC.getReference() == v.getReference() &&
                eMapVC.getAlternateAllele(0) == v.getAlternateAllele(0) &&
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
        long genomicStartPosition = refHap.getStartPosition();
        long refOffset = genomicStartPosition;

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

            int refStartIndex = (int) (refOffset - genomicStartPosition);
            int refEndIndex = Math.toIntExact(vc.getStart() - genomicStartPosition);
            if (refEndIndex - refStartIndex > 0) {
                runningCigar.add(new CigarElement(refEndIndex - refStartIndex, CigarOperator.M));
            }
            // Include the ref base for indel if the base immediately proceeding this event is not already tracked
            boolean includeRefBaseForIndel = vc.isIndel() && (refEndIndex - refStartIndex > 0 || refOffset == genomicStartPosition);

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

            if (refEndIndex - refStartIndex > 0) {
                newRefBases = ArrayUtils.addAll(newRefBases, ArrayUtils.subarray(refbases, refStartIndex, (int) (vc.getStart() - genomicStartPosition))); // bases before the variant
            }
            // Handle the ref base for indels that exlcude their ref bases
            if (refAllele.length() != altAllele.length() && !includeRefBaseForIndel) {
                newRefBases = ArrayUtils.addAll(newRefBases, Arrays.copyOfRange(altAllele.getBases(),1, altAllele.length()));
            // else add the snp
            } else {
                newRefBases = ArrayUtils.addAll(newRefBases, altAllele.getBases()); // refbases added
            }
            refOffset = vc.getEnd() + 1; //TODO this is probably not set for future reference
        }

        // Finish off the haplotype with the final bases
        int refStartIndex = (int) (refOffset - genomicStartPosition);
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
        long refOffset = genomicStartPosition;

        byte[] refBasesToAddTo = base.getBases();
        CigarBuilder runningCigar = new CigarBuilder();
        byte[] newHaplotypeBasees = {};
        byte[] pdBytes = {};

        //ASSUME sorted for now
        // use the reverse list to save myself figuring out cigars for right now
        for (VariantContext vc : constituentEvents) {
            int refStartIndex = (int) (refOffset - genomicStartPosition);
            int refEndIndex = Math.toIntExact(vc.getStart() - genomicStartPosition);
            //Check that we are allowed to add this event (and if not we are)
            if ((vc.isIndel() && refEndIndex - refStartIndex < -1) || (!vc.isIndel() && refEndIndex - refStartIndex < 0)) {
                throw new RuntimeException("Variant "+vc+" is out of order in the PD event list: "+constituentEvents);
            }

            // Add the cigar for bases we skip over
            if (refEndIndex - refStartIndex > 0) {
                runningCigar.add(new CigarElement(refEndIndex - refStartIndex, CigarOperator.M));
            }

            // Include the ref base for indel if the base immediately proceeding this event is not already tracked
            boolean includeRefBaseForIndel = vc.isIndel() && (refEndIndex - refStartIndex > 0 || refOffset == genomicStartPosition);

            // Determine the alleles to add
            Allele refAllele = vc.getReference();
            Allele altAllele = vc.getAlternateAllele(0);
            boolean shouldFlip = altAllele.length() > refAllele.length();
            boolean isEvent = false;
            // WILL ACCEPT MULTIPLE ALT ALLELES FOR NOW
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
            if (refEndIndex - refStartIndex > 0) {
                newHaplotypeBasees = ArrayUtils.addAll(newHaplotypeBasees, ArrayUtils.subarray(refBasesToAddTo, refStartIndex, (int) (vc.getStart() - genomicStartPosition))); // bases before the variant
                pdBytes = ArrayUtils.addAll(pdBytes, new byte[vc.getStart() - (int)refOffset]); // bases before the variant
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
            refOffset = vc.getEnd() + 1; //TODO this is probably not set for future reference
        }

        // Finish off the haplotype with the final bases
        int refStartIndex = (int) (refOffset - genomicStartPosition);
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
}
