package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.api.client.util.Lists;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;

import java.util.*;
import java.util.stream.Collectors;

/**
 * TODO
 */
public class PartiallyDeterminedHaplotypeComputationEngine {

    /**
     * TODO
     */
    public static AssemblyResultSet generatePDHaplotypes(AssemblyResultSet sourceSet,
                                                         Haplotype referenceHaplotype,
                                                         SortedSet<VariantContext> assemblyVariants,
                                                         List<VariantContext> pileupAllelesFoundShouldFilter,
                                                         List<VariantContext> pileupAllelesPassingFilters,
                                                         int snpAdjacentToIndelLimit) {
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
                System.err.println("Removing assembly variants due to columnwise heurisits: " + variantsToRemove);
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

        System.out.println("Variants to PDHapDetermination:\n"+
                vcs.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("\n")));

        // TODO this is where we filter out if indels > 32

        // TODO this is where we do eventTestSet work

        // TODO lets break everything into event sets
        Map<Integer, List<VariantContext>> eventGroups = new LinkedHashMap<>();
        int lastEventEnd = -1;
        int lastEventStart = -1;

        for (VariantContext vc : vcs) {
            // Break everything into independent groups (don't worry about transitivitiy right now)
            if (vc.getStart() > lastEventEnd) {
                eventGroups.put(vc.getStart(), new ArrayList<VariantContext>(Collections.singletonList(vc)));
                lastEventStart = vc.getStart();
                lastEventEnd = vc.getEnd();//TODO check this

            } else {
                eventGroups.get(lastEventStart).add(vc);
                lastEventEnd = Math.max(lastEventEnd, vc.getEnd());
            }
        }

        // now iterate over the linked hash map...
        //TODO get branches of compatable variants to treat in tandem
        List<List<VariantContext>> branchedHaplotypeGroups = getBranchedHaplotypeGroups(eventGroups);

        System.out.println("Branches of PD Haps to construct:\n"+
                branchedHaplotypeGroups.stream().map(variantContexts ->
                variantContexts.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")))
                .collect(Collectors.joining("\n")));

        Set<Haplotype> outputHaplotypes = new HashSet<>();
        outputHaplotypes.add(referenceHaplotype);
        //Generate the PDHaplotypes based on each of the branches
        for (List<VariantContext> branch : branchedHaplotypeGroups) {
            System.out.println("Handling Branch " + branch.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
            List<PartiallyDeterminedHaplotype> branchHaps = new ArrayList<>();
            for (VariantContext variantContext : branch) {
                // Add both the variant and the event
                branchHaps.add(createNewPDHaplotypeFromEvents(referenceHaplotype, variantContext, true, branch));
                branchHaps.add(createNewPDHaplotypeFromEvents(referenceHaplotype, variantContext, false, branch));
            }
            System.out.println("Constructed PD Haps:" + branchHaps.stream().map(PartiallyDeterminedHaplotype::toString).collect(Collectors.joining("\n")));
            outputHaplotypes.addAll(branchHaps);
        }
        sourceSet.replaceAllHaplotypes(outputHaplotypes);

        return sourceSet;
    }

    //Helper method that does the work of computing what braches need to exist for later PDHap construction
    //TODO this will EVENTUALLY take as input the results fro smith waterman
    //TODO this approach will not eventually work whne we need to split groups
    private static List<List<VariantContext>> getBranchedHaplotypeGroups(Map<Integer, List<VariantContext>> eventGroups) {
        return recursiveHapGroupHelper(new ArrayList<>(eventGroups.entrySet()), -1, new ArrayList<>());
    }


    private static List<List<VariantContext>> recursiveHapGroupHelper(List<Map.Entry<Integer, List<VariantContext>>> listOfGroups, int positionOfLastAddedElement, List<VariantContext> continueGroup) {
        List<List<VariantContext>> output = new ArrayList<>();

        for (int i = positionOfLastAddedElement + 1; i < listOfGroups.size(); i++) {

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

    @VisibleForTesting
    public static Haplotype constructHaplotypeFromVariants(final Haplotype refHap, List<VariantContext> variantContexts, boolean longestIndelAlelle) {
        //ASSERT that the base is ref and cool
        if (!refHap.isReference() || refHap.getCigar().numCigarElements() > 1) {
            throw new RuntimeException("This is not a valid base haplotype for construction");
        }
        long startPosition = refHap.getStartPosition();
        long refOffset = startPosition;

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

            int refStartIndex = (int) (refOffset - startPosition);
            int refEndIndex = Math.toIntExact(vc.getStart() - startPosition);
            runningCigar.add(new CigarElement(refEndIndex - refStartIndex, CigarOperator.M));
            CigarElement newCigarElement;
            if (refAllele.length() == altAllele.length()) {
                newCigarElement = new CigarElement(refAllele.length(), CigarOperator.X);
            } else {
                runningCigar.add(new CigarElement(1, CigarOperator.M)); //When we add an indel we end up inserting a matching base
                newCigarElement = new CigarElement(Math.abs(altAllele.length() - refAllele.length()),
                        refAllele.length() > altAllele.length() ?
                                CigarOperator.D : CigarOperator.I);
            }
            runningCigar.add(newCigarElement);

            newRefBases = ArrayUtils.addAll(newRefBases, ArrayUtils.subarray(refbases, refStartIndex, (int) (vc.getStart() - startPosition))); // bases before the variant
            newRefBases = ArrayUtils.addAll(newRefBases, altAllele.getBases()); // refbases added
            refOffset = vc.getEnd() + 1; //TODO this is probably not set for future reference
        }

        // Finish off the haplotype with the final bases
        int refStartIndex = (int) (refOffset - startPosition);
        newRefBases = ArrayUtils.addAll(newRefBases, ArrayUtils.subarray(refbases, refStartIndex, refbases.length));
        runningCigar.add(new CigarElement(refbases.length - refStartIndex, CigarOperator.M));

        return new Haplotype(newRefBases, false, refHap.getGenomeLocation(), runningCigar.make());
    }

    /**
     * Construct a PD haplotype from scratch
     *
     *
     * Generally we are constructing a new haplotype with all of the reference bases for SNP events and with the longest possible allele for INDEL events.
     * For deleteions, we extend the haplotype by the ref lenght
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
        long startPosition = base.getStartPosition();
        long refOffset = startPosition;

        byte[] refBasesToAddTo = base.getBases();
        CigarBuilder runningCigar = new CigarBuilder();
        byte[] newHaplotypeBasees = {};
        byte[] pdBytes = {};

        //ASSUME sorted for now
        // use the reverse list to save myself figuring out cigars for right now
        for (VariantContext vc : constituentEvents) {
            int refStartIndex = (int) (refOffset - startPosition);
            int refEndIndex = Math.toIntExact(vc.getStart() - startPosition);
            // Add the cigar for bases we skip over
            runningCigar.add(new CigarElement(refEndIndex - refStartIndex, CigarOperator.M));

            // Determine the alleles to add
            Allele refAllele = vc.getReference();
            List<Allele> altAlleles = vc.getAlternateAlleles();
            Allele longestAltAllele = altAlleles.stream().max(Comparator.comparing(Allele::length)).get();
            boolean shouldFlip = longestAltAllele.length() > refAllele.length();
            // We want to always add the longest allele
            byte[] basesToAdd = refAllele.length() >= longestAltAllele.length() ? refAllele.getBases() : longestAltAllele.getBases();
            boolean isEvent = false;
            // WILL ACCEPT MULTIPLE ALT ALLELES FOR NOW
            if (vc.getStart()==eventWithVariant.getStart()) {
                if (!useRef && eventWithVariant.getAlternateAlleles().size() > 1) {
                    throw new RuntimeException("the Blessed variant must be monoallelic");
                }
                isEvent = true;
                basesToAdd = useRef? refAllele.getBases() : longestAltAllele.getBases();

            }

            // Figure out the cigar to add:
            CigarElement newCigarElement;
            if (isEvent && useRef || refAllele.length() == longestAltAllele.length()) {
                newCigarElement = new CigarElement(refAllele.length(), CigarOperator.M);
            } else {
                if (shouldFlip) {
                    runningCigar.add(new CigarElement(1, CigarOperator.M)); //When we add an indel we end up inserting a matching base
                    newCigarElement = new CigarElement(longestAltAllele.length() - refAllele.length(), CigarOperator.I);
                } else {
                    newCigarElement = new CigarElement(basesToAdd.length, CigarOperator.M);
                }
            }
            runningCigar.add(newCigarElement);

            newHaplotypeBasees = ArrayUtils.addAll(newHaplotypeBasees, ArrayUtils.subarray(refBasesToAddTo, refStartIndex, (int) (vc.getStart() - startPosition))); // bases before the variant
            newHaplotypeBasees = ArrayUtils.addAll(newHaplotypeBasees, basesToAdd); // refbases added
            pdBytes = ArrayUtils.addAll(pdBytes, new byte[vc.getStart() - (int)refOffset]); // bases before the variant
            pdBytes = ArrayUtils.addAll(pdBytes, isEvent? new byte[basesToAdd.length] : PartiallyDeterminedHaplotype.getPDBytesForHaplotypes(shouldFlip? longestAltAllele : refAllele, shouldFlip? Collections.singletonList(refAllele) : altAlleles)); // refbases added
            refOffset = vc.getEnd() + 1; //TODO this is probably not set for future reference
        }

        // Finish off the haplotype with the final bases
        int refStartIndex = (int) (refOffset - startPosition);
        newHaplotypeBasees = ArrayUtils.addAll(newHaplotypeBasees, ArrayUtils.subarray(refBasesToAddTo, refStartIndex, refBasesToAddTo.length));
        pdBytes = ArrayUtils.addAll(pdBytes, new byte[refBasesToAddTo.length - refStartIndex]);
        runningCigar.add(new CigarElement(refBasesToAddTo.length - refStartIndex, CigarOperator.M));

        return new PartiallyDeterminedHaplotype(
                //TODO that genomeloc is bad
                new Haplotype(newHaplotypeBasees, false, base.getGenomeLocation(), runningCigar.make()),
                false,
                pdBytes,
                constituentEvents,
                eventWithVariant,
                runningCigar.make());
    }
}
