package org.broadinstitute.hellbender.tools.walkers;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

/**
 * Variant context utilities related to merging variant-context instances.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@SuppressWarnings({"rawtypes","unchecked"}) //TODO fix uses of untyped Comparable.
final class ReferenceConfidenceVariantContextMerger {

    /**
     * Combine annotation values by computing the medianish middle element
     * @param values a list of integers
     * @return the median element
     */
    private static <T extends Comparable<? super T>> T combineAnnotationValues(final List<T> values) {
        Utils.nonEmpty(values);
        final int size = values.size();
        if ( size == 1 ) { return values.get(0); }
        else {
            final List<T> sorted = values.stream().sorted().collect(Collectors.toList());
            return sorted.get(size / 2);
        }
    }

    /**
     * Merges VariantContexts from gVCFs into a single hybrid.
     * Assumes that none of the input records are filtered.
     *
     * @param vcs     collection of unsorted genomic vcs
     * @param loc     the current location
     * @param refBase the reference allele to use if all contexts in the VC are spanning (i.e. don't start at the location in loc); if null, we'll return null in this case
     * @param removeNonRefSymbolicAllele if true, remove the <NON_REF> allele from the merged VC
     * @param samplesAreUniquified  if true, sample names have been uniquified
     * @return new VariantContext representing the merge of all vcs or null if it not relevant
     */
    public static VariantContext merge(final List<VariantContext> vcs, final Locatable loc, final Byte refBase,
                                       final boolean removeNonRefSymbolicAllele, final boolean samplesAreUniquified) {
        Utils.nonEmpty(vcs);

        // establish the baseline info (sometimes from the first VC)
        final String name = vcs.get(0).getSource();

        // ref allele
        final Allele refAllele = determineReferenceAlleleGivenReferenceBase(vcs, loc, refBase);
        if ( refAllele == null ) {
            return null;
        }

        // In this list we hold the mapping of each variant context alleles.
        final List<VCWithNewAlleles> vcAndNewAllelePairs = new ArrayList<>(vcs.size());

        // cycle through and add info from the other vcs
        for ( final VariantContext vc : vcs ) {
            // if this context doesn't start at the current location then it must be a spanning event (deletion or ref block)
            final boolean isSpanningEvent = loc.getStart() != vc.getStart();
            vcAndNewAllelePairs.add(new VCWithNewAlleles(vc, isSpanningEvent ? replaceWithNoCallsAndDels(vc) : remapAlleles(vc, refAllele), isSpanningEvent));
        }

        final List<Allele> allelesList = collectTargetAlleles(vcAndNewAllelePairs, refAllele, removeNonRefSymbolicAllele);

        final Set<String> rsIDs = new LinkedHashSet<>(1); // most of the time there's one id
        int depth = 0;
        final Map<String, List<Comparable>> annotationMap = new LinkedHashMap<>();

        final GenotypesContext genotypes = GenotypesContext.create();

        for ( final VCWithNewAlleles vcWithNewAlleles : vcAndNewAllelePairs ) {
            final VariantContext vc = vcWithNewAlleles.getVc();
            final List<Allele> remappedAlleles = vcWithNewAlleles.getNewAlleles();

            genotypes.addAll(mergeRefConfidenceGenotypes(vc, remappedAlleles, allelesList, samplesAreUniquified));
            depth += calculateVCDepth(vc);

            if ( loc.getStart() != vc.getStart() ) {
                continue;
            }

            // special case ID (just preserve it)
            if ( vc.hasID() ) {
                rsIDs.add(vc.getID());
            }

            // add attributes
            addReferenceConfidenceAttributes(vc.getAttributes(), annotationMap);
        }

        final Map<String, Object> attributes = mergeAttributes(depth, annotationMap);

        final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : String.join(",", rsIDs);

        // note that in order to calculate the end position, we need a list of alleles that doesn't include anything symbolic
        final VariantContextBuilder builder = new VariantContextBuilder().source(name).id(ID).alleles(allelesList)
                .chr(loc.getContig()).start(loc.getStart()).computeEndFromAlleles(nonSymbolicAlleles(allelesList), loc.getStart(), loc.getStart())
                .genotypes(genotypes).unfiltered().attributes(new TreeMap<>(attributes)).log10PError(CommonInfo.NO_LOG10_PERROR);  // we will need to re-genotype later

        return builder.make();
    }

    /**
     * Replaces any alleles in the VariantContext with NO CALLS or the symbolic deletion allele as appropriate, except for the generic ALT allele
     *
     * @param vc   VariantContext with the alleles to replace
     * @return non-null list of alleles
     */
    private static List<Allele> replaceWithNoCallsAndDels(final VariantContext vc) {
        Utils.nonNull(vc);

        final List<Allele> result = new ArrayList<>(vc.getNAlleles());

        // no-call the reference allele
        result.add(Allele.NO_CALL);

        // handle the alternate alleles
        for ( final Allele allele : vc.getAlternateAlleles() ) {
            final Allele replacement;
            if ( allele.equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) ) {
                replacement = allele;
            } else if ( allele.length() < vc.getReference().length() ) {
                replacement = GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED;
            } else {
                replacement = Allele.NO_CALL;
            }

            result.add(replacement);
        }
        return result;
    }

    /**
     * This method does a couple of things:
     * <ul><li>
     *     remaps the vc alleles considering the differences between the final reference allele and its own reference,</li>
     * <li>
     *     collects alternative alleles present in variant context and add them to the {@code finalAlleles} set.
     * </li></ul>
     *
     * @param vc           the variant context.
     * @param refAllele    final reference allele.
     * @return never {@code null}
     */
    //TODO as part of a larger refactoring effort {@link #remapAlleles} can be merged with {@link GATKVariantContextUtils#remapAlleles}.
    @VisibleForTesting
    static List<Allele> remapAlleles(final VariantContext vc, final Allele refAllele) {
        final Allele vcRef = vc.getReference();
        final byte[] refBases = refAllele.getBases();
        final int extraBaseCount = refBases.length - vcRef.getBases().length;
        if (extraBaseCount < 0) {
            throw new IllegalStateException("the wrong reference was selected");
        }

        final List<Allele> result = new ArrayList<>(vc.getNAlleles());
        result.add(refAllele);

        for (final Allele a : vc.getAlternateAlleles()) {
            if (a.isSymbolic()) {
                result.add(a);
            } else if (a.isCalled()) {
                result.add(extendAllele(a, extraBaseCount, refBases));
            } else { // NO_CALL and strange miscellanea
                result.add(a);
            }
        }
        return result;
    }

    private static class VCWithNewAlleles {
        private final VariantContext vc;
        private final List<Allele> newAlleles;
        private final boolean isSpanningEvent;

        private VCWithNewAlleles(final VariantContext vc, final List<Allele> newAlleles, final boolean isSpanningEvent) {
            this.vc = vc;
            this.newAlleles = newAlleles;
            this.isSpanningEvent = isSpanningEvent;
        }

        public List<Allele> filterAllelesForFinalSet() {
            return newAlleles.stream().filter(a -> !a.equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE))
                    .filter(a -> !a.isReference())
                    .filter(a -> !(a.isSymbolic() && vc.isSymbolic())) // skip <*DEL> if there isn't a real alternate allele.
                    .filter(Allele::isCalled) // skip NO_CALL
                    .collect(Collectors.toList());
        }

        // record whether it's also a spanning deletion/event (we know this because the VariantContext type is no
        // longer "symbolic" but "mixed" because there are real alleles mixed in with the symbolic non-ref allele)
        public boolean isSpanningDeletion() {
            return  (isSpanningEvent && vc.isMixed()) || vc.getAlleles().contains(GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED);
        }

        public boolean isNonSpanningEvent() {
            return !isSpanningEvent && vc.isMixed();

        }
        public VariantContext getVc() {
            return vc;
        }

        public List<Allele> getNewAlleles() {
            return newAlleles;
        }
    }

    private static List<Allele> collectTargetAlleles(final List<VCWithNewAlleles> vcAndNewAllelePairs, final Allele refAllele, final boolean removeNonRefSymbolicAllele) {
        // FinalAlleleSet contains the alleles of the new resulting VC
        // Using linked set in order to guarantee a stable order
        final Set<Allele> finalAlleleSet = new LinkedHashSet<>(10);
        // Reference goes first
        finalAlleleSet.add(refAllele);

        final List<Allele> collected = vcAndNewAllelePairs.stream()
                .flatMap(vcWithNewAlleles -> vcWithNewAlleles.filterAllelesForFinalSet().stream())
                .collect(Collectors.toList());
        finalAlleleSet.addAll(collected);

        final boolean sawSpanningDeletion = vcAndNewAllelePairs.stream().anyMatch(VCWithNewAlleles::isSpanningDeletion);
        final boolean sawNonSpanningEvent = vcAndNewAllelePairs.stream().anyMatch(VCWithNewAlleles::isNonSpanningEvent);

        // Add <DEL> and <NON_REF> to the end if at all required in in the output.
        if ( sawSpanningDeletion && (sawNonSpanningEvent || !removeNonRefSymbolicAllele) ) {
            finalAlleleSet.add(GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED);
        }
        if (!removeNonRefSymbolicAllele) {
            finalAlleleSet.add(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        }

        return new ArrayList<>(finalAlleleSet);
    }

    /**
     * lookup the depth from the VC DP field or calculate by summing the depths of the genotypes
     */
    private static int calculateVCDepth(VariantContext vc) {

        if ( vc.hasAttribute(VCFConstants.DEPTH_KEY) ) {
            return vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
        } else { // handle the gVCF case from the HaplotypeCaller
            return vc.getGenotypes().stream()
                    .mapToInt(ReferenceConfidenceVariantContextMerger::getBestDepthValue)
                    .sum();
        }
    }

    private static Map<String, Object> mergeAttributes(int depth, Map<String, List<Comparable>> annotationMap) {
        final Map<String, Object> attributes = new LinkedHashMap<>();

        // when combining annotations use the median value from all input vcs which had annotations provided
        annotationMap.entrySet().stream()
                .filter(p -> !p.getValue().isEmpty())
                .forEachOrdered(p -> attributes.put(p.getKey(), combineAnnotationValues(p.getValue()))
        );

        if ( depth > 0 ) {
            attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));
        }

        // remove stale AC and AF based attributes
        removeStaleAttributesAfterMerge(attributes);
        return attributes;
    }

    private static int getBestDepthValue(final Genotype gt) {
        if (gt.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY)) {
            return Integer.parseInt((String) gt.getAnyAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY));
        } else {
            return gt.hasDP() ? gt.getDP() : 0;
        }
    }

    /**
     * @param alleles  the original alleles alleles
     * @return a non-null alleles of non-symbolic alleles
     */
    private static List<Allele> nonSymbolicAlleles(final List<Allele> alleles) {
        return alleles.stream()
                .filter(a -> !a.isSymbolic())
                .collect(Collectors.toList());
    }

    /**
     * Determines the ref allele given the provided reference base at this position
     *
     * @param VCs     collection of unsorted genomic VCs
     * @param loc     the current location
     * @param refBase the reference allele to use if all contexts in the VC are spanning
     * @return new Allele or null if no reference allele/base is available
     */
    private static Allele determineReferenceAlleleGivenReferenceBase(final List<VariantContext> VCs, final Locatable loc, final Byte refBase) {
        final Allele refAllele = GATKVariantContextUtils.determineReferenceAllele(VCs, loc);
        if ( refAllele == null ) {
            if (refBase == null) {
                return null;
            } else {
                return Allele.create(refBase, true);
            }
        } else {
            return refAllele;
        }
    }

    /**
     * Remove the stale attributes from the merged set
     *
     * @param attributes the attribute map
     */
    private static void removeStaleAttributesAfterMerge(final Map<String, Object> attributes) {
        attributes.remove(VCFConstants.ALLELE_COUNT_KEY);
        attributes.remove(VCFConstants.ALLELE_FREQUENCY_KEY);
        attributes.remove(VCFConstants.ALLELE_NUMBER_KEY);
        attributes.remove(GATKVCFConstants.MLE_ALLELE_COUNT_KEY);
        attributes.remove(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY);
        attributes.remove(VCFConstants.END_KEY);
    }

    /**
     * Adds attributes to the global map from the new context in a sophisticated manner
     *
     * @param myAttributes               attributes to add from
     * @param annotationMap              map of annotations for combining later
     */
    @VisibleForTesting
    static <T extends Comparable<? super T>> void addReferenceConfidenceAttributes(final Map<String, Object> myAttributes,
                                                         final Map<String, List<Comparable>> annotationMap) {
        for ( final Map.Entry<String, Object> p : myAttributes.entrySet() ) {
            final String key = p.getKey();
            final Object value = p.getValue();

            // add the annotation values to a list for combining later
            List<Comparable> values = annotationMap.get(key);
            if( values == null ) {
                values = new ArrayList<>();
                annotationMap.put(key, values);
            }
            try {
                final String stringValue = value.toString();
                // Branch to avoid unintentional, implicit type conversions that occur with the ? operator.
                if (stringValue.contains(".")) {
                    values.add(Double.parseDouble(stringValue));
                } else {
                    values.add(Integer.parseInt(stringValue));
                }
            } catch (final NumberFormatException e) {
                // nothing to do
                //TODO This comes directly from gatk3, we should determine if we should log this or if it happens so much it is not useful to log
            }
        }
    }

    /**
     * prefix an allele with additional reference bases if extraBaseCount > 0
     */
    private static Allele extendAllele(Allele allele, int extraBaseCount, byte[] refBases) {
        if (extraBaseCount > 0) {
            final byte[] oldBases = allele.getBases();
            final byte[] newBases = Arrays.copyOf(oldBases, oldBases.length + extraBaseCount);
            System.arraycopy(refBases, refBases.length - extraBaseCount, newBases, oldBases.length, extraBaseCount);
            return Allele.create(newBases,false);
        } else {
            return allele;
        }
    }


    /**
     * Merge into the context a new genotype represented by the given VariantContext for the provided list of target alleles.
     * This method assumes that none of the alleles in the VC overlaps with any of the alleles in the set.
     *  @param vc                    the Variant Context for the sample
     * @param remappedAlleles       the list of remapped alleles for the sample
     * @param targetAlleles         the list of target alleles
     * @param samplesAreUniquified  true if sample names have been uniquified
     */
    private static GenotypesContext mergeRefConfidenceGenotypes(final VariantContext vc,
                                                    final List<Allele> remappedAlleles,
                                                    final List<Allele> targetAlleles,
                                                    final boolean samplesAreUniquified) {
        final GenotypesContext mergedGenotypes = GenotypesContext.create();
        final int maximumPloidy = vc.getMaxPloidy(GATKVariantContextUtils.DEFAULT_PLOIDY);
        // the map is different depending on the ploidy, so in order to keep this method flexible (mixed ploidies)
        // we need to get a map done (lazily inside the loop) for each ploidy, up to the maximum possible.
        final int[][] genotypeIndexMapsByPloidy = new int[maximumPloidy + 1][];
        final int maximumAlleleCount = Math.max(remappedAlleles.size(),targetAlleles.size());
        int[] perSampleIndexesOfRelevantAlleles;

        for ( final Genotype g : vc.getGenotypes() ) {
            final String name;
            if (samplesAreUniquified) {
                name = g.getSampleName() + "." + vc.getSource();
            } else {
                name = g.getSampleName();
            }
            final int ploidy = g.getPloidy();
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(g).alleles(GATKVariantContextUtils.noCallAlleles(g.getPloidy()));
            genotypeBuilder.name(name);
            if (g.hasPL()) {
                // lazy initialization of the genotype index map by ploidy.
                perSampleIndexesOfRelevantAlleles = getIndexesOfRelevantAlleles(remappedAlleles, targetAlleles, vc.getStart(), g);
                final GenotypeLikelihoodCalculators calculators = new GenotypeLikelihoodCalculators();
                final int[] genotypeIndexMapByPloidy = genotypeIndexMapsByPloidy[ploidy] == null
                            ? calculators.getInstance(ploidy, maximumAlleleCount).genotypeIndexMap(perSampleIndexesOfRelevantAlleles, calculators) //probably horribly slow
                            : genotypeIndexMapsByPloidy[ploidy];
                final int[] PLs = generatePL(g, genotypeIndexMapByPloidy);
                final int[] AD = g.hasAD() ? generateAD(g.getAD(), perSampleIndexesOfRelevantAlleles) : null;
                genotypeBuilder.PL(PLs).AD(AD);
            }
            mergedGenotypes.add(genotypeBuilder.make());
        }

        return mergedGenotypes;
    }

    /**
     * Composes a new likelihood array given the original genotype and the genotype index map.
     *
     * @param g the original genotype.
     * @param genotypeIndexMapByPloidy genotype index map. The ith element indicates what genotype in {@code g} corresponds
     *                                 to the ith genotype in the return likelihoods array.
     *
     * @throws NullPointerException if {@code g} or {@code genotypeIndexMapByPloidy} is {@code null}, or if {@code g}
     *    does not contain likelihoods.
     * @throws IndexOutOfBoundsException if {@code genotypeIndexMapByPloidy} contain non valid
     *  genotype indices given the likelihood array in {@code g}.
     *
     * @return never {@code null} but an array of exactly {@code genotypeIndexMapByPloidy.length} positions.
     */
    private static int[] generatePL(final Genotype g, final int[] genotypeIndexMapByPloidy) {
        final int[] PLs = new int[genotypeIndexMapByPloidy.length];
        final int[] oldPLs = g.getPL();
        for (int i = 0; i < PLs.length; i++) {
            PLs[i] = oldPLs[genotypeIndexMapByPloidy[i]];
        }
        return PLs;
    }

    /**
     * Determines the allele mapping from myAlleles to the targetAlleles, substituting the generic "<ALT>" as appropriate.
     * If the myAlleles set does not contain "<ALT>" as an allele, it throws an exception.
     *
     * @param remappedAlleles   the list of alleles to evaluate
     * @param targetAlleles     the target list of alleles
     * @param position          position to output error info
     * @param g                 genotype from which targetAlleles are derived
     * @return non-null array of ints representing indexes
     */
    @VisibleForTesting
    static int[] getIndexesOfRelevantAlleles(final List<Allele> remappedAlleles, final List<Allele> targetAlleles, final int position, final Genotype g) {

        Utils.nonEmpty(remappedAlleles);
        Utils.nonEmpty(targetAlleles);

        if ( !remappedAlleles.contains(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) ) {
            throw new UserException("The list of input alleles must contain " + GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE + " as an allele but that is not the case at position " + position + "; please use the Haplotype Caller with gVCF output to generate appropriate records");
        }

        final int indexOfNonRef = remappedAlleles.indexOf(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        final int[] indexMapping = new int[targetAlleles.size()];

        // the reference likelihoods should always map to each other (even if the alleles don't)
        indexMapping[0] = 0;

        // create the index mapping, using the <NON-REF> allele whenever such a mapping doesn't exist
        for ( int i = 1; i < targetAlleles.size(); i++ ) {
            final int indexOfRemappedAllele = remappedAlleles.indexOf(targetAlleles.get(i));
            indexMapping[i] = indexOfRemappedAllele == -1 ? indexOfNonRef : indexOfRemappedAllele;
        }

        return indexMapping;
    }

    /**
     * Generates a new AD array by adding zeros for missing alleles given the set of indexes of the Genotype's current
     * alleles from the original AD.
     *
     * @param originalAD    the original AD to extend
     * @param indexesOfRelevantAlleles the indexes of the original alleles corresponding to the new alleles
     * @return non-null array of new AD values
     */
    @VisibleForTesting
    static int[] generateAD(final int[] originalAD, final int[] indexesOfRelevantAlleles) {
        Utils.nonNull(originalAD);
        Utils.nonNull(indexesOfRelevantAlleles);

        final int numADs = indexesOfRelevantAlleles.length;
        final int[] newAD = new int[numADs];

        for ( int i = 0; i < numADs; i++ ) {
            final int oldIndex = indexesOfRelevantAlleles[i];
            if ( oldIndex >= originalAD.length ) {
                newAD[i] = 0;
            } else {
                newAD[i] = originalAD[oldIndex];
            }
        }

        return newAD;
    }


}
