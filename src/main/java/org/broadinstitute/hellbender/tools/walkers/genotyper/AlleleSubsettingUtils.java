package org.broadinstitute.hellbender.tools.walkers.genotyper;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.genotyper.GenotypePriorCalculator;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.Permutation;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Utilities class containing methods for restricting {@link VariantContext} and {@link GenotypesContext} objects to a
 * reduced set of alleles, as well as for choosing the best set of alleles to keep and for cleaning up annotations and
 * genotypes after subsetting.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleSubsettingUtils {
    private AlleleSubsettingUtils() {}  // prevent instantiation
    private static final int PL_INDEX_OF_HOM_REF = 0;
    public static final int NUM_OF_STRANDS = 2; // forward and reverse strands

    private static final GenotypeLikelihoodCalculators GL_CALCS = new GenotypeLikelihoodCalculators();


    /**
     * Create the new GenotypesContext with the subsetted PLs and ADs
     *
     * Will reorder subsetted alleles according to the ordering provided by the list allelesToKeep
     *
     * @param originalGs               the original GenotypesContext
     * @param originalAlleles          the original alleles
     * @param allelesToKeep            the subset of alleles to use with the new Genotypes
     * @param assignmentMethod         assignment strategy for the (subsetted) PLs
     * @param depth                    the original variant DP or 0 if there was no DP
     * @return                         a new non-null GenotypesContext
     */
    public static GenotypesContext subsetAlleles(final GenotypesContext originalGs, final int defaultPloidy,
                                                 final List<Allele> originalAlleles,
                                                 final List<Allele> allelesToKeep,
                                                 final GenotypePriorCalculator gpc,
                                                 final GenotypeAssignmentMethod assignmentMethod,
                                                 final int depth) {
        Utils.nonNull(originalGs, "original GenotypesContext must not be null");
        Utils.nonNull(allelesToKeep, "allelesToKeep is null");
        Utils.nonEmpty(allelesToKeep, "must keep at least one allele");
        Utils.validateArg(allelesToKeep.get(0).isReference(), "First allele must be the reference allele");

        final GenotypesContext newGTs = GenotypesContext.create(originalGs.size());
        final Permutation<Allele> allelePermutation = new IndexedAlleleList<>(originalAlleles).permutation(new IndexedAlleleList<>(allelesToKeep));

        final Map<Integer, int[]> subsettedLikelihoodIndicesByPloidy = new TreeMap<>();
        for (final Genotype g : originalGs) {
            final int ploidy = g.getPloidy() > 0 ? g.getPloidy() : defaultPloidy;
            if (!subsettedLikelihoodIndicesByPloidy.containsKey(ploidy)) {
                subsettedLikelihoodIndicesByPloidy.put(ploidy, subsettedPLIndices(ploidy, originalAlleles, allelesToKeep));
            }
            final int[] subsettedLikelihoodIndices = subsettedLikelihoodIndicesByPloidy.get(ploidy);

            final int expectedNumLikelihoods = GenotypeLikelihoods.numLikelihoods(originalAlleles.size(), ploidy);
            // create the new likelihoods array from the alleles we are allowed to use
            double[] newLikelihoods = null;
            double newLog10GQ = -1;
            if (g.hasLikelihoods()) {
                final double[] originalLikelihoods = g.getLikelihoods().getAsVector();
                newLikelihoods = originalLikelihoods.length == expectedNumLikelihoods ?
                        MathUtils.scaleLogSpaceArrayForNumericalStability(Arrays.stream(subsettedLikelihoodIndices)
                                .mapToDouble(idx -> originalLikelihoods[idx]).toArray()) : null;
                if (newLikelihoods != null) {
                    final int PLindex = MathUtils.maxElementIndex(newLikelihoods);
                    newLog10GQ = GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, newLikelihoods);
                }

            } else if (g.hasGQ()) {
                newLog10GQ = -0.1*g.getGQ();
            }

            final boolean useNewLikelihoods = newLikelihoods != null && (depth != 0 || GATKVariantContextUtils.isInformative(newLikelihoods));
            final GenotypeBuilder gb = new GenotypeBuilder(g);
            if (newLog10GQ != -1) {
                final Map<String, Object> attributes = new HashMap<>(g.getExtendedAttributes());
                gb.log10PError(newLog10GQ);
                attributes.remove(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY);
                gb.noAttributes().attributes(attributes);
            }
            if (useNewLikelihoods) {
                gb.PL(newLikelihoods);
            }
            GATKVariantContextUtils.makeGenotypeCall(g.getPloidy(), gb, assignmentMethod, newLikelihoods, allelesToKeep, g.getAlleles(), gpc);

            // restrict SAC to the new allele subset
            if (g.hasExtendedAttribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY)) {
                final int[] newSACs = subsetSACAlleles(g,  originalAlleles, allelesToKeep);
                gb.attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, newSACs);
            }

            // restrict AD to the new allele subset
            if(g.hasAD()) {
                final int[] oldAD = g.getAD();
                final int[] newAD = IntStream.range(0, allelesToKeep.size()).map(n -> oldAD[allelePermutation.fromIndex(n)]).toArray();
                final int nonRefIndex = allelesToKeep.indexOf(Allele.NON_REF_ALLELE);
                if (nonRefIndex != -1 && nonRefIndex < newAD.length) {
                    newAD[nonRefIndex] = 0;  //we will "lose" coverage here, but otherwise merging NON_REF AD counts with other alleles "creates" reads
                }
                gb.AD(newAD);
            }
            newGTs.add(gb.make());
        }
        return newGTs;
    }


    /**
     *  Remove alternate alleles from a set of genotypes turning removed alleles to no-call and dropping other per-allele attributes
     *
     * @param outputHeader  header for the final output VCF, used to validate annotation counts and types
     * @param originalGs    genotypes with full set of alleles
     * @param allelesToKeep contains the reference allele and may contain the NON_REF
     * @param relevantIndices   indices of allelesToKeep w.r.t. the original VC (including ref and possibly NON_REF)
     * @return
     */
    public static GenotypesContext subsetSomaticAlleles(final VCFHeader outputHeader, final GenotypesContext originalGs,
                                                        final List<Allele> allelesToKeep, final int[] relevantIndices) {
        final GenotypesContext newGTs = GenotypesContext.create(originalGs.size());
        GenotypeBuilder gb;
        for (final Genotype g : originalGs) {
            gb = new GenotypeBuilder(g);
            gb.noAttributes();
            List<Allele> keepGTAlleles = new ArrayList<>(g.getAlleles());
            //keep the "ploidy", (i.e. number of different called alleles) the same, but no-call the ones we drop
            for (Allele a : keepGTAlleles) {
                if (!allelesToKeep.contains(a)) {
                    keepGTAlleles.set(keepGTAlleles.indexOf(a), Allele.NO_CALL);
                }
            }
            gb.alleles(keepGTAlleles);
            gb.AD(generateAD(g.getAD(), relevantIndices));
            Set<String> keys = g.getExtendedAttributes().keySet();
            for (final String key : keys) {
                final VCFFormatHeaderLine headerLine = outputHeader.getFormatHeaderLine(key);
                gb.attribute(key, ReferenceConfidenceVariantContextMerger.generateAnnotationValueVector(headerLine.getCountType(),
                            VariantContextGetters.attributeToList(g.getAnyAttribute(key)), relevantIndices));
            }
            newGTs.add(gb.make());
            }
        return newGTs;
    }

    /**
     * Add the VCF INFO field annotations for the used alleles when subsetting
     *
     * @param vc                    original variant context
     * @param builder               variant context builder with subset of original variant context's alleles
     * @param keepOriginalChrCounts keep the original chromosome counts before subsetting
     * @return variant context builder with updated INFO field attribute values
     */
    public static void addInfoFieldAnnotations(final VariantContext vc, final VariantContextBuilder builder,
                                               final boolean keepOriginalChrCounts) {
        Utils.nonNull(vc);
        Utils.nonNull(builder);
        Utils.nonNull(builder.getAlleles());

        final List<Allele> alleles = builder.getAlleles();
        if (alleles.size() < 2)
            throw new IllegalArgumentException("the variant context builder must contain at least 2 alleles");

        // don't have to subset, the original vc has the same number and hence, the same alleles
        boolean keepOriginal = (vc.getAlleles().size() == alleles.size());

        List<Integer> alleleIndices = builder.getAlleles().stream().map(vc::getAlleleIndex).collect(Collectors.toList());
        if (keepOriginalChrCounts) {
            if (vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY))
                builder.attribute(GATKVCFConstants.ORIGINAL_AC_KEY, keepOriginal ?
                        vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY) : alleleIndices.stream().filter(i -> i > 0).map(j -> vc.getAttributeAsList(VCFConstants.ALLELE_COUNT_KEY).get(j - 1)).collect(Collectors.toList()).get(0));
            if (vc.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY))
                builder.attribute(GATKVCFConstants.ORIGINAL_AF_KEY, keepOriginal ?
                        vc.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) : alleleIndices.stream().filter(i -> i > 0).map(j -> vc.getAttributeAsList(VCFConstants.ALLELE_FREQUENCY_KEY).get(j - 1)).collect(Collectors.toList()).get(0));
            if (vc.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY)) {
                builder.attribute(GATKVCFConstants.ORIGINAL_AN_KEY, vc.getAttribute(VCFConstants.ALLELE_NUMBER_KEY));
            }
        }

        VariantContextUtils.calculateChromosomeCounts(builder, true);
    }

    /**
     * From a given genotype, extract a given subset of alleles and return the new SACs
     *
     * @param g                             genotype to subset
     * @param originalAlleles               the original alleles before subsetting
     * @param allelesToUse                  alleles to use in subset
     * @return                              the subsetted SACs
     */
    private static int[] subsetSACAlleles(final Genotype g, final List<Allele> originalAlleles, final List<Allele> allelesToUse) {

        // Keep original SACs if using all of the alleles
        if ( originalAlleles.size() == allelesToUse.size() ) {
            return getSACs(g);
        } else {
            return makeNewSACs(g, originalAlleles, allelesToUse);
        }
    }

    /**
     * Make a new SAC array from the a subset of the genotype's original SAC
     *
     * @param g               the genotype
     * @param originalAlleles the original alleles before subsetting
     * @param allelesToUse    alleles to use in subset
     * @return subset of SACs from the original genotype, the original SACs if sacIndicesToUse is null
     */
    private static int[] makeNewSACs(final Genotype g, final List<Allele> originalAlleles, final List<Allele> allelesToUse) {

        final int[] oldSACs = getSACs(g);
        final int[] newSACs = new int[NUM_OF_STRANDS * allelesToUse.size()];

        int newIndex = 0;
        for (int alleleIndex = 0; alleleIndex < originalAlleles.size(); alleleIndex++) {
            if (allelesToUse.contains(originalAlleles.get(alleleIndex))) {
                newSACs[NUM_OF_STRANDS * newIndex] = oldSACs[NUM_OF_STRANDS * alleleIndex];
                newSACs[NUM_OF_STRANDS * newIndex + 1] = oldSACs[NUM_OF_STRANDS * alleleIndex + 1];
                newIndex++;
            }
        }

        return newSACs;
    }

    /**
     * Get the genotype SACs
     *
     * @param g the genotype
     * @return an arrays of SACs
     * @throws IllegalArgumentException if the genotype does not have an SAC attribute
     * @throws GATKException if the type of the SACs is unexpected
     */
    private static int[] getSACs(final Genotype g) {

        if ( !g.hasExtendedAttribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY) ) {
            throw new IllegalArgumentException("Genotype must have SAC");
        }

        Class<?> clazz = g.getExtendedAttributes().get(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY).getClass();

        if ( clazz.equals(String.class) ) {
            final String SACsString = (String) g.getExtendedAttributes().get(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY);
            String[] stringSACs = SACsString.split(",");
            final int[] intSACs = new int[stringSACs.length];
            int i = 0;
            for (String sac : stringSACs) {
                intSACs[i++] = Integer.parseInt(sac);
            }

            return intSACs;
        }
        else if ( clazz.equals(int[].class) ) {
            return (int[]) g.getExtendedAttributes().get(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY);
        }
        else {
            throw new GATKException("Unexpected SAC type");
        }
    }

    /**
     * Returns the new set of alleles to use based on a likelihood score: alleles' scores are the sum of their counts in
     * sample genotypes, weighted by the confidence in the genotype calls.
     *
     * In the case of ties, the alleles will be chosen from lowest index to highest index.
     *
     * @param vc target variant context.
     * @param numAltAllelesToKeep number of alt alleles to keep.
     * @return the list of alleles to keep, including the reference and {@link Allele#NON_REF_ALLELE} if present
     *
     */
    public static List<Allele> calculateMostLikelyAlleles(final VariantContext vc, final int defaultPloidy,
                                                          final int numAltAllelesToKeep) {
        Utils.nonNull(vc, "vc is null");
        Utils.validateArg(defaultPloidy > 0, () -> "default ploidy must be > 0 but defaultPloidy=" + defaultPloidy);
        Utils.validateArg(numAltAllelesToKeep > 0, () -> "numAltAllelesToKeep must be > 0, but numAltAllelesToKeep=" + numAltAllelesToKeep);

        final boolean hasSymbolicNonRef = vc.hasAllele(Allele.NON_REF_ALLELE);
        final int numberOfAllelesThatArentProperAlts = hasSymbolicNonRef ? 2 : 1; 
        final int numberOfProperAltAlleles = vc.getNAlleles() - numberOfAllelesThatArentProperAlts;

        if (numAltAllelesToKeep >= numberOfProperAltAlleles) {
            return vc.getAlleles();
        }

        final double[] likelihoodSums = calculateLikelihoodSums(vc, defaultPloidy);
        return filterToMaxNumberOfAltAllelesBasedOnScores(numAltAllelesToKeep, vc.getAlleles(), likelihoodSums);
    }


    /**
     * @param alleles a list of alleles including the reference and possible the NON_REF
     * @return a list of the best proper alt alleles based on the likelihood sums, keeping the reference allele and {@link Allele#NON_REF_ALLELE}
     * the list will include no more than {@code numAltAllelesToKeep + 2} alleles and will maintain the order of the original alleles in {@code vc}
     *
     */
    public static List<Allele> filterToMaxNumberOfAltAllelesBasedOnScores(int numAltAllelesToKeep, List<Allele> alleles, double[] likelihoodSums) {
        final int nonRefAltAlleleIndex = alleles.indexOf(Allele.NON_REF_ALLELE);
        final int numAlleles = alleles.size();
        final Set<Integer> properAltIndexesToKeep = IntStream.range(1, numAlleles).filter(n -> n != nonRefAltAlleleIndex).boxed()
                                                             .sorted(Comparator.comparingDouble((Integer n) -> likelihoodSums[n]).reversed())
                                                             .limit(numAltAllelesToKeep)
                                                             .collect(Collectors.toSet());

        return IntStream.range(0, numAlleles)
                .filter( i ->  i == 0 || i == nonRefAltAlleleIndex || properAltIndexesToKeep.contains(i)  )
                .mapToObj(alleles::get)
                .collect(Collectors.toList());
    }

    /** the likelihood sum for an alt allele is the sum over all samples whose likeliest genotype contains that allele of
     * the GL difference between the most likely genotype and the hom ref genotype
     *
     * Since GLs are log likelihoods, this quantity is thus
     * SUM_{samples whose likeliest genotype contains this alt allele} log(likelihood alt / likelihood hom ref)
     */
    @VisibleForTesting
    static double[] calculateLikelihoodSums(final VariantContext vc, final int defaultPloidy) {
        final double[] likelihoodSums = new double[vc.getNAlleles()];
        for ( final Genotype genotype : vc.getGenotypes().iterateInSampleNameOrder() ) {
            final GenotypeLikelihoods gls = genotype.getLikelihoods();
            if (gls == null) {
                continue;
            }
            final double[] glsVector = gls.getAsVector();
            final int indexOfMostLikelyGenotype = MathUtils.maxElementIndex(glsVector);
            final double GLDiffBetweenRefAndBest = glsVector[indexOfMostLikelyGenotype] - glsVector[PL_INDEX_OF_HOM_REF];
            final int ploidy = genotype.getPloidy() > 0 ? genotype.getPloidy() : defaultPloidy;

            final int[] alleleCounts = new GenotypeLikelihoodCalculators()
                    .getInstance(ploidy, vc.getNAlleles()).genotypeAlleleCountsAt(indexOfMostLikelyGenotype)
                    .alleleCountsByIndex(vc.getNAlleles() - 1);

            for (int allele = 1; allele < alleleCounts.length; allele++) {
                if (alleleCounts[allele] > 0) {
                    likelihoodSums[allele] += GLDiffBetweenRefAndBest;
                }
            }
        }
        return likelihoodSums;
    }

    /**
     * Given a list of original alleles and a subset of new alleles to retain, find the array of old PL indices that correspond
     * to new PL indices i.e. result[7] = old PL index of genotype containing same alleles as the new genotype with PL index 7.
     *
     * This method is written in terms f indices rather than subsetting PLs directly in order to produce output that can be
     * recycled from sample to sample, provided that the ploidy is the same.
     *
     * @param ploidy                Ploidy (number of chromosomes describing PL's)
     * @param originalAlleles       List of original alleles
     * @param newAlleles            New alleles -- must be a subset of {@code originalAlleles}
     * @return                      old PL indices of new genotypes
     */
    public static int[] subsettedPLIndices(final int ploidy, final List<Allele> originalAlleles, final List<Allele> newAlleles) {
        final int[] result = new int[GenotypeLikelihoods.numLikelihoods(newAlleles.size(), ploidy)];
        final Permutation<Allele> allelePermutation = new IndexedAlleleList<>(originalAlleles).permutation(new IndexedAlleleList<>(newAlleles));

        final GenotypeLikelihoodCalculator glCalc = GL_CALCS.getInstance(ploidy, originalAlleles.size());
        for (int oldPLIndex = 0; oldPLIndex < glCalc.genotypeCount(); oldPLIndex++) {
            final GenotypeAlleleCounts oldAlleleCounts = glCalc.genotypeAlleleCountsAt(oldPLIndex);

            final boolean containsOnlyNewAlleles = IntStream.range(0, oldAlleleCounts.distinctAlleleCount())
                    .map(oldAlleleCounts::alleleIndexAt).allMatch(allelePermutation::isKept);

            if (containsOnlyNewAlleles) {
                // make an array in the format described in {@link GenotypeAlleleCounts}:
                // [(new) index of first allele, count of first allele, (new) index of second allele, count of second allele. . .]
                final int[] newAlleleCounts = IntStream.range(0, newAlleles.size()).flatMap(newAlleleIndex ->
                        IntStream.of(newAlleleIndex, oldAlleleCounts.alleleCountFor(allelePermutation.fromIndex(newAlleleIndex)))).toArray();

                final int newPLIndex = glCalc.alleleCountsToIndex(newAlleleCounts);
                result[newPLIndex] = oldPLIndex;
            }
        }
        return  result;
    }

    /**
     * Determines the allele mapping from myAlleles to the targetAlleles, substituting the generic "<NON_REF>" as appropriate.
     * If the remappedAlleles set does not contain "<NON_REF>" as an allele, it throws an exception.
     *
     * @param remappedAlleles   the list of alleles to evaluate
     * @param targetAlleles     the target list of alleles
     * @param position          position to output error info
     * @param g                 genotype from which targetAlleles are derived
     * @return non-null array of ints representing indexes
     */
     public static int[] getIndexesOfRelevantAllelesForGVCF(final List<Allele> remappedAlleles, final List<Allele> targetAlleles, final int position, final Genotype g, final boolean doSomaticMerge) {

        Utils.nonEmpty(remappedAlleles);
        Utils.nonEmpty(targetAlleles);

        if ( !remappedAlleles.contains(Allele.NON_REF_ALLELE) ) {
            throw new UserException("The list of input alleles must contain " + Allele.NON_REF_ALLELE + " as an allele but that is not the case at position " + position + "; please use the Haplotype Caller with gVCF output to generate appropriate records");
        }

        final int indexOfNonRef = remappedAlleles.indexOf(Allele.NON_REF_ALLELE);
        final int[] indexMapping = new int[targetAlleles.size()];

        // the reference likelihoods should always map to each other (even if the alleles don't)
        indexMapping[0] = 0;

        // create the index mapping, using the <NON-REF> allele whenever such a mapping doesn't exist
        for ( int i = 1; i < targetAlleles.size(); i++ ) {
            // if there's more than 1 spanning deletion (*) allele then we need to use the best one
            if (targetAlleles.get(i) == Allele.SPAN_DEL && !doSomaticMerge && g.hasPL()) {
                final int occurrences = Collections.frequency(remappedAlleles, Allele.SPAN_DEL);
                if (occurrences > 1) {
                    final int indexOfBestDel = indexOfBestDel(remappedAlleles, g.getPL(), g.getPloidy());
                    indexMapping[i] = (indexOfBestDel == -1 ? indexOfNonRef : indexOfBestDel);
                    continue;
                }
            }

            final int indexOfRemappedAllele = remappedAlleles.indexOf(targetAlleles.get(i));
            indexMapping[i] = indexOfRemappedAllele == -1 ? indexOfNonRef : indexOfRemappedAllele;
        }

        return indexMapping;
    }

    public static int[] getIndexesOfRelevantAlleles(final List<Allele> remappedAlleles, final List<Allele> targetAlleles, final int position, final Genotype g) {
        Utils.nonEmpty(remappedAlleles);
        Utils.nonEmpty(targetAlleles);

        final int[] indexMapping = new int[targetAlleles.size()];

        // the reference likelihoods should always map to each other (even if the alleles don't)
        indexMapping[0] = 0;

        for ( int i = 1; i < targetAlleles.size(); i++ ) {
            // if there's more than 1 spanning deletion (*) allele then we need to use the best one
            if (targetAlleles.get(i) == Allele.SPAN_DEL && g.hasPL()) {
                final int occurrences = Collections.frequency(remappedAlleles, Allele.SPAN_DEL);
                if (occurrences > 1) {
                    final int indexOfBestDel = indexOfBestDel(remappedAlleles, g.getPL(), g.getPloidy());
                    if (indexOfBestDel == -1) {
                        throw new IllegalArgumentException("At position " + position + " targetAlleles contains a spanning deletion, but remappedAlleles does not.");
                    }
                    indexMapping[i] = indexOfBestDel;
                    continue;
                }
            }

            final int indexOfRemappedAllele = remappedAlleles.indexOf(targetAlleles.get(i));
            if (indexOfRemappedAllele == -1) {
                throw new IllegalArgumentException("At position " + position + " targetAlleles contains a " + targetAlleles.get(i) + " allele, but remappedAlleles does not.");
            }
            indexMapping[i] = indexOfRemappedAllele;
        }

        return indexMapping;
    }

    /**
     * Returns the index of the best spanning deletion allele based on AD counts
     *
     * @param alleles   the list of alleles
     * @param PLs       the list of corresponding PL values
     * @param ploidy    the ploidy of the sample
     * @return the best index or -1 if not found
     */
    private static int indexOfBestDel(final List<Allele> alleles, final int[] PLs, final int ploidy) {
        int bestIndex = -1;
        int bestPL = Integer.MAX_VALUE;

        for ( int i = 0; i < alleles.size(); i++ ) {
            if ( alleles.get(i) == Allele.SPAN_DEL ) {
                final int homAltIndex = findHomIndex(GL_CALCS.getInstance(ploidy, alleles.size()), i, ploidy);
                final int PL = PLs[homAltIndex];
                if ( PL < bestPL ) {
                    bestIndex = i;
                    bestPL = PL;
                }
            }
        }

        return bestIndex;
    }

    /** //TODO simplify these methods
     * Returns the index of the PL that represents the homozygous genotype of the given i'th allele
     *
     * @param i           the index of the allele with the list of alleles
     * @param ploidy      the ploidy of the sample
     * @return the hom index
     */
    private static int findHomIndex(final GenotypeLikelihoodCalculator calculator, final int i, final int ploidy) {
        // some quick optimizations for the common case
        if ( ploidy == 2 )
            return GenotypeLikelihoods.calculatePLindex(i, i);
        if ( ploidy == 1 )
            return i;

        final int[] alleleIndexes = new int[ploidy];
        Arrays.fill(alleleIndexes, i);
        return calculator.allelesToIndex(alleleIndexes);
    }

    /**
     * Generates a new AD array by adding zeros for missing alleles given the set of indexes of the Genotype's current
     * alleles from the original AD.
     *
     * @param originalAD    the original AD to extend
     * @param indexesOfRelevantAlleles the indexes of the original alleles corresponding to the new alleles
     * @return non-null array of new AD values
     */
    public static int[] generateAD(final int[] originalAD, final int[] indexesOfRelevantAlleles) {
        final List<Integer> adList = remapRLengthList(Arrays.stream(originalAD).boxed().collect(Collectors.toList()), indexesOfRelevantAlleles, 0);
        return Ints.toArray(adList);
    }

    /**
     * Generates a new AF (allele fraction) array
     * @param originalAF
     * @param indexesOfRelevantAlleles
     * @return non-null array of new AFs
     */
    public static double[] generateAF(final double[] originalAF, final int[] indexesOfRelevantAlleles) {
        final List<Double> afList = remapALengthList(Arrays.stream(originalAF).boxed().collect(Collectors.toList()), indexesOfRelevantAlleles, 0.0);
        return Doubles.toArray(afList);
    }

    /**
     * Given a list of per-allele attributes including the reference allele, subset to relevant alleles
     * @param originalList
     * @param indexesOfRelevantAlleles
     * @return
     */
    public static <T> List<T> remapRLengthList(final List<T> originalList, final int[] indexesOfRelevantAlleles, T filler) {
        Utils.nonNull(originalList);
        Utils.nonNull(indexesOfRelevantAlleles);

        return remapList(originalList, indexesOfRelevantAlleles, 0, filler);
    }

    /**
     * Given a list of per-alt-allele attributes, subset to relevant alt alleles
     * @param originalList
     * @param indexesOfRelevantAlleles
     * @return
     */
    public static <T> List<T> remapALengthList(final List<T> originalList, final int[] indexesOfRelevantAlleles, T filler) {
        Utils.nonNull(originalList);
        Utils.nonNull(indexesOfRelevantAlleles);

        return remapList(originalList, indexesOfRelevantAlleles, 1, filler);
    }

    /**
     * Subset a list of per-allele attributes
     *
     * @param originalList  input per-allele attributes
     * @param indexesOfRelevantAlleles  indexes of alleles to keep, including the reference
     * @param offset    used to indicate whether to include the ref allele values in the output or not
     * @param filler default value to use if no value is mapped
     * @return  a non-null List
     */
    private static <T> List<T> remapList(final List<T> originalList, final int[] indexesOfRelevantAlleles,
                                            final int offset, T filler) {
        final int numValues = indexesOfRelevantAlleles.length - offset; //since these are log odds, this should just be alts
        final List<T> newValues = new ArrayList<>();

        for ( int i = offset; i < numValues + offset; i++ ) {
            final int oldIndex = indexesOfRelevantAlleles[i];
            if ( oldIndex >= originalList.size() + offset ) {
                newValues.add(i-offset, filler);
            } else {
                newValues.add(i-offset, originalList.get(oldIndex-offset));
            }
        }
        return newValues;
    }
}
