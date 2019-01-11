package org.broadinstitute.hellbender.tools.walkers.genotyper;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.Permutation;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

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

            }

            final boolean useNewLikelihoods = newLikelihoods != null && (depth != 0 || GATKVariantContextUtils.isInformative(newLikelihoods));
            final GenotypeBuilder gb = new GenotypeBuilder(g);
            if (useNewLikelihoods) {
                final Map<String, Object> attributes = new HashMap<>(g.getExtendedAttributes());
                gb.PL(newLikelihoods).log10PError(newLog10GQ);
                attributes.remove(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY);
                gb.noAttributes().attributes(attributes);
            }
            else {
                gb.noPL().noGQ();
            }
            GATKVariantContextUtils.makeGenotypeCall(g.getPloidy(), gb, assignmentMethod, newLikelihoods, allelesToKeep, g.getAlleles());

            // restrict SAC to the new allele subset
            if (g.hasExtendedAttribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY)) {
                final int[] newSACs = subsetSACAlleles(g,  originalAlleles, allelesToKeep);
                gb.attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, newSACs);
            }

            // restrict AD to the new allele subset
            if(g.hasAD()) {
                final int[] oldAD = g.getAD();
                final int[] newAD = IntStream.range(0, allelesToKeep.size()).map(n -> oldAD[allelePermutation.fromIndex(n)]).toArray();
                gb.AD(newAD);
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
     * @param keepOriginalChrCounts keep the orignal chromosome counts before subsetting
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

        List<Integer> alleleIndecies = builder.getAlleles().stream().map(a -> vc.getAlleleIndex(a)).collect(Collectors.toList());
        if (keepOriginalChrCounts) {
            if (vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY))
                builder.attribute(GATKVCFConstants.ORIGINAL_AC_KEY, keepOriginal ?
                        vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY) : alleleIndecies.stream().filter(i -> i > 0).map(j -> vc.getAttributeAsList(VCFConstants.ALLELE_COUNT_KEY).get(j - 1)).collect(Collectors.toList()).get(0));
            if (vc.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY))
                builder.attribute(GATKVCFConstants.ORIGINAL_AF_KEY, keepOriginal ?
                        vc.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) : alleleIndecies.stream().filter(i -> i > 0).map(j -> vc.getAttributeAsList(VCFConstants.ALLELE_FREQUENCY_KEY).get(j - 1)).collect(Collectors.toList()).get(0));
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
     * @return a list of the best proper alt alleles based on the likelihood sums, keeping the reference allele and {@link Allele#NON_REF_ALLELE}
     * the list will include no more than {@code numAltAllelesToKeep + 2} alleles and will maintain the order of the original alleles in {@code vc}
     *
     */
    @VisibleForTesting
    static List<Allele> filterToMaxNumberOfAltAllelesBasedOnScores(int numAltAllelesToKeep, List<Allele> alleles, double[] likelihoodSums) {
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
}
