package org.broadinstitute.hellbender.utils.genotyper;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.apache.commons.collections.ListUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;

import java.util.*;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Evidence-likelihoods container implementation based on integer indexed arrays.
 *
 * @param <A> the type of the allele the likelihood makes reference to.
 *
 * Note: this class uses FastUtil collections for speed.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class AlleleLikelihoods<EVIDENCE extends Locatable, A extends Allele> implements SampleList, AlleleList<A> {

    public static final double LOG_10_INFORMATIVE_THRESHOLD = 0.2;
    public static final double NATURAL_LOG_INFORMATIVE_THRESHOLD = MathUtils.log10ToLog(LOG_10_INFORMATIVE_THRESHOLD);

    protected boolean isNaturalLog = false;

    private double getInformativeThreshold() {
        return isNaturalLog ? NATURAL_LOG_INFORMATIVE_THRESHOLD : LOG_10_INFORMATIVE_THRESHOLD;
    }

    /**
     * Index indicating that the reference allele is missing.
     */
    private static final int MISSING_REF = -1;

    /**
     * Evidence by sample index. Each sub array contains reference to the evidence of the ith sample.
     */
    protected final List<List<EVIDENCE>> evidenceBySampleIndex;

    /**
     * Indexed per sample, allele and finally evidence (within sample).
     * <p>
     *     valuesBySampleIndex[s][a][r] == lnLk(R_r | A_a) where R_r comes from Sample s.
     * </p>
     */
    protected final double[][][] valuesBySampleIndex;

    /**
     * Sample list
     */
    protected final SampleList samples;

    /**
     * Allele list
     */
    protected AlleleList<A> alleles;

    /**
     * Maps from each unit of evidence to its index within the sample.
     *
     * <p>In order to save CPU time the indices contained in this array (not the array itself) is
     * lazily initialized by invoking {@link #evidenceIndexBySampleIndex(int)}.</p>
     */
    protected final List<Object2IntMap<EVIDENCE>> evidenceIndexBySampleIndex;

    /**
     * Index of the reference allele if any, otherwise {@link #MISSING_REF}.
     */
    private int referenceAlleleIndex = MISSING_REF;

    /**
     * Sample matrices lazily initialized (the elements not the array) by invoking {@link #sampleMatrix(int)}.
     */
    private final LikelihoodMatrix<EVIDENCE,A>[] sampleMatrices;

    /**
     * Is this container expected to have the per-allele liklihoods calculations filled in.
     */
    public boolean hasFilledLikelihoods() {
        return true;
    }

    /**
     * Constructs a new evidence-likelihood collection.
     *
     * <p>
     *     The initial likelihoods for all allele-evidence combinations are
     *     0.
     * </p>
     *
     * @param samples all supported samples in the collection.
     * @param alleles all supported alleles in the collection.
     * @param evidenceBySample evidence stratified per sample.
     *
     * @throws IllegalArgumentException if any of {@code allele}, {@code samples}
     * or {@code evidenceBySample} is {@code null},
     *  or if they contain null values.
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    public AlleleLikelihoods(final SampleList samples,
                             final AlleleList<A> alleles,
                             final Map<String, List<EVIDENCE>> evidenceBySample) {
        Utils.nonNull(alleles);
        Utils.nonNull(samples);
        Utils.nonNull(evidenceBySample);

        this.samples = samples;
        this.alleles = alleles;

        final int sampleCount = samples.numberOfSamples();
        final int alleleCount = alleles.numberOfAlleles();

        evidenceBySampleIndex = new ArrayList<>(sampleCount);
        valuesBySampleIndex = new double[sampleCount][][];
        referenceAlleleIndex = findReferenceAllele(alleles);

        evidenceIndexBySampleIndex = new ArrayList<>(Collections.nCopies(sampleCount, null));

        setupIndexes(evidenceBySample, sampleCount, alleleCount);

        sampleMatrices = (LikelihoodMatrix<EVIDENCE, A>[]) new LikelihoodMatrix[sampleCount];
    }


    // Internally used constructor.
    @SuppressWarnings({"unchecked", "rawtypes"})
    AlleleLikelihoods(final AlleleList alleles,
                      final SampleList samples,
                      final List<List<EVIDENCE>> evidenceBySampleIndex,
                      final double[][][] values) {
        this.samples = samples;
        this.alleles = alleles;
        this.evidenceBySampleIndex = evidenceBySampleIndex;
        this.valuesBySampleIndex = values;
        final int sampleCount = samples.numberOfSamples();
        evidenceIndexBySampleIndex = new ArrayList<>(Collections.nCopies(sampleCount, null));

        referenceAlleleIndex = findReferenceAllele(alleles);
        sampleMatrices = (LikelihoodMatrix<EVIDENCE,A>[]) new LikelihoodMatrix[sampleCount];
    }

    // Add all the indices to alleles, sample and evidence in the look-up maps.
    private void setupIndexes(final Map<String, List<EVIDENCE>> evidenceBySample, final int sampleCount, final int alleleCount) {
        for (int s = 0; s < sampleCount; s++) {
            final String sample = samples.getSample(s);

            evidenceBySampleIndex.add(evidenceBySample.getOrDefault(sample, new ArrayList<>()));
            final int sampleEvidenceCount = evidenceBySampleIndex.get(s).size();

            final double[][] sampleValues = new double[alleleCount][sampleEvidenceCount];
            valuesBySampleIndex[s] = sampleValues;
        }
    }

    // Search for the reference allele, if not found the index is {@link MISSING_REF}.
    private static int findReferenceAllele(final AlleleList<?> alleles) {
        return IntStream.range(0, alleles.numberOfAlleles()).filter(i -> alleles.getAllele(i).isReference()).findAny().orElse(MISSING_REF);
    }

    /**
     * Returns the index of a sample within the likelihood collection.
     *
     * @param sample the query sample.
     *
     * @throws IllegalArgumentException if {@code sample} is {@code null}.
     * @return -1 if the allele is not included, 0 or greater otherwise.
     */
    @Override
    public int indexOfSample(final String sample) {
        return samples.indexOfSample(sample);
    }

    /**
     * Number of samples included in the likelihood collection.
     * @return 0 or greater.
     */
    @Override
    public int numberOfSamples() {
        return samples.numberOfSamples();
    }

    /**
     * Returns sample name given its index.
     *
     * @param sampleIndex query index.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is negative.
     *
     * @return never {@code null}.
     */
    @Override
    public String getSample(final int sampleIndex) {
        return samples.getSample(sampleIndex);
    }

    /**
     * Returns the index of an allele within the likelihood collection.
     *
     * @param allele the query allele.
     *
     * @throws IllegalArgumentException if {@code allele} is {@code null}.
     *
     * @return -1 if the allele is not included, 0 or greater otherwise.
     */
    @Override
    public int indexOfAllele(final A allele) {
        return alleles.indexOfAllele(allele);
    }

    /**
     * Returns number of alleles in the collection.
     * @return 0 or greater.
     */
    @Override
    public int numberOfAlleles() {
        return alleles.numberOfAlleles();
    }

    /**
     * Returns the allele given its index.
     *
     * @param alleleIndex the allele index.
     *
     * @throws IllegalArgumentException the allele index is {@code null}.
     *
     * @return never {@code null}.
     */
    @Override
    public A getAllele(final int alleleIndex) {
        return alleles.getAllele(alleleIndex);
    }

    /**
     * Returns the units of evidence that belong to a sample sorted by their index (within that sample).
     *
     * @param sampleIndex the requested sample.
     * @return never {@code null} but perhaps a zero-length array if there is no evidence in sample. No element in
     *   the array will be null.
     */
    public List<EVIDENCE> sampleEvidence(final int sampleIndex) {
        return Collections.unmodifiableList(evidenceBySampleIndex.get(sampleIndex));
    }

    /**
     * Returns an evidence vs allele likelihood matrix corresponding to a sample.
     *
     * @param sampleIndex target sample.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is not null.
     *
     * @return never {@code null}
     */
    public LikelihoodMatrix<EVIDENCE,A> sampleMatrix(final int sampleIndex) {
        Utils.validIndex(sampleIndex, samples.numberOfSamples());
        final LikelihoodMatrix<EVIDENCE,A> extantResult = sampleMatrices[sampleIndex];
        if (extantResult == null) {
            return sampleMatrices[sampleIndex] = new SampleMatrix(sampleIndex);
        } else {
            return extantResult;
        }
    }

    public void switchToNaturalLog() {
        Utils.validate(!isNaturalLog, "Likelihoods have already been switched to natural log");
        final int sampleCount = samples.numberOfSamples();
        final int alleleCount = alleles.numberOfAlleles();

        for (int s = 0; s < sampleCount; s++) {
            final int evidenceCount = sampleEvidenceCount(s);
            for (int a = 0; a < alleleCount; a++) {
                for (int e = 0; e < evidenceCount; e++) {
                    valuesBySampleIndex[s][a][e] = MathUtils.log10ToLog(valuesBySampleIndex[s][a][e]);
                }
            }
        }
        isNaturalLog = true;
    }

    /**
     * Downsamples reads based on contamination fractions making sure that all alleles are affected proportionally.
     *
     * @param perSampleDownsamplingFraction contamination sample map where the sample name are the keys and the
     *                                       fractions are the values.
     *
     * @throws IllegalArgumentException if {@code perSampleDownsamplingFraction} is {@code null}.
     */
    public void contaminationDownsampling(final Map<String, Double> perSampleDownsamplingFraction) {
        Utils.nonNull(perSampleDownsamplingFraction);

        final int alleleCount = alleles.numberOfAlleles();
        for (int s = 0; s < samples.numberOfSamples(); s++) {
            final double fraction = perSampleDownsamplingFraction.getOrDefault(samples.getSample(s), 0.0);
            if (Double.isNaN(fraction) || fraction <= 0.0) {
                continue;
            }

            final Map<A, List<EVIDENCE>> alleleEvidenceMap = evidenceByBestAlleleMap(s);
            final Collection<EVIDENCE> evidenceToRemove = fraction >= 1 ? evidenceBySampleIndex.get(s) :
                    AlleleBiasedDownsamplingUtils.selectAlleleBiasedEvidence(alleleEvidenceMap, fraction);
            removeSampleEvidence(s, evidenceToRemove, alleleCount);
        }
    }

    /**
     * Adjusts likelihoods so that for each unit of evidence, the best allele likelihood is 0 and caps the minimum likelihood
     * of any allele for each unit of evidence based on the maximum alternative allele likelihood.
     *
     * @param maximumLikelihoodDifferenceCap maximum difference between the best alternative allele likelihood
     *                                           and any other likelihood.
     *
     * @throws IllegalArgumentException if {@code maximumDifferenceWithBestAlternative} is not 0 or less.
     */
    public void normalizeLikelihoods(final double maximumLikelihoodDifferenceCap) {
        Utils.validateArg(maximumLikelihoodDifferenceCap < 0.0 && !Double.isNaN(maximumLikelihoodDifferenceCap),
                "the minimum reference likelihood fall must be negative");

        if (maximumLikelihoodDifferenceCap == Double.NEGATIVE_INFINITY) {
            return;
        }

        final int alleleCount = alleles.numberOfAlleles();
        if (alleleCount == 0){ // trivial case there is no alleles.
            return;
        } else if (alleleCount == 1) {
            return;
        }

        for (int s = 0; s < valuesBySampleIndex.length; s++) {
            final double[][] sampleValues = valuesBySampleIndex[s];
            final int evidenceCount = evidenceBySampleIndex.get(s).size();
            for (int r = 0; r < evidenceCount; r++) {
                normalizeLikelihoodsPerEvidence(maximumLikelihoodDifferenceCap, sampleValues, s, r);
            }
        }
    }

    // Does the normalizeLikelihoods job for each piece of evidence.
    private void normalizeLikelihoodsPerEvidence(final double maximumBestAltLikelihoodDifference,
                                                 final double[][] sampleValues, final int sampleIndex, final int evidenceIndex) {

        //allow the best "alternative" allele to be the reference because asymmetry leads to strange artifacts like het calls with >90% alt reads
        final BestAllele bestAlternativeAllele = searchBestAllele(sampleIndex,evidenceIndex,true);

        final double worstLikelihoodCap = bestAlternativeAllele.likelihood + maximumBestAltLikelihoodDifference;

        final int alleleCount = alleles.numberOfAlleles();

        // Guarantee to be the case by enclosing code.
        for (int a = 0; a < alleleCount; a++) {
            if (sampleValues[a][evidenceIndex] < worstLikelihoodCap) {
                sampleValues[a][evidenceIndex] = worstLikelihoodCap;
            }
        }

    }

    /**
     * Returns the samples in this evidence-likelihood collection.
     * <p>
     *     Samples are sorted by their index in the collection.
     * </p>
     *
     * <p>
     *     The returned list is an unmodifiable view on the evidence-likelihoods sample list.
     * </p>
     *
     * @return never {@code null}.
     */
    public List<String> samples() {
        return samples.asListOfSamples();
    }

    /**
     * Returns the samples in this evidence-likelihood collection.
     * <p>
     *     Samples are sorted by their index in the collection.
     * </p>
     *
     * <p>
     *     The returned list is an unmodifiable. It will not be updated if the collection
     *     allele list changes.
     * </p>
     *
     * @return never {@code null}.
     */
    public List<A> alleles() {
        return alleles.asListOfAlleles();
    }


    /**
     * Search the best allele for a unit of evidence.
     *
     * @param sampleIndex including sample index.
     * @param evidenceIndex  target evidence index.
     *
     * @param priorities An array of allele priorities (higher values have higher priority) to be used, if present, to break ties for
     *                   uninformative likelihoods, in which case the evidence is assigned to the allele with the higher score.
     * @return never {@code null}, but with {@link BestAllele#allele allele} == {@code null}
     * if non-could be found.
     */
    private BestAllele searchBestAllele(final int sampleIndex, final int evidenceIndex, final boolean canBeReference, Optional<double[]> priorities) {
        final int alleleCount = alleles.numberOfAlleles();
        if (alleleCount == 0 || (alleleCount == 1 && referenceAlleleIndex == 0 && !canBeReference)) {
            return new BestAllele(sampleIndex, evidenceIndex, -1, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        }

        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        int bestAlleleIndex = canBeReference || referenceAlleleIndex != 0 ? 0 : 1;

        int secondBestIndex = 0;
        double bestLikelihood = sampleValues[bestAlleleIndex][evidenceIndex];
        double secondBestLikelihood = Double.NEGATIVE_INFINITY;

        for (int a = bestAlleleIndex + 1; a < alleleCount; a++) {
            if (!canBeReference && referenceAlleleIndex == a) {
                continue;
            }
            final double candidateLikelihood = sampleValues[a][evidenceIndex];
            if (candidateLikelihood > bestLikelihood) {
                secondBestIndex = bestAlleleIndex;
                bestAlleleIndex = a;
                secondBestLikelihood = bestLikelihood;
                bestLikelihood = candidateLikelihood;
            } else if (candidateLikelihood > secondBestLikelihood) {
                secondBestIndex = a;
                secondBestLikelihood = candidateLikelihood;
            }
        }

        if (priorities.isPresent() && bestLikelihood - secondBestLikelihood < getInformativeThreshold()) {
            double bestPriority = priorities.get()[bestAlleleIndex];
            double secondBestPriority = priorities.get()[secondBestIndex];
            for (int a = 0; a < alleleCount; a++) {
                final double candidateLikelihood = sampleValues[a][evidenceIndex];
                if (a == bestAlleleIndex || (!canBeReference && a == referenceAlleleIndex) || bestLikelihood - candidateLikelihood > getInformativeThreshold()) {
                    continue;
                }
                final double candidatePriority = priorities.get()[a];

                if (candidatePriority > bestPriority) {
                    secondBestIndex = bestAlleleIndex;
                    bestAlleleIndex = a;
                    secondBestPriority = bestPriority;
                    bestPriority = candidatePriority;
                } else if (candidatePriority > secondBestPriority) {
                    secondBestIndex = a;
                    secondBestPriority = candidatePriority;
                }
            }
        }

        bestLikelihood = sampleValues[bestAlleleIndex][evidenceIndex];
        secondBestLikelihood = secondBestIndex != bestAlleleIndex ? sampleValues[secondBestIndex][evidenceIndex] : Double.NEGATIVE_INFINITY;

        return new BestAllele(sampleIndex, evidenceIndex, bestAlleleIndex, bestLikelihood, secondBestLikelihood);
    }

    private BestAllele searchBestAllele(final int sampleIndex, final int evidenceIndex, final boolean canBeReference) {
        return searchBestAllele(sampleIndex, evidenceIndex, canBeReference, Optional.empty());
    }

    public void changeEvidence(final Map<EVIDENCE, EVIDENCE> evidenceReplacements) {
        final int sampleCount = samples.numberOfSamples();
        for (int s = 0; s < sampleCount; s++) {
            final List<EVIDENCE> sampleEvidence = evidenceBySampleIndex.get(s);
            final Object2IntMap<EVIDENCE> evidenceIndex = evidenceIndexBySampleIndex.get(s);
            final int sampleEvidenceCount = sampleEvidence.size();
            for (int r = 0; r < sampleEvidenceCount; r++) {
                final EVIDENCE evidence = sampleEvidence.get(r);
                final EVIDENCE replacement = evidenceReplacements.get(evidence);
                if (replacement == null) {
                    continue;
                }
                sampleEvidence.set(r, replacement);
                if (evidenceIndex != null) {
                    evidenceIndex.remove(evidence);
                    evidenceIndex.put(replacement, r);
                }
            }
        }
    }

    /**
     * Add alleles that are missing in the evidence-likelihoods collection giving all evidence a default
     * likelihood value.
     * @param candidateAlleles the potentially missing alleles.
     * @param defaultLikelihood the default evidence likelihood value for that allele.
     *
     * @return {@code true} iff the the evidence-likelihood collection was modified by the addition of the input alleles.
     *  So if all the alleles in the input collection were already present in the evidence-likelihood collection this method
     *  will return {@code false}.
     *
     * @throws IllegalArgumentException if {@code candidateAlleles} is {@code null} or there is more than
     * one missing allele that is a reference or there is one but the collection already has
     * a reference allele.
     */
    public boolean addMissingAlleles(final Collection<A> candidateAlleles, final double defaultLikelihood) {
        Utils.nonNull(candidateAlleles, "the candidateAlleles list cannot be null");
        if (candidateAlleles.isEmpty()) {
            return false;
        }
        final List<A> allelesToAdd = candidateAlleles.stream().filter(allele -> !alleles.containsAllele(allele)).collect(Collectors.toList());

        if (allelesToAdd.isEmpty()) {
            return false;
        }

        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = alleles.numberOfAlleles() + allelesToAdd.size();

        int referenceIndex = this.referenceAlleleIndex;

        @SuppressWarnings("unchecked")
        final List<A> newAlleles = ListUtils.union(alleles.asListOfAlleles(), allelesToAdd);
        alleles = new IndexedAlleleList<>(newAlleles);

        // if we previously had no reference allele, update the reference index if a reference allele is added
        // if we previously had a reference and try to add another, throw an exception
        final OptionalInt indexOfReferenceInAllelesToAdd = IntStream.range(0, allelesToAdd.size())
                .filter(n -> allelesToAdd.get(n).isReference()).findFirst();
        if (referenceIndex != MISSING_REF) {
            Utils.validateArg(!indexOfReferenceInAllelesToAdd.isPresent(), "there can only be one reference allele");
        } else if (indexOfReferenceInAllelesToAdd.isPresent()){
            referenceAlleleIndex = oldAlleleCount + indexOfReferenceInAllelesToAdd.getAsInt();
        }

        //copy old allele likelihoods and set new allele likelihoods to the default value
        for (int s = 0; s < samples.numberOfSamples(); s++) {
            final int sampleEvidenceCount = evidenceBySampleIndex.get(s).size();
            final double[][] newValuesBySampleIndex = Arrays.copyOf(valuesBySampleIndex[s], newAlleleCount);
            for (int a = oldAlleleCount; a < newAlleleCount; a++) {
                newValuesBySampleIndex[a] = new double[sampleEvidenceCount];
                if (defaultLikelihood != 0.0) {
                    Arrays.fill(newValuesBySampleIndex[a], defaultLikelihood);
                }
            }
            valuesBySampleIndex[s] = newValuesBySampleIndex;
        }
        return true;
    }

    /**
     * Group evidence into lists of evidence -- for example group by read name to force read pairs to support a single haplotype.
     *
     * Log Likelihoods are summed over all evidence
     * in a group, corresponding to an independent evidence assumption.  Since this container's likelihoods generally pertain to
     * sequencing only (and not sample prep etc) this is usually a good assumption.
     *
     * @param groupingFunction Attribute function for grouping evidence, for example GATKRead::getName
     * @param gather Transformation applied to collections of evidence with same value of groupingFunction.  For example, Fragment::new
     *               to construct a fragment out of a pair of reads with the same name
     *
     * @return a new AlleleLikelihoods based on the grouped, transformed evidence.
     */
    public <U, NEW_EVIDENCE_TYPE extends Locatable> AlleleLikelihoods<NEW_EVIDENCE_TYPE, A> groupEvidence(final Function<EVIDENCE, U> groupingFunction, final Function<List<EVIDENCE>, NEW_EVIDENCE_TYPE> gather) {
        final int sampleCount = samples.numberOfSamples();
        final double[][][] newLikelihoodValues = new double[sampleCount][][];
        final int alleleCount = alleles.numberOfAlleles();

        final List<List<NEW_EVIDENCE_TYPE>> newEvidenceBySampleIndex = new ArrayList<>(sampleCount);

        for (int s = 0; s < sampleCount; s++) {
            final List<List<EVIDENCE>> evidenceGroups = new ArrayList<>(sampleEvidence(s).stream()
                    .collect(Collectors.groupingBy(groupingFunction)).values());


            final int newEvidenceCount = evidenceGroups.size();

            final double[][] oldSampleValues = valuesBySampleIndex[s];
            newLikelihoodValues[s] = new double[alleleCount][newEvidenceCount];

            // For each old allele and read we update the new table keeping the maximum likelihood.
            for (int newEvidenceIndex = 0; newEvidenceIndex < newEvidenceCount; newEvidenceIndex++) {
                for (int a = 0; a < alleleCount; a++) {
                    for (final EVIDENCE evidence : evidenceGroups.get(newEvidenceIndex)) {
                        final int oldEvidenceIndex = evidenceIndex(s, evidence);
                        newLikelihoodValues[s][a][newEvidenceIndex] += oldSampleValues[a][oldEvidenceIndex];
                    }
                }
            }
            newEvidenceBySampleIndex.add(evidenceGroups.stream().map(gather).collect(Collectors.toList()));
        }

        // Finally we create the new read-likelihood
        final AlleleLikelihoods<NEW_EVIDENCE_TYPE, A> result = new AlleleLikelihoods<NEW_EVIDENCE_TYPE, A>(
                alleles,
                samples,
                newEvidenceBySampleIndex,
                newLikelihoodValues);

        result.isNaturalLog = this.isNaturalLog;
        return result;
    }

    /**
     * Perform marginalization from an allele set to another (smaller one) taking the maximum value
     * for each evidence in the original allele subset.
     *
     * @param newToOldAlleleMap map where the keys are the new alleles and the value list the original
     *                          alleles that correspond to the new one.
     * @return never {@code null}. The result will have the requested set of new alleles (keys in {@code newToOldAlleleMap}, and
     * the same set of samples and evidence as the original.
     *
     * @throws IllegalArgumentException is {@code newToOldAlleleMap} is {@code null} or contains {@code null} values,
     *  or its values contain reference to non-existing alleles in this evidence-likelihood collection. Also no new allele
     *  can have zero old alleles mapping nor two new alleles can make reference to the same old allele.
     */
    @SuppressWarnings({"unchecked", "rawtypes"})
    public <B extends Allele> AlleleLikelihoods<EVIDENCE, B> marginalize(final Map<B, List<A>> newToOldAlleleMap) {
        Utils.nonNull(newToOldAlleleMap);

        final B[] newAlleles = newToOldAlleleMap.keySet().toArray((B[]) new Allele[newToOldAlleleMap.size()]);
        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = newAlleles.length;

        // we get the index correspondence between new old -> new allele, -1 entries mean that the old
        // allele does not map to any new; supported but typically not the case.
        final int[] oldToNewAlleleIndexMap = oldToNewAlleleIndexMap(newToOldAlleleMap, oldAlleleCount, newAlleles);

        // We calculate the marginal likelihoods.
        final double[][][] newLikelihoodValues = marginalLikelihoods(oldAlleleCount, newAlleleCount, oldToNewAlleleIndexMap, null);

        final int sampleCount = samples.numberOfSamples();

        final List<List<EVIDENCE>> newEvidenceBySampleIndex = new ArrayList<>(sampleCount);

        for (int s = 0; s < sampleCount; s++) {
            newEvidenceBySampleIndex.add(new ArrayList<>(evidenceBySampleIndex.get(s)));
        }

        // Finally we create the new evidence-likelihood
        final AlleleLikelihoods<EVIDENCE, B> result = new AlleleLikelihoods<EVIDENCE, B>(
                new IndexedAlleleList(newAlleles),
                samples,
                newEvidenceBySampleIndex,
                newLikelihoodValues);
        result.isNaturalLog = isNaturalLog;
        return result;
    }


    /**
     * Perform marginalization from an allele set to another (smaller one) taking the maximum value
     * for each unit of evidence in the original allele subset.
     *
     * @param newToOldAlleleMap map where the keys are the new alleles and the value list the original
     *                          alleles that correspond to the new one.
     * @return never {@code null}. The result will have the requested set of new alleles (keys in {@code newToOldAlleleMap}, and
     * the same set of samples and evidence as the original.
     *
     * @param overlap if not {@code null}, only units of evidence that overlap the location (with unclipping) will be present in
     *                        the output evidence-collection.
     *
     * @throws IllegalArgumentException is {@code newToOldAlleleMap} is {@code null} or contains {@code null} values,
     *  or its values contain reference to non-existing alleles in this evidence-likelihood collection. Also no new allele
     *  can have zero old alleles mapping nor two new alleles can make reference to the same old allele.
     */
    public <B extends Allele> AlleleLikelihoods<EVIDENCE, B> marginalize(final Map<B, List<A>> newToOldAlleleMap, final Locatable overlap) {
        Utils.nonNull(newToOldAlleleMap, "the input allele mapping cannot be null");
        if (overlap == null) {
            return marginalize(newToOldAlleleMap);
        }

        @SuppressWarnings("unchecked")
        final B[] newAlleles = newToOldAlleleMap.keySet().toArray((B[]) new Allele[newToOldAlleleMap.size()]);
        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = newAlleles.length;

        // we get the index correspondence between new old -> new allele, -1 entries mean that the old
        // allele does not map to any new; supported but typically not the case.
        final int[] oldToNewAlleleIndexMap = oldToNewAlleleIndexMap(newToOldAlleleMap, oldAlleleCount, newAlleles);

        final int[][] evidenceToKeep = overlappingEvidenceIndicesBySampleIndex(overlap);

        final double[][][] newLikelihoodValues = marginalLikelihoods(oldAlleleCount, newAlleleCount, oldToNewAlleleIndexMap, evidenceToKeep);

        final int sampleCount = samples.numberOfSamples();

        @SuppressWarnings({"rawtypes","unchecked"})
        final List<List<EVIDENCE>> newEvidenceBySampleIndex = new ArrayList<>(sampleCount);

        for (int s = 0; s < sampleCount; s++) {
            final int[] sampleEvidenceToKeep = evidenceToKeep[s];
            final List<EVIDENCE> oldSampleEvidence = evidenceBySampleIndex.get(s);
            final int oldSampleEvidenceCount = oldSampleEvidence.size();
            final int newSampleEvidenceCount = sampleEvidenceToKeep.length;

            final List<EVIDENCE> newSampleEvidence = newSampleEvidenceCount == oldSampleEvidenceCount ? new ArrayList<>(oldSampleEvidence) :
                    Arrays.stream(sampleEvidenceToKeep).mapToObj(oldSampleEvidence::get).collect(Collectors.toList());
            newEvidenceBySampleIndex.add(newSampleEvidence);
        }

        // Finally we create the new evidence-likelihood
        final AlleleLikelihoods<EVIDENCE, B> result = new AlleleLikelihoods<>(new IndexedAlleleList<>(newAlleles), samples,
                newEvidenceBySampleIndex,
                newLikelihoodValues);
        result.isNaturalLog = isNaturalLog;
        return result;
    }

    private int[][] overlappingEvidenceIndicesBySampleIndex(final Locatable overlap) {
        if (overlap == null) {
            return null;
        }
        final int sampleCount = samples.numberOfSamples();
        final int[][] result = new int[sampleCount][];
        final IntArrayList buffer = new IntArrayList(200);

        for (int s = 0; s < sampleCount; s++) {
            buffer.clear();
            final List<EVIDENCE> sampleEvidence = evidenceBySampleIndex.get(s);
            final int sampleEvidenceCount = sampleEvidence.size();
            buffer.ensureCapacity(sampleEvidenceCount);
            for (int r = 0; r < sampleEvidenceCount; r++) {
                if (sampleEvidence.get(r).overlaps(overlap)) {
                    buffer.add(r);
                }
            }
            result[s] = buffer.toIntArray();
        }
        return result;
    }

    // Calculate the marginal likelihoods considering the old -> new allele index mapping.
    private double[][][] marginalLikelihoods(final int oldAlleleCount, final int newAlleleCount, final int[] oldToNewAlleleIndexMap, final int[][] evidenceToKeep) {

        final int sampleCount = samples.numberOfSamples();
        final double[][][] result = new double[sampleCount][][];

        for (int s = 0; s < sampleCount; s++) {
            final int sampleEvidenceCount = evidenceBySampleIndex.get(s).size();
            final double[][] oldSampleValues = valuesBySampleIndex[s];
            final int[] sampleEvidenceToKeep = evidenceToKeep == null || evidenceToKeep[s].length == sampleEvidenceCount ? null : evidenceToKeep[s];
            final int newSampleEvidenceCount = sampleEvidenceToKeep == null ? sampleEvidenceCount : sampleEvidenceToKeep.length;
            final double[][] newSampleValues = result[s] = new double[newAlleleCount][newSampleEvidenceCount];
            // We initiate all likelihoods to -Inf.
            for (int a = 0; a < newAlleleCount; a++) {
                Arrays.fill(newSampleValues[a], Double.NEGATIVE_INFINITY);
            }
            // For each old allele and unit of evidence we update the new table keeping the maximum likelihood.
            for (int r = 0; r < newSampleEvidenceCount; r++) {
                for (int a = 0; a < oldAlleleCount; a++) {
                    final int oldEvidenceIndex = newSampleEvidenceCount == sampleEvidenceCount ? r : sampleEvidenceToKeep[r];
                    final int newAlleleIndex = oldToNewAlleleIndexMap[a];
                    if (newAlleleIndex == -1) {
                        continue;
                    }
                    final double likelihood = oldSampleValues[a][oldEvidenceIndex];
                    if (likelihood > newSampleValues[newAlleleIndex][r]) {
                        newSampleValues[newAlleleIndex][r] = likelihood;
                    }
                }
            }
        }
        return result;
    }

    // calculates an old to new allele index map array.
    private <B extends Allele> int[] oldToNewAlleleIndexMap(final Map<B, List<A>> newToOldAlleleMap, final int oldAlleleCount, final B[] newAlleles) {
        Arrays.stream(newAlleles).forEach(Utils::nonNull);
        Utils.containsNoNull(newToOldAlleleMap.values(), "no new allele list can be null");
        newToOldAlleleMap.values().stream().forEach(oldList -> Utils.containsNoNull(oldList,"old alleles cannot be null"));

        final int[] oldToNewAlleleIndexMap = new int[oldAlleleCount];
        Arrays.fill(oldToNewAlleleIndexMap, -1);  // -1 indicate that there is no new allele that make reference to that old one.

        for (int newIndex = 0; newIndex < newAlleles.length; newIndex++) {
            final B newAllele = newAlleles[newIndex];
            for (final A oldAllele : newToOldAlleleMap.get(newAllele)) {
                final int oldAlleleIndex = indexOfAllele(oldAllele);
                if (oldAlleleIndex == -1) {
                    throw new IllegalArgumentException("missing old allele " + oldAllele + " in likelihood collection ");
                }
                if (oldToNewAlleleIndexMap[oldAlleleIndex] != -1) {
                    throw new IllegalArgumentException("collision: two new alleles make reference to the same old allele");
                }
                oldToNewAlleleIndexMap[oldAlleleIndex] = newIndex;
            }
        }
        return oldToNewAlleleIndexMap;
    }

    /**
     * Add more evidence to the collection.
     *
     * @param evidenceBySample evidence to add.
     * @param initialLikelihood the likelihood for the new entries.
     *
     * @throws IllegalArgumentException if {@code evidenceBySample} is {@code null} or {@code evidenceBySample} contains
     *  {@code null} evidence, or {@code evidenceBySample} contains evidence already present in the evidence-likelihood
     *  collection.
     */
    public void addEvidence(final Map<String,List<EVIDENCE>> evidenceBySample, final double initialLikelihood) {
        for (final Map.Entry<String,List<EVIDENCE>> entry : evidenceBySample.entrySet()) {
            final String sample = entry.getKey();
            final List<EVIDENCE> newSampleEvidence = entry.getValue();
            final int sampleIndex = samples.indexOfSample(sample);

            if (sampleIndex == -1) {
                throw new IllegalArgumentException("input sample " + sample + " is not part of the evidence-likelihoods collection");
            }

            if (newSampleEvidence == null || newSampleEvidence.isEmpty()) {
                continue;
            }

            final int sampleEvidenceCount = evidenceBySampleIndex.get(sampleIndex).size();
            final int newSampleEvidenceCount = sampleEvidenceCount + newSampleEvidence.size();

            appendEvidence(newSampleEvidence, sampleIndex);
            extendsLikelihoodArrays(initialLikelihood, sampleIndex, sampleEvidenceCount, newSampleEvidenceCount);
        }
    }

    // Extends the likelihood arrays-matrices.
    private void extendsLikelihoodArrays(final double initialLikelihood, final int sampleIndex, final int sampleEvidenceCount, final int newSampleEvidenceCount) {
        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        final int alleleCount = alleles.numberOfAlleles();
        for (int a = 0; a < alleleCount; a++) {
            sampleValues[a] = Arrays.copyOf(sampleValues[a], newSampleEvidenceCount);
        }
        if (initialLikelihood != 0.0) // the default array new value.
        {
            for (int a = 0; a < alleleCount; a++) {
                Arrays.fill(sampleValues[a], sampleEvidenceCount, newSampleEvidenceCount, initialLikelihood);
            }
        }
    }

    // Append the new evidence reference into the structure per-sample.
    private void appendEvidence(final List<EVIDENCE> newSampleEvidence, final int sampleIndex) {

        final int nextEvidenceIndex = evidenceBySampleIndex.get(sampleIndex).size();

         evidenceBySampleIndex.get(sampleIndex).addAll(newSampleEvidence);

        final Object2IntMap<EVIDENCE> sampleEvidenceIndex = evidenceIndexBySampleIndex.get(sampleIndex);
        for (final EVIDENCE newEvidence : newSampleEvidence) {
            //    if (sampleEvidenceIndex.containsKey(newEvidence)) // might be worth handle this without exception (ignore the evidence?) but in practice should never be the case.
            //        throw new IllegalArgumentException("you cannot add evidence already in evidence-likelihood collection");
            if (sampleEvidenceIndex != null ) {
                sampleEvidenceIndex.put(newEvidence, nextEvidenceIndex);
            }
        }
    }

    /**
     * Adds the non-reference allele to the evidence-likelihood collection setting each evidence likelihood to the second
     * best found (or best one if only one allele has likelihood).
     *
     * <p>Nothing will happen if the evidence-likelihoods collection already includes the non-ref allele</p>
     *
     * <p>
     *     <i>Implementation note: even when strictly speaking we do not need to demand the calling code to pass
     *     the reference the non-ref allele, we still demand it in order to lead the
     *     the calling code to use the right generic type for this likelihoods
     *     collection {@link Allele}.</i>
     * </p>
     *
     * @param nonRefAllele the non-ref allele.
     *
     * @throws IllegalArgumentException if {@code nonRefAllele} is anything but the designated &lt;NON_REF&gt;
     * symbolic allele {@link Allele#NON_REF_ALLELE}.
     */
    public void addNonReferenceAllele(final A nonRefAllele) {
        Utils.nonNull(nonRefAllele, "non-ref allele cannot be null");
        if (!nonRefAllele.equals(Allele.NON_REF_ALLELE)) {
            throw new IllegalArgumentException("the non-ref allele is not valid");
        } else if (alleles.containsAllele(nonRefAllele)) {
            return;
        } else if (addMissingAlleles(Collections.singleton(nonRefAllele), Double.NEGATIVE_INFINITY)) {
            updateNonRefAlleleLikelihoods();
        }
    }

    /**
     * Updates the likelihoods of the non-ref allele, if present, considering all non-symbolic alleles avaialble.
     */
    public void updateNonRefAlleleLikelihoods() {
        updateNonRefAlleleLikelihoods(alleles);
    }

    /**
     * Updates the likelihood of the NonRef allele (if present) based on the likelihoods of a set of non-symbolic
     * <p>
     *     This method does
     * </p>
     *
     *
     * @param allelesToConsider
     */
    @SuppressWarnings("unchecked")  // for the cast (A) Allele.NON_REF_ALLELE below
    public void updateNonRefAlleleLikelihoods(final AlleleList<A> allelesToConsider) {
        final int nonRefAlleleIndex = indexOfAllele((A) Allele.NON_REF_ALLELE);
        if ( nonRefAlleleIndex < 0) {
            return;
        }
        final int alleleCount = alleles.numberOfAlleles();
        final int nonSymbolicAlleleCount = alleleCount - 1;
        // likelihood buffer reused across evidence:
        final double[] qualifiedAlleleLikelihoods = new double[nonSymbolicAlleleCount];
        final Median medianCalculator = new Median();
        for (int s = 0; s < samples.numberOfSamples(); s++) {
            final double[][] sampleValues = valuesBySampleIndex[s];
            final int evidenceCount = sampleValues[0].length;
            for (int r = 0; r < evidenceCount; r++) {
                final BestAllele bestAllele = searchBestAllele(s, r, true);
                int numberOfQualifiedAlleleLikelihoods = 0;
                for (int i = 0; i < alleleCount; i++) {
                    final double alleleLikelihood = sampleValues[i][r];
                    if (i != nonRefAlleleIndex && alleleLikelihood < bestAllele.likelihood
                            && !Double.isNaN(alleleLikelihood) && allelesToConsider.indexOfAllele(alleles.getAllele(i)) != -1) {
                        qualifiedAlleleLikelihoods[numberOfQualifiedAlleleLikelihoods++] = alleleLikelihood;
                    }
                }
                final double nonRefLikelihood = medianCalculator.evaluate(qualifiedAlleleLikelihoods, 0, numberOfQualifiedAlleleLikelihoods);
                // when the median is NaN that means that all applicable likekihoods are the same as the best
                // so the evidence is not informative at all given the existing alleles. Unless there is only one (or zero) concrete
                // alleles with give the same (the best) likelihood to the NON-REF. When there is only one (or zero) concrete
                // alleles we set the NON-REF likelihood to NaN.
                sampleValues[nonRefAlleleIndex][r] = !Double.isNaN(nonRefLikelihood) ? nonRefLikelihood
                        : nonSymbolicAlleleCount <= 1 ? Double.NaN : bestAllele.likelihood;
            }
        }
    }

    /**
     * Returns the collection of best allele estimates for the evidence based on the evidence-likelihoods.
     * "Ties" where the ref likelihood is within {@code AlleleLikelihoods.INFORMATIVE_THRESHOLD} of the greatest likelihood
     * are broken by the {@code tieBreakingPriority} function.
     *
     * @throws IllegalStateException if there is no alleles.
     *
     * @return never {@code null}, one element per unit fo evidence in the evidence-likelihoods collection.
     */
    public Collection<BestAllele> bestAllelesBreakingTies(final Function<A, Double> tieBreakingPriority) {
        return IntStream.range(0, numberOfSamples()).boxed().flatMap(n -> bestAllelesBreakingTies(n, tieBreakingPriority).stream()).collect(Collectors.toList());
    }

    /**
     * Default version where ties are broken in favor of the reference allele
     */
    public Collection<BestAllele> bestAllelesBreakingTies() {
        return IntStream.range(0, numberOfSamples()).boxed().flatMap(n -> bestAllelesBreakingTies(n).stream()).collect(Collectors.toList());
    }

    /**
     * Returns the collection of best allele estimates for one sample's evidence based on the evidence-likelihoods.
     * "Ties" where the ref likelihood is within {@code AlleleLikelihoods.INFORMATIVE_THRESHOLD} of the greatest likelihood
     * are broken by the {@code tieBreakingPriority} function.
     *
     * @throws IllegalStateException if there is no alleles.
     *
     * @return never {@code null}, one element per unit of evidence in the evidence-likelihoods collection.
     */
    public Collection<BestAllele> bestAllelesBreakingTies(final String sample, final Function<A, Double> tieBreakingPriority) {
        final int sampleIndex = indexOfSample(sample);
        return bestAllelesBreakingTies(sampleIndex, tieBreakingPriority);
    }

    /**
     * Default version where ties are broken in favor of the reference allele
     */
    public Collection<BestAllele> bestAllelesBreakingTies(final String sample) {
        final int sampleIndex = indexOfSample(sample);
        return bestAllelesBreakingTies(sampleIndex);
    }

    /**
     * Returns the collection of best allele estimates for one sample's evidence based on the evidence-likelihoods.
     * "Ties" where the ref likelihood is within {@code AlleleLikelihoods.INFORMATIVE_THRESHOLD} of the greatest likelihood
     * are broken by the {@code tieBreakingPriority} function.
     *
     * @throws IllegalStateException if there is no alleles.
     *
     * @return never {@code null}, one element per unit of evidence in the evidence-likelihoods collection.
     */
    private Collection<BestAllele> bestAllelesBreakingTies(final int sampleIndex, final Function<A, Double> tieBreakingPriority) {
        Utils.validIndex(sampleIndex, numberOfSamples());

        //TODO: this currently just does ref vs alt.  Really we want CIGAR complexity.
        final Optional<double[]> priorities = alleles == null ? Optional.empty() :
                Optional.of(alleles.asListOfAlleles().stream().mapToDouble(tieBreakingPriority::apply).toArray());

        final int evidenceCount = evidenceBySampleIndex.get(sampleIndex).size();
        final List<BestAllele> result = new ArrayList<>(evidenceCount);
        for (int r = 0; r < evidenceCount; r++) {
            result.add(searchBestAllele(sampleIndex, r, true, priorities));
        }

        return result;
    }

    /**
     * Default version where ties are broken in favor of the reference allele
     */
    private Collection<BestAllele> bestAllelesBreakingTies(final int sampleIndex) {
        return bestAllelesBreakingTies(sampleIndex, a -> a.isReference() ? 1.0 : 0);
    }

    /**
     * Returns evidence stratified by best allele.
     * @param sampleIndex the target sample.
     * @return never {@code null}, perhaps empty.
     */
    protected Map<A,List<EVIDENCE>> evidenceByBestAlleleMap(final int sampleIndex) {
        Utils.validIndex(sampleIndex, numberOfSamples());
        final int alleleCount = alleles.numberOfAlleles();
        final int sampleEvidenceCount = evidenceBySampleIndex.get(sampleIndex).size();
        final Map<A,List<EVIDENCE>> result = new LinkedHashMap<>(alleleCount);
        for (int a = 0; a < alleleCount; a++) {
            result.put(alleles.getAllele(a), new ArrayList<>(sampleEvidenceCount));
        }
        evidenceByBestAlleleMap(sampleIndex,result);
        return result;
    }

    /**
     * Returns evidence stratified by best allele.
     * @return never {@code null}, perhaps empty.
     */
    @VisibleForTesting
    Map<A,List<EVIDENCE>> evidenceByBestAlleleMap() {
        final int alleleCount = alleles.numberOfAlleles();
        final Map<A,List<EVIDENCE>> result = new LinkedHashMap<>(alleleCount);
        final int totalEvidenceCount = evidenceCount();
        for (int a = 0; a < alleleCount; a++) {
            result.put(alleles.getAllele(a), new ArrayList<>(totalEvidenceCount));
        }
        final int sampleCount = samples.numberOfSamples();
        for (int s = 0; s < sampleCount; s++) {
            evidenceByBestAlleleMap(s, result);
        }
        return result;
    }

    private void evidenceByBestAlleleMap(final int sampleIndex, final Map<A, List<EVIDENCE>> result) {
        final int evidenceCount = evidenceBySampleIndex.get(sampleIndex).size();

        for (int r = 0; r < evidenceCount; r++) {
            final BestAllele bestAllele = searchBestAllele(sampleIndex,r,true);
            if (!bestAllele.isInformative()) {
                continue;
            }
            result.get(bestAllele.allele).add(bestAllele.evidence);
        }
    }

    /**
     * Returns the index of evidence within a sample evidence-likelihood sub collection.
     * @param sampleIndex the sample index.
     * @param evidence the query evidence.
     * @return -1 if there is no such evidence in that sample, 0 or greater otherwise.
     */
    @VisibleForTesting
    int evidenceIndex(final int sampleIndex, final EVIDENCE evidence) {
        return evidenceIndexBySampleIndex(sampleIndex).getOrDefault(evidence, -1);
    }

    /**
     * Returns the total count of evidence in the evidence-likelihood collection.
     */
    public int evidenceCount() {
        return evidenceBySampleIndex.stream().mapToInt(List::size).sum();
    }

    /**
     * Returns the quantity of evidence that belongs to a sample in the evidence-likelihood collection.
     * @param sampleIndex the query sample index.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is not a valid sample index.
     * @return 0 or greater.
     */
    public int sampleEvidenceCount(final int sampleIndex) {
        Utils.validIndex(sampleIndex, samples.numberOfSamples());
        return evidenceBySampleIndex.get(sampleIndex).size();
    }

    /**
     * Removes evidence that does not overlap certain genomic location.
     *
     * <p>
     *     This method modifies the current evidence-likelihoods collection.
     * </p>
     *
     * @param location the target location.
     *
     * @throws IllegalArgumentException the location cannot be {@code null} nor unmapped.
     */
    public void filterToOnlyOverlappingEvidence(final SimpleInterval location) {
        Utils.nonNull(location, "the location cannot be null");

        final int sampleCount = samples.numberOfSamples();

        final int alleleCount = alleles.numberOfAlleles();
        for (int s = 0; s < sampleCount; s++) {
            final List<EVIDENCE> evidenceToRemove = evidenceBySampleIndex.get(s).stream()
                    .filter(r -> !r.overlaps(location))
                    .collect(Collectors.toList());
            removeSampleEvidence(s, evidenceToRemove, alleleCount);
        }
    }

    protected double maximumLikelihoodOverAllAlleles(final int sampleIndex, final int evidenceIndex) {
        double result = Double.NEGATIVE_INFINITY;
        final int alleleCount = alleles.numberOfAlleles();
        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        for (int a = 0; a < alleleCount; a++) {
            if (sampleValues[a][evidenceIndex] > result) {
                result = sampleValues[a][evidenceIndex];
            }
        }
        return result;
    }

    /**
     * Contains information about the best allele for a unit of evidence.
     */
    public final class BestAllele {

        /**
         * Null if there is no possible match (no allele?).
         */
        public final A allele;

        /**
         * The containing sample.
         */
        public final String sample;

        /**
         * The query evidence.
         */
        public final EVIDENCE evidence;

        /**
         * If allele != null, the indicates the likelihood of the evidence.
         */
        public final double likelihood;

        /**
         * Confidence that the evidence actually was generated under that likelihood.
         * This is equal to the difference between this and the second best allele match.
         */
        public final double confidence;

        private BestAllele(final int sampleIndex, final int evidenceIndex, final int bestAlleleIndex,
                           final double likelihood, final double secondBestLikelihood) {
            allele = bestAlleleIndex == -1 ? null : alleles.getAllele(bestAlleleIndex);
            this.likelihood = likelihood;
            sample = samples.getSample(sampleIndex);
            evidence = evidenceBySampleIndex.get(sampleIndex).get(evidenceIndex);
            confidence = likelihood == secondBestLikelihood ? 0 : likelihood - secondBestLikelihood;
        }

        public boolean isInformative() {
            return confidence > LOG_10_INFORMATIVE_THRESHOLD;
        }
    }

    // Requires that the collection passed iterator can remove elements, and it can be modified.
    public void removeSampleEvidence(final int sampleIndex, final Collection<EVIDENCE> evidenceToRemove, final int alleleCount) {
        if (evidenceToRemove.isEmpty()) {
            return;
        }

        final List<EVIDENCE> sampleEvidence = evidenceBySampleIndex.get(sampleIndex);

        final Set<EVIDENCE> evidenceToRemoveSet = new HashSet<>(evidenceToRemove);
        final int[] evidenceIndicesToKeep = IntStream.range(0, sampleEvidence.size()).filter(r -> !evidenceToRemoveSet.contains(sampleEvidence.get(r))).toArray();

        // Nothing to remove we just finish here.
        if (evidenceIndicesToKeep.length == sampleEvidence.size()) {
            return;
        }

        final List<EVIDENCE> newSampleEvidence = Arrays.stream(evidenceIndicesToKeep).mapToObj(sampleEvidence::get).collect(Collectors.toList());
        final Object2IntMap<EVIDENCE> indexByEvidence = evidenceIndexBySampleIndex(sampleIndex);
        indexByEvidence.clear();
        new IndexRange(0, newSampleEvidence.size()).forEach(r -> indexByEvidence.put(newSampleEvidence.get(r), r));

        // Then we skim out the likelihoods of the removed evidence.
        final double[][] oldSampleValues = valuesBySampleIndex[sampleIndex];
        final double[][] newSampleValues = new double[alleleCount][newSampleEvidence.size()];

        for (int a = 0; a < alleleCount; a++) {
            for (int r = 0; r < evidenceIndicesToKeep.length; r++) {
                newSampleValues[a][r] = oldSampleValues[a][evidenceIndicesToKeep[r]];
            }
        }

        valuesBySampleIndex[sampleIndex] = newSampleValues;
        evidenceBySampleIndex.set(sampleIndex, newSampleEvidence);
    }

    /**
     * Removes those read that the best possible likelihood given any allele is just too low.
     *
     * <p>
     *     This is determined by a maximum error per read-base against the best likelihood possible.
     * </p>
     *
     * @param log10MinTrueLikelihood Function that returns the minimum likelihood that the best allele for a unit of evidence must have
     * @throws IllegalStateException is not supported for read-likelihood that do not contain alleles.
     *
     * @throws IllegalArgumentException if {@code maximumErrorPerBase} is negative.
     */
    public void filterPoorlyModeledEvidence(final ToDoubleFunction<EVIDENCE> log10MinTrueLikelihood) {
        Utils.validateArg(alleles.numberOfAlleles() > 0, "unsupported for read-likelihood collections with no alleles");

        final int numberOfSamples = samples.numberOfSamples();
        for (int s = 0; s < numberOfSamples; s++) {
            final List<EVIDENCE> sampleEvidence = evidenceBySampleIndex.get(s);
            final List<EVIDENCE> evidenceToRemove = new ArrayList<EVIDENCE>(sampleEvidence.size());

            final int numberOfEvidence = sampleEvidence.size();
            for (int e = 0; e < numberOfEvidence; e++) {
                if (maximumLikelihoodOverAllAlleles(s, e) < log10MinTrueLikelihood.applyAsDouble(sampleEvidence.get(e))) {
                    evidenceToRemove.add(sampleEvidence.get(e));
                }
            }
            removeSampleEvidence(s, evidenceToRemove, alleles.numberOfAlleles());
        }
    }


    private Object2IntMap<EVIDENCE> evidenceIndexBySampleIndex(final int sampleIndex) {
        if (evidenceIndexBySampleIndex.get(sampleIndex) == null) {
            final List<EVIDENCE> sampleEvidence = evidenceBySampleIndex.get(sampleIndex);
            final int sampleEvidenceCount = sampleEvidence.size();
            evidenceIndexBySampleIndex.set(sampleIndex, new Object2IntOpenHashMap<>(sampleEvidenceCount));
            for (int r = 0; r < sampleEvidenceCount; r++) {
                evidenceIndexBySampleIndex.get(sampleIndex).put(sampleEvidence.get(r), r);
            }
        }
        return evidenceIndexBySampleIndex.get(sampleIndex);
    }


    /**
     * Collect a map stratified per-sample of the base pileups at the provided Location
     * NOTE: Since we shouldn't need to use the pileup if we have more reliable liklihoods, we want to discourage their use
     *
     * @param loc reference location to construct pileups for
     * @return
     */
    public Map<String, List<PileupElement>> getStratifiedPileups(final Locatable loc) {
        return null;
    }

    /**
     * Implements a likelihood matrix per sample given its index.
     */
    private final class SampleMatrix implements LikelihoodMatrix<EVIDENCE, A> {

        private final int sampleIndex;

        private SampleMatrix(final int sampleIndex) {
            this.sampleIndex = sampleIndex;
        }

        public List<EVIDENCE> evidence() {
            return sampleEvidence(sampleIndex);
        }

        @Override
        public List<A> alleles() {
            return AlleleLikelihoods.this.alleles();
        }

        @Override
        public void set(final int alleleIndex, final int evidenceIndex, final double value) {
            Utils.validIndex(alleleIndex, valuesBySampleIndex[sampleIndex].length);
            Utils.validIndex(evidenceIndex, valuesBySampleIndex[sampleIndex][alleleIndex].length);
            valuesBySampleIndex[sampleIndex][alleleIndex][evidenceIndex] = value;
        }

        @Override
        public double get(final int alleleIndex, final int evidenceIndex) {
            Utils.validIndex(alleleIndex, valuesBySampleIndex[sampleIndex].length);
            Utils.validIndex(evidenceIndex, valuesBySampleIndex[sampleIndex][alleleIndex].length);
            return valuesBySampleIndex[sampleIndex][alleleIndex][evidenceIndex];
        }

        @Override
        public int indexOfAllele(final A allele) {
            Utils.nonNull(allele);
            return AlleleLikelihoods.this.indexOfAllele(allele);
        }

        @Override
        public int indexOfEvidence(final EVIDENCE evidence) {
            Utils.nonNull(evidence);
            return AlleleLikelihoods.this.evidenceIndex(sampleIndex, evidence);
        }

        @Override
        public int numberOfAlleles() {
            return alleles.numberOfAlleles();
        }

        @Override
        public int evidenceCount() {
            return evidenceBySampleIndex.get(sampleIndex).size();
        }

        @Override
        public A getAllele(final int alleleIndex) {
            return AlleleLikelihoods.this.getAllele(alleleIndex);
        }

        @Override
        public EVIDENCE getEvidence(final int evidenceIndex) {
            final List<EVIDENCE> sampleEvidence = evidenceBySampleIndex.get(sampleIndex);
            Utils.validIndex(evidenceIndex, sampleEvidence.size());
            return sampleEvidence.get(evidenceIndex);
        }

        @Override
        public void copyAlleleLikelihoods(final int alleleIndex, final double[] dest, final int offset) {
            Utils.nonNull(dest);
            Utils.validIndex(alleleIndex, valuesBySampleIndex[sampleIndex].length);
            System.arraycopy(valuesBySampleIndex[sampleIndex][alleleIndex], 0, dest, offset, evidenceCount());
        }
    }
}
