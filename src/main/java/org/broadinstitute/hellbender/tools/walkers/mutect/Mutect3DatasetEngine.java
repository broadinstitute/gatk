package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.AssemblyComplexity;
import org.broadinstitute.hellbender.tools.walkers.annotator.FeaturizedReadSets;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReferenceBases;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Mutect3DatasetEngine implements AutoCloseable {

    public static final int CAPACITY = 100000;

    private enum VariantType {
        SNV, INSERTION, DELETION
    }

    private enum Label {
        ARTIFACT, VARIANT, UNLABELED, IGNORE
    }

    // number of features for each vectorized read
    private static final int FEATURES_PER_READ = FeaturizedReadSets.FEATURES_PER_READ;

    // number of additional variant features for assembly complexity, reference context
    private static final int NUM_EXTRA_FEATURES = 9;

    // threshold of negative log-10 population allele frequency to consider something an artifact for the purposes of training data
    // we want to be really sure we don't get germline variants
    // TODO: is this really necessary?
    private static final double RARE_POPAF_THRESHOLD = 5.9;

    // very cautious threshold of negative log-10 population allele frequency to consider something germline for training data.
    // There are so many germline variants we can be wasteful!
    private static final double COMMON_POPAF_THRESHOLD = 1;

    // below this tumor log odds we don't consider it an artifact, just a sequencing error
    private static final double TLOD_THRESHOLD = 6.0;

    private final int maxRefCount;
    private final int maxAltCount;

    // TODO: is this necessary?
    private static final int MIN_REF = 5;

    private final PrintWriter printWriter;

    // number of nonartifact data to keep for each artifact datum
    private final int nonArtifactPerArtifact;

    // are we generating dataset for training a model or for filtering calls with a pre-trained model?
    private final boolean trainingMode;

    private final Set<String> normalSamples;

    // simple method to balance data: for each k-alt-read artifact there are
    // nonArtifactPerArtifact (downsampled) k-alt-read non-artifacts.
    private final EnumMap<VariantType, ArrayBlockingQueue<Integer>> unmatchedArtifactAltCounts;


    public Mutect3DatasetEngine(final File datasetFile, final boolean trainingMode, final int maxRefCount, final int maxAltCount, final int nonArtifactPerArtifact, final Set<String> normalSamples) {
        try {
            printWriter = new PrintWriter(new FileWriter(Utils.nonNull(datasetFile)));
        } catch (IOException ex) {
            throw new UserException.BadInput("Could not create dataset file writer");
        }

        this.normalSamples = Utils.nonNull(normalSamples);
        this.trainingMode = trainingMode;
        this.nonArtifactPerArtifact = nonArtifactPerArtifact;
        this.maxRefCount = maxRefCount;
        this.maxAltCount = maxAltCount;

        unmatchedArtifactAltCounts = new EnumMap<>(VariantType.class);
        for (final VariantType type : VariantType.values()) {
            unmatchedArtifactAltCounts.put(type, new ArrayBlockingQueue<>(CAPACITY));
        }
    }

    // add one datum per alt allele
    public void addData(final ReferenceContext ref, final VariantContext vc, final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                         final AlleleLikelihoods<Fragment, Haplotype> logFragmentLikelihoods) {
        final String refBases = ReferenceBases.annotate(ref, vc);
        final String refAllele = vc.getReference().getBaseString();
        final String contig = vc.getContig();
        final int position = vc.getStart();
        final Set<String> tumorSamples = likelihoods.samples().stream().filter(sample -> !normalSamples.contains(sample)).collect(Collectors.toSet());
        final int numAlt = vc.getNAlleles() - 1;


        // the variant has already been annotated, so we have POPAF and AD
        final double[] popafs = VariantContextGetters.getAttributeAsDoubleArray(vc, GATKVCFConstants.POPULATION_AF_KEY);
        //final double[] altPopulationAFs = MathUtils.applyToArray(popafs, x -> Math.pow(10, -x ));
        final double[] tumorLods = Mutect2FilteringEngine.getTumorLogOdds(vc);
        final int[] tumorADs = sumADsOverSamples(vc, tumorSamples);
        final int[] normalADs = sumADsOverSamples(vc, normalSamples);
        final int tumorDepth = (int) MathUtils.sum(tumorADs);
        final int normalDepth = (int) MathUtils.sum(normalADs);
        final boolean hasNormal = normalDepth > 0;

        final List<Label> labels = new ArrayList<>(numAlt);
        final Map<Allele, Integer> altDownsampleMap= new HashMap<>();

        for (int n = 0; n < numAlt; n++) {
            final double tumorAF = tumorADs[n+1] / ((double) tumorDepth);
            final double normalAF = hasNormal ? normalADs[n+1] / ((double) normalDepth) : 0.0;
            final String altAllele = vc.getAlternateAllele(n).getBaseString();
            final int diff = altAllele.length() - refAllele.length();
            final VariantType type = diff == 0 ? VariantType.SNV : ( diff > 0 ? VariantType.INSERTION : VariantType.DELETION);

            if (trainingMode) {
                final ArrayBlockingQueue<Integer> unmatchedQueue = unmatchedArtifactAltCounts.get(type);
                final boolean likelySeqError = tumorLods[n] < TLOD_THRESHOLD;
                final boolean likelyGermline = hasNormal && normalAF > 0.2;

                // extremely strict criteria because there are so many germline variants we can afford to waste a lot
                final boolean definiteGermline = !likelySeqError && popafs[n] < COMMON_POPAF_THRESHOLD &&
                        tumorAF > 0.35 && (!hasNormal || normalAF > 0.35);

                // low AF in tumor and normal, rare in population implies artifact
                if  (!(likelySeqError || likelyGermline) && tumorAF < 0.2 && popafs[n] > RARE_POPAF_THRESHOLD) {

                    if (unmatchedQueue.size() > 0.9*CAPACITY) { // this should rarely come up
                        labels.add(Label.IGNORE);
                    } else {
                        labels.add(Label.ARTIFACT);
                        unmatchedQueue.addAll(Collections.nCopies(nonArtifactPerArtifact, tumorADs[n + 1]));
                    }
                } else if (definiteGermline && !unmatchedQueue.isEmpty()) {
                    // high AF in tumor and normal, common in population implies germline, which we downsample
                    labels.add(Label.VARIANT);
                    altDownsampleMap.put(vc.getAlternateAllele(n), unmatchedQueue.poll());
                } else if (tumorLods[n] > 4.0 && tumorAF < 0.3) {
                    labels.add(Label.UNLABELED);
                } else {
                    labels.add(Label.IGNORE);
                }
            } else {
                labels.add(Label.UNLABELED);
            }
        }

        Utils.validate(labels.size() == numAlt, "We have not labeled every alt, or have labeled too much");
        // we check this later allele by allele, but we can save a lot of compute here if we realize no alt allele yields training data
        if (trainingMode && labels.stream().allMatch(label -> label == Label.IGNORE)) {
            return;
        }

        // haplotype equivalence counts, haplotype complexity, haplotype dominance
        final Triple<int[], int[], double[]> assemblyComplexity = AssemblyComplexity.annotate(vc, logFragmentLikelihoods);

        // TODO: for now we don't really need normal reads
        final List<List<List<Integer>>> normalReadVectorsByAllele =  FeaturizedReadSets.getReadVectors(vc, normalSamples, likelihoods, logFragmentLikelihoods, maxRefCount, maxAltCount);
        final List<List<List<Integer>>> tumorReadVectorsByAllele =  FeaturizedReadSets.getReadVectors(vc, tumorSamples, likelihoods, logFragmentLikelihoods, maxRefCount, maxAltCount, altDownsampleMap);

        // ref and alt reads have already been downsampled by the read featurizer
        final List<List<Integer>> tumorRefReads = tumorReadVectorsByAllele.get(0);
        final List<List<Integer>> normalRefReads = normalReadVectorsByAllele.get(0);

        for (int n = 0; n < numAlt; n++) {
            if (labels.get(n) ==  Label.IGNORE) {
                continue;
            }

            final String altAllele = vc.getAlternateAllele(n).getBaseString();
            final List<Double> variantFeatureVector = variantFeatures(n, assemblyComplexity, refBases);
            final List<List<Integer>> tumorAltReads = tumorReadVectorsByAllele.get(n+1);
            final List<List<Integer>> normalAltReads = normalReadVectorsByAllele.get(n+1);

            printWriter.println(labels.get(n).toString());
            printWriter.printf("%s:%d,%s->%s%n", contig, position, refAllele, altAllele);
            printWriter.println(refBases);
            printWriter.println(numberString(variantFeatureVector, "%.2f", " "));
            //printWriter.printf("%d %d %d %d%n", tumorRefReads.size(), tumorAltReads.size(), normalRefReads.size(), normalAltReads.size());
            // zeroes because we currently don't use normal read vectors
            printWriter.printf("%d %d %d %d%n", tumorRefReads.size(), tumorAltReads.size(), 0, 0);

            tumorRefReads.forEach(r -> printWriter.println(integerString(r)));
            tumorAltReads.forEach(r -> printWriter.println(integerString(r)));
            //normalRefReads.forEach(r -> printWriter.print(numberString(r)));
            //normalAltReads.forEach(r -> printWriter.print(numberString(r)));
            printWriter.printf("%d %d %d %d%n", tumorDepth, tumorADs[n+1], normalDepth, normalADs[n+1]);  // pre-downsampling counts for normal artifact model
            }
    }

    private String integerString(final List<Integer> numbers) {
        return numberString(numbers, "%d", " ");
    }

    private String numberString(final List<? extends Number> numbers, final String formatString, final String separator) {
        final boolean decimal = formatString.endsWith("f");
        return numbers.stream().map(x -> String.format(formatString, decimal ? x.floatValue() : x)).collect(Collectors.joining(separator));
    }

    private List<Double> variantFeatures(final int altAlleleIndex, Triple<int[], int[], double[]> assemblyComplexity, final String refBases) {
        final int[] haplotypeEquivalenceCounts = assemblyComplexity.getLeft();
        final int haplotypeComplexity = assemblyComplexity.getMiddle()[altAlleleIndex];
        final double haplotypeDominance = assemblyComplexity.getRight()[altAlleleIndex];

        final List<Double> result = new ArrayList<>(NUM_EXTRA_FEATURES);

        // take integer haplotype equivalence counts (already in order from greatest to least from Mutect)
        // and calculate the fractional share of the 2nd and 3rd, or 0 if none exist
        final double total = MathUtils.sum(haplotypeEquivalenceCounts);
        result.add(haplotypeEquivalenceCounts.length < 2 ? 0.0 : haplotypeEquivalenceCounts[1] / total);
        result.add(haplotypeEquivalenceCounts.length < 3 ? 0.0 : haplotypeEquivalenceCounts[2] / total);
        result.add((double) haplotypeComplexity);
        result.add(haplotypeDominance);

        IntStream.range(1, 6).forEach(repeatLength -> result.add((double) countRepeats(refBases.getBytes(), repeatLength)));

        Utils.validate(result.size() == NUM_EXTRA_FEATURES, "produced a variant feature vector of wrong size");
        return result;
    }

    // count how many repeats of length k surround the middle base
    // example: countRepeats(GACTACTACTG,3) = 3
    private int countRepeats(final byte[] refBases, final int k) {
        final int N = refBases.length;
        final int n = (N - 1) / 2;
        Utils.validateArg(k <= n, "Too few ref bases for given repeat length");

        // extend a repeat forward, to front and then to back(both exclusive)
        // note that this only extends backward if they match the bases going forward
        // that is AAAAAGTGTCC(first G in the middle) will get extended forward through
        // both GTs but won 't be extended back
        int front = n + k;
        while (front < N && refBases[front] == refBases[front - k]) {
            front++;
        }
        int back = n - 1;
        while (back >= 0 && refBases[back] == refBases[back + k]){
            back--;
        }
        final int forwardRepeats = (front - back - 1) / k;

        // same idea but extending backwards first (now back is exclusive and front is inclusive)
        back = n - k;
        while (back >= 0 && refBases[back] == refBases[back + k]) {
            back--;
        }
        front = n + 1;
        while (front < N && refBases[front] == refBases[front - k]) {
            front++;
        }
        final int backwardRepeats = (front - back - 1) / k;

        return FastMath.max(forwardRepeats, backwardRepeats);
    }

    private int[] sumADsOverSamples(final VariantContext vc, final Set<String> samples) {
        final int[] ADs = new int[vc.getNAlleles()];
        vc.getGenotypes(samples).stream().map(Genotype::getAD).forEach(ad -> new IndexRange(0, vc.getNAlleles()).forEach(n -> ADs[n] += ad[n]));
        return ADs;
    }

    @Override
    public void close() {
        printWriter.close();
    }
}
