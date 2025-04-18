package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.AssemblyComplexity;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReferenceBases;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
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

public class PermutectDatasetEngine implements AutoCloseable {

    public static final int CAPACITY = 100000;

    private enum VariantType {
        SNV, INSERTION, DELETION
    }

    private enum Label {
        ARTIFACT, VARIANT, UNLABELED, IGNORE
    }

    private final SAMSequenceDictionary sequenceDictionary;

    private final Map<String, Integer> readGroupIndices = new HashMap<>();

    // number of additional variant features for assembly complexity (4), TLOD / tumor depth (1)
    private static final int NUM_EXTRA_FEATURES = 5;

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
    private final PrintWriter contigPrintWriter;
    private final PrintWriter readGroupPrintWriter;

    // number of nonartifact data to keep for each artifact datum
    private final int nonArtifactPerArtifact;

    // are we generating dataset for training a model or for filtering calls with a pre-trained model?
    private final boolean trainingMode;

    private final Set<String> normalSamples;

    // simple method to balance data: for each k-alt-read artifact there are
    // nonArtifactPerArtifact (downsampled) k-alt-read non-artifacts.
    private final EnumMap<VariantType, ArrayBlockingQueue<Integer>> unmatchedArtifactAltCounts;

    private final Random random = new Random(1);


    public PermutectDatasetEngine(final File datasetFile, final boolean trainingMode, final int maxRefCount,
                                  final int maxAltCount, final int nonArtifactPerArtifact, final Set<String> normalSamples,
                                  final SAMFileHeader header, final SAMSequenceDictionary sequenceDictionary) {
        try {
            printWriter = new PrintWriter(new FileWriter(Utils.nonNull(datasetFile)));
            final File contigTableFile = datasetFile.toPath().resolveSibling("contigs.table").toFile();
            final File readGroupTableFile = datasetFile.toPath().resolveSibling("read-groups.table").toFile();
            contigPrintWriter = new PrintWriter(new FileWriter(contigTableFile));
            readGroupPrintWriter = new PrintWriter(new FileWriter(readGroupTableFile));
        } catch (IOException ex) {
            throw new UserException.BadInput("Could not create dataset file writer");
        }

        this.normalSamples = Utils.nonNull(normalSamples);
        this.trainingMode = trainingMode;
        this.nonArtifactPerArtifact = nonArtifactPerArtifact;
        this.maxRefCount = maxRefCount;
        this.maxAltCount = maxAltCount;

        this.sequenceDictionary = sequenceDictionary;
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        for (int n = 0; n < readGroups.size(); n++) {
            readGroupIndices.put(readGroups.get(n).getReadGroupId(), n);
        }

        unmatchedArtifactAltCounts = new EnumMap<>(VariantType.class);
        for (final VariantType type : VariantType.values()) {
            unmatchedArtifactAltCounts.put(type, new ArrayBlockingQueue<>(CAPACITY));
        }
    }

    // add one datum per alt allele
    public void addData(final ReferenceContext ref, final VariantContext vc, Optional<List<VariantContext>> truthVCs,
                        final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                        final AlleleLikelihoods<Fragment, Haplotype> logFragmentLikelihoods,
                        final AlleleLikelihoods<Fragment, Allele> logFragmentAlleleLikelihoods,
                        final M2ArgumentCollection.PermutectDatasetMode permutectDatasetMode) {
        final String refBases = ReferenceBases.annotate(ref, vc);
        final String refAllele = vc.getReference().getBaseString();
        final int contigIndex = sequenceDictionary.getSequenceIndex(vc.getContig());
        final int position = vc.getStart();
        final Set<String> tumorSamples = likelihoods.samples().stream().filter(sample -> !normalSamples.contains(sample)).collect(Collectors.toSet());
        final int numAlt = vc.getNAlleles() - 1;


        // the variant has already been annotated, so we have POPAF, AD, and in-PON status
        final double[] popafs = VariantContextGetters.getAttributeAsDoubleArray(vc, GATKVCFConstants.POPULATION_AF_KEY);
        final double[] tumorLods = Mutect2FilteringEngine.getTumorLogOdds(vc);
        // the PoN is NOT for calling! We only use it to assign some unlabeled data as artifacts in training
        final boolean isBlacklisted = vc.hasAttribute(GATKVCFConstants.IN_PON_KEY);

        // These ADs, which later make up the pre-downsampling depths, come from the genotype AD field applied by Mutect2.
        // This means that uninformative reads are not discarded; rather the expected non-integral ADs are rounded to
        // the nearest integer.  Also, these ADs are *fragment* ADs from the log fragment-allele likelihoods in the
        // SomaticGenotypingEngine.
        final int[] tumorADs = sumADsOverSamples(vc, tumorSamples);
        final int[] normalADs = sumADsOverSamples(vc, normalSamples);
        final int tumorDepth = (int) MathUtils.sum(tumorADs);
        final int normalDepth = (int) MathUtils.sum(normalADs);
        final boolean hasNormal = normalDepth > 0;

        final List<Allele> allRefAlleles = new ArrayList<>();
        allRefAlleles.add(vc.getReference());
        truthVCs.ifPresent(vcs -> vcs.forEach(var -> allRefAlleles.add(var.getReference())));
        final Allele longestRef = allRefAlleles.stream().sorted(Comparator.comparingInt(Allele::length).reversed()).findFirst().get();

        // skip(1) excludes the reference allele
        final List<Allele> remappedAltAlleles = ReferenceConfidenceVariantContextMerger.remapAlleles(vc, longestRef).stream()
                .skip(1).toList();

        final Set<Allele> truthAlleles = !truthVCs.isPresent() ? Collections.emptySet() : truthVCs.get().stream()
                .filter(truthVC -> ! truthVC.isFiltered())
                .flatMap(truthVC -> ReferenceConfidenceVariantContextMerger.remapAlleles(truthVC, longestRef).stream())
                .collect(Collectors.toSet());

        final List<Label> labels = new ArrayList<>(numAlt);
        final Map<Allele, Integer> altDownsampleMap= new HashMap<>();

        for (int n = 0; n < numAlt; n++) {
            final double tumorAF = tumorADs[n+1] / ((double) tumorDepth);
            final double normalAF = hasNormal ? normalADs[n+1] / ((double) normalDepth) : 0.0;
            Allele altAllele = vc.getAlternateAllele(n);
            final Allele remappedAltAlelle = remappedAltAlleles.get(n);
            final String altAlleleString = altAllele.getBaseString();
            final int diff = altAlleleString.length() - refAllele.length();
            final VariantType type = diff == 0 ? VariantType.SNV : ( diff > 0 ? VariantType.INSERTION : VariantType.DELETION);

            if (trainingMode) { // training mode -- collecting tensors to train the Permutect artifact model
                final ArrayBlockingQueue<Integer> unmatchedQueue = unmatchedArtifactAltCounts.get(type);
                final boolean likelySeqError = tumorLods[n] < TLOD_THRESHOLD;

                // extremely strict criteria because there are so many germline variants we can afford to waste a lot
                final boolean definiteGermline = !likelySeqError && popafs[n] < COMMON_POPAF_THRESHOLD &&
                        tumorAF > 0.35 && (!hasNormal || normalAF > 0.35);
                final boolean trueVariant = truthVCs.isPresent() ? truthAlleles.contains(remappedAltAlelle) : definiteGermline;

                // low AF in tumor and normal, rare in population implies artifact
                boolean probableArtifact = !likelySeqError && (truthVCs.isPresent() ? !truthAlleles.contains(remappedAltAlelle) :
                        (tumorAF < 0.2 && popafs[n] > RARE_POPAF_THRESHOLD));

                if  (probableArtifact) {
                    if (unmatchedQueue.size() > 0.9 * CAPACITY) { // this should rarely come up
                        labels.add(Label.IGNORE);
                    } else if (tumorAF < 0.3) {
                        labels.add(Label.ARTIFACT);
                        unmatchedQueue.addAll(Collections.nCopies(nonArtifactPerArtifact, tumorADs[n + 1]));
                    } else {    // very high AF could just be a germline variant missing from the truth VCF.  Leave it unlabeled for semisupervised learning to deal with
                        labels.add(Label.UNLABELED);
                    }
                } else if (trueVariant && !unmatchedQueue.isEmpty()) {
                    // high AF in tumor and normal, common in population implies germline, which we downsample
                    labels.add(Label.VARIANT);
                    // Note: we used to have "altDownsampleMap.put(altAllele, unmatchedQueue.poll())" here, which forces the alt count distribution of
                    // variants to match that of artifacts.  Due to Permutect's count-agnostic architecture this is
                    // unnecessary, and it gives the model too few high-AF variant examples on which to learn calibration.
                    unmatchedQueue.poll(); // just to pop off a value
                    altDownsampleMap.put(altAllele, random.nextInt(5, 20));
                } else if (tumorLods[n] > 5.0 && isBlacklisted) {   // being in pon turns what would otherwise be unlabeled into artifact
                    labels.add(Label.ARTIFACT);
                } else if (tumorLods[n] > 4.0 && tumorAF < 0.3) {
                    labels.add(Label.UNLABELED);
                } else {
                    labels.add(Label.IGNORE);
                }
            } else {    // not training mode -- we are generating tensors in order to apply the Permutect artifact model to a callset
                if (truthVCs.isPresent()) {
                    // here, for the purposes of test data, both sequencing errors and technical artifacts get the "ARTIFACT" label
                    labels.add(truthAlleles.contains(remappedAltAlelle) ? Label.VARIANT : Label.ARTIFACT);
                } else {
                    labels.add(Label.UNLABELED);
                }
            }
        }

        Utils.validate(labels.size() == numAlt, "We have not labeled every alt, or have labeled too much");
        // we check this later allele by allele, but we can save a lot of compute here if we realize no alt allele yields training data
        if (trainingMode && labels.stream().allMatch(label -> label == Label.IGNORE)) {
            return;
        }

        // haplotype equivalence counts, haplotype complexity, haplotype dominance
        final Triple<int[], int[], double[]> assemblyComplexity = AssemblyComplexity.annotate(vc, logFragmentLikelihoods, false);

        // TODO: for now we don't really need normal reads
        // note that the following use the VC's allele order, not necessarily the likelihoods' allele order
        final List<List<List<Integer>>> normalReadVectorsByAllele =  FeaturizedReadSets.getReadVectors(vc, normalSamples,
                likelihoods, logFragmentLikelihoods, maxRefCount, maxAltCount, permutectDatasetMode, readGroupIndices);
        final List<List<List<Integer>>> tumorReadVectorsByAllele =  FeaturizedReadSets.getReadVectors(vc, tumorSamples,
                likelihoods, logFragmentLikelihoods, maxRefCount, maxAltCount, altDownsampleMap, permutectDatasetMode, readGroupIndices);

        // ref and alt reads have already been downsampled by the read featurizer
        final List<List<Integer>> tumorRefReads = tumorReadVectorsByAllele.get(0);
        final List<List<Integer>> normalRefReads = normalReadVectorsByAllele.get(0);

        final List<LikelihoodMatrix<Fragment,Allele>> tumorMatrices = tumorSamples.stream()
                .map(s -> logFragmentAlleleLikelihoods.sampleMatrix(logFragmentAlleleLikelihoods.indexOfSample(s)))
                .collect(Collectors.toList());

        for (int n = 0; n < numAlt; n++) {
            if (labels.get(n) ==  Label.IGNORE) {
                continue;
            }
            final double tlod = vc.getAttributeAsDoubleList("TLOD", 0).get(n);

            final String altAllele = vc.getAlternateAllele(n).getBaseString();
            final List<Double> variantFeatureVector = variantFeatures(n, assemblyComplexity, tlod, tumorADs[n+1]);
            final List<List<Integer>> tumorAltReads = tumorReadVectorsByAllele.get(n+1);
            final List<List<Integer>> normalAltReads = normalReadVectorsByAllele.get(n+1);

            printWriter.println(labels.get(n).toString());
            printWriter.printf("%d:%d,%s->%s%n", contigIndex, position, refAllele, altAllele);
            printWriter.println(refBases);
            printWriter.println(numberString(variantFeatureVector, "%.2f", " "));
            //printWriter.printf("%d %d %d %d%n", tumorRefReads.size(), tumorAltReads.size(), normalRefReads.size(), normalAltReads.size());
            // zeroes because we currently don't use normal read vectors
            printWriter.printf("%d %d %d %d%n", tumorRefReads.size(), tumorAltReads.size(), 0, 0);

            tumorRefReads.forEach(r -> printWriter.println(integerString(r)));
            tumorAltReads.forEach(r -> printWriter.println(integerString(r)));
            //normalRefReads.forEach(r -> printWriter.print(numberString(r)));
            //normalAltReads.forEach(r -> printWriter.print(numberString(r)));
            printWriter.printf("%d %d %d %d%n", tumorDepth, tumorADs[n+1], normalDepth, normalADs[n+1]);  // pre-downsampling counts
            // this is approximately the likelihood that these particular reads are alt given sequencing error, excluding
            // the depth C N_alt combinatorial factor that is common to all likelihoods in Permutect
            // basically, it's the TLOD with a correction for the marginalized flat prior from M2
            final double seqErrorLogLikelihood = -MathUtils.log10ToLog(tlod) - Math.log(tumorDepth + 1);
            printWriter.printf("%.3f%n", seqErrorLogLikelihood);

            // and do the same for the normal
            final double nalod = normalSamples.isEmpty() ? 0 : vc.getAttributeAsDoubleList("NALOD", 0).get(n);
            final double normalSeqErrorLogLikelihood = -MathUtils.log10ToLog(nalod) - Math.log(normalDepth + 1);
            printWriter.printf("%.3f%n", normalSeqErrorLogLikelihood);
        }
    }

    private String integerString(final List<Integer> numbers) {
        return numberString(numbers, "%d", " ");
    }

    private String numberString(final List<? extends Number> numbers, final String formatString, final String separator) {
        final boolean decimal = formatString.endsWith("f");
        return numbers.stream().map(x -> String.format(formatString, decimal ? x.floatValue() : x)).collect(Collectors.joining(separator));
    }

    private List<Double> variantFeatures(final int altAlleleIndex, Triple<int[], int[], double[]> assemblyComplexity, final double tlod, final int tumorAltCount) {
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
        result.add(tumorAltCount == 0 ? 0.0 : tlod / tumorAltCount);

        Utils.validate(result.size() == NUM_EXTRA_FEATURES, "produced a variant feature vector of wrong size");
        return result;
    }

    private int[] sumADsOverSamples(final VariantContext vc, final Set<String> samples) {
        final int[] ADs = new int[vc.getNAlleles()];
        vc.getGenotypes(samples).stream().map(Genotype::getAD).forEach(ad -> new IndexRange(0, vc.getNAlleles()).forEach(n -> ADs[n] += ad[n]));
        return ADs;
    }

    @Override
    public void close() {
        printWriter.close();

        for (final SAMSequenceRecord contigRecord : sequenceDictionary.getSequences()) {
            contigPrintWriter.println(String.format("%s\t%d", contigRecord.getContig(), contigRecord.getSequenceIndex()));
        }

        for (final Map.Entry<String, Integer> entry : readGroupIndices.entrySet()) {
            readGroupPrintWriter.println(String.format("%s\t%d", entry.getKey(), entry.getValue()));
        }

        contigPrintWriter.close();
        readGroupPrintWriter.close();
    }
}
