package org.broadinstitute.hellbender.tools.walkers.annotator;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;

/**
 * For each sample and for each allele a list feature vectors of supporting reads
 * In order to reduce the number of delimiter characters, we flatten featurized reads.  For example, suppose allele 1 has featurized reads
 * [1,2] and [3,4] and allele 2 has featurized reads [5,6] and [7,8], the annotation is
 * 1,2,3,4|5,6,7,8
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY,
        summary="Featurized read sets for Mutect3 training data")
public class FeaturizedReadSets implements JumboGenotypeAnnotation {
    public static final int DEFAULT_BASE_QUALITY = 25;

    private static final int DEFAULT_MAX_REF_COUNT = Integer.MAX_VALUE;

    public static final int FEATURES_PER_READ = 11;

    private static final SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);

    // downsample ref reads to this count if needed
    private final int maxRefCount;

    public FeaturizedReadSets(final int maxRefCount) {
        this.maxRefCount = maxRefCount;
    }

    public FeaturizedReadSets() {
        this(DEFAULT_MAX_REF_COUNT);
    }

    @Override
    public void annotate(final ReferenceContext ref,
                         final FeatureContext features,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                         final AlleleLikelihoods<Fragment, Allele> fragmentLikelihoods,
                         final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods) {
        Utils.nonNull(vc);
        Utils.nonNull(g);
        Utils.nonNull(gb);

        if ( likelihoods == null) {
            return;
        }

        final List<List<List<Integer>>> readVectorsByAllele = getReadVectors(vc, Collections.singletonList(g.getSampleName()),
                likelihoods, haplotypeLikelihoods, maxRefCount, Integer.MAX_VALUE);

        // flatten twice: over all reads supporting each allele and over all alleles
        // we can partition by allele with the countsInAlleleOrder annotation
        // and by read using the constant feature vector size
        final int[] flattenedTensorInAlleleOrder = readVectorsByAllele.stream()
                .flatMap(listOfLists -> listOfLists.stream().flatMap(List::stream))
                .mapToInt(n -> n)
                .toArray();

        final int[] countsInAlleleOrder = readVectorsByAllele.stream().mapToInt(List::size).toArray();

        gb.attribute(GATKVCFConstants.FEATURIZED_READ_SETS_KEY, flattenedTensorInAlleleOrder);
        gb.attribute(GATKVCFConstants.FEATURIZED_READ_SETS_COUNTS_KEY, countsInAlleleOrder);
    }

    public static List<List<List<Integer>>> getReadVectors(final VariantContext vc,
                                                           final Collection<String> samples,
                                                           final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                                           final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods,
                                                           final int refDownsample,
                                                           final int altDownsample) {
        return getReadVectors(vc, samples, likelihoods, haplotypeLikelihoods, refDownsample, altDownsample, Collections.emptyMap());
    }

    // returns Lists (in allele order) of lists of read vectors supporting each allele
    public static List<List<List<Integer>>> getReadVectors(final VariantContext vc,
                                                           final Collection<String> samples,
                                                           final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                                           final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods,
                                                           final int refDownsample,
                                                           final int altDownsample,
                                                           final Map<Allele, Integer> altDownsampleMap) {
        final Map<Allele, List<GATKRead>> readsByAllele = likelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new ArrayList<>()));

        samples.stream().flatMap(s -> likelihoods.bestAllelesBreakingTies(s).stream())
                .filter(ba -> ba.isInformative())
                .forEach(ba -> readsByAllele.get(ba.allele).add(ba.evidence));

        // downsample if necessary
        final Allele refAllele = likelihoods.alleles().stream().filter(Allele::isReference).findFirst().get();
        for (final Allele allele : likelihoods.alleles()) {
            final int downsample = allele.isReference() ? refDownsample : altDownsampleMap.getOrDefault(allele, altDownsample);
            if (readsByAllele.get(allele).size() > downsample) {
                Collections.shuffle(readsByAllele.get(allele));
                readsByAllele.put(allele, readsByAllele.get(allele).subList(0, downsample));
            }
        }

        final Map<GATKRead, Haplotype> bestHaplotypes = new HashMap<>();
        samples.stream().flatMap(s -> haplotypeLikelihoods.bestAllelesBreakingTies(s).stream())
                .forEach(ba -> ba.evidence.getReads().forEach(read -> bestHaplotypes.put(read, ba.allele)));

        return vc.getAlleles().stream()
                .map(allele -> readsByAllele.get(allele).stream().map(read -> featurize(read, vc, bestHaplotypes)).collect(Collectors.toList()))
                .collect(Collectors.toList());
    }


    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.FEATURIZED_READ_SETS_KEY, GATKVCFConstants.FEATURIZED_READ_SETS_COUNTS_KEY);
    }

    private static List<Integer> featurize(final GATKRead read, final VariantContext vc, final Map<GATKRead, Haplotype> bestHaplotypes) {
        final List<Integer> result = new ArrayList<>();
        result.add(read.getMappingQuality());
        result.add(BaseQuality.getBaseQuality(read, vc).orElse(DEFAULT_BASE_QUALITY));
        result.add(read.isFirstOfPair() ? 1 : 0);
        result.add(read.isReverseStrand() ? 1 : 0);

        // distances from ends of read
        final int readPosition = ReadPosition.getPosition(read, vc).orElse(0);
        result.add(readPosition);
        result.add(read.getLength() - readPosition);


        result.add(Math.abs(read.getFragmentLength()));

        // distances from ends of fragment
        final int fragmentStart = Math.min(read.getMateStart(), read.getUnclippedStart());
        final int fragmentEnd = fragmentStart + Math.abs(read.getFragmentLength());
        result.add(vc.getStart() - fragmentStart);
        result.add(fragmentEnd - vc.getEnd());

        // mismatches versus best haplotype
        final Haplotype haplotype = bestHaplotypes.get(read);
        final SmithWatermanAlignment readToHaplotypeAlignment = aligner.align(haplotype.getBases(), read.getBases(), SmithWatermanAlignmentConstants.ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, SWOverhangStrategy.SOFTCLIP);
        final GATKRead copy = read.copy();
        copy.setCigar(readToHaplotypeAlignment.getCigar());
        final int mismatchCount = AlignmentUtils.getMismatchCount(copy, haplotype.getBases(), readToHaplotypeAlignment.getAlignmentOffset()).numMismatches;
        result.add(mismatchCount);

        final long indelsVsBestHaplotype = readToHaplotypeAlignment.getCigar().getCigarElements().stream().filter(el -> el.getOperator().isIndel()).count();
        result.add((int) indelsVsBestHaplotype);
        Utils.validate(result.size() == FEATURES_PER_READ, "Wrong number of features");

        return result;
    }

}
