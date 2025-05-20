package org.broadinstitute.hellbender.tools.walkers.mutect;


import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.tools.walkers.annotator.BaseQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReadPosition;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * For each sample and for each allele a list feature vectors of supporting reads
 * In order to reduce the number of delimiter characters, we flatten featurized reads.  For example, suppose allele 1 has featurized reads
 * [1,2] and [3,4] and allele 2 has featurized reads [5,6] and [7,8], the annotation is
 * 1,2,3,4|5,6,7,8
 */
public class FeaturizedReadSets {
    private static final Logger logger = LogManager.getLogger(FeaturizedReadSets.class);

    public static final int DEFAULT_BASE_QUALITY = 25;

    private static final SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
    private static final int FEATURES_PER_RANGE = 5;
    private static final List<Integer> RANGES = List.of(5, 10, 20);
    public static final int NUM_RANGED_FEATURES = FEATURES_PER_RANGE * RANGES.size();
    private static final int VERY_BAD_QUAL_THRESHOLD = 10;
    private static final int BAD_QUAL_THRESHOLD = 20;

    private FeaturizedReadSets() { }

    public static List<List<List<Integer>>> getReadVectors(final VariantContext vc,
                                                           final Collection<String> samples,
                                                           final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                                           final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods,
                                                           final int refDownsample,
                                                           final int altDownsample,
                                                           final M2ArgumentCollection.PermutectDatasetMode permutectDatasetMode,
                                                           final Map<String, Integer> readGroupIndices) {
        return getReadVectors(vc, samples, likelihoods, haplotypeLikelihoods, refDownsample, altDownsample, Collections.emptyMap(), permutectDatasetMode, readGroupIndices);
    }

    // returns Lists (in allele order) of lists of read vectors supporting each allele
    public static List<List<List<Integer>>> getReadVectors(final VariantContext vc,
                                                           final Collection<String> samples,
                                                           final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                                           final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods,
                                                           final int refDownsample,
                                                           final int altDownsample,
                                                           final Map<Allele, Integer> altDownsampleMap,
                                                           final M2ArgumentCollection.PermutectDatasetMode permutectDatasetMode,
                                                           final Map<String, Integer> readGroupIndices) {
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

        // we don't want too many potential "best" haplotypes because otherwise errors may have
        // their own personal best haplotypes.  That said, we don't want a genuine variant to be missing a haplotype.
        // In most cases -- diploid germline when generating training data, diploid germline plus one somatic variant for
        // most tumor calling -- there is a strong prior of at most two germline haplotypes and at most one other
        // variant haplotype, and thus this isn't actually that delicate of a balance.  In some cases, however, such as
        // mitochondria where there are many valid haplotypes with low allele fraction, we may eventually need to be
        // more careful.  For now, here is our heuristic for which haplotypes to keep:
        // 1) Any allele in the variant context at this locus gets at least one haplotype, the one with the most read
        //  support among haplotypes containing that allele.
        // 2) Keep any other haplotype that is best-supported by at least 5 reads and at least 10% of all reads.
        // TODO: replace this by a simple likelihood calculation so that the percentage required decreases with total
        // TODO: read support
        // TODO: distinguish between normal and tumor sample read support

        // count the number of reads that best support each haplotype
        final Map<Haplotype, MutableInt> haplotypeSupportCounts = new HashMap<>();
        for (final String sample : samples) {
            for (final AlleleLikelihoods<Fragment, Haplotype>.BestAllele fragmentBestHaplotype : haplotypeLikelihoods.bestAllelesBreakingTies(sample)) {
                final Haplotype bestHaplotype = fragmentBestHaplotype.allele;
                haplotypeSupportCounts.putIfAbsent(bestHaplotype, new MutableInt(0));
                haplotypeSupportCounts.get(bestHaplotype).increment();
            }
        }

        final Map<Allele, Set<Haplotype>> haplotypesByAllele = new HashMap<>();

        // determine which haplotypes contain which alleles.  Note that the haplotypes' event maps are in trimmed form
        // while the VC's alleles might not be if there are multiple alt alleles.  For example:
        // VC ref = GT, alt = G, GTT
        // the event for the second alt allele is G -> GT
        // since events' alt alleles are in minimal form and left-aligned, we only have to pad them by whatever difference
        // there exists at the end of the ref allele.
        // for example, in the G -> GT event, the ref has one base while the VC ref has two bases, thus we must append
        // the (2-1) = 1 last base of the VC ref to get GT + T = GTT
        // There does exist an edge case where the event ref may be longer than the VC if it corresponds to a deletion
        // that didn't make it into the VC (not enough evidence) that is longer than anything in the VC.  This case is definitely
        // an event that's not in the VC, though, so we can filter it.
        final Allele vcRef = vc.getReference();
        for (final Haplotype haplotype : haplotypeLikelihoods.alleles()) {
            final Allele alleleAtThisSite = haplotype.getEventMap().getOverlappingEvents(vc.getStart()).stream()
                    .filter(e -> e.getStart() == vc.getStart()) // there can be at most one event at the VC start
                    .filter(e -> e.refAllele().length() <= vcRef.length())
                    .map(event -> {
                        final Allele eventRef = event.refAllele();
                        final int excess = vcRef.length() - eventRef.length();  // we have ensured that this is non-negative
                        final String basesToAppend = vcRef.getBaseString().substring(vcRef.length()-excess);
                        final String altString = event.altAllele().getBaseString() + basesToAppend;
                        final Allele paddedAllele = Allele.create(altString, false);
                        return paddedAllele;
                    })
                    .findFirst()
                    .orElse(vc.getReference());

            haplotypesByAllele.putIfAbsent(alleleAtThisSite, new HashSet<>());
            haplotypesByAllele.get(alleleAtThisSite).add(haplotype);
        }

        final Set<Haplotype> haplotypesToKeep = new HashSet<>();

        // add the best haplotype for each allele, including reference, at this locus
        for (final Map.Entry<Allele, Set<Haplotype>> entry : haplotypesByAllele.entrySet()) {
            final Optional<Haplotype> mostSupportedHaplotype = entry.getValue().stream()
                    .max(Comparator.comparingInt(h -> haplotypeSupportCounts.getOrDefault(h, new MutableInt(0)).intValue()));
            mostSupportedHaplotype.ifPresent(h -> haplotypesToKeep.add(h));
        }

        // TODO: HEURISTIC -- IMPROVE THIS!!!
        // add any haplotype supported by at least 5 reads AND 10% of reads
        for (final Map.Entry<Haplotype, MutableInt> haplotypeAndSupportCount : haplotypeSupportCounts.entrySet()) {
            if (haplotypeAndSupportCount.getValue().intValue() >= 5) {
                haplotypesToKeep.add(haplotypeAndSupportCount.getKey());
            }
        }

        // create likelihoods of just the haplotypes we believe in
        final AlleleLikelihoods<Fragment, Haplotype> restrictedHaplotypeLikelihoods =
                haplotypeLikelihoods.removeAllelesToSubset(haplotypesToKeep);


        // TODO: left off with need to explore this breakpoint
        final Map<GATKRead, Haplotype> bestHaplotypes = new HashMap<>();
        samples.stream().flatMap(s -> restrictedHaplotypeLikelihoods.bestAllelesBreakingTies(s).stream())
                .forEach(ba -> ba.evidence.getReads().forEach(read -> bestHaplotypes.put(read, ba.allele)));

        // Step 1: find the most-supported haplotype containing each alt allele at this locus, i.e. if it's multiallelic
        // with A->C and A->G substitutions, get the most-supported C and the most supported G haplotpyes.  Do this regardless
        // of support.
        // Step 2: find all haplotypes with sufficient support to think they could be germline.
        // Step 3: restrict haplotypeLikelihoods to just these haplotypes
        // Step 4: featurize with respect to best haplotypes as before

        // In the VC's allele order (careful, this might be different from the likelihoods' order!), for each allele
        // take all reads
        return vc.getAlleles().stream()
                .map(allele -> readsByAllele.get(allele).stream().map(read -> featurize(read, vc, bestHaplotypes, permutectDatasetMode, readGroupIndices)).collect(Collectors.toList()))
                .collect(Collectors.toList());
    }


    private static List<Integer> featurize(final GATKRead read, final VariantContext vc,
                                           final Map<GATKRead, Haplotype> bestHaplotypes,
                                           final M2ArgumentCollection.PermutectDatasetMode permutectDatasetMode,
                                           final Map<String, Integer> readGroupIndices) {
        final List<Integer> result = new ArrayList<>();
        result.add(readGroupIndices.get(read.getReadGroup()));  // this is read group metadata rather than part of the tensor
        result.add(read.getMappingQuality());

        // TODO: why not add BQ before and after as well -- especially useful for indels
        result.add(BaseQuality.getBaseQuality(read, vc).orElse(DEFAULT_BASE_QUALITY));
        result.add(read.isFirstOfPair() ? 1 : 0);
        result.add(read.isReverseStrand() ? 1 : 0);

        // distances from ends of read
        // this DOES account for the hard clips due to fitting the read inside the assembly window!!
        final int readPositionOfVariantStart = ReadPosition.getPosition(read, vc).orElse(0);
        result.add(readPositionOfVariantStart);

        // read.getLength(), however, does not account for hard clips and we need to add the length of any leading or trailing hard clips in the read's CIGAR
        final int totalHardClips = read.getCigarElements().stream().filter(el -> el.getOperator() == CigarOperator.HARD_CLIP).mapToInt(el -> el.getLength()).sum();
        final int actualReadLength = read.getLength() + totalHardClips;
        result.add(actualReadLength - readPositionOfVariantStart);


        result.add(Math.abs(read.getFragmentLength()));

        if (permutectDatasetMode == M2ArgumentCollection.PermutectDatasetMode.ILLUMINA) {
            // distances from ends of fragment
            final int fragmentStart = Math.min(read.getMateStart(), read.getUnclippedStart());
            final int fragmentEnd = fragmentStart + Math.abs(read.getFragmentLength());
            result.add(vc.getStart() - fragmentStart);
            result.add(fragmentEnd - vc.getEnd());
        }

        // Ultima-specific read tags
        if (permutectDatasetMode == M2ArgumentCollection.PermutectDatasetMode.ULTIMA) {
            result.add(read.getAttributeAsInteger("si"));   // si is an integer on the order of 100s or 1000s
            result.add((int) (1000*read.getAttributeAsFloat("rq")));    // rq is a float from 0 and 1, so we multiply by 1000 and round
        }

        // mismatches versus best haplotype
        final Haplotype haplotype = bestHaplotypes.get(read);

        // TODO: fix this
        // I have no idea why this edge case occurs in Ultima data and maybe sometimes in Illumina
        if (!bestHaplotypes.containsKey(read)) {
            logger.warn(String.format("Best haplotypes don't contain read %s at position %s:%d.", read.getName(),
                    vc.getContig(), vc.getStart()));
            //result.add(3);
            //result.add(2);

            for (int n = 0; n < NUM_RANGED_FEATURES; n++) {
                result.add(0);
            }
        } else {
            byte[] haplotypeBases = haplotype.getBases();
            final SmithWatermanAlignment readToHaplotypeAlignment = aligner.align(haplotypeBases, read.getBases(), SmithWatermanAlignmentConstants.ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, SWOverhangStrategy.SOFTCLIP);
            final GATKRead copy = read.copy();
            copy.setCigar(readToHaplotypeAlignment.getCigar());
            //final int mismatchCount = AlignmentUtils.getMismatchCount(copy, haplotypeBases, readToHaplotypeAlignment.getAlignmentOffset()).numMismatches;
            //result.add(mismatchCount);

            //final long indelsVsBestHaplotype = readToHaplotypeAlignment.getCigar().getCigarElements().stream().filter(el -> el.getOperator().isIndel()).count();
            //result.add((int) indelsVsBestHaplotype);

            final int readStartInHaplotype = readToHaplotypeAlignment.getAlignmentOffset();
            final AlignmentStateMachine asm = new AlignmentStateMachine(copy);
            asm.stepForwardOnGenome();
            final List<int[]> rangedFeatures = RANGES.stream().map(range -> new int[FEATURES_PER_RANGE]).toList();

            while (!asm.isRightEdge()) {
                final PileupElement pe = asm.makePileupElement();
                final int distanceFromVariant = Math.abs(asm.getReadOffset() - readPositionOfVariantStart);

                // pick which array's features we are accounting.  If the ranges are 5, 10, 20, and the distance is, say 8, then the '<= 10' range is relevant
                final OptionalInt relevantRange = IntStream.range(0, RANGES.size()).filter(n -> distanceFromVariant <= RANGES.get(n)).findFirst();
                if (relevantRange.isPresent()) {
                    final int[] featuresToAddTo = rangedFeatures.get(relevantRange.getAsInt());
                    if (pe.isBeforeInsertion()) {
                        featuresToAddTo[0] += pe.getLengthOfImmediatelyFollowingIndel();
                    }

                    if (pe.isDeletion()) {
                        featuresToAddTo[1]++;
                    } else {
                        final byte base = pe.getBase();
                        final byte qual = pe.getQual();
                        final byte haplotypeBase = haplotypeBases[asm.getGenomeOffset() + readStartInHaplotype];

                        if (base != haplotypeBase) {
                            featuresToAddTo[2]++;
                        }

                        if (qual < VERY_BAD_QUAL_THRESHOLD) {
                            featuresToAddTo[3]++;
                        } else if (qual < BAD_QUAL_THRESHOLD) {
                            featuresToAddTo[4]++;
                        }
                    }
                }
                asm.stepForwardOnGenome();
            }

            for (final int[] featuresToAdd : rangedFeatures) {
                for (final int val : featuresToAdd) {
                    result.add(val);
                }
            }
        }
        // the +1 is for the read group index that comes before the features
        Utils.validate(result.size() == permutectDatasetMode.getNumReadFeatures() + 1, "Wrong number of features");

        return result;
    }

}
