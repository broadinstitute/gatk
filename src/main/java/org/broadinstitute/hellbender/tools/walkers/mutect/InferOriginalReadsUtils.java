package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class InferOriginalReadsUtils {

    // TODO: can this be refactored to share code with AssemblyBasedCallerUtils.getVariantContextsFromActiveHaplotypes?
    // TODO: write tests
    public static Map<LocationAndAlleles, List<Haplotype>> groupHaplotypesByAllelesAtThisLoc(final int loc,
                                                                                             final List<Haplotype> haplotypes,
                                                                                             final boolean includeSpanningEvents) {
        final Map<LocationAndAlleles, List<Haplotype>> results = new HashMap<>();

        haplotypes.stream()
                .map(h -> new ImmutablePair<>(h, h.getEventMap().getVariantContexts()))
                .filter(pair -> Objects.nonNull(pair.right))
                .forEach(pair -> {
                    for (VariantContext vc : pair.right){
                        if (!includeSpanningEvents && vc.getStart() != loc) {
                            continue;
                        }

                        final LocationAndAlleles locationAndAlleles = new LocationAndAlleles(vc.getStart(), vc.getAlleles());
                        if (results.containsKey(locationAndAlleles)) {
                            results.get(locationAndAlleles).add(pair.left);
                        } else {
                            final List<Haplotype> haplotypeList = new ArrayList<>();
                            haplotypeList.add(pair.left);
                            results.put(locationAndAlleles, haplotypeList);
                        }
                    }
                });

        return results;
    }

    // TODO: write tests
    public static double computeVarianceAroundMostCommon(final Map<LocationAndAlleles, List<Haplotype>> haplotypesByAllele,
                                                         final LocationAndAlleles mostCommonAllele) {
        final int alleleLengthOfMostCommon = getAlleleLength(mostCommonAllele.getAlleles().get(0), mostCommonAllele.getAlleles().get(1));

        double variance = 0.0;
        for (Map.Entry<LocationAndAlleles, List<Haplotype>> pair : haplotypesByAllele.entrySet()){
            final LocationAndAlleles allele = pair.getKey();
            final int alleleLength = getAlleleLength(allele.getAlleles().get(0), allele.getAlleles().get(1));
            final int count = pair.getValue().size();
            variance += count * Math.pow(alleleLength - alleleLengthOfMostCommon, 2);
        }

        return variance;
    }

    /** *signed* length of an alt allele (deletion is negative) **/
    private static int getAlleleLength(final Allele ref, final Allele alt){
        return ref.length() - alt.length();
    }

    public byte computeHaplotypePosteriorWithPairHMM(final ReadLikelihoodCalculationEngine likelihoodCalculationEngine,
                                                     final List<GATKRead> reads,
                                                     final IndexedSampleList indexedSampleList,
                                                     final AssemblyResultSet assemblyResult,
                                                     final SAMFileHeader header,
                                                     final String umi,
                                                     final InferOriginalReadEngine.ReadNum readNum,
                                                     final InferOriginalReadEngine.Strand strand){
        final Map<String, List<GATKRead>> readsBySample = AssemblyBasedCallerUtils.splitReadsBySample(indexedSampleList, header, reads);
        final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult, indexedSampleList, readsBySample);

        /** Now call genotypes. Should I just do assignGenotypeLikelihoods? {@link HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods}
         * Can think of a few approaches here.
         * 1. I can get the events. Filter to indels. Then call the indel length.
         * 2. Call haplotype as a whole, not a single event.
         * Is there a way to get a weighted set of haplotypes? That is set with counts of elements that went into it.
         *
         * The way haplotype caller does genotyping is as follows;
         * - Find events in the haplotypes. Split them by start positions.
         * - If deletions line up at the same position, for instance, we call them against ref, which...we don't want.
         *
         * Marginalization is more like ... homomorphism?
         * When do the ref bases become important? For now, give it
         *
         * Also see: {@link likelihoodToGenotypesHC}
         */

        final Pair<List<Haplotype>, double[]> haplotype2ProbabilityMap = calculateHaplotypePosterior(readLikelihoods, assemblyResult.getPaddedReferenceLoc());
        final int indexOfMaxPosteriorHaplotype = MathUtils.maxElementIndex(haplotype2ProbabilityMap.getRight());
        final double maximumPosterior = MathUtils.arrayMax(haplotype2ProbabilityMap.getRight());
        final byte indelQuality = posteriorToGapOpeningQualityMap(maximumPosterior);

        final Haplotype consensusBases = haplotype2ProbabilityMap.getLeft().get(indexOfMaxPosteriorHaplotype);
        final GATKRead consensusRead = haplotype2Read(consensusBases, reads.get(0), umi, strand, readNum, header);
        final EventMap eventMap = new EventMap(consensusBases, assemblyResult.getFullReferenceWithPadding(), assemblyResult.getPaddedReferenceLoc(), "c", 1);
        final List<VariantContext> indelEvents = eventMap.getVariantContexts().stream().filter(vc -> vc.isIndel()).collect(Collectors.toList());
        // Safe guard against getting weird haplotypes/list assumptions here
        final List<Integer> readLengths = reads.stream().mapToInt(r -> r.getLength()).distinct().boxed().collect(Collectors.toList());
        if (readLengths.size() > 2){
            int d = 3;
        }

        return indelQuality;
    }

    /**  Code skeleton borrowed from HaplotypeBAMWriter.writeHaplotype **/
    public static GATKRead haplotype2Read(final Haplotype haplotype,
                                          final GATKRead sampleRead,
                                          final String umi,
                                          final InferOriginalReadEngine.Strand strand,
                                          final InferOriginalReadEngine.ReadNum readNumber,
                                          final SAMFileHeader header) {
        // TOOD: perhaps benefitial to create a class for haplotype/GATKRead pair
        Utils.nonNull(haplotype, "haplotype cannot be null");

        /** Placeholders for read attributes that I will eventually have to assign properly **/
        final int MQ_PLACEHOLDER = 60;
        final String READ_NAME_PLACEHOLDER = umi + "_" + strand + "_" + readNumber;
        final SimpleInterval MATE_LOCATION_PLACEHOLDER = new SimpleInterval("17", 1, 151);
        final String READ_GROUP_PLACEHOLDER = "read_group1";

        final GATKRead read = new SAMRecordToGATKReadAdapter(new SAMRecord(header));
        read.setAttribute(InferOriginalReadEngine.CONSENSUS_READ_TAG, "CR"); // Is there a flag?
        read.setBases(haplotype.getBases());
        // read.setAlignmentStart(paddedRefLoc.getStart() + haplotype.getAlignmentStartHapwrtRef());
        read.setPosition(new SimpleInterval(haplotype.getLocation()));
        // Use a base quality value "!" for it's display value (quality value is not meaningful)
        read.setCigar(AlignmentUtils.consolidateCigar(haplotype.getCigar()));
        read.setMappingQuality(MQ_PLACEHOLDER);
        read.setName(READ_NAME_PLACEHOLDER);
        read.setReadGroup(READ_GROUP_PLACEHOLDER);

        // Mate info
        read.setMatePosition(MATE_LOCATION_PLACEHOLDER);

        if (sampleRead.isSecondOfPair()){
            read.isSecondOfPair();
        }

        if (sampleRead.isReverseStrand()){
            read.isReverseStrand();
        }

        return read;
    }

    private byte posteriorToGapOpeningQualityMap(final double posterirConsensusProbability){
        // New error probability = alpha*posteriorConsensusProb
        // Default error prob = 10^-4.5
        // Example: posteriorConsensusProb = p = 0.5
        // We assume that the mapping is a cubic function. Then the new gap opening probability is
        // f(p) = p^3 * 10^-4.5 = 0.125 * 10^-4.5 ~ 10^-5.5
        // In phred space, this is -10*log[p^3*(-4.5)]
        //
        final byte log10CorrectionFactor = (byte) ( - 3 * Math.log10(posterirConsensusProbability));
        return (byte) (log10CorrectionFactor + ReadUtils.DEFAULT_INSERTION_DELETION_QUAL);
    }

    /** This is the simpler implementation of the ploidy one genotype likelihood **/
    private Pair<List<Haplotype>, double[]> calculateHaplotypePosterior(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                        final SimpleInterval referenceLoc){
        final LikelihoodMatrix<GATKRead,Haplotype> matrix = readLikelihoods.sampleMatrix(0);
        final List<Haplotype> haplotypes = matrix.alleles();
        // Could use eventMaps and stuff here. But for now why not just use the entire haplotype?
        final int numHaplotypes = haplotypes.size();

        final boolean DEBUG = true;
// temporarily block
//        if (DEBUG){
//            final Path bam = Paths.get("/dsde/working/tsato/consensus/tp53/test/haplotypes.bam");
//            final HaplotypeBAMWriter hbw = new HaplotypeBAMWriter(HaplotypeBAMWriter.WriterType.CALLED_HAPLOTYPES, bam, true, false, header);
//            hbw.writeHaplotypesAsReads(haplotypes, new TreeSet<>(haplotypes), referenceLoc, referenceLoc);
//            hbw.close();
//        }

        // multiply \prod_i p(read_i|haplotype_k) for each of k haplotypes
        final double[] log10HaplotypeLikelihoods = new double[numHaplotypes];
        for (int k = 0; k < numHaplotypes; k++){
            for (int i = 0; i < matrix.evidenceCount(); i++){
                log10HaplotypeLikelihoods[k] += matrix.get(k, i);
            }
        }

        final double[] haplotypePosterior = MathUtils.normalizeLog10(log10HaplotypeLikelihoods, false, false);
//        final Map<Haplotype, Double> result = new HashMap<>();
//        IntStream.range(0, haplotypes.size()).forEach(k -> result.put(haplotypes.get(k), haplotypePosterior[k]));
        // (12/18/19) This is just an enigma. What i need is like a plot of default gap opening vs posterior of each (and of each read...) yea, let's do that
        return new ImmutablePair<>(haplotypes, haplotypePosterior);
    }

}
