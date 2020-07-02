package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.fragments.FragmentCollection;
import org.broadinstitute.hellbender.utils.fragments.FragmentUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * Created by davidben on 9/8/16.
 */
public final class AssemblyBasedCallerUtils {

    static final int REFERENCE_PADDING_FOR_ASSEMBLY = 500;
    public static final int NUM_HAPLOTYPES_TO_INJECT_FORCE_CALLING_ALLELES_INTO = 5;
    public static final String SUPPORTED_ALLELES_TAG="XA";
    public static final String CALLABLE_REGION_TAG = "CR";
    public static final String ALIGNMENT_REGION_TAG = "AR";
    public static final Function<Haplotype, Double> HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY = h -> {
        final Cigar cigar = h.getCigar();
        final int referenceTerm = (h.isReference() ? 1 : 0);
        final int cigarTerm = cigar == null ? 0 : (1 - cigar.numCigarElements());
        return (double) referenceTerm + cigarTerm;
    };

    // After trimming to fit the assembly window, throw away read stubs shorter than this length
    // if we don't, the several bases left of reads that end just within the assembly window can
    // get realigned incorrectly.  See https://github.com/broadinstitute/gatk/issues/5060
    public static final int MINIMUM_READ_LENGTH_AFTER_TRIMMING = 10;

    private static final String phase01 = "0|1";
    private static final String phase10 = "1|0";

    /**
     * Returns a map with the original read as a key and the realigned read as the value.
     * <p>
     *     Missing keys or equivalent key and value pairs mean that the read was not realigned.
     * </p>
     * @return never {@code null}
     */
    public static Map<GATKRead, GATKRead> realignReadsToTheirBestHaplotype(final AlleleLikelihoods<GATKRead, Haplotype> originalReadLikelihoods, final Haplotype refHaplotype, final Locatable paddedReferenceLoc, final SmithWatermanAligner aligner) {
        final Collection<AlleleLikelihoods<GATKRead, Haplotype>.BestAllele> bestAlleles = originalReadLikelihoods.bestAllelesBreakingTies(HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY);
        final Map<GATKRead, GATKRead> result = new HashMap<>(bestAlleles.size());

        for (final AlleleLikelihoods<GATKRead, Haplotype>.BestAllele bestAllele : bestAlleles) {
            final GATKRead originalRead = bestAllele.evidence;
            final Haplotype bestHaplotype = bestAllele.allele;
            final boolean isInformative = bestAllele.isInformative();
            final GATKRead realignedRead = AlignmentUtils.createReadAlignedToRef(originalRead, bestHaplotype, refHaplotype, paddedReferenceLoc.getStart(), isInformative, aligner);
            result.put(originalRead, realignedRead);
        }
        return result;
    }

    public static void finalizeRegion(final AssemblyRegion region,
                                      final boolean errorCorrectReads,
                                      final boolean skipSoftClips,
                                      final byte minTailQuality,
                                      final SAMFileHeader readsHeader,
                                      final SampleList samplesList,
                                      final boolean correctOverlappingBaseQualities) {
        if ( region.isFinalized() ) {
            return;
        }

        final byte minTailQualityToUse = errorCorrectReads ? HaplotypeCallerEngine.MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION : minTailQuality;
        final List<GATKRead> readsToUse = region.getReads().stream()
                // TODO unclipping soft clips may introduce bases that aren't in the extended region if the unclipped bases
                // TODO include a deletion w.r.t. the reference.  We must remove kmers that occur before the reference haplotype start
                .map(read -> skipSoftClips || ! ReadUtils.hasWellDefinedFragmentSize(read) ?
                    ReadClipper.hardClipSoftClippedBases(read) : ReadClipper.revertSoftClippedBases(read))
                .map(read -> ReadClipper.hardClipLowQualEnds(read, minTailQualityToUse))
                .filter(read -> read.getStart() <= read.getEnd())
                .map(read -> read.isUnmapped() ? read : ReadClipper.hardClipAdaptorSequence(read))
                .filter(read ->  !read.isEmpty() && read.getCigar().getReadLength() > 0)
                .map(read -> ReadClipper.hardClipToRegion(read, region.getPaddedSpan().getStart(), region.getPaddedSpan().getEnd() ))
                .filter(read -> read.getStart() <= read.getEnd() && read.getLength() > 0 && read.overlaps(region.getPaddedSpan()))
                .sorted(new ReadCoordinateComparator(readsHeader)) // TODO: sort may be unnecessary here
                .collect(Collectors.toList());

        // handle overlapping read pairs from the same fragment
        if (correctOverlappingBaseQualities) {
            cleanOverlappingReadPairs(readsToUse, samplesList, readsHeader, true, OptionalInt.empty(), OptionalInt.empty());
        }

        region.clearReads();
        region.addAll(readsToUse);
        region.setFinalized(true);
    }

    /**
     *  Modify base qualities when paired reads overlap to account for the possibility of PCR error.
     *
     *  Overlapping mates provded independent evidence as far as sequencing error is concerned, but their PCR errors
     *  are correlated.  The base qualities are thus limited by the sequencing base quality as well as half of the PCR
     *  quality.  We use half of the PCR quality because downstream we treat read pairs as independent, and summing two halves
     *  effectively gives the PCR quality of the pairs when taken together.
     *
     * @param reads the list of reads to consider
     * @param samplesList   list of samples
     * @param readsHeader   bam header of reads' source
     * @param setConflictingToZero if true, set base qualities to zero when mates have different base at overlapping position
     * @param halfOfPcrSnvQual half of phred-scaled quality of substitution errors from PCR
     * @param halfOfPcrIndelQual half of phred-scaled quality of indel errors from PCR
     */
    public static void cleanOverlappingReadPairs(final List<GATKRead> reads, final SampleList samplesList, final SAMFileHeader readsHeader,
                                                 final boolean setConflictingToZero, final OptionalInt halfOfPcrSnvQual, final OptionalInt halfOfPcrIndelQual) {
        Utils.nonNull(reads);
        Utils.nonNull(samplesList);
        Utils.nonNull(halfOfPcrSnvQual);
        Utils.nonNull(halfOfPcrSnvQual);
        for ( final List<GATKRead> perSampleReadList : splitReadsBySample(samplesList, readsHeader, reads).values() ) {
            final FragmentCollection<GATKRead> fragmentCollection = FragmentCollection.create(perSampleReadList);
            for ( final Pair<GATKRead, GATKRead> overlappingPair : fragmentCollection.getOverlappingPairs() ) {
                FragmentUtils.adjustQualsOfOverlappingPairedFragments(overlappingPair, setConflictingToZero, halfOfPcrSnvQual, halfOfPcrIndelQual);
            }
        }
    }

    public static Map<String, List<GATKRead>> splitReadsBySample( final SampleList samplesList, final SAMFileHeader header, final Collection<GATKRead> reads ) {
        final Map<String, List<GATKRead>> returnMap = new HashMap<>();
        for (final String sample : samplesList.asListOfSamples()) {
            returnMap.put(sample, new ArrayList<>());
        }

        for ( final GATKRead read : reads ) {
            returnMap.get(ReadUtils.getSampleName(read, header)).add(read);
        }

        return returnMap;
    }

    /**
     * Helper function to create the reference haplotype out of the active region and a padded loc
     * @param region the active region from which to generate the reference haplotype
     * @param paddedReferenceLoc the interval which includes padding and shows how big the reference haplotype should be
     * @return a non-null haplotype
     */
    public static Haplotype createReferenceHaplotype(final AssemblyRegion region, final SimpleInterval paddedReferenceLoc, final ReferenceSequenceFile referenceReader) {
        return ReferenceConfidenceModel.createReferenceHaplotype(region, region.getAssemblyRegionReference(referenceReader), paddedReferenceLoc);
    }

    public static SimpleInterval getPaddedReferenceLoc(final AssemblyRegion region, final int referencePadding, final ReferenceSequenceFile referenceReader) {
        final int padLeft = Math.max(region.getPaddedSpan().getStart() - referencePadding, 1);
        final int padRight = Math.min(region.getPaddedSpan().getEnd() + referencePadding, referenceReader.getSequenceDictionary().getSequence(region.getPaddedSpan().getContig()).getSequenceLength());
        return new SimpleInterval(region.getPaddedSpan().getContig(), padLeft, padRight);
    }

    public static CachingIndexedFastaSequenceFile createReferenceReader(final GATKPath referenceInput) {
        // fasta reference reader to supplement the edges of the reference sequence
        return new CachingIndexedFastaSequenceFile(referenceInput.toPath());
    }

    /**
     * Instantiates the appropriate likelihood calculation engine.
     *
     * @return never {@code null}.
     */
    public static ReadLikelihoodCalculationEngine createLikelihoodCalculationEngine(final LikelihoodEngineArgumentCollection likelihoodArgs) {
        //AlleleLikelihoods::normalizeLikelihoods uses Double.NEGATIVE_INFINITY as a flag to disable capping
        final double log10GlobalReadMismappingRate = likelihoodArgs.phredScaledGlobalReadMismappingRate < 0 ? Double.NEGATIVE_INFINITY
                : QualityUtils.qualToErrorProbLog10(likelihoodArgs.phredScaledGlobalReadMismappingRate);

        return new PairHMMLikelihoodCalculationEngine((byte) likelihoodArgs.gcpHMM, likelihoodArgs.pairHMMNativeArgs.getPairHMMArgs(), likelihoodArgs.pairHMM, log10GlobalReadMismappingRate, likelihoodArgs.pcrErrorModel, likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD);
    }

    public static Optional<HaplotypeBAMWriter> createBamWriter(final AssemblyBasedCallerArgumentCollection args,
                                                               final boolean createBamOutIndex,
                                                               final boolean createBamOutMD5,
                                                               final SAMFileHeader header) {
        return args.bamOutputPath != null ?
                Optional.of(new HaplotypeBAMWriter(args.bamWriterType, IOUtils.getPath(args.bamOutputPath), createBamOutIndex, createBamOutMD5, header)) :
                Optional.empty();
    }

    // Contract: the List<Allele> alleles of the resulting VariantContext is the ref allele followed by alt alleles in the
    // same order as in the input vcs
    public static VariantContext makeMergedVariantContext(final List<VariantContext> vcs) {
        if (vcs.isEmpty()) {
            return null;
        }
        final List<String> haplotypeSources = vcs.stream().map(VariantContext::getSource).collect(Collectors.toList());
        return GATKVariantContextUtils.simpleMerge(vcs, haplotypeSources,
                GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, false);
    }


    /**
     * High-level function that runs the assembler on the given region's reads,
     * returning a data structure with the resulting information needed
     * for further HC steps
     */
    public static AssemblyResultSet assembleReads(final AssemblyRegion region,
                                                  final List<VariantContext> givenAlleles,
                                                  final AssemblyBasedCallerArgumentCollection argumentCollection,
                                                  final SAMFileHeader header,
                                                  final SampleList sampleList,
                                                  final Logger logger,
                                                  final ReferenceSequenceFile referenceReader,
                                                  final ReadThreadingAssembler assemblyEngine,
                                                  final SmithWatermanAligner aligner,
                                                  final boolean correctOverlappingBaseQualities){
        finalizeRegion(region, argumentCollection.assemblerArgs.errorCorrectReads, argumentCollection.dontUseSoftClippedBases, (byte)(argumentCollection.minBaseQualityScore - 1), header, sampleList, correctOverlappingBaseQualities);
        if( argumentCollection.assemblerArgs.debugAssembly) {
            logger.info("Assembling " + region.getSpan() + " with " + region.size() + " reads:    (with overlap region = " + region.getPaddedSpan() + ")");
        }

        final byte[] fullReferenceWithPadding = region.getAssemblyRegionReference(referenceReader, REFERENCE_PADDING_FOR_ASSEMBLY);
        final SimpleInterval paddedReferenceLoc = getPaddedReferenceLoc(region, REFERENCE_PADDING_FOR_ASSEMBLY, referenceReader);
        final Haplotype refHaplotype = createReferenceHaplotype(region, paddedReferenceLoc, referenceReader);

        final ReadErrorCorrector readErrorCorrector = argumentCollection.assemblerArgs.pileupErrorCorrectionLogOdds == Double.NEGATIVE_INFINITY ?
                (argumentCollection.assemblerArgs.errorCorrectReads ?
                        new NearbyKmerErrorCorrector(argumentCollection.assemblerArgs.kmerLengthForReadErrorCorrection,
                                HaplotypeCallerEngine.MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION,
                                argumentCollection.assemblerArgs.minObservationsForKmerToBeSolid,
                                argumentCollection.assemblerArgs.debugAssembly,
                                fullReferenceWithPadding) :
                        null)
                : new PileupReadErrorCorrector(argumentCollection.assemblerArgs.pileupErrorCorrectionLogOdds, header);
        try {
            final AssemblyResultSet assemblyResultSet = assemblyEngine.runLocalAssembly(region, refHaplotype, fullReferenceWithPadding,
                    paddedReferenceLoc, readErrorCorrector, header, aligner);
            if (!givenAlleles.isEmpty()) {
                addGivenAlleles(region.getPaddedSpan().getStart(), givenAlleles, argumentCollection.maxMnpDistance, aligner, refHaplotype, assemblyResultSet);
            }

            assemblyResultSet.setDebug(argumentCollection.assemblerArgs.debugAssembly);
            assemblyResultSet.debugDump(logger);
            return assemblyResultSet;
        } catch (final Exception e){
            // Capture any exception that might be thrown, and write out the assembly failure BAM if requested
            if (argumentCollection.assemblerArgs.captureAssemblyFailureBAM){
                try (final SAMFileWriter writer = ReadUtils.createCommonSAMWriter(new File("assemblyFailure.bam"), null, header, false, false, false)){
                    for (final GATKRead read : region.getReads()) {
                        writer.addAlignment(read.convertToSAMRecord(header));
                    }
                }
            }
            throw e;
        }
    }

    @VisibleForTesting
    static void addGivenAlleles(final int assemblyRegionStart, final List<VariantContext> givenAlleles, final int maxMnpDistance,
                                final SmithWatermanAligner aligner, final Haplotype refHaplotype, final AssemblyResultSet assemblyResultSet) {
        final int activeRegionStart = refHaplotype.getAlignmentStartHapwrtRef();
        final Map<Integer, VariantContext> assembledVariants = assemblyResultSet.getVariationEvents(maxMnpDistance).stream()
                .collect(Collectors.groupingBy(VariantContext::getStart, Collectors.collectingAndThen(Collectors.toList(), AssemblyBasedCallerUtils::makeMergedVariantContext)));

        final List<Haplotype> assembledHaplotypes = assemblyResultSet.getHaplotypeList();
        for (final VariantContext givenVC : givenAlleles) {
            final VariantContext assembledVC = assembledVariants.get(givenVC.getStart());
            final int givenVCRefLength = givenVC.getReference().length();
            final Allele longerRef = (assembledVC == null || givenVCRefLength > assembledVC.getReference().length()) ? givenVC.getReference() : assembledVC.getReference();
            final List<Allele> unassembledGivenAlleles;
            if (assembledVC == null) {
                unassembledGivenAlleles = givenVC.getAlternateAlleles();
            } else {
                // map all alleles to the longest common reference
                final Set<Allele> assembledAlleleSet = new HashSet<>(longerRef.length() == assembledVC.getReference().length() ? assembledVC.getAlternateAlleles() :
                        ReferenceConfidenceVariantContextMerger.remapAlleles(assembledVC, longerRef));
                final Set<Allele> givenAlleleSet = new HashSet<>(longerRef.length() == givenVCRefLength ? givenVC.getAlternateAlleles() :
                        ReferenceConfidenceVariantContextMerger.remapAlleles(givenVC, longerRef));
                unassembledGivenAlleles = givenAlleleSet.stream().filter(a -> !assembledAlleleSet.contains(a)).collect(Collectors.toList());
            }

            final List<Allele> unassembledNonSymbolicAlleles = unassembledGivenAlleles.stream().filter(a -> {
                final byte[] bases = a.getBases();
                return !(Allele.wouldBeNoCallAllele(bases) || Allele.wouldBeNullAllele(bases) || Allele.wouldBeStarAllele(bases) || Allele.wouldBeSymbolicAllele(bases));
            }).collect(Collectors.toList());

            // choose the highest-scoring haplotypes along with the reference for building force-calling haplotypes
            final List<Haplotype> baseHaplotypes = unassembledNonSymbolicAlleles.isEmpty() ? Collections.emptyList() : assembledHaplotypes.stream()
                    .sorted(Comparator.comparingInt((Haplotype hap) -> hap.isReference() ? 1 : 0).thenComparingDouble(hap -> hap.getScore()).reversed())
                    .limit(NUM_HAPLOTYPES_TO_INJECT_FORCE_CALLING_ALLELES_INTO)
                    .collect(Collectors.toList());

            for (final Allele givenAllele : unassembledNonSymbolicAlleles) {
                for (final Haplotype baseHaplotype : baseHaplotypes) {
                    // make sure this allele doesn't collide with a variant on the haplotype
                    if (baseHaplotype.getEventMap()!= null && baseHaplotype.getEventMap().getVariantContexts().stream().anyMatch(vc -> vc.overlaps(givenVC))) {
                        continue;
                    }

                    final Haplotype insertedHaplotype = baseHaplotype.insertAllele(longerRef, givenAllele, activeRegionStart + givenVC.getStart() - assemblyRegionStart, givenVC.getStart());
                    if (insertedHaplotype != null) { // can be null if the requested allele can't be inserted into the haplotype
                        final Cigar cigar = CigarUtils.calculateCigar(refHaplotype.getBases(), insertedHaplotype.getBases(), aligner, SWOverhangStrategy.INDEL);
                        insertedHaplotype.setCigar(cigar);
                        insertedHaplotype.setGenomeLocation(refHaplotype.getGenomeLocation());
                        insertedHaplotype.setAlignmentStartHapwrtRef(activeRegionStart);
                        assemblyResultSet.add(insertedHaplotype);
                    }
                }
            }
        }
        assemblyResultSet.regenerateVariationEvents(maxMnpDistance);
    }

    /**
     * Annotates reads in AlleleLikelihoods with alignment region (the ref region spanned by the haplotype the read is aligned to) and
     * callable region (the ref region over which a caller is using these AlleleLikelihoods to call variants)
     *
     * @param likelihoods AlleleLikelihoods containing reads to be annotated along with haplotypes to which these reads have been aligned
     * @param callableRegion ref region over which caller is using these AlleleLikelihoods to call variants
     */
    public static void annotateReadLikelihoodsWithRegions(final AlleleLikelihoods<GATKRead, Haplotype> likelihoods,
                                                          final Locatable callableRegion) {
        //assign alignment regions to each read
        final Collection<AlleleLikelihoods<GATKRead, Haplotype>.BestAllele> bestHaplotypes = likelihoods.bestAllelesBreakingTies(HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY);
        bestHaplotypes.forEach(bh -> bh.evidence.setAttribute(ALIGNMENT_REGION_TAG, bh.allele.getGenomeLocation().toString()));
        for (final AlleleLikelihoods<GATKRead, Haplotype>.BestAllele bestHaplotype : bestHaplotypes) {
            final GATKRead read = bestHaplotype.evidence;
            final Haplotype haplotype = bestHaplotype.allele;
            read.setAttribute(ALIGNMENT_REGION_TAG, haplotype.getGenomeLocation().toString());
        }

        //assign callable region to each read
        final int sampleCount = likelihoods.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            for (final GATKRead read : likelihoods.sampleEvidence(i)) {
                read.setAttribute(CALLABLE_REGION_TAG, callableRegion.toString());
            }
        }
    }

    /**
     * For the given variant, reads are annotated with which alleles they support, if any.  If a read already has a
     * supported alleles annotation this additional annotation is appended to the previous annotation, it does not replace it.
     * @param vc The variant for which to annotate the reads
     * @param likelihoodsAllele ReadLiklihoods containing reads to be annotated along with alleles of the variant vc
     */
    public static void annotateReadLikelihoodsWithSupportedAlleles(final VariantContext vc,
                                                                     final AlleleLikelihoods<GATKRead, Allele> likelihoodsAllele) {
        annotateReadLikelihoodsWithSupportedAlleles(vc, likelihoodsAllele, Collections::singletonList);
    }

    /**
     * For the given variant, reads are annotated with which alleles they support, if any.  If a read already has a
     * supported alleles annotation this additional annotation is appended to the previous annotation, it does not replace it.
     * @param vc The variant for which to annotate the reads
     * @param likelihoodsAllele ReadLiklihoods containing reads to be annotated along with alleles of the variant vc
     */
    public static <U extends Locatable> void annotateReadLikelihoodsWithSupportedAlleles(final VariantContext vc, final AlleleLikelihoods<U, Allele> likelihoodsAllele,
                                                                       final Function<U, Collection<GATKRead>> readCollectionFunc) {
        //assign supported alleles to each read
        final Map<Allele, List<Allele>> alleleSubset = vc.getAlleles().stream().collect(Collectors.toMap(a -> a, Arrays::asList));
        final AlleleLikelihoods<U, Allele> subsettedLikelihoods = likelihoodsAllele.marginalize(alleleSubset);
        final Collection<AlleleLikelihoods<U, Allele>.BestAllele> bestAlleles = subsettedLikelihoods.bestAllelesBreakingTies().stream()
                .filter(ba -> ba.isInformative()).collect(Collectors.toList());
        for (AlleleLikelihoods<U, Allele>.BestAllele bestAllele : bestAlleles) {
            final Allele allele = bestAllele.allele;
            for (final GATKRead read : readCollectionFunc.apply(bestAllele.evidence)) {
                final String prevAllelesString = read.hasAttribute(SUPPORTED_ALLELES_TAG) ? read.getAttributeAsString(SUPPORTED_ALLELES_TAG) + ", " : "";
                final String newAllelesString = vc.getContig() + ":" + vc.getStart() + "=" + vc.getAlleleIndex(allele);
                read.setAttribute(SUPPORTED_ALLELES_TAG, prevAllelesString + newAllelesString);
            }
        }
    }

    /*
     * Create a context that maps each read to the reference haplotype with log10 L of 0
     * @param refHaplotype a non-null reference haplotype
     * @param samples a list of all samples
     * @param readsHeader SAM header to use for querying sample name from read
     * @param region the assembly region containing reads
     * @return a placeholder AlleleLikelihoods data structure with likelihoods all set to zero
     */
    public static AlleleLikelihoods<GATKRead, Haplotype> createDummyStratifiedReadMap(final Haplotype refHaplotype,
                                                                          final SampleList samples,
                                                                          final SAMFileHeader readsHeader,
                                                                          final AssemblyRegion region) {
        return new AlleleLikelihoods<>(samples, new IndexedAlleleList<>(refHaplotype),
                splitReadsBySample(samples, readsHeader, region.getReads()));
    }

    /**
     * Get a list of pileups that span the entire active region span, in order, one for each position
     */
    public static List<ReadPileup> getPileupsOverReference(final SAMFileHeader readsHeader,
                                                     final SimpleInterval activeRegionSpan,
                                                     final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                     final SampleList samples) {
        final List<GATKRead> reads = new ArrayList<>(readLikelihoods.sampleEvidence(0));
        reads.sort(new ReadCoordinateComparator(readsHeader));  //because we updated the reads based on the local realignments we have to re-sort or the pileups will be... unpredictable

        final LocusIteratorByState libs = new LocusIteratorByState(reads.iterator(), LocusIteratorByState.NO_DOWNSAMPLING,
                false, samples.asSetOfSamples(), readsHeader, true);

        final int startPos = activeRegionSpan.getStart();
        final List<ReadPileup> pileups = new ArrayList<>(activeRegionSpan.getEnd() - startPos);
        AlignmentContext next = libs.advanceToLocus(startPos, true);
        for ( int curPos = startPos; curPos <= activeRegionSpan.getEnd(); curPos++ ) {
            if ( next != null && next.getLocation().getStart() == curPos ) {
                pileups.add(next.getBasePileup());
                next = libs.hasNext() ? libs.next() : null;
            } else {
                // no data, so we create empty pileups
                pileups.add(new ReadPileup(new SimpleInterval(activeRegionSpan.getContig(), curPos, curPos)));
            }
        }

        return pileups;
    }

    /**
     * Returns the list of given alleles active at this location. This method will include events that span the current
     * location if includeSpanningEvents is set to true; otherwise it will only include events that have loc as their \
     * start position.
     * @param loc The start position we are genotyping
     * @param activeAllelesToGenotype The list of given alleles for the current active region, empty unless we are in GGA mode
     * @param includeSpanningEvents If true, will also return events that span loc
     */
    public static List<VariantContext> getVariantContextsFromGivenAlleles(final int loc,
                                                                          final List<VariantContext> activeAllelesToGenotype,
                                                                          final boolean includeSpanningEvents) {
        final Set<LocationAndAlleles> uniqueLocationsAndAlleles = new HashSet<>();
        final List<VariantContext> results = new ArrayList<>();

        int givenAlleleSourceCount = 0;
        for( final VariantContext givenAlleleVC : activeAllelesToGenotype ) {
            if( givenAlleleVC.getStart() <= loc && givenAlleleVC.getEnd() >= loc) {
                if (! (includeSpanningEvents || givenAlleleVC.getStart() == loc)) {
                    continue;
                }
                int alleleCount = 0;
                for( final Allele givenAltAllele : givenAlleleVC.getAlternateAlleles() ) {
                    final List<Allele> alleleSet = Arrays.asList(givenAlleleVC.getReference(), givenAltAllele);

                    //TODO: this source name seems arbitrary and probably just has to be unique
                    //TODO: how about replace it by vcSourceName = String.parseInt(nameCounter++)?
                    final String vcSourceName = "Comp" + givenAlleleSourceCount + "Allele" + alleleCount;
                    // check if this event is already in the list of events due to a repeat in the input alleles track
                    final VariantContext candidateEventToAdd = new VariantContextBuilder(givenAlleleVC).alleles(alleleSet)
                            .genotypes(GenotypesContext.NO_GENOTYPES).source(vcSourceName).make();

                    final LocationAndAlleles locationAndAlleles = new LocationAndAlleles(candidateEventToAdd.getStart(), candidateEventToAdd.getAlleles());
                    if (! uniqueLocationsAndAlleles.contains(locationAndAlleles)) {
                        uniqueLocationsAndAlleles.add(locationAndAlleles);
                        results.add(candidateEventToAdd);
                    }

                    alleleCount++;
                }
            }
            givenAlleleSourceCount++;
        }
        return results;
    }

    /**
     * Returns the list of events discovered in assembled haplotypes that are active at this location. The results will
     * include events that span the current location if includeSpanningEvents is set to true; otherwise it will only
     * include events that have loc as their start position.
     * @param loc The start position we are genotyping
     * @param haplotypes list of active haplotypes at the current location
     * @param includeSpanningEvents If true, will also return events that span loc
     */
    public static List<VariantContext> getVariantContextsFromActiveHaplotypes(final int loc,
                                                                                 final List<Haplotype> haplotypes,
                                                                                 final boolean includeSpanningEvents) {
        final List<VariantContext> results = new ArrayList<>();
        final Set<LocationAndAlleles> uniqueLocationsAndAlleles = new HashSet<>();

        haplotypes.stream()
                .flatMap(h -> Utils.stream(h.getEventMap().getOverlappingEvents(loc)))
                .filter(Objects::nonNull)
                .filter(v -> (includeSpanningEvents || v.getStart() == loc))
                .forEach(v -> {
                    final LocationAndAlleles locationAndAlleles = new LocationAndAlleles(v.getStart(), v.getAlleles());
                    if (! uniqueLocationsAndAlleles.contains(locationAndAlleles)) {
                        uniqueLocationsAndAlleles.add(locationAndAlleles);
                        results.add(v);
                    }
                });
        return results;
    }

    /**
     * Returns a mapping from Allele in the mergedVC, which represents all of the alleles being genotyped at loc,
     * to a list of Haplotypes that support that allele. If the mergedVC includes a spanning deletion allele, all haplotypes
     * that support spanning deletions will be assigned to that allele in the map.
     *
     * @param mergedVC The merged variant context for the locus, which includes all active alternate alleles merged to a single reference allele
     * @param loc The active locus being genotyped
     * @param haplotypes Haplotypes for the current active region
     * @return
     */
    public static Map<Allele, List<Haplotype>> createAlleleMapper(final VariantContext mergedVC, final int loc, final List<Haplotype> haplotypes) {

        final Map<Allele, List<Haplotype>> result = new LinkedHashMap<>();

        final Allele ref = mergedVC.getReference();
        result.put(ref, new ArrayList<>());

        //Note: we can't use the alleles implied by eventsAtThisLoc because they are not yet merged to a common reference
        //For example, a homopolymer AAAAA reference with a single and double deletion would yield (i) AA* A and (ii) AAA* A
        //in eventsAtThisLoc, when in mergedVC it would yield AAA* AA A
        mergedVC.getAlternateAlleles().stream().filter(a -> !a.isSymbolic()).forEach(a -> result.put(a, new ArrayList<>()));

        for (final Haplotype h : haplotypes) {

            final List<VariantContext> spanningEvents = h.getEventMap().getOverlappingEvents(loc);

            if (spanningEvents.isEmpty()) {    //no events --> this haplotype supports the reference at this locus
                result.get(ref).add(h);
                continue;
            }

            for (VariantContext spanningEvent : spanningEvents) {
                if (spanningEvent.getStart() == loc) {
                    // the event starts at the current location

                    if (spanningEvent.getReference().length() == mergedVC.getReference().length()) {
                        // reference allele lengths are equal; we can just use the spanning event's alt allele
                        // in the case of GGA mode the spanning event might not match an allele in the mergedVC
                        if (result.containsKey(spanningEvent.getAlternateAllele(0))) {
                            // variant contexts derived from the event map have only one alt allele each, so we can just
                            // grab the first one (we're not assuming that the sample is biallelic)
                            result.get(spanningEvent.getAlternateAllele(0)).add(h);
                        }
                    } else if (spanningEvent.getReference().length() < mergedVC.getReference().length()) {
                        // spanning event has shorter ref allele than merged VC; we need to pad out its alt allele
                        final Map<Allele, Allele> spanningEventAlleleMappingToMergedVc
                                = GATKVariantContextUtils.createAlleleMapping(mergedVC.getReference(), spanningEvent);
                        final Allele remappedSpanningEventAltAllele = spanningEventAlleleMappingToMergedVc.get(spanningEvent.getAlternateAllele(0));
                        // in the case of GGA mode the spanning event might not match an allele in the mergedVC
                        if (result.containsKey(remappedSpanningEventAltAllele)) {
                            result.get(remappedSpanningEventAltAllele).add(h);
                        }
                    } else {
                        // the process of creating the merged VC in AssemblyBasedCallerUtils::makeMergedVariantContext should have
                        // already padded out the reference allele, therefore this spanning VC must not be in events at this site
                        // because we're in GGA mode and it's not an allele we want
                        continue;
                    }

                } else {
                    // the event starts prior to the current location, so it's a spanning deletion
                    if (! result.containsKey(Allele.SPAN_DEL)) {
                        result.put(Allele.SPAN_DEL, new ArrayList<>());
                    }
                    result.get(Allele.SPAN_DEL).add(h);
                    break;
                }
            }

        }
        return result;
    }

    /**
     * Tries to phase the individual alleles based on pairwise comparisons to the other alleles based on all called haplotypes
     *
     * @param calls             the list of called alleles
     * @param calledHaplotypes  the set of haplotypes used for calling
     * @return a non-null list which represents the possibly phased version of the calls
     */
    public static List<VariantContext> phaseCalls(final List<VariantContext> calls, final Set<Haplotype> calledHaplotypes) {

        // construct a mapping from alternate allele to the set of haplotypes that contain that allele
        final Map<VariantContext, Set<Haplotype>> haplotypeMap = constructHaplotypeMapping(calls, calledHaplotypes);

        // construct a mapping from call to phase set ID
        final Map<VariantContext, Pair<Integer, String>> phaseSetMapping = new HashMap<>();
        final int uniqueCounterEndValue = constructPhaseSetMapping(calls, haplotypeMap, calledHaplotypes.size() - 1, phaseSetMapping);

        // we want to establish (potential) *groups* of phased variants, so we need to be smart when looking at pairwise phasing partners
        return constructPhaseGroups(calls, phaseSetMapping, uniqueCounterEndValue);
    }

    /**
     * Construct the mapping from alternate allele to the set of haplotypes that contain that allele
     *
     * @param originalCalls    the original unphased calls
     * @param calledHaplotypes  the set of haplotypes used for calling
     * @return non-null Map
     */
    @VisibleForTesting
    static Map<VariantContext, Set<Haplotype>> constructHaplotypeMapping(final List<VariantContext> originalCalls,
                                                                                   final Set<Haplotype> calledHaplotypes) {
        final Map<VariantContext, Set<Haplotype>> haplotypeMap = new HashMap<>(originalCalls.size());
        for ( final VariantContext call : originalCalls ) {
            // don't try to phase if there is not exactly 1 alternate allele
            if ( ! isBiallelic(call) ) {
                haplotypeMap.put(call, Collections.<Haplotype>emptySet());
                continue;
            }

            // keep track of the haplotypes that contain this particular alternate allele
            final Allele alt = call.getAlternateAllele(0);
            final Predicate<VariantContext> hasThisAlt = vc -> (vc.getStart() == call.getStart() && vc.getAlternateAlleles().contains(alt)) ||
                    (Allele.SPAN_DEL.equals(alt) && vc.getStart() < call.getStart() && vc.getEnd() >= call.getStart());
            final Set<Haplotype> hapsWithAllele = calledHaplotypes.stream()
                    .filter(h -> h.getEventMap().getVariantContexts().stream().anyMatch(hasThisAlt))
                    .collect(Collectors.toCollection(HashSet<Haplotype>::new));

            haplotypeMap.put(call, hapsWithAllele);
        }

        return haplotypeMap;
    }

    /**
     * Construct the mapping from call (variant context) to phase set ID
     *
     * @param originalCalls    the original unphased calls
     * @param haplotypeMap     mapping from alternate allele to the set of haplotypes that contain that allele
     * @param totalAvailableHaplotypes the total number of possible haplotypes used in calling
     * @param phaseSetMapping  the map to populate in this method;
     *                         note that it is okay for this method NOT to populate the phaseSetMapping at all (e.g. in an impossible-to-phase situation)
     * @return the next incremental unique index
     */
    @VisibleForTesting
    static int constructPhaseSetMapping(final List<VariantContext> originalCalls,
                                                  final Map<VariantContext, Set<Haplotype>> haplotypeMap,
                                                  final int totalAvailableHaplotypes,
                                                  final Map<VariantContext, Pair<Integer, String>> phaseSetMapping) {

        final int numCalls = originalCalls.size();
        int uniqueCounter = 0;

        // use the haplotype mapping to connect variants that are always/never present on the same haplotypes
        for ( int i = 0; i < numCalls - 1; i++ ) {
            final VariantContext call = originalCalls.get(i);
            final Set<Haplotype> haplotypesWithCall = haplotypeMap.get(call);
            if ( haplotypesWithCall.isEmpty() ) {
                continue;
            }

            final boolean callIsOnAllHaps = haplotypesWithCall.size() == totalAvailableHaplotypes;

            for ( int j = i+1; j < numCalls; j++ ) {
                final VariantContext comp = originalCalls.get(j);
                final Set<Haplotype> haplotypesWithComp = haplotypeMap.get(comp);
                if ( haplotypesWithComp.isEmpty() ) {
                    continue;
                }

                // if the variants are together on all haplotypes, record that fact.
                // another possibility is that one of the variants is on all possible haplotypes (i.e. it is homozygous).
                final boolean compIsOnAllHaps = haplotypesWithComp.size() == totalAvailableHaplotypes;
                if ( (haplotypesWithCall.size() == haplotypesWithComp.size() && haplotypesWithCall.containsAll(haplotypesWithComp)) || callIsOnAllHaps || compIsOnAllHaps ) {

                    // create a new group if these are the first entries
                    if ( ! phaseSetMapping.containsKey(call) ) {
                        // note that if the comp is already in the map then that is very bad because it means that there is
                        // another variant that is in phase with the comp but not with the call.  Since that's an un-phasable
                        // situation, we should abort if we encounter it.
                        if ( phaseSetMapping.containsKey(comp) ) {
                            phaseSetMapping.clear();
                            return 0;
                        }

                        // An important note: even for homozygous variants we are setting the phase as "0|1" here.
                        // We do this because we cannot possibly know for sure at this time that the genotype for this
                        // sample will actually be homozygous downstream: there are steps in the pipeline that are liable
                        // to change the genotypes.  Because we can't make those assumptions here, we have decided to output
                        // the phase as if the call is heterozygous and then "fix" it downstream as needed.
                        phaseSetMapping.put(call, Pair.of(uniqueCounter, phase01));
                        phaseSetMapping.put(comp, Pair.of(uniqueCounter, phase01));
                        uniqueCounter++;
                    }
                    // otherwise it's part of an existing group so use that group's unique ID
                    else if ( ! phaseSetMapping.containsKey(comp) ) {
                        final Pair<Integer, String> callPhase = phaseSetMapping.get(call);
                        phaseSetMapping.put(comp, Pair.of(callPhase.getLeft(), callPhase.getRight()));
                    }
                }
                // if the variants are apart on *all* haplotypes, record that fact
                else if ( haplotypesWithCall.size() + haplotypesWithComp.size() == totalAvailableHaplotypes ) {

                    final Set<Haplotype> intersection = new HashSet<>();
                    intersection.addAll(haplotypesWithCall);
                    intersection.retainAll(haplotypesWithComp);
                    if ( intersection.isEmpty() ) {
                        // create a new group if these are the first entries
                        if ( ! phaseSetMapping.containsKey(call) ) {
                            // note that if the comp is already in the map then that is very bad because it means that there is
                            // another variant that is in phase with the comp but not with the call.  Since that's an un-phasable
                            // situation, we should abort if we encounter it.
                            if ( phaseSetMapping.containsKey(comp) ) {
                                phaseSetMapping.clear();
                                return 0;
                            }

                            phaseSetMapping.put(call, Pair.of(uniqueCounter, phase01));
                            phaseSetMapping.put(comp, Pair.of(uniqueCounter, phase10));
                            uniqueCounter++;
                        }
                        // otherwise it's part of an existing group so use that group's unique ID
                        else if ( ! phaseSetMapping.containsKey(comp) ){
                            final Pair<Integer, String> callPhase = phaseSetMapping.get(call);
                            phaseSetMapping.put(comp, Pair.of(callPhase.getLeft(), callPhase.getRight().equals(phase01) ? phase10 : phase01));
                        }
                    }
                }
            }
        }

        return uniqueCounter;
    }

    /**
     * Assemble the phase groups together and update the original calls accordingly
     *
     * @param originalCalls    the original unphased calls
     * @param phaseSetMapping  mapping from call (variant context) to phase group ID
     * @param indexTo          last index (exclusive) of phase group IDs
     * @return a non-null list which represents the possibly phased version of the calls
     */
    @VisibleForTesting
    static List<VariantContext> constructPhaseGroups(final List<VariantContext> originalCalls,
                                                               final Map<VariantContext, Pair<Integer, String>> phaseSetMapping,
                                                               final int indexTo) {
        final List<VariantContext> phasedCalls = new ArrayList<>(originalCalls);

        // if we managed to find any phased groups, update the VariantContexts
        for ( int count = 0; count < indexTo; count++ ) {
            // get all of the (indexes of the) calls that belong in this group (keeping them in the original order)
            final List<Integer> indexes = new ArrayList<>();
            for ( int index = 0; index < originalCalls.size(); index++ ) {
                final VariantContext call = originalCalls.get(index);
                if ( phaseSetMapping.containsKey(call) && phaseSetMapping.get(call).getLeft() == count ) {
                    indexes.add(index);
                }
            }
            if ( indexes.size() < 2 ) {
                throw new IllegalStateException("Somehow we have a group of phased variants that has fewer than 2 members");
            }

            // create a unique ID based on the leftmost one
            final String uniqueID = createUniqueID(originalCalls.get(indexes.get(0)));

            // create the phase set identifier, which is the position of the first variant in the set
            final int phaseSetID = originalCalls.get(indexes.get(0)).getStart();

            // update the VCs
            for ( final int index : indexes ) {
                final VariantContext originalCall = originalCalls.get(index);
                final VariantContext phasedCall = phaseVC(originalCall, uniqueID, phaseSetMapping.get(originalCall).getRight(), phaseSetID);
                phasedCalls.set(index, phasedCall);
            }
        }

        return phasedCalls;
    }

    /**
     * Is this variant bi-allelic?  This implementation is very much specific to this class so shouldn't be pulled out into a generalized place.
     *
     * @param vc the variant context
     * @return true if this variant context is bi-allelic, ignoring the NON-REF symbolic allele, false otherwise
     */
    private static boolean isBiallelic(final VariantContext vc) {
        return vc.isBiallelic() || (vc.getNAlleles() == 3 && vc.getAlternateAlleles().contains(Allele.NON_REF_ALLELE));
    }

    /**
     * Create a unique identifier given the variant context
     *
     * @param vc   the variant context
     * @return non-null String
     */
    private static String createUniqueID(final VariantContext vc) {
        return String.format("%d_%s_%s", vc.getStart(), vc.getReference().getDisplayString(), vc.getAlternateAllele(0).getDisplayString());
    }

    /**
     * Add physical phase information to the provided variant context
     *
     * @param vc   the variant context
     * @param ID   the ID to use
     * @param phaseGT the phase GT string to use
     * @return phased non-null variant context
     */
    private static VariantContext phaseVC(final VariantContext vc, final String ID, final String phaseGT, final int phaseSetID) {
        final List<Genotype> phasedGenotypes = new ArrayList<>();
        for ( final Genotype g : vc.getGenotypes() ) {
            final List<Allele> alleles = g.getAlleles();
            if (phaseGT.equals(phase10) && g.isHet()) Collections.reverse(alleles); // swap the alleles if heterozygous
            final Genotype genotype = new GenotypeBuilder(g)
                .alleles(alleles)
                .phased(true)
                .attribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, ID)
                .attribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, phaseGT)
                .attribute(VCFConstants.PHASE_SET_KEY, phaseSetID)
                .make();
            phasedGenotypes.add(genotype);
        }
        return new VariantContextBuilder(vc).genotypes(phasedGenotypes).make();
    }

    public static Set<Allele> getAllelesConsistentWithGivenAlleles(final List<VariantContext> givenAlleles, final VariantContext mergedVC) {
        if (givenAlleles.isEmpty()) {
            return Collections.emptySet();
        }

        final List<Pair<Allele, Allele>> givenAltAndRefAllelesInOriginalContext =  getVariantContextsFromGivenAlleles(mergedVC.getStart(), givenAlleles, false).stream()
                .flatMap(vc -> vc.getAlternateAlleles().stream().map(allele -> ImmutablePair.of(allele, vc.getReference()))).collect(Collectors.toList());

        return mergedVC.getAlternateAlleles().stream()
                .map(allele -> ImmutablePair.of(allele, mergedVC.getReference()))
                .filter(altAndRef -> givenAltAndRefAllelesInOriginalContext.stream().anyMatch(givenAltAndRef -> allelesAreConsistent(givenAltAndRef, altAndRef)))
                .map(altAndRefPair -> altAndRefPair.getLeft())
                .collect(Collectors.toSet());
    }

    // check whether two alleles coming from different variant contexts and with possibly different reference alleles
    // could in fact be the same.  The condition is that one is a prefix of the other
    private static boolean allelesAreConsistent(final Pair<Allele,Allele> altAndRef1, final Pair<Allele,Allele> altAndRef2) {
        final Allele alt1 = altAndRef1.getLeft();
        final Allele alt2 = altAndRef2.getLeft();
        if (alt1.isSymbolic() || alt2.isSymbolic()) {
            return false;
        } else {
            final int sizeDiff1 = alt1.length() - altAndRef1.getRight().length();
            final int sizeDiff2 = alt2.length() - altAndRef2.getRight().length();
            return (sizeDiff1 == sizeDiff2) && (alt1.length() < alt2.length() ?
                    alt1.basesMatch(Arrays.copyOf(alt2.getBases(), alt1.length())) :
                    alt2.basesMatch(Arrays.copyOf(alt1.getBases(), alt2.length())));
        }
    }
}
