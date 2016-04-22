package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

/**
 * Utility class that error-corrects reads.
 * Main idea: An error in a read will appear as a bubble in a k-mer (de Bruijn) graph and such bubble will have very low multiplicity.
 * Hence, read errors will appear as "sparse" kmers with very little support.
 * Historically, the most common approach to error-correct reads before assembly has been to first compute the kmer spectrum of the reads,
 * defined as the kmer composition of a set of reads along with the multiplicity of each kmer.
 * First-generation correctors like the Euler corrector (Pevzner 2001) mapped low frequency kmers (kmers appearing say below N times)
 * into high frequency ones that lied within a certain Hamming or edit distance.
 * This is doable, but has some drawbacks:
 * - Kmers used for error correction become tied to kmers used for graph building.
 * - Hence, large kmers (desireable for graph building because they can resolve repeats better) are a hindrance for error correction,
 * because they are seen less often.
 * - After error correction, there is no guarantee that a sequence of kmers corresponds to an "actual" read.
 *
 * An error-corrected set of reads also makes a much smoother graph without the need to resolving so many bubbles.
 *
 * Idea hence is to correct reads based on their kmer content, but in a context independent from graph building.
 * In order to do this, the following steps are taken:
 * - The k-mer spectrum of a set of reads in computed. However, we are at freedom to choose the most convenient k-mer size (typicially around
 * read length /2).
 * - We partition the set of observed k-mers into "solid" kmers which have multiplicity > M, and "insolid" ones otherwise (Pevzner 2001).
 *
 * - Main idea of the algorithm is to try to substitute a sequence of bases in a read by a sequence better supported by kmers.
 * - For each "unsolid" kmer observed in reads, we try to find a "solid" kmer within a maximum Hamming distance.
 * - If such solid kmer exists, then this unsolid kmer is "correctable", otherwise, uncorrectable.
 * - For each read, then:
 * -- Walk through  read and visit all kmers.
 * -- If kmer is solid, continue to next kmer.
 * -- If not, and if it's correctable (i.e. there exists a mapping from an unsolid kmer to a solid kmer within a given Hamming distance),
 *    add the bases and offsets corresponding to differing positions between unsolid and solid kmer to correction list.
 * -- At the end, each base in read will have a list of corrections associated with it. We can then choose to correct or not.
 *    If read has only consistent corrections, then we can correct base to common base in corrections.
 *
 *    TODO:
 *    todo Q: WHAT QUALITY TO USE??
 *    todo how do we deal with mate pairs?
 *
 *
*/
public final class ReadErrorCorrector {
    private static final Logger logger = LogManager.getLogger(ReadErrorCorrector.class);
    /**
     * A map of for each kmer to its num occurrences in addKmers
     */
    final KMerCounter countsByKMer;

    private final Map<Kmer,Kmer> kmerCorrectionMap = new HashMap<>();
    private final Map<Kmer,Pair<int[],byte[]>> kmerDifferingBases = new HashMap<>();
    private final int kmerLength;
    private final boolean debug;
    private final boolean trimLowQualityBases;
    private final byte minTailQuality;
    private final int maxMismatchesToCorrect;
    private final byte qualityOfCorrectedBases;
    private final int maxObservationsForKmerToBeCorrectable;
    private final int maxHomopolymerLengthInRegion;
    private final int minObservationsForKmerToBeSolid;

    // default values, for debugging
    private static final boolean doInplaceErrorCorrection = false;    // currently not used, since we want corrected reads to be used only for assembly
    private static final int MAX_MISMATCHES_TO_CORRECT = 2;
    private static final byte QUALITY_OF_CORRECTED_BASES = 30; // what's a reasonable value here?
    private static final int MAX_OBSERVATIONS_FOR_KMER_TO_BE_CORRECTABLE = 1;
    private static final boolean TRIM_LOW_QUAL_TAILS = false;
    private static final boolean DONT_CORRECT_IN_LONG_HOMOPOLYMERS = false;
    private static final int MAX_HOMOPOLYMER_THRESHOLD = 12;

    // debug counter structure
    private final ReadErrorCorrectionStats readErrorCorrectionStats = new ReadErrorCorrectionStats();

    /**
     * Create a new kmer corrector
     *
     * @param kmerLength the length of kmers we'll be counting to error correct, must be >= 1
     * @param maxMismatchesToCorrect e >= 0
     * @param qualityOfCorrectedBases  Bases to be corrected will be assigned this quality
     */
    public ReadErrorCorrector(final int kmerLength,
                              final int maxMismatchesToCorrect,
                              final int maxObservationsForKmerToBeCorrectable,
                              final byte qualityOfCorrectedBases,
                              final int minObservationsForKmerToBeSolid,
                              final boolean trimLowQualityBases,
                              final byte minTailQuality,
                              final boolean debug,
                              final byte[] fullReferenceWithPadding) {
        if ( kmerLength < 1 ) {
            throw new IllegalArgumentException("kmerLength must be > 0 but got " + kmerLength);
        }
        if ( maxMismatchesToCorrect < 1 ) {
            throw new IllegalArgumentException("maxMismatchesToCorrect must be >= 1 but got " + maxMismatchesToCorrect);
        }
        if ( qualityOfCorrectedBases < 2 || qualityOfCorrectedBases > QualityUtils.MAX_REASONABLE_Q_SCORE) {
            throw new IllegalArgumentException("qualityOfCorrectedBases must be >= 2 and <= MAX_REASONABLE_Q_SCORE but got " + qualityOfCorrectedBases);
        }

        countsByKMer = new KMerCounter(kmerLength);
        this.kmerLength = kmerLength;
        this.maxMismatchesToCorrect = maxMismatchesToCorrect;
        this.qualityOfCorrectedBases = qualityOfCorrectedBases;
        this.minObservationsForKmerToBeSolid = minObservationsForKmerToBeSolid;
        this.trimLowQualityBases = trimLowQualityBases;
        this.minTailQuality = minTailQuality;
        this.debug = debug;
        this.maxObservationsForKmerToBeCorrectable = maxObservationsForKmerToBeCorrectable;

        // when region has long homopolymers, we may want not to correct reads, since assessment is complicated,
        // so we may decide to skip error correction in these regions
        maxHomopolymerLengthInRegion = computeMaxHLen(fullReferenceWithPadding);
    }

    /**
     * Simple constructor with sensible defaults
     * @param kmerLength            K-mer length for error correction (not necessarily the same as for assembly graph)
     * @param minTailQuality        Minimum tail quality: remaining bases with Q's below this value are hard-clipped after correction
     * @param debug                 Output debug information
     */
    public ReadErrorCorrector(final int kmerLength, final byte minTailQuality, final int minObservationsForKmerToBeSolid, final boolean debug, final byte[] fullReferenceWithPadding) {
        this(kmerLength, MAX_MISMATCHES_TO_CORRECT, MAX_OBSERVATIONS_FOR_KMER_TO_BE_CORRECTABLE, QUALITY_OF_CORRECTED_BASES, minObservationsForKmerToBeSolid, TRIM_LOW_QUAL_TAILS, minTailQuality, debug,fullReferenceWithPadding);
    }

    /**
    * Main entry routine to add all kmers in a read to the read map counter
    * @param read                        Read to add bases
    */
    protected void addReadKmers(final GATKRead read) {
        Utils.nonNull(read);
        if (DONT_CORRECT_IN_LONG_HOMOPOLYMERS && maxHomopolymerLengthInRegion > MAX_HOMOPOLYMER_THRESHOLD) {
            return;
        }

        final byte[] readBases = read.getBases();
        for (int offset = 0; offset <= readBases.length-kmerLength; offset++ )  {
            countsByKMer.addKmer(new Kmer(readBases,offset,kmerLength),1);
        }
    }

    /**
     * Correct a collection of reads based on stored k-mer counts
     * @param reads
     */
    public final List<GATKRead> correctReads(final Collection<GATKRead> reads) {

        final List<GATKRead> correctedReads = new ArrayList<>(reads.size());
        if (DONT_CORRECT_IN_LONG_HOMOPOLYMERS && maxHomopolymerLengthInRegion > MAX_HOMOPOLYMER_THRESHOLD) {
            // just copy reads into output and exit
            correctedReads.addAll(reads);
        }
        else {
            computeKmerCorrectionMap();
            for (final GATKRead read: reads) {
                final GATKRead correctedRead = correctRead(read);
                if (trimLowQualityBases) {
                    correctedReads.add(ReadClipper.hardClipLowQualEnds(correctedRead, minTailQuality));
                } else {
                    correctedReads.add(correctedRead);
                }
            }
            if (debug) {
                logger.info("Number of corrected bases:"+readErrorCorrectionStats.numBasesCorrected);
                logger.info("Number of corrected reads:"+readErrorCorrectionStats.numReadsCorrected);
                logger.info("Number of skipped reads:"+readErrorCorrectionStats.numReadsUncorrected);
                logger.info("Number of solid kmers:"+readErrorCorrectionStats.numSolidKmers);
                logger.info("Number of corrected kmers:"+readErrorCorrectionStats.numCorrectedKmers);
                logger.info("Number of uncorrectable kmers:"+readErrorCorrectionStats.numUncorrectableKmers);
            }
        }
        return correctedReads;
    }


    /**
     * Do actual read correction based on k-mer map. First, loop through stored k-mers to get a list of possible corrections
     * for each position in the read. Then correct read based on all possible consistent corrections.
     * @param inputRead                               Read to correct
     * @return                                        Corrected read (can be same reference as input if doInplaceErrorCorrection is set)
     */
    private GATKRead correctRead(final GATKRead inputRead) {
        Utils.nonNull(inputRead);
        // do actual correction
        boolean corrected = false;
        final byte[] correctedBases = inputRead.getBases();
        final byte[] correctedQuals = inputRead.getBaseQualities();

        // array to store list of possible corrections for read
        final CorrectionSet correctionSet = buildCorrectionMap(correctedBases);

        for (int offset = 0; offset < correctedBases.length; offset++) {
            final Byte b = correctionSet.getConsensusCorrection(offset);
            if (b != null && b != correctedBases[offset]) {
                correctedBases[offset] = b;
                correctedQuals[offset] = qualityOfCorrectedBases;
                corrected = true;
            }
            readErrorCorrectionStats.numBasesCorrected++;
        }

        if (corrected) {
            readErrorCorrectionStats.numReadsCorrected++;
            if (doInplaceErrorCorrection) {
                inputRead.setBases(correctedBases);
                inputRead.setBaseQualities(correctedQuals);
                return inputRead;
            }
            else {
                final GATKRead correctedRead = inputRead.deepCopy();

                //  do the actual correction
                // todo - do we need to clone anything else from read?
                correctedRead.setBaseQualities(inputRead.getBaseQualities());
                correctedRead.setBases(inputRead.getBases());
                correctedRead.setReadGroup(inputRead.getReadGroup());
                return correctedRead;
            }
        }
        else {
            readErrorCorrectionStats.numReadsUncorrected++;
            return inputRead;
        }
    }

    /**
     * Build correction map for each of the bases in read.
     * For each of the constituent kmers in read:
     * a) See whether the kmer has been mapped to a corrected kmer.
     * b) If so, get list of differing positions and corresponding bases.
     * c) Add then list of new bases to index in correction list.
     * Correction list is of read size, and holds a list of bases to correct.
     * @param correctedBases                        Bases to attempt to correct
     * @return                                      CorrectionSet object.
     */
    private CorrectionSet buildCorrectionMap(final byte[] correctedBases) {
        Utils.nonNull(correctedBases);
        // array to store list of possible corrections for read
        final CorrectionSet correctionSet = new CorrectionSet(correctedBases.length);

        for (int offset = 0; offset <= correctedBases.length-kmerLength; offset++ )  {
            final Kmer kmer = new Kmer(correctedBases,offset,kmerLength);
            final Kmer newKmer = kmerCorrectionMap.get(kmer);
            if (newKmer != null && !newKmer.equals(kmer)){
                final Pair<int[],byte[]> differingPositions = kmerDifferingBases.get(kmer);
                final int[] differingIndeces = differingPositions.getLeft();
                final byte[] differingBases = differingPositions.getRight();

                for (int k=0; k < differingIndeces.length; k++) {
                    // get list of differing positions for corrected kmer
                    // for each of these, add correction candidate to correction set
                    correctionSet.add(offset + differingIndeces[k],differingBases[k]);
                }
            }
        }
        return correctionSet;
    }


    /**
     * Top-level entry point that adds a collection of reads to our kmer list.
     * For each read in list, its constituent kmers will be logged in our kmer table.
     * @param reads
     */
    public void addReadsToKmers(final Collection<GATKRead> reads) {
        Utils.nonNull(reads);
        for (final GATKRead read: reads) {
            addReadKmers(read);
        }

        if (debug) {
            for (final KMerCounter.CountedKmer countedKmer : countsByKMer.getCountedKmers()) {
                logger.info(String.format("%s\t%d\n", countedKmer.kmer, countedKmer.count));
            }
        }
    }


    /**
     * For each kmer we've seen, do the following:
     * a) If kmer count > threshold1, this kmer is good, so correction map will be to itself.
     * b) If kmer count <= threshold2, this kmer is bad.
     *    In that case, loop through all other kmers. If kmer is good, compute distance, and get minimal distance.
     *    If such distance is < some threshold, map to this kmer, and record differing positions and bases.
     *
     */
    private void computeKmerCorrectionMap() {
        for (final KMerCounter.CountedKmer storedKmer : countsByKMer.getCountedKmers()) {
            if (storedKmer.getCount() >= minObservationsForKmerToBeSolid) {
                // this kmer is good: map to itself
                kmerCorrectionMap.put(storedKmer.getKmer(),storedKmer.getKmer());
                kmerDifferingBases.put(storedKmer.getKmer(), Pair.of(new int[0],new byte[0])); // dummy empty array
                readErrorCorrectionStats.numSolidKmers++;
            }
            else if (storedKmer.getCount() <= maxObservationsForKmerToBeCorrectable) {
                // loop now thru all other kmers to find nearest neighbor
                final Pair<Kmer,Pair<int[],byte[]>> nearestNeighbor = findNearestNeighbor(storedKmer.getKmer(),countsByKMer,maxMismatchesToCorrect);

                // check if nearest neighbor lies in a close vicinity. If so, log the new bases and the correction map
                if (nearestNeighbor != null) { // ok, found close neighbor
                    kmerCorrectionMap.put(storedKmer.getKmer(), nearestNeighbor.getLeft());
                    kmerDifferingBases.put(storedKmer.getKmer(), nearestNeighbor.getRight());
                    readErrorCorrectionStats.numCorrectedKmers++;
//                    if (debug)
//                        logger.info("Original kmer:"+storedKmer + "\tCorrected kmer:"+nearestNeighbor.first+"\tDistance:"+dist);
                }
                else {
                    readErrorCorrectionStats.numUncorrectableKmers++;
                }

            }
         }
    }

    /**
     * Finds nearest neighbor of a given k-mer, among a list of counted K-mers, up to a given distance.
     * If many k-mers share same closest distance, an arbitrary k-mer is picked
     * @param kmer                        K-mer of interest
     * @param countsByKMer                KMerCounter storing set of counted k-mers (may include kmer of interest)
     * @param maxDistance                 Maximum distance to search
     * @return                            Pair of values: closest K-mer in Hamming distance and list of differing bases.
     *                                      If no neighbor can be found up to given distance, returns null
     */
    private Pair<Kmer,Pair<int[],byte[]>> findNearestNeighbor(final Kmer kmer,
                                                             final KMerCounter  countsByKMer,
                                                             final int maxDistance) {
        Utils.nonNull(kmer, "KMER");
        Utils.nonNull(countsByKMer, "countsByKMer");
        Utils.validateArg(maxDistance >= 1, "countsByKMer");

        int minimumDistance = Integer.MAX_VALUE;
        Kmer closestKmer = null;

        final int[] differingIndeces = new int[maxDistance+1];
        final byte[] differingBases = new byte[maxDistance+1];

        final int[] closestDifferingIndices = new int[maxDistance+1];
        final byte[] closestDifferingBases = new byte[maxDistance+1];

        for (final KMerCounter.CountedKmer candidateKmer : countsByKMer.getCountedKmers()) {
            // skip if candidate set includes test kmer
            if (candidateKmer.getKmer().equals(kmer)) {
                continue;
            }

            final int hammingDistance  = kmer.getDifferingPositions(candidateKmer.getKmer(), maxDistance, differingIndeces, differingBases);
            if (hammingDistance < 0) // can't compare kmer? skip
            {
                continue;
            }

            if (hammingDistance < minimumDistance)  {
                minimumDistance = hammingDistance;
                closestKmer = candidateKmer.getKmer();
                System.arraycopy(differingBases,0,closestDifferingBases,0,differingBases.length);
                System.arraycopy(differingIndeces,0,closestDifferingIndices,0,differingIndeces.length);
            }
        }
        return Pair.of(closestKmer, Pair.of(closestDifferingIndices,closestDifferingBases));
    }


    /**
     * experimental function to compute max homopolymer length in a given reference context
     * @param fullReferenceWithPadding                Reference context of interest
     * @return                                        Max homopolymer length in region
     */
    private static int computeMaxHLen(final byte[] fullReferenceWithPadding) {
        Utils.nonNull(fullReferenceWithPadding);
        int leftRun = 1;
        int maxRun = 1;
        for ( int i = 1; i < fullReferenceWithPadding.length; i++) {
            if ( fullReferenceWithPadding[i] == fullReferenceWithPadding[i-1] ) {
                leftRun++;
            } else {
                leftRun = 1;
            }
            }
            if (leftRun > maxRun) {
                maxRun = leftRun;
            }


        return maxRun;
    }

    private static final class ReadErrorCorrectionStats {
        public int numReadsCorrected;
        public int numReadsUncorrected;
        public int numBasesCorrected;
        public int numSolidKmers;
        public int numUncorrectableKmers;
        public int numCorrectedKmers;
    }

    /**
     * Wrapper utility class that holds, for each position in read, a list of bytes representing candidate corrections.
     * So, a read ACAGT where the middle A has found to be errorful might look like:
     * 0: {}
     * 1: {}
     * 2: {'C','C','C'}
     * 3: {}
     * 4: {}
     *
     * It's up to the method getConsensusCorrection()  to decide how to use the correction sets for each position.
     * By default, only strict consensus is allowed right now.
     *
     */
    protected static class CorrectionSet {
        private final int size;
        private ArrayList<List<Byte>> corrections;

        /**
         * Main class constructor.
         * @param size      Size of correction set, needs to be set equal to the read being corrected
         */
        public CorrectionSet(final int size) {
            this.size = size;
            corrections = new ArrayList<>(size);
            for (int k=0; k < size; k++) {
                corrections.add(k, new ArrayList<>());
            }
        }

        /**
         * Add a base to this correction set at a particular offset, measured from the start of the read
         * @param offset                          Offset from start of read
         * @param base                            base to be added to list of corrections at this offset
         */
        public void add(final int offset, final byte base) {
            if (offset >= size || offset < 0) {
                throw new IllegalStateException("Bad entry into CorrectionSet: offset > size");
            }
            if (!BaseUtils.isRegularBase(base)) {
                return; // no irregular base correction
            }

            final List<Byte> storedBytes = corrections.get(offset);
            storedBytes.add(base);
        }

        /**
         * Get list of corrections for a particular offset
         * @param offset                            Offset of interest
         * @return                                  List of bases representing possible corrections at this offset
         */
        public List<Byte> get(final int offset) {
            if (offset >= size || offset < 0) {
                throw new IllegalArgumentException("Illegal call of CorrectionSet.get(): offset must be < size");
            }
            return corrections.get(offset);
        }

        /**
         * Get consensus correction for a particular offset. In this implementation, it just boils down to seeing if
         * byte list associated with offset has identical values. If so, return this base, otherwise return null.
         * @param offset
         * @return                                 Consensus base, or null if no consensus possible.
         */
        public Byte getConsensusCorrection(final int offset) {
            if (offset >= size || offset < 0) {
                throw new IllegalArgumentException("Illegal call of CorrectionSet.getConsensusCorrection(): offset must be < size");
            }
            final List<Byte> storedBytes = corrections.get(offset);
            if (storedBytes.isEmpty()) {
                return null;
            }

            // todo - is there a cheaper/nicer way to compare if all elements in list are identical??
            final byte lastBase = storedBytes.remove(storedBytes.size()-1);
            for (final Byte b: storedBytes) {
                // strict correction rule: all bases must match
                if (b != lastBase) {
                    return null;
                }
            }

            // all bytes then are equal:
            return lastBase;

        }



    }
}
