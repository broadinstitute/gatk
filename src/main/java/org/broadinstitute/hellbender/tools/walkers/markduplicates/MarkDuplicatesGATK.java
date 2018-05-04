package org.broadinstitute.hellbender.tools.walkers.markduplicates;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.SortingLongCollection;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.*;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * A better duplication marking algorithm that handles all cases including clipped
 * and gapped alignments.
 *
 * @author Tim Fennell
 */
@ExperimentalFeature
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "(Experimental) Examines aligned records in the supplied SAM/BAM/CRAM file to locate duplicate molecules. " +
                "All records are then written to the output file with the duplicate records flagged.",
        oneLineSummary = "Examines aligned records in the supplied SAM/BAM/CRAM file to locate duplicate molecules.",
        programGroup = ReadDataManipulationProgramGroup.class
)
public final class MarkDuplicatesGATK extends AbstractMarkDuplicatesCommandLineProgram {
    /**
     * If more than this many sequences in SAM file, don't spill to disk because there will not
     * be enough file handles.
     */

    @Argument(shortName = "MAX_SEQS",
            doc = "This option is obsolete. ReadEnds will always be spilled to disk.")
    public int MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP = 50000;

    @Argument(shortName = "MAX_FILE_HANDLES",
            doc = "Maximum number of file handles to keep open when spilling read ends to disk. " +
                    "Set this number a little lower than the per-process maximum number of file that may be open. " +
                    "This number can be found by executing the 'ulimit -n' command on a Unix system.")
    public int MAX_FILE_HANDLES_FOR_READ_ENDS_MAP = 8000;

    @Argument(doc = "This number, plus the maximum RAM available to the JVM, determine the memory footprint used by " +
            "some of the sorting collections.  If you are running out of memory, try reducing this number.")
    public double SORTING_COLLECTION_SIZE_RATIO = 0.25;

    @Argument(doc = "Report Memory Stats at various times during the run")
    public boolean reportMemoryStats = false;


    private SortingCollection<ReadEndsForMarkDuplicates> pairSort;
    private SortingCollection<ReadEndsForMarkDuplicates> fragSort;
    private SortingLongCollection duplicateIndexes;
    private int numDuplicateIndices = 0;

    private LibraryIdGenerator libraryIdGenerator = null; // this is initialized in buildSortedReadEndLists

    public MarkDuplicatesGATK() {
        DUPLICATE_SCORING_STRATEGY = DuplicateScoringStrategy.ScoringStrategy.SUM_OF_BASE_QUALITIES;
    }

    /**
     * Main work method.  Reads the BAM file once and collects sorted information about
     * the 5' ends of both ends of each read (or just one end in the case of pairs).
     * Then makes a pass through those determining duplicates before re-reading the
     * input file and writing it out with duplication flags set correctly.
     */
    @Override
    protected Object doWork() {
        IOUtil.assertFilesAreReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(METRICS_FILE);

        reportMemoryStats("Start of doWork");
        logger.info("Reading input file and constructing read end information.");
        buildSortedReadEndLists();
        reportMemoryStats("After buildSortedReadEndLists");
        generateDuplicateIndexes();
        reportMemoryStats("After generateDuplicateIndexes");
        logger.info("Marking " + this.numDuplicateIndices + " records as duplicates.");

        if (this.opticalDuplicatesArgumentCollection.READ_NAME_REGEX == null) {
            logger.warn("Skipped optical duplicate cluster discovery; library size estimation may be inaccurate!");
        } else {
            logger.info("Found " + (this.libraryIdGenerator.getNumberOfOpticalDuplicateClusters()) + " optical duplicate clusters.");
        }

        try( final SamHeaderAndIterator headerAndIterator = openInputs()) {
            final SAMFileHeader header = headerAndIterator.header;

            final SAMFileHeader outputHeader = ReadUtils.cloneSAMFileHeader(header);
            outputHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            for (final String comment : COMMENT) outputHeader.addComment(comment);

            // Key: previous PG ID on a SAM Record (or null).  Value: New PG ID to replace it.
            final Map<String, String> chainedPgIds = getChainedPgIds(outputHeader);

            try (final SAMFileWriter out = createSAMWriter(OUTPUT, REFERENCE_SEQUENCE, outputHeader, true)) {

                // Now copy over the file while marking all the necessary indexes as duplicates
                long recordInFileIndex = 0;
                long nextDuplicateIndex = (this.duplicateIndexes.hasNext() ? this.duplicateIndexes.next() : -1);

                final ProgressLogger progress = new ProgressLogger(logger, (int) 1e7, "Written");
                try (final CloseableIterator<SAMRecord> iterator = headerAndIterator.iterator) {
                    while (iterator.hasNext()) {
                        final SAMRecord rec = iterator.next();
                        if (!rec.isSecondaryOrSupplementary()) {
                            final String library = LibraryIdGenerator.getLibraryName(header, rec);
                            DuplicationMetrics metrics = libraryIdGenerator.getMetricsByLibrary(library);
                            if (metrics == null) {
                                metrics = new DuplicationMetrics();
                                metrics.LIBRARY = library;
                                libraryIdGenerator.addMetricsByLibrary(library, metrics);
                            }

                            // First bring the simple metrics up to date
                            if (rec.getReadUnmappedFlag()) {
                                ++metrics.UNMAPPED_READS;
                            } else if (!rec.getReadPairedFlag() || rec.getMateUnmappedFlag()) {
                                ++metrics.UNPAIRED_READS_EXAMINED;
                            } else {
                                ++metrics.READ_PAIRS_EXAMINED; // will need to be divided by 2 at the end
                            }


                            if (recordInFileIndex == nextDuplicateIndex) {
                                rec.setDuplicateReadFlag(true);

                                // Update the duplication metrics
                                if (!rec.getReadPairedFlag() || rec.getMateUnmappedFlag()) {
                                    ++metrics.UNPAIRED_READ_DUPLICATES;
                                } else {
                                    ++metrics.READ_PAIR_DUPLICATES;// will need to be divided by 2 at the end
                                }

                                // Now try and figure out the next duplicate index
                                if (this.duplicateIndexes.hasNext()) {
                                    nextDuplicateIndex = this.duplicateIndexes.next();
                                } else {
                                    // Only happens once we've marked all the duplicates
                                    nextDuplicateIndex = -1;
                                }
                            } else {
                                rec.setDuplicateReadFlag(false);
                            }
                        }
                        recordInFileIndex++;

                        if (!this.REMOVE_DUPLICATES || !rec.getDuplicateReadFlag()) {
                            if (PROGRAM_RECORD_ID != null) {
                                rec.setAttribute(SAMTag.PG.name(), chainedPgIds.get(rec.getStringAttribute(SAMTag.PG.name())));
                            }
                            out.addAlignment(rec);
                            progress.record(rec);
                        }
                    }
                }
                this.duplicateIndexes.cleanup();

                reportMemoryStats("Before output close");
            }
            reportMemoryStats("After output close");
        }

        // Write out the metrics
        finalizeAndWriteMetrics(libraryIdGenerator);

        return null;
    }

    @VisibleForTesting
    long numOpticalDuplicates() { return ((long) this.libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap().getSumOfValues()); } // cast as long due to returning a double

    /** Print out some quick JVM memory stats. */
    private void reportMemoryStats(final String stage) {
        if(reportMemoryStats) {
            System.gc();
            final Runtime runtime = Runtime.getRuntime();
            logger.info(stage + " freeMemory: " + runtime.freeMemory() + "; totalMemory: " + runtime.totalMemory() +
                    "; maxMemory: " + runtime.maxMemory());
        }
    }

    /**
     * Goes through all the records in a file and generates a set of ReadEndsForMarkDuplicates objects that
     * hold the necessary information (reference sequence, 5' read coordinate) to do
     * duplication, caching to disk as necessary to sort them.
     */
    private void buildSortedReadEndLists() {
        final int maxInMemory = (int) ((Runtime.getRuntime().maxMemory() * SORTING_COLLECTION_SIZE_RATIO) / ReadEndsForMarkDuplicates.SIZE_OF);
        logger.info("Will retain up to " + maxInMemory + " data points before spilling to disk.");

        this.pairSort = SortingCollection.newInstanceFromPaths(ReadEndsForMarkDuplicates.class,
                new ReadEndsForMarkDuplicatesCodec(),
                new ReadEndsMDComparator(),
                maxInMemory,
                TMP_DIR.stream().map(File::toPath).collect(Collectors.toList()));

        this.fragSort = SortingCollection.newInstanceFromPaths(ReadEndsForMarkDuplicates.class,
                new ReadEndsForMarkDuplicatesCodec(),
                new ReadEndsMDComparator(),
                maxInMemory,
                TMP_DIR.stream().map(File::toPath).collect(Collectors.toList()));

        try(final SamHeaderAndIterator headerAndIterator = openInputs()) {
            final SAMFileHeader header = headerAndIterator.header;
            final ReadEndsForMarkDuplicatesMap tmp = new DiskBasedReadEndsForMarkDuplicatesMap(MAX_FILE_HANDLES_FOR_READ_ENDS_MAP);
            long index = 0;
            final ProgressLogger progress = new ProgressLogger(logger, (int) 1e6, "Read");
            final CloseableIterator<SAMRecord> iterator = headerAndIterator.iterator;

            if (null == this.libraryIdGenerator) {
                this.libraryIdGenerator = new LibraryIdGenerator(header);
            }

            while (iterator.hasNext()) {
                final SAMRecord rec = iterator.next();

                // This doesn't have anything to do with building sorted ReadEnd lists, but it can be done in the same pass
                // over the input
                if (PROGRAM_RECORD_ID != null) {
                    // Gather all PG IDs seen in merged input files in first pass.  These are gathered for two reasons:
                    // - to know how many different PG records to create to represent this program invocation.
                    // - to know what PG IDs are already used to avoid collisions when creating new ones.
                    // Note that if there are one or more records that do not have a PG tag, then a null value
                    // will be stored in this set.
                    pgIdsSeen.add(rec.getStringAttribute(SAMTag.PG.name()));
                }

                if (rec.getReadUnmappedFlag()) {
                    if (rec.getReferenceIndex() == -1) {
                        // When we hit the unmapped reads with no coordinate, no reason to continue.
                        break;
                    }
                    // If this read is unmapped but sorted with the mapped reads, just skip it.
                } else if (!rec.isSecondaryOrSupplementary()) {
                    final ReadEndsForMarkDuplicates fragmentEnd = buildReadEnds(header, index, rec);
                    this.fragSort.add(fragmentEnd);

                    if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
                        final String key = rec.getAttribute(ReservedTagConstants.READ_GROUP_ID) + ":" + rec.getReadName();
                        ReadEndsForMarkDuplicates pairedEnds = tmp.remove(rec.getReferenceIndex(), key);

                        // See if we've already seen the first end or not
                        if (pairedEnds == null) {
                            pairedEnds = buildReadEnds(header, index, rec);
                            tmp.put(pairedEnds.read2ReferenceIndex, key, pairedEnds);
                        } else {
                            final int sequence = fragmentEnd.read1ReferenceIndex;
                            final int coordinate = fragmentEnd.read1Coordinate;

                            // Set orientationForOpticalDuplicates, which always goes by the first then the second end for the strands.  NB: must do this
                            // before updating the orientation later.
                            if (rec.getFirstOfPairFlag()) {
                                pairedEnds.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(rec.getReadNegativeStrandFlag(), pairedEnds.orientation == ReadEnds.R);
                            } else {
                                pairedEnds.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(pairedEnds.orientation == ReadEnds.R, rec.getReadNegativeStrandFlag());
                            }

                            // If the second read is actually later, just add the second read data, else flip the reads
                            if (sequence > pairedEnds.read1ReferenceIndex ||
                                    (sequence == pairedEnds.read1ReferenceIndex && coordinate >= pairedEnds.read1Coordinate)) {
                                pairedEnds.read2ReferenceIndex = sequence;
                                pairedEnds.read2Coordinate = coordinate;
                                pairedEnds.read2IndexInFile = index;
                                pairedEnds.orientation = ReadEnds.getOrientationByte(pairedEnds.orientation == ReadEnds.R,
                                        rec.getReadNegativeStrandFlag());
                            } else {
                                pairedEnds.read2ReferenceIndex = pairedEnds.read1ReferenceIndex;
                                pairedEnds.read2Coordinate = pairedEnds.read1Coordinate;
                                pairedEnds.read2IndexInFile = pairedEnds.read1IndexInFile;
                                pairedEnds.read1ReferenceIndex = sequence;
                                pairedEnds.read1Coordinate = coordinate;
                                pairedEnds.read1IndexInFile = index;
                                pairedEnds.orientation = ReadEnds.getOrientationByte(rec.getReadNegativeStrandFlag(),
                                        pairedEnds.orientation == ReadEnds.R);
                            }

                            pairedEnds.score += DuplicateScoringStrategy.computeDuplicateScore(rec, this.DUPLICATE_SCORING_STRATEGY);
                            this.pairSort.add(pairedEnds);
                        }
                    }
                }

                // Print out some stats every 1m reads
                ++index;
                if (progress.record(rec)) {
                    logger.info("Tracking " + tmp.size() + " as yet unmatched pairs. " + tmp.sizeInRam() + " records in RAM.");
                }
            }

            logger.info("Read " + index + " records. " + tmp.size() + " pairs never matched.");
            iterator.close();
        }

        // Tell these collections to free up memory if possible.
        this.pairSort.doneAdding();
        this.fragSort.doneAdding();
    }

    /** Builds a read ends object that represents a single read. */
    private ReadEndsForMarkDuplicates buildReadEnds(final SAMFileHeader header, final long index, final SAMRecord rec) {
        final ReadEndsForMarkDuplicates ends = new ReadEndsForMarkDuplicates();
        ends.read1ReferenceIndex = rec.getReferenceIndex();
        ends.read1Coordinate = rec.getReadNegativeStrandFlag() ? rec.getUnclippedEnd() : rec.getUnclippedStart();
        ends.orientation = rec.getReadNegativeStrandFlag() ? ReadEnds.R : ReadEnds.F;
        ends.read1IndexInFile = index;
        ends.score = DuplicateScoringStrategy.computeDuplicateScore(rec, this.DUPLICATE_SCORING_STRATEGY);

        // Doing this lets the ends object know that it's part of a pair
        if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
            ends.read2ReferenceIndex = rec.getMateReferenceIndex();
        }

        // Fill in the library ID
        ends.libraryId = libraryIdGenerator.getLibraryId(rec);

        // Fill in the location information for optical duplicates
        if (this.opticalDuplicateFinder.addLocationInformation(rec.getReadName(), ends)) {
            // calculate the RG number (nth in list)
            ends.readGroup = 0;
            final String rg = (String) rec.getAttribute("RG");
            final List<SAMReadGroupRecord> readGroups = header.getReadGroups();

            if (rg != null && readGroups != null) {
                for (final SAMReadGroupRecord readGroup : readGroups) {
                    if (readGroup.getReadGroupId().equals(rg)) break;
                    else ends.readGroup++;
                }
            }
        }

        return ends;
    }

    /**
     * Goes through the accumulated ReadEndsForMarkDuplicates objects and determines which of them are
     * to be marked as duplicates.
     *
     * @return an array with an ordered list of indexes into the source file
     */
    private void generateDuplicateIndexes() {
        // Keep this number from getting too large even if there is a huge heap.
        final int maxInMemory = (int) Math.min((Runtime.getRuntime().maxMemory() * 0.25) / SortingLongCollection.SIZEOF,
                (double) (Integer.MAX_VALUE - 5));
        logger.info("Will retain up to " + maxInMemory + " duplicate indices before spilling to disk.");
        this.duplicateIndexes = new SortingLongCollection(maxInMemory, TMP_DIR.toArray(new File[TMP_DIR.size()]));

        ReadEndsForMarkDuplicates firstOfNextChunk = null;
        final List<ReadEndsForMarkDuplicates> nextChunk = new ArrayList<>(200);

        // First just do the pairs
        logger.info("Traversing read pair information and detecting duplicates.");
        for (final ReadEndsForMarkDuplicates next : this.pairSort) {
            if (firstOfNextChunk == null) {
                firstOfNextChunk = next;
                nextChunk.add(firstOfNextChunk);
            } else if (areComparableForDuplicates(firstOfNextChunk, next, true)) {
                nextChunk.add(next);
            } else {
                if (nextChunk.size() > 1) {
                    markDuplicatePairs(nextChunk);
                }

                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
            }
        }
        if (nextChunk.size() > 1) markDuplicatePairs(nextChunk);
        this.pairSort.cleanup();
        this.pairSort = null;

        // Now deal with the fragments
        logger.info("Traversing fragment information and detecting duplicates.");
        boolean containsPairs = false;
        boolean containsFrags = false;

        for (final ReadEndsForMarkDuplicates next : this.fragSort) {
            if (firstOfNextChunk != null && areComparableForDuplicates(firstOfNextChunk, next, false)) {
                nextChunk.add(next);
                containsPairs = containsPairs || next.isPaired();
                containsFrags = containsFrags || !next.isPaired();
            } else {
                if (nextChunk.size() > 1 && containsFrags) {
                    markDuplicateFragments(nextChunk, containsPairs);
                }

                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
                containsPairs = next.isPaired();
                containsFrags = !next.isPaired();
            }
        }
        markDuplicateFragments(nextChunk, containsPairs);
        this.fragSort.cleanup();
        this.fragSort = null;

        logger.info("Sorting list of duplicate records.");
        this.duplicateIndexes.doneAddingStartIteration();
    }

    private static boolean areComparableForDuplicates(final ReadEndsForMarkDuplicates lhs, final ReadEndsForMarkDuplicates rhs, final boolean compareRead2) {
        boolean retval = lhs.libraryId == rhs.libraryId &&
                lhs.read1ReferenceIndex == rhs.read1ReferenceIndex &&
                lhs.read1Coordinate == rhs.read1Coordinate &&
                lhs.orientation == rhs.orientation;

        if (retval && compareRead2) {
            retval = lhs.read2ReferenceIndex == rhs.read2ReferenceIndex &&
                    lhs.read2Coordinate == rhs.read2Coordinate;
        }

        return retval;
    }

    private void addIndexAsDuplicate(final long bamIndex) {
        this.duplicateIndexes.add(bamIndex);
        ++this.numDuplicateIndices;
    }

    /**
     * Takes a list of ReadEndsForMarkDuplicates objects and removes from it all objects that should
     * not be marked as duplicates.  This assumes that the list contains objects representing pairs.
     *
     * @param list
     */
    private void markDuplicatePairs(final List<ReadEndsForMarkDuplicates> list) {
        short maxScore = 0;
        ReadEndsForMarkDuplicates best = null;

        /** All read ends should have orientation FF, FR, RF, or RR **/
        for (final ReadEndsForMarkDuplicates end : list) {
            if (end.score > maxScore || best == null) {
                maxScore = end.score;
                best = end;
            }
        }

        for (final ReadEndsForMarkDuplicates end : list) {
            if (end != best) {
                addIndexAsDuplicate(end.read1IndexInFile);
                addIndexAsDuplicate(end.read2IndexInFile);
            }
        }

        if (this.opticalDuplicatesArgumentCollection.READ_NAME_REGEX != null) {
            AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates(list, best, opticalDuplicateFinder, libraryIdGenerator);
        }
    }

    /**
     * Takes a list of ReadEndsForMarkDuplicates objects and removes from it all objects that should
     * not be marked as duplicates.  This will set the duplicate index for only list items are fragments.
     *
     * @param list
     * @param containsPairs true if the list also contains objects containing pairs, false otherwise.
     */
    private void markDuplicateFragments(final List<ReadEndsForMarkDuplicates> list, final boolean containsPairs) {
        if (containsPairs) {
            for (final ReadEndsForMarkDuplicates end : list) {
                if (!end.isPaired()) addIndexAsDuplicate(end.read1IndexInFile);
            }
        } else {
            short maxScore = 0;
            ReadEndsForMarkDuplicates best = null;
            for (final ReadEndsForMarkDuplicates end : list) {
                if (end.score > maxScore || best == null) {
                    maxScore = end.score;
                    best = end;
                }
            }

            for (final ReadEndsForMarkDuplicates end : list) {
                if (end != best) {
                    addIndexAsDuplicate(end.read1IndexInFile);
                }
            }
        }
    }

    /** Comparator for ReadEndsForMarkDuplicates that orders by read1 position then pair orientation then read2 position. */
    static class ReadEndsMDComparator implements Comparator<ReadEndsForMarkDuplicates>, Serializable {
        private static final long serialVersionUID = 1L;
        @Override
        public int compare(final ReadEndsForMarkDuplicates lhs, final ReadEndsForMarkDuplicates rhs) {
            int retval = lhs.libraryId - rhs.libraryId;
            if (retval == 0) retval = lhs.read1ReferenceIndex - rhs.read1ReferenceIndex;
            if (retval == 0) retval = lhs.read1Coordinate - rhs.read1Coordinate;
            if (retval == 0) retval = lhs.orientation - rhs.orientation;
            if (retval == 0) retval = lhs.read2ReferenceIndex - rhs.read2ReferenceIndex;
            if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
            if (retval == 0) retval = (int) (lhs.read1IndexInFile - rhs.read1IndexInFile);
            if (retval == 0) retval = (int) (lhs.read2IndexInFile - rhs.read2IndexInFile);
            return retval;
        }
    }
}
