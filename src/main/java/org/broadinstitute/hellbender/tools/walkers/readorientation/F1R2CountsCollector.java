package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Collects the F1R2 counts of each sample from which {@link LearnReadOrientationModel} learns, into a tar.gz of, for each sample:
 * <ul>
 *     <li>a ref histogram table and an alt histogram table of the sites without alt reads and with one alt read, by block
 *     of the genome, reference context and depth (see {@link F1R2HistogramTable})</li>
 *     <li>an alt table of the sites with more alt reads of each reference context (see {@link AltSiteRecord})</li>
 * </ul>
 */
public class F1R2CountsCollector {

    public static final String ALT_TABLE_EXTENSION = ".alt_table";
    public static final String ALT_HIST_EXTENSION = ".alt_histogram";
    public static final String REF_HIST_EXTENSION = ".ref_histogram";

    private final Set<String> samples;

    private final SAMFileHeader header;

    private final CollectF1R2CountsArgumentCollection CF1R2Args;

    // the statistics of the current sample's pileup at the current locus, reused for every sample and locus
    private final PileupStatistics pileupStatistics;

    // the reference bases around the current locus
    private final ReferenceContextWindow referenceContextWindow = new ReferenceContextWindow(ReferenceContextWindow.DEFAULT_WINDOW_SIZE);

    // For each sample, the counts of sites without alt reads and with one alt read of the current block of the genome
    private final Map<String, BlockSiteCounts> blockSiteCounts;

    // For each sample, the writers of its ref and alt histogram tables, to which the counts of each block are written
    // when the next block begins, since loci come in order
    private final Map<String, TableWriter<F1R2HistogramTable.Row>> refHistogramWriters;
    private final Map<String, TableWriter<F1R2HistogramTable.Row>> altHistogramWriters;

    // the contig of the current block and its index along the contig, or null and -1 before the first locus
    private String blockContig = null;
    private int blockIndex = -1;

    // the blocks that have begun, to make sure that loci come in order
    private final Set<String> blocks = new HashSet<>();

    // alt table writer for each sample
    // This is by far the largest quantity of data.  In a genome of 3 billion bases where every hundredth site or so has
    // more than 1 alt base the alt records could come out to a gigabyte or more.  Thus it makes sense to split the alt
    // tables by sample *and* by 3-base context so that downstream tools (LearnReadOrientationModel) need only load the
    // records for one context at a time.
    // Each sample's table is split after collection, one sample at a time, so that one file per sample is open during
    // collection rather than one per sample and context.  The unsplit tables have a column of contexts, which the split
    // tables have in their names instead.
    private Map<String, AltSiteRecord.AltSiteRecordTableWriter> altTableWriters;

    // the alt table of each sample, before it is split by 3-base context
    private final File unsplitAltTableDir = IOUtils.createTempDir("unsplit");

    // We put all pre-tarred output files in a temporary directory.  These include, for each sample, a ref and an alt
    // histogram table, each with information for all 3-base contexts, and an alt table for each 3-base context.
    private final File tmpDir = IOUtils.createTempDir("untarred");

    // the final output file is produced by archiving every file in the temporary directory
    private final File outputTarGzFile;

    public F1R2CountsCollector(final CollectF1R2CountsArgumentCollection CF1R2Args, final SAMFileHeader header, final File outputTarGzFile, final Collection<String> samples) {
        this.CF1R2Args = CF1R2Args;
        this.samples = samples.size() == 1 ? Collections.singleton(samples.iterator().next()) : new HashSet<>(samples);
        this.header = header;
        this.outputTarGzFile = outputTarGzFile;
        pileupStatistics = new PileupStatistics(CF1R2Args.minMedianMapQual);

        blockSiteCounts = new HashMap<>(samples.size());
        refHistogramWriters = new HashMap<>(samples.size());
        altHistogramWriters = new HashMap<>(samples.size());
        altTableWriters = new HashMap<>(samples.size());

        for (final String sample : samples) {
            final BlockSiteCounts counts = new BlockSiteCounts(CF1R2Args.maxDepth);
            blockSiteCounts.put(sample, counts);

            // Intentionally not use try-with-resources so that the writers stay open outside of the try block
            final File altTableFile = getUnsplitAltTable(sample);
            try {
                refHistogramWriters.put(sample, F1R2HistogramTable.writer(new File(tmpDir, IOUtils.urlEncode(sample) + REF_HIST_EXTENSION).toPath(),
                        sample, false, counts.getBinDepths()));
                altHistogramWriters.put(sample, F1R2HistogramTable.writer(new File(tmpDir, IOUtils.urlEncode(sample) + ALT_HIST_EXTENSION).toPath(),
                        sample, true, counts.getBinDepths()));
                altTableWriters.put(sample, new AltSiteRecord.AltSiteRecordTableWriter(IOUtils.fileToPath(altTableFile), sample, true));
            } catch (IOException e) {
                throw new UserException(String.format("Encountered an IO exception creating the writers of the F1R2 counts of %s", sample), e);
            }
        }
    }

    public void process(final ReadPileup pileup, final ReferenceContext referenceContext) {
        final String contig = referenceContext.getInterval().getContig();
        final int position = referenceContext.getInterval().getStart();
        updateBlock(contig, position);

        // optimize the common case of a single tumor sample
        final String onlySample = samples.size() == 1 ? samples.iterator().next() : null;
        final Map<String, ReadPileup> splitPileup = samples.size() == 1 ?
                Collections.singletonMap(onlySample, pileup.getPileupForSample(onlySample, header)) : pileup.splitBySample(header, null);

        for (final Map.Entry<String, ReadPileup> entry : splitPileup.entrySet()) {
            // one pass over the elements whose base qualities are above the minimum
            pileupStatistics.reset();
            for (final PileupElement pe : entry.getValue()) {
                if (pe.getQual() > CF1R2Args.minBaseQuality) {
                    pileupStatistics.add(pe);
                }
            }
            process(entry.getKey(), pileupStatistics, contig, position, referenceContext);
        }
    }

    /**
     * Collect the F1R2 counts of a sample at a locus from the statistics of its pileup there, which the caller accumulated
     * over the elements whose base qualities exceed {@link #getMinBaseQuality}, as Mutect2 does one read at a time.  Loci
     * must come in order.
     *
     * @param referenceContext a reference context on the contig, whose data source supplies the reference bases
     */
    public void process(final String sample, final PileupStatistics statistics, final String contig, final int position,
                        final ReferenceContext referenceContext) {
        updateBlock(contig, position);
        // null at the ends of the contig and at Ns
        final String refContext = referenceContextWindow.contextAt(referenceContext, contig, position);
        if (refContext == null || !statistics.isGood()) {
            return;
        }
        final Nucleotide refBase = Nucleotide.decode((byte) refContext.charAt(F1R2FilterConstants.MIDDLE_INDEX));

        final int[] baseCounts = statistics.baseCounts;
        final int depth = statistics.depth;

        // Make a copy of base counts and update the counts of ref to -1. Now the maxElementIndex of the array gives us
        // the alt base.
        final int[] baseCountsCopy = Arrays.copyOf(baseCounts, baseCounts.length);
        baseCountsCopy[refBase.ordinal()] = -1;
        final int altBaseIndex = MathUtils.maxElementIndex(baseCountsCopy);
        final boolean referenceSite = baseCounts[altBaseIndex] == 0;

        // If the site is ref, we simply update the coverage histogram
        if (referenceSite) {
            blockSiteCounts.get(sample).addRefSite(refContext, depth);
            return;
        }

        // If we got here, we have an alt site with a single alt base
        final Nucleotide altBase = Nucleotide.decode(BaseUtils.baseIndexToSimpleBase(altBaseIndex));

        final int refCount = baseCounts[refBase.ordinal()];
        final int altCount = baseCounts[altBaseIndex];
        Utils.validate(altCount > 0, "We must have a nonzero alt read but got " + altCount);
        final int refF1R2 = statistics.f1r2BaseCounts[refBase.ordinal()];
        final int altF1R2 = statistics.f1r2BaseCounts[altBase.ordinal()];

        if (altCount == 1) {
            final ReadOrientation type = altF1R2 == 1 ? ReadOrientation.F1R2 : ReadOrientation.F2R1;
            blockSiteCounts.get(sample).addOneAltReadSite(refContext, altBase, type, depth);
            return;
        }

        try {
            altTableWriters.get(sample).writeRecord(new AltSiteRecord(refContext, refCount, altCount, refF1R2, altF1R2, altBase, contig, position));
        } catch (IOException e) {
            throw new UserException("Encountered an IO Exception writing to the alt data table", e);
        }
    }

    /**
     * @return the minimum base quality, which the base qualities of the pileup elements whose statistics we use must exceed
     */
    public int getMinBaseQuality() {
        return CF1R2Args.minBaseQuality;
    }

    /**
     * @return empty statistics of a pileup, to which the caller adds elements
     */
    public PileupStatistics makePileupStatistics() {
        return new PileupStatistics(CF1R2Args.minMedianMapQual);
    }

    /**
     * If a locus is in a different block of the genome than the previous one, write the histograms of the previous block
     */
    private void updateBlock(final String contig, final int position) {
        final int index = (position - 1) / F1R2FilterConstants.BLOCK_SIZE;
        if (index == blockIndex && contig.equals(blockContig)) {
            return;
        }

        writeHistograms();
        blockContig = contig;
        blockIndex = index;
        Utils.validate(blocks.add(contig + ":" + index), () -> String.format("F1R2 counts must be collected in order of " +
                "position, but the locus %s:%d is in a block of the genome that came before other blocks", contig, position));
    }

    /**
     * Write each sample's counts of the current block of the genome to its histogram tables, and reset them for the next
     * block.  Before the first locus there is nothing to write.
     */
    public void writeHistograms() {
        if (blockContig == null) {
            return;
        }
        for (final String sample : samples) {
            try {
                blockSiteCounts.get(sample).writeRowsAndReset(refHistogramWriters.get(sample), altHistogramWriters.get(sample),
                        blockContig, blockIndex * F1R2FilterConstants.BLOCK_SIZE + 1);
            } catch (IOException e) {
                throw new UserException(String.format("Encountered an IO exception writing the histograms of %s", sample), e);
            }
        }
    }

    public void closeAndArchiveFiles() {
        for (final Map<String, ? extends TableWriter<?>> writers : List.of(refHistogramWriters, altHistogramWriters, altTableWriters)) {
            for (final TableWriter<?> writer : writers.values()) {
                try {
                    writer.close();
                } catch (IOException e) {
                    throw new UserException("Encountered an IO exception while closing the F1R2 count tables", e);
                }
            }
        }

        for (final String sample : samples) {
            try {
                splitAltTable(sample);
            } catch (IOException e) {
                throw new UserException(String.format("Encountered an IO exception while splitting the alt table of %s by reference context", sample), e);
            }
        }

        try {
            IOUtils.writeTarGz(outputTarGzFile.getAbsolutePath(), tmpDir.listFiles());
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile("Could not create tar.gz file", ex);
        }
    }

    /**
     * Split a sample's alt table into a table for each 3-base context, keeping the order of the records.  Only this
     * sample's tables are open at once.
     */
    private void splitAltTable(final String sample) throws IOException {
        final Map<String, AltSiteRecord.AltSiteRecordTableWriter> contextWriters = new HashMap<>(F1R2FilterConstants.ALL_KMERS.size());
        for (final String context : F1R2FilterConstants.ALL_KMERS) {
            final File altTableFile = new File(tmpDir, getAltTableName(sample, context));
            contextWriters.put(context, new AltSiteRecord.AltSiteRecordTableWriter(IOUtils.fileToPath(altTableFile), sample));
        }

        final File unsplitAltTable = getUnsplitAltTable(sample);
        try (AltSiteRecord.AltSiteRecordTableReader reader = new AltSiteRecord.AltSiteRecordTableReader(IOUtils.fileToPath(unsplitAltTable))) {
            for (final AltSiteRecord record : reader) {
                contextWriters.get(record.getReferenceContext()).writeRecord(record);
            }
        }

        for (final AltSiteRecord.AltSiteRecordTableWriter writer : contextWriters.values()) {
            writer.close();
        }
        Files.delete(unsplitAltTable.toPath());
    }

    private File getUnsplitAltTable(final String sample) {
        return new File(unsplitAltTableDir, IOUtils.urlEncode(sample) + ALT_TABLE_EXTENSION);
    }

    /**
     * The statistics of a sample's pileup at a locus, accumulated in one pass over the elements whose base qualities are
     * above the minimum: the counts of their A, C, G and T bases, overall and in F1R2 reads, and what decides whether the
     * pileup is good, ie whether it has such bases, few elements next to indels, and a high median mapping quality.
     */
    public static final class PileupStatistics {
        private final int minMedianMappingQuality;

        // the counts of A, C, G and T, indexed as by BaseUtils.simpleBaseToBaseIndex, and their sum
        final int[] baseCounts = new int[4];
        int depth;

        // the counts of A, C, G and T in F1R2 reads, indexed by Nucleotide.decode(base).ordinal()
        final int[] f1r2BaseCounts = new int[4];

        private int numElements;
        private int numIndelElements;

        // the number of mapping qualities at least the minimum median, the largest below it and the smallest at least it
        private int numAtLeastMinimum;
        private int maxBelowMinimum;
        private int minAtLeastMinimum;

        PileupStatistics(final int minMedianMappingQuality) {
            this.minMedianMappingQuality = minMedianMappingQuality;
            reset();
        }

        public void reset() {
            Arrays.fill(baseCounts, 0);
            Arrays.fill(f1r2BaseCounts, 0);
            depth = 0;
            numElements = 0;
            numIndelElements = 0;
            numAtLeastMinimum = 0;
            maxBelowMinimum = Integer.MIN_VALUE;
            minAtLeastMinimum = Integer.MAX_VALUE;
        }

        void add(final PileupElement pe) {
            add(pe.isDeletion(), pe.getBase(), ReadUtils.isF1R2(pe.getRead()),
                    pe.isDeletion() || pe.isAfterInsertion() || pe.isBeforeDeletionStart(), pe.getMappingQual());
        }

        /**
         * Add a pileup element, whose base quality must exceed the minimum
         *
         * @param isNextToIndel whether the element is a deletion, the first after an insertion or the last before a deletion
         */
        public void add(final boolean isDeletion, final byte base, final boolean isF1R2, final boolean isNextToIndel,
                        final int mappingQuality) {
            numElements++;
            if (!isDeletion) {
                final int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);
                if (baseIndex != -1) {
                    baseCounts[baseIndex]++;
                    depth++;
                }
                final int nucleotide = Nucleotide.decode(base).ordinal();
                if (isF1R2 && nucleotide < f1r2BaseCounts.length) {
                    f1r2BaseCounts[nucleotide]++;
                }
            }
            if (isNextToIndel) {
                numIndelElements++;
            }
            if (mappingQuality >= minMedianMappingQuality) {
                numAtLeastMinimum++;
                minAtLeastMinimum = Math.min(minAtLeastMinimum, mappingQuality);
            } else {
                maxBelowMinimum = Math.max(maxBelowMinimum, mappingQuality);
            }
        }

        /**
         * Use a series of heuristics to detect a bad pileup.  Without A, C, G or T bases, eg if all the reads have deleted
         * bases here, or if more than 1% as many elements as those bases are next to indels, consider the site an indel.
         */
        boolean isGood() {
            return depth > 0 && numIndelElements <= depth / 100 &&
                    medianIsAtLeast(numElements, numAtLeastMinimum, maxBelowMinimum, minAtLeastMinimum, minMedianMappingQuality);
        }

        @Override
        public String toString() {
            return String.format("bases %s, F1R2 bases %s, %d elements, %d next to indels, MQ: %d >= minimum, max below %d, min at least %d",
                    Arrays.toString(baseCounts), Arrays.toString(f1r2BaseCounts), numElements, numIndelElements,
                    numAtLeastMinimum, maxBelowMinimum, minAtLeastMinimum);
        }
    }

    /**
     * Whether the median of some integers, as {@link MathUtils#median} computes it, is at least a threshold, given the
     * number of integers, the number at least the threshold, the largest below it and the smallest at least it.  The median
     * of an odd number of values is the middle one, and of an even number the mean of the middle two.
     */
    @VisibleForTesting
    static boolean medianIsAtLeast(final int numValues, final int numAtLeast, final int maxBelow, final int minAtLeast,
                                   final int threshold) {
        Utils.validateArg(numValues > 0, "the median of no values is undefined");
        if (2 * numAtLeast != numValues) {
            // the middle value, or both middle values, are on the side of the threshold of the majority of values
            return 2 * numAtLeast > numValues;
        }
        // an even number of values whose middle two are the largest below the threshold and the smallest at least it
        return maxBelow + minAtLeast >= 2 * threshold;
    }

    public static List<File> getRefHistogramsFromExtractedTar(final File extractedTarDir) {
        return Arrays.stream(extractedTarDir.listFiles()).filter(file -> file.getAbsolutePath().endsWith(REF_HIST_EXTENSION)).collect(Collectors.toList());
    }

    public static List<File> getAltHistogramsFromExtractedTar(final File extractedTarDir) {
        return Arrays.stream(extractedTarDir.listFiles()).filter(file -> file.getAbsolutePath().endsWith(ALT_HIST_EXTENSION)).collect(Collectors.toList());
    }

    public static List<File> getAltTablesFromExtractedTar(final File extractedTarDir) {
        return Arrays.stream(extractedTarDir.listFiles()).filter(file -> file.getAbsolutePath().endsWith(ALT_TABLE_EXTENSION)).collect(Collectors.toList());
    }

    // alt tables are named <sample>.<reference context>.alt_table
    public static String getAltTableName(final String sample, final String refContext) {
        return IOUtils.urlEncode(sample) + "." + refContext + ALT_TABLE_EXTENSION;
    }

    /**
     * The reference context of an alt table, which is in its name, <sample>.<reference context>.alt_table
     */
    public static String getAltTableContext(final File altTable) {
        final String name = altTable.getName();
        final String nameWithoutExtension = name.substring(0, Math.max(name.length() - ALT_TABLE_EXTENSION.length(), 0));
        final String context = nameWithoutExtension.substring(nameWithoutExtension.lastIndexOf('.') + 1);
        if (!name.endsWith(ALT_TABLE_EXTENSION) || !F1R2FilterConstants.ALL_KMERS.contains(context)) {
            throw new UserException.BadInput(String.format("%s is not named <sample>.<reference context>%s", altTable, ALT_TABLE_EXTENSION));
        }
        return context;
    }
}
