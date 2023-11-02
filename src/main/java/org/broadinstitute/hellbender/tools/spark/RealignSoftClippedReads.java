package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.MergingIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.SequenceDictionaryValidationArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaReadAligner;
import org.broadinstitute.hellbender.utils.read.*;
import org.codehaus.plexus.util.FileUtils;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.File;
import java.util.*;

/**
 * Realigns soft-clipped reads. Intended for use with short-read Dragen v3.7.8 BAMs/CRAMs.
 */
@CommandLineProgramProperties(
        summary = "Realigns soft-clipped reads to a given reference.",
        oneLineSummary = "Realigns soft-clipped reads to a given reference.",
        programGroup = OtherProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public final class RealignSoftClippedReads extends MultiplePassReadWalker {
    private static final long serialVersionUID = 1L;

    @Argument(doc="Output bam file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public GATKPath output;

    @Argument(doc="BWA index image path.",
            fullName = BwaArgumentCollection.BWA_MEM_INDEX_IMAGE_FULL_NAME)
    public GATKPath indexImage;

    @Argument(doc="Minimum length of soft clips for realignment.",
            fullName = "min-clipped-length",
            optional = true,
            minValue = 0)
    public int minSoftClipLength = 0;

    @Argument(doc="Aligner buffer size.",
            fullName = "buffer-size",
            optional = true,
            minValue = 100)
    public int bufferSize = 10000;

    @Argument(doc="Number of bwa threads.",
            fullName = "bwa-threads",
            optional = true,
            minValue = 1)
    public int bwaThreads = 12;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public SequenceDictionaryValidationArgumentCollection getSequenceDictionaryValidationArgumentCollection(){
        return new SequenceDictionaryValidationArgumentCollection.NoValidationCollection();
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    private Set<String> softclippedReadNames;
    private ValidSoftClipReadFilter softClipReadFilter;
    private GATKRead lastRead = null;
    private List<GATKRead> readBuffer;

    @Override
    public void onTraversalStart() {
        softclippedReadNames = new HashSet<>();
        softClipReadFilter = new ValidSoftClipReadFilter(minSoftClipLength);
        readBuffer = new ArrayList<>(bufferSize);
    }

    @Override
    public void traverseReads() {
        forEachRead((GATKRead read, ReferenceContext reference, FeatureContext features) ->
                firstPass(read)
        );
        logger.info("Found " + softclippedReadNames.size() + " soft-clipped read pairs");
        final SAMFileHeader tempHeader = getHeaderForReads().clone();
        tempHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);

        final File clippedReadBam = FileUtils.createTempFile("realign_1", ".bam", null);
        final File unclippedReadBam = FileUtils.createTempFile("realign_2", ".bam", null);
        try (final SAMFileWriter clippedWriter = new SAMFileWriterFactory().makeWriter(tempHeader, false,
                clippedReadBam, null);
             final SAMFileGATKReadWriter nonClippedWriter = createSAMWriter(new GATKPath(unclippedReadBam.getPath()), true)) {

            logger.info("Splitting clipped and unclipped reads...");
            forEachRead((GATKRead read, ReferenceContext reference, FeatureContext features) ->
                    secondPass(read, clippedWriter, nonClippedWriter)
            );
        }

        final File alignmentBam = FileUtils.createTempFile("realign_3", ".bam", null);
        try (final ReadsPathDataSource clippedSource = new ReadsPathDataSource(clippedReadBam.toPath());
             final SAMFileGATKReadWriter alignmentWriter = createSAMWriter(new GATKPath(alignmentBam.getPath()), false)) {
            final BwaReadAligner aligner = new BwaReadAligner(indexImage.toString(), clippedSource.getHeader(), true, true);
            aligner.setNThreads(bwaThreads);
            logger.info("Aligning reads...");
            clippedSource.forEach(read -> thirdPass(read, aligner, alignmentWriter));
            flushBuffer(aligner, alignmentWriter);
        }

        try (final ReadsPathDataSource unclippedSource = new ReadsPathDataSource(unclippedReadBam.toPath());
             final ReadsPathDataSource alignmentSource = new ReadsPathDataSource(alignmentBam.toPath());
             final SAMFileGATKReadWriter writer = createSAMWriter(output, true)) {

            final CloseableIterator<SAMRecord> unclippedRecordsIter = new ClosableReadSourceIterator(unclippedSource);
            final CloseableIterator<SAMRecord> alignedRecordsIter = new ClosableReadSourceIterator(alignmentSource);
            final Comparator<SAMRecord> weakComparator = new Comparator<>() {

                final SAMRecordComparator strongComparator = getHeaderForReads().getSortOrder().getComparatorInstance();
                @Override
                public int compare(SAMRecord o1, SAMRecord o2) {
                    return strongComparator.fileOrderCompare(o1, o2);
                }
            };
            final MergingIterator<SAMRecord> iter = new MergingIterator<>(weakComparator, List.of(unclippedRecordsIter, alignedRecordsIter));
            iter.stream().map(SAMRecordToGATKReadAdapter::new).forEach(writer::addRead);
        }
    }

    /**
     * Helper class to map a {@link ReadsDataSource} iterator to a {@link CloseableIterator} and maps to
     * {@link SAMRecord} as well.
     */
    private static final class ClosableReadSourceIterator implements CloseableIterator<SAMRecord> {

        final ReadsDataSource source;
        final Iterator<GATKRead> iterator;

        public ClosableReadSourceIterator(final ReadsDataSource source) {
            this.source = source;
            this.iterator = source.iterator();
        }

        @Override
        public void close() {
            source.close();
        }

        @Override
        public boolean hasNext() {
            if (!iterator.hasNext()) {
                close();
                return false;
            }
            return iterator.hasNext();
        }

        @Override
        public SAMRecord next() {
            return iterator.next().convertToSAMRecord(source.getHeader());
        }
    }

    /**
     * Gather names of soft-clipped reads
     */
    private void firstPass(final GATKRead read) {
        if (softClipReadFilter.test(read)) {
            softclippedReadNames.add(read.getName());
        }
    }

    /**
     * Divide reads into temporary clipped and unclipped read files
     */
    private void secondPass(final GATKRead read, final SAMFileWriter clippedWriter, final GATKReadWriter nonClippedWriter) {
        if (softclippedReadNames.contains(read.getName())) {
            if (ReadFilterLibrary.PRIMARY_LINE.test(read)) {
                clippedWriter.addAlignment(read.convertToSAMRecord(clippedWriter.getFileHeader()));
            }
        } else {
            nonClippedWriter.addRead(read);
        }
    }

    /**
     * Realign the clipped reads
     */
    private void thirdPass(final GATKRead read, final BwaReadAligner aligner, final GATKReadWriter writer) {
        if (lastRead == null) {
            lastRead = read;
        } else if (lastRead.getName().equals(read.getName())) {
            if (lastRead.isReverseStrand()) {
                lastRead.reverseComplement();
            }
            if (read.isReverseStrand()) {
                read.reverseComplement();
            }
            readBuffer.add(lastRead);
            readBuffer.add(read);
            if (readBuffer.size() == bufferSize || readBuffer.size() + 1 == bufferSize) {
                flushBuffer(aligner, writer);
            }
            lastRead = null;
        } else {
            lastRead = read;
        }
    }

    private void flushBuffer(final BwaReadAligner aligner, final GATKReadWriter writer) {
        logger.info("  Running batch of " + readBuffer.size() + " reads");
        aligner.apply(readBuffer.iterator()).forEachRemaining(writer::addRead);
        readBuffer.clear();
    }

    static GATKRead unmapRead(final GATKRead read, final SAMFileHeader header) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(read.getName());
        record.setReadBases(read.getBases());
        record.setBaseQualities(read.getBaseQualities());
        record.setAttribute(SAMTag.RG, read.getReadGroup());
        record.setDuplicateReadFlag(read.isDuplicate());
        if (read.isFirstOfPair()) {
            read.setIsFirstOfPair();
        } else if (read.isSecondOfPair()) {
            read.setIsSecondOfPair();
        }
        if (read.isReverseStrand()) {
            record.reverseComplement();
        }
        return SAMRecordToGATKReadAdapter.headerlessReadAdapter(record);
    }

    static GATKRead setRealigned(GATKRead read) {
        final GATKRead copy = read.copy();
        copy.setAttribute("RA", 1);
        return copy;
    }

    private static final class ValidSoftClipReadFilter extends ReadFilter {
        private static final long serialVersionUID = 1L;
        private final int minSoftClipLength;

        public ValidSoftClipReadFilter(final int minSoftClipLength) {
            this.minSoftClipLength = minSoftClipLength;
        }

        @Override
        public boolean test(final GATKRead read) {
            return ReadFilterLibrary.PRIMARY_LINE.test(read)
                    && read.getCigarElements().stream().anyMatch(c -> c.getOperator() == CigarOperator.SOFT_CLIP
                    && c.getLength() >= this.minSoftClipLength);
        }

    }
}
