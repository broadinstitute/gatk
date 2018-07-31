package org.broadinstitute.hellbender.tools.spark.bwa;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.SparkFiles;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.*;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.*;

/**
 * The BwaSparkEngine provides a simple interface for transforming a JavaRDD<GATKRead> in which the reads are paired
 * and unaligned, into a JavaRDD<GATKRead> of aligned reads, and does so lazily.
 * Use it like this:
 *     Make one, call the {@link #align} method for each of your input RDDs in a pipeline that runs some action, close it.
 *
 * The reason that the pipeline must culminate in some action, is because this class implements a lazy
 * transform, and nothing will happen otherwise.
 *
 * See {@link BwaSpark#runTool runTool} for an example.
 */
public final class BwaSparkEngine implements AutoCloseable {
    private static final String REFERENCE_INDEX_IMAGE_FILE_SUFFIX = ".img";
    private final JavaSparkContext ctx;
    private final String indexFileName;
    private final boolean resolveIndexFileName;
    private final Broadcast<SAMFileHeader> broadcastHeader;

    /**
     * @param ctx           the Spark context
     * @param referenceFile the path to the reference file named <i>_prefix_.fa</i>, which is used to find the image file with name <i>_prefix_.fa.img</i>.
     *                      Can be <code>null</code> if the indexFileName is provided.
     * @param indexFileName the index image file name that already exists, or <code>null</code> to have the image file automatically distributed.
     * @param inputHeader   the SAM file header to use for reads
     * @param refDictionary the sequence dictionary to use for reads if the SAM file header doesn't have one (or it's empty)
     */
    public BwaSparkEngine(final JavaSparkContext ctx,
                          final String referenceFile,
                          final String indexFileName,
                          SAMFileHeader inputHeader,
                          final SAMSequenceDictionary refDictionary) {
        Utils.nonNull(referenceFile);
        Utils.nonNull(inputHeader);
        this.ctx = ctx;
        if (indexFileName != null) {
            this.indexFileName = indexFileName;
            this.resolveIndexFileName = false;
        } else {
            String indexFile = referenceFile + REFERENCE_INDEX_IMAGE_FILE_SUFFIX;
            ctx.addFile(indexFile); // distribute index file to all executors
            this.indexFileName = IOUtils.getPath(indexFile).getFileName().toString();
            this.resolveIndexFileName = true;
        }

        if (inputHeader.getSequenceDictionary() == null || inputHeader.getSequenceDictionary().isEmpty()) {
            Utils.nonNull(refDictionary);
            inputHeader = inputHeader.clone();
            inputHeader.setSequenceDictionary(refDictionary);
        }
        broadcastHeader = ctx.broadcast(inputHeader);
    }

    public SAMFileHeader getHeader() { return broadcastHeader.getValue(); }

    /**
     * Performs pair-end alignment on a RDD.
     * @param unalignedReads the read-pairs to align.
     * @return never {@code null}.
     */
    public JavaRDD<GATKRead> alignPaired(final JavaRDD<GATKRead> unalignedReads) {
        return align(unalignedReads, true);
    }

    /**
     * Performs single-end alignment on a RDD.
     *
     * @param unalignedReads the reads to align.
     * @return never {@code null}.
     */
    public JavaRDD<GATKRead> alignUnpaired(final JavaRDD<GATKRead> unalignedReads) {
        return align(unalignedReads, false);
    }

    /**
     * Performs read alignment on a RDD.
     * @param unalignedReads the reads to align.
     * @param pairedAlignment whether it should perform pair-end alignment ({@code true}) or single-end alignment ({@code false}).
     * @return never {@code null}.
     */
    public JavaRDD<GATKRead> align(final JavaRDD<GATKRead> unalignedReads, final boolean pairedAlignment) {
        final Broadcast<SAMFileHeader> broadcastHeader = this.broadcastHeader;
        final String indexFileName = this.indexFileName;
        final boolean resolveIndexFileName = this.resolveIndexFileName;
        return unalignedReads.mapPartitions(itr ->
                new ReadAligner(resolveIndexFileName ? SparkFiles.get(indexFileName) : indexFileName, broadcastHeader.value(), pairedAlignment).apply(itr));
    }

    @Override
    public void close() {
        broadcastHeader.destroy();
        BwaMemIndexCache.closeAllDistributedInstances(ctx);
    }

    private static final class ReadAligner {
        private final BwaMemIndex bwaMemIndex;
        private final SAMFileHeader readsHeader;
        private final boolean alignsPairs;

        // assumes 128Mb partitions, with reads needing about 100bytes each when BAM compressed
        private static final int READS_PER_PARTITION_GUESS = 1500000;

        ReadAligner( final String indexFileName, final SAMFileHeader readsHeader, final boolean alignsPairs) {
            this.bwaMemIndex = BwaMemIndexCache.getInstance(indexFileName);
            this.readsHeader = readsHeader;
            this.alignsPairs = alignsPairs;
            if ( alignsPairs && readsHeader.getSortOrder() != SAMFileHeader.SortOrder.queryname ) {
                throw new UserException("Input must be queryname sorted unless you use single-ended alignment mode.");
            }
        }

        Iterator<GATKRead> apply( final Iterator<GATKRead> readItr ) {
            final List<GATKRead> inputReads = new ArrayList<>(READS_PER_PARTITION_GUESS);
            while ( readItr.hasNext() ) {
                inputReads.add(readItr.next());
            }
            final int nReads = inputReads.size();
            if ( alignsPairs ) {
                if ( (nReads & 1) != 0 ) {
                    throw new GATKException("We're supposed to be aligning paired reads, but there are an odd number of them.");
                }
                for ( int idx = 0; idx != nReads; idx += 2 ) {
                    final String readName1 = inputReads.get(idx).getName();
                    final String readName2 = inputReads.get(idx+1).getName();
                    if ( !Objects.equals(readName1,readName2) ) {
                        throw new GATKException("Read pair has varying template name: "+readName1+" .vs "+readName2);
                    }
                }
            }
            final List<List<BwaMemAlignment>> allAlignments;
            if ( nReads == 0 ) allAlignments = Collections.emptyList();
            else {
                final List<byte[]> seqs = new ArrayList<>(nReads);
                for (final GATKRead read : inputReads) {
                    seqs.add(read.getBases());
                }
                final BwaMemAligner aligner = new BwaMemAligner(bwaMemIndex);
                // we are dealing with interleaved, paired reads.  tell BWA that they're paired.
                if (alignsPairs) {
                    aligner.alignPairs();
                }
                allAlignments = aligner.alignSeqs(seqs);
            }
            final List<String> refNames = bwaMemIndex.getReferenceContigNames();
            final List<GATKRead> outputReads = new ArrayList<>(allAlignments.stream().mapToInt(List::size).sum());
            for ( int idx = 0; idx != nReads; ++idx ) {
                final GATKRead originalRead = inputReads.get(idx);
                final String readName = originalRead.getName();
                final byte[] bases = originalRead.getBases();
                final byte[] quals = originalRead.getBaseQualities();
                final String readGroup = originalRead.getReadGroup();
                final List<BwaMemAlignment> alignments = allAlignments.get(idx);
                final Map<BwaMemAlignment,String> saTagMap = BwaMemAlignmentUtils.createSATags(alignments,refNames);
                for ( final BwaMemAlignment alignment : alignments ) {
                    final SAMRecord samRecord =
                            BwaMemAlignmentUtils.applyAlignment(readName, bases, quals, readGroup,
                                                                alignment, refNames, readsHeader, false, true);
                    final GATKRead rec = SAMRecordToGATKReadAdapter.headerlessReadAdapter(samRecord);
                    final String saTag = saTagMap.get(alignment);
                    if ( saTag != null ) rec.setAttribute("SA", saTag);
                    outputReads.add(rec);
                }
            }
            return outputReads.iterator();
        }
    }
}
