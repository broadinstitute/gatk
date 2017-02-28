package org.broadinstitute.hellbender.tools.spark.bwa;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.seqdoop.hadoop_bam.BAMInputFormat;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * The BwaSparkEngine provides a simple interface for transforming a JavaRDD<GATKRead> in which the reads are paired
 * and unaligned, into a JavaRDD<GATKRead> of aligned reads, and does so lazily.
 * Use it like this:
 *     Make one, call the align method for each of your input RDDs in a pipeline that runs some action, close it.
 *
 * The reason that the pipeline must culminate in some action, is because this class implements a lazy
 * transform, and nothing will happen otherwise.
 *
 * See {@link BwaSpark#runTool runTool} for an example.
 */
public final class BwaSparkEngine implements AutoCloseable {
    private final JavaSparkContext ctx;
    private final String indexFileName;
    private final Broadcast<SAMFileHeader> broadcastHeader;

    public BwaSparkEngine(final JavaSparkContext ctx,
                          final String indexFileName,
                          SAMFileHeader inputHeader,
                          final SAMSequenceDictionary refDictionary) {
        Utils.nonNull(indexFileName);
        Utils.nonNull(inputHeader);
        this.ctx = ctx;
        this.indexFileName = indexFileName;
        if (inputHeader.getSequenceDictionary() == null || inputHeader.getSequenceDictionary().isEmpty()) {
            Utils.nonNull(refDictionary);
            inputHeader = inputHeader.clone();
            inputHeader.setSequenceDictionary(refDictionary);
        }
        broadcastHeader = ctx.broadcast(inputHeader);
    }

    public SAMFileHeader getHeader() { return broadcastHeader.getValue(); }

    public JavaRDD<GATKRead> align(final JavaRDD<GATKRead> unalignedReads) {
        final Broadcast<SAMFileHeader> broadcastHeader = this.broadcastHeader;
        final String indexFileName = this.indexFileName;
        return unalignedReads.mapPartitions(itr -> new ReadAligner(indexFileName, broadcastHeader.value()).apply(itr));
    }

    @Override
    public void close() {
        broadcastHeader.destroy();
        BwaMemIndexSingleton.closeAllDistributedInstances(ctx);
    }

    private static final class ReadAligner {
        private final BwaMemIndex bwaMemIndex;
        private final SAMFileHeader readsHeader;

        // assumes 128Mb partitions, with reads needing about 100bytes each when BAM compressed
        private static final int READS_PER_PARTITION_GUESS = 1500000;

        ReadAligner( final String indexFileName, final SAMFileHeader readsHeader ) {
            this.bwaMemIndex = BwaMemIndexSingleton.getInstance(indexFileName);
            this.readsHeader = readsHeader;
        }

        Iterator<GATKRead> apply( final Iterator<GATKRead> readItr ) {
            final List<GATKRead> inputReads = new ArrayList<>(READS_PER_PARTITION_GUESS);
            while ( readItr.hasNext() ) {
                inputReads.add(readItr.next());
            }
            final int nReads = inputReads.size();
            if ( (nReads & 1) != 0 ) {
                throw new GATKException("We're supposed to be aligning paired reads, but there are an odd number of them.");
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
                aligner.alignPairs();
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
                final long nPrimaries = alignments.stream().filter(aln->(aln.getSamFlag()&SAMFlag.NOT_PRIMARY_ALIGNMENT.intValue())==0).count();
                for ( final BwaMemAlignment alignment : alignments ) {
                    final SAMRecord samRecord =
                            BwaMemAlignmentUtils.applyAlignment(readName, bases, quals, readGroup,
                                                                alignment, refNames, readsHeader, false, true);
                    final GATKRead rec = SAMRecordToGATKReadAdapter.headerlessReadAdapter(samRecord);
                    if ( nPrimaries > 1 ) {
                        final StringBuilder saTag = new StringBuilder();
                        alignments.stream()
                                .filter(aln->aln != alignment)
                                .filter(aln->(aln.getSamFlag()&SAMFlag.NOT_PRIMARY_ALIGNMENT.intValue())==0)
                                .forEach(aln->saTag.append(BwaMemAlignmentUtils.asTag(aln, refNames)));
                        rec.setAttribute("SA", saTag.toString());
                    }
                    outputReads.add(rec);
                }
            }
            return outputReads.iterator();
        }
    }
}
