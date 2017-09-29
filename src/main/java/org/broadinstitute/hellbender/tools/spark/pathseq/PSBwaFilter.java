package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Aligns using BWA and filters out reads above the minimum coverage and identity.
 * Reads are NOT otherwise modified, i.e. the alignments are not copied to the reads.
 */
public class PSBwaFilter {

    private final BwaMemIndex bwaIndex;
    private final int minIdent;
    private final int minSeedLength, numThreads;
    private final boolean bPaired;

    public PSBwaFilter(final String indexFileName, final int minIdent,
                       final int minSeedLength, final int numThreads, final boolean bPaired) {
        this.bwaIndex = BwaMemIndexCache.getInstance(indexFileName);
        this.minIdent = minIdent;
        this.minSeedLength = minSeedLength;
        this.numThreads = numThreads;
        this.bPaired = bPaired;
    }

    public Iterator<GATKRead> apply(final Iterator<GATKRead> itr) {

        //Initialize aligner
        final BwaMemAligner aligner = new BwaMemAligner(bwaIndex);
        if (bPaired) {
            aligner.alignPairs();
        }
        aligner.setMaxXAHitsAltOption(0);
        aligner.setMaxXAHitsOption(0);
        aligner.setMinSeedLengthOption(minSeedLength);
        aligner.setOutputScoreThresholdOption(0);
        aligner.setNThreadsOption(numThreads);

        //Collect reads
        final List<GATKRead> reads = new ArrayList<>();
        while (itr.hasNext()) {
            reads.add(itr.next());
        }

        final int numReads = reads.size();
        if (bPaired && (numReads & 1) != 0) {
            throw new GATKException("Cannot do paired alignment with an odd number of reads");
        }
        if (numReads == 0) {
            return Collections.emptyIterator();
        }

        //Do alignment
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(reads, GATKRead::getBases);

        //Filter reads if they map sufficiently well to the reference
        final HostAlignmentReadFilter hostFilter = new HostAlignmentReadFilter(minIdent);
        final ArrayList<GATKRead> results = new ArrayList<>(numReads);
        for (int i = 0; i < numReads; i++) {
            if (testReadAlignments(hostFilter, alignments.get(i))) {
                results.add(reads.get(i));
            }
        }
        results.trimToSize();
        return results.iterator();
    }

    private static boolean testReadAlignments(final HostAlignmentReadFilter hostFilter, final List<BwaMemAlignment> alignmentList) {
        for (final BwaMemAlignment alignment : alignmentList) {
            if (alignment.getCigar() != null && !alignment.getCigar().isEmpty()
                    && !hostFilter.test(TextCigarCodec.decode(alignment.getCigar()), alignment.getNMismatches())) {
                return false;
            }
        }
        return true;
    }
}