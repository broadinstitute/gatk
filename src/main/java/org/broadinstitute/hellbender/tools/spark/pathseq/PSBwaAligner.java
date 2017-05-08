package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Aligns reads using BWA. If bPreserveAlignments is set, existing alignment data in each read is not overwritten, and
 * any new alignments are added to the SA tag.
 */
public final class PSBwaAligner {

    private final BwaMemIndex bwaIndex;
    private final SAMFileHeader readsheader;
    private final int minSeedLength, numThreads, maxAltHits, scoreThreshold;
    private final boolean bPaired, bPreserveAlignments;

    public PSBwaAligner(final String indexFileName, final SAMFileHeader readsHeader, final int minSeedLength,
                        final int numThreads, final int maxAltHits, final int scoreThreshold, final boolean bPaired,
                        final boolean bPreserveAlignments) {

        this.bwaIndex = BwaMemIndexSingleton.getInstance(indexFileName);
        this.readsheader = readsHeader;
        this.minSeedLength = minSeedLength;
        this.numThreads = numThreads;
        this.maxAltHits = maxAltHits;
        this.scoreThreshold = scoreThreshold;
        this.bPaired = bPaired;
        this.bPreserveAlignments = bPreserveAlignments;
    }

    public Iterator<GATKRead> apply(final Iterator<GATKRead> itr) {
        final BwaMemAligner aligner = new BwaMemAligner(bwaIndex);
        if (bPaired) {
            aligner.alignPairs();
        }
        aligner.setMaxXAHitsAltOption(maxAltHits);
        aligner.setMaxXAHitsOption(maxAltHits);
        aligner.setMinSeedLengthOption(minSeedLength);
        aligner.setOutputScoreThresholdOption(scoreThreshold);
        aligner.setNThreadsOption(numThreads);

        final List<GATKRead> reads = new ArrayList<>();
        while (itr.hasNext()) {
            reads.add(itr.next());
        }
        final int numReads = reads.size();
        final List<byte[]> seqs = reads.stream().map(GATKRead::getBases).collect(Collectors.toList());
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(seqs);
        final List<Tuple2<GATKRead, List<BwaMemAlignment>>> results = new ArrayList<>(numReads);
        for (int i = 0; i < numReads; i++) {
            results.add(new Tuple2<>(reads.get(i), alignments.get(i)));
        }
        final List<String> refNames = bwaIndex.getReferenceContigNames();
        return Utils.stream(results).map(item -> {
            final GATKRead read = item._1;
            if (!bPreserveAlignments) {
                read.setIsUnmapped();
                read.clearAttribute("SA");
                read.clearAttribute("XA");
                read.clearAttribute("NM");
            }
            SAMRecord rec = read.convertToSAMRecord(readsheader);
            final List<BwaMemAlignment> alignmentList = item._2;
            for (final BwaMemAlignment alignment : alignmentList) {
                if (SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(alignment.getSamFlag())
                        && SAMFlag.SUPPLEMENTARY_ALIGNMENT.isUnset(alignment.getSamFlag())) {
                    if (rec.getReadUnmappedFlag()) {
                        rec = BwaMemAlignmentUtils.applyAlignment(rec.getReadName(), rec.getReadBases(), rec.getBaseQualities(),
                                rec.getReadGroup().getId(), alignment, refNames, readsheader, false, false);
                    } else if (alignment.getRefId() != -1) {
                        final String alignmentTag = BwaMemAlignmentUtils.asTag(alignment, refNames);
                        final String currentTag = rec.getStringAttribute("SA") != null ? rec.getStringAttribute("SA") : "";
                        rec.setAttribute("SA", currentTag + alignmentTag);
                    }
                }
            }
            return SAMRecordToGATKReadAdapter.headerlessReadAdapter(rec).copy();
        }).iterator();
    }
}