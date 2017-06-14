package org.broadinstitute.hellbender.tools.spark.pathseq;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Loads Bwa index and aligns reads. This is not a general Bwa aligner because it only retrieves primary alignments.
 * If preserveAlignments is set, existing alignment data in each read is not overwritten, and any new alignments
 * are added to the SA tag.
 */
public final class PSBwaAligner {

    private final BwaMemIndex bwaIndex;
    private final int minSeedLength, numThreads, maxAltHits, scoreThreshold;
    private final boolean pairedAlignment;

    public PSBwaAligner(final PSBwaArgumentCollection bwaArgs, final boolean pairedAlignment) {
        this.bwaIndex = BwaMemIndexSingleton.getInstance(bwaArgs.bwaImage);
        this.minSeedLength = bwaArgs.seedLength;
        this.numThreads = bwaArgs.bwaThreads;
        this.maxAltHits = bwaArgs.maxAlternateHits;
        this.scoreThreshold = bwaArgs.scoreThreshold;
        this.pairedAlignment = pairedAlignment;
    }

    public Iterator<GATKRead> apply(final Iterator<GATKRead> itr, final SAMFileHeader header) {
        //Create aligner and set options
        final BwaMemAligner aligner = new BwaMemAligner(bwaIndex);
        if (pairedAlignment) {
            aligner.alignPairs();
        }
        aligner.setMaxXAHitsAltOption(maxAltHits);
        aligner.setMaxXAHitsOption(maxAltHits);
        aligner.setMinSeedLengthOption(minSeedLength);
        aligner.setOutputScoreThresholdOption(scoreThreshold);
        aligner.setNThreadsOption(numThreads);

        //Get list of reads on the partition
        final List<GATKRead> reads = new ArrayList<>();
        while (itr.hasNext()) {
            reads.add(itr.next());
        }
        final int numReads = reads.size();
        if (pairedAlignment && numReads % 2 != 0) {
            throw new UserException.BadInput("Expected paired reads but there are an odd number");
        }

        //Align read sequences
        final List<Tuple2<GATKRead, List<BwaMemAlignment>>> results;
        if (numReads == 0) {
            results = Collections.emptyList();
        } else {
            final List<byte[]> seqs = reads.stream().map(GATKRead::getBases).collect(Collectors.toList());
            final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(seqs);
            results = new ArrayList<>(numReads);
            for (int i = 0; i < numReads; i++) {
                results.add(new Tuple2<>(reads.get(i), alignments.get(i)));
            }
        }

        //Map alignments back to the reads
        final List<String> refNames = bwaIndex.getReferenceContigNames();
        return Utils.stream(results).map(item -> {
            final GATKRead read = item._1;
            SAMRecord rec = read.convertToSAMRecord(header);
            final List<BwaMemAlignment> alignmentList = item._2;
            for (final BwaMemAlignment alignment : alignmentList) {
                //Only get primary alignments
                if (SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(alignment.getSamFlag())
                        && SAMFlag.SUPPLEMENTARY_ALIGNMENT.isUnset(alignment.getSamFlag())) {
                    if (rec.getReadUnmappedFlag()) {
                        //Record is currently unmapped, so apply first alignment
                        rec = BwaMemAlignmentUtils.applyAlignment(rec.getReadName(), rec.getReadBases(), rec.getBaseQualities(),
                                rec.getReadGroup().getId(), alignment, refNames, header, false, false);
                    } else if (alignment.getRefId() != -1) {
                        //Record is currently mapped, but if the alignment is not empty then add as an alternate alignment
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