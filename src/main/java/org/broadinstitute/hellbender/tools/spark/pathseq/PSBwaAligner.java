package org.broadinstitute.hellbender.tools.spark.pathseq;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.bwa.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Loads Bwa index and aligns reads. This is not a general Bwa aligner because it only retrieves primary alignments.
 * If preserveAlignments is set, existing alignment data in each read is not overwritten, and any new alignments
 * are added to the SA tag.
 */
public final class PSBwaAligner {

    private final BwaMemIndex bwaIndex;
    private final PSBwaArgumentCollection bwaArgs;
    private final boolean pairedAlignment;

    public PSBwaAligner(final PSBwaArgumentCollection bwaArgs, final boolean pairedAlignment) {
        this.bwaIndex = BwaMemIndexCache.getInstance(bwaArgs.bwaImage);
        this.bwaArgs = bwaArgs;
        this.pairedAlignment = pairedAlignment;
    }

    private static GATKRead applyAlignments(GATKRead read, final List<BwaMemAlignment> alignmentList,
                                            final List<String> refNames, final SAMFileHeader header) {
        final List<String> saTags = new ArrayList<>(alignmentList.size());
        for (final BwaMemAlignment alignment : alignmentList) {
            //Only get primary alignments
            if (SAMFlag.SECONDARY_ALIGNMENT.isUnset(alignment.getSamFlag())
                    && SAMFlag.SUPPLEMENTARY_ALIGNMENT.isUnset(alignment.getSamFlag())) {
                if (read.isUnmapped()) {
                    //Record is currently unmapped, so apply first alignment
                    read = new SAMRecordToGATKReadAdapter(BwaMemAlignmentUtils.applyAlignment(read.getName(),
                            read.getBases(), read.getBaseQualities(), read.getReadGroup(), alignment, refNames,
                            header, false, false));
                } else if (alignment.getRefId() != -1) {
                    //Record is currently mapped, but if the alignment is not empty then add as an alternate alignment
                    //This can happen if the input read is already mapped
                    saTags.add(BwaMemAlignmentUtils.asTag(alignment, refNames));
                }
            }
        }
        if (!saTags.isEmpty()) {
            final String currentTag = read.getAttributeAsString("SA") != null ? read.getAttributeAsString("SA") : "";
            read.setAttribute("SA", currentTag + String.join(";", saTags));
        }
        return read;
    }

    public Iterator<GATKRead> apply(final Iterator<GATKRead> itr, final SAMFileHeader header) {
        //Create aligner and set options
        final BwaMemAligner aligner = new BwaMemAligner(bwaIndex);
        if (pairedAlignment) {
            aligner.alignPairs();
        }
        aligner.setMaxXAHitsAltOption(bwaArgs.maxAlternateHits);
        aligner.setMaxXAHitsOption(bwaArgs.maxAlternateHits);
        aligner.setMinSeedLengthOption(bwaArgs.seedLength);
        aligner.setOutputScoreThresholdOption(bwaArgs.scoreThreshold);
        aligner.setNThreadsOption(bwaArgs.bwaThreads);

        //Get list of reads on the partition
        final List<GATKRead> reads = new ArrayList<>();
        while (itr.hasNext()) {
            reads.add(itr.next());
        }

        final int numReads = reads.size();
        if (pairedAlignment && numReads % 2 != 0) {
            throw new UserException.BadInput("Expected paired reads but there are an odd number");
        }
        if (numReads == 0) {
            return new ArrayList<GATKRead>(0).iterator();
        }

        //Align read sequences
        final List<String> refNames = bwaIndex.getReferenceContigNames();
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(reads, GATKRead::getBases);
        for (int i = 0; i < reads.size(); i++) {
            reads.set(i, applyAlignments(reads.get(i), alignments.get(i), refNames, header));
        }
        return reads.iterator();
    }
}