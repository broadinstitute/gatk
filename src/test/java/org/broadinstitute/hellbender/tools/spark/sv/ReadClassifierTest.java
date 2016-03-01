package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ReadClassifierTest extends BaseTest {
    @Test(groups = "spark")
    void restOfFragmentSizeTest() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(2, 1, 10000000, 1);
        final String groupName = header.getReadGroups().get(0).getReadGroupId();
        final int readSize = 151;
        final int fragmentLen = 475;
        final ReadMetadata.ReadGroupFragmentStatistics groupStats = new ReadMetadata.ReadGroupFragmentStatistics(fragmentLen, 25.f);
        final ReadMetadata readMetadata = new ReadMetadata(header, Collections.singletonList(groupStats), groupStats);
        final String templateName = "xyzzy";
        final int leftStart = 1010101;
        final int rightStart = leftStart + fragmentLen - readSize;
        final List<GATKRead> readPair = ArtificialReadUtils.createPair(header, templateName, readSize, leftStart, rightStart, true, false);
        final GATKRead read = readPair.get(0);
        read.setReadGroup(groupName);
        final ReadClassifier classifier = new ReadClassifier(readMetadata);
        checkClassification(classifier, read, Collections.emptyList());
        read.setCigar(ReadClassifier.MIN_SOFT_CLIP_LEN+"S"+(readSize-ReadClassifier.MIN_SOFT_CLIP_LEN)+"M");
        checkClassification(classifier, read, Collections.singletonList(new BreakpointEvidence.SplitRead(read, readMetadata, true)));
        read.setCigar((readSize-ReadClassifier.MIN_SOFT_CLIP_LEN)+"M"+ReadClassifier.MIN_SOFT_CLIP_LEN+"S");
        checkClassification(classifier, read, Collections.singletonList(new BreakpointEvidence.SplitRead(read, readMetadata, false)));
        read.setCigar((readSize/2)+"M"+ReadClassifier.MIN_INDEL_LEN+"D"+((readSize+1)/2)+"M");
        final int locus = leftStart + readSize/2 + ReadClassifier.MIN_INDEL_LEN/2;
        checkClassification(classifier, read, Collections.singletonList(new BreakpointEvidence.LargeIndel(read, readMetadata, locus)));
        read.setCigar(readSize+"M");
        read.setMateIsUnmapped();
        checkClassification(classifier, read, Collections.singletonList(new BreakpointEvidence.MateUnmapped(read, readMetadata)));
        read.setMatePosition(read.getContig(), rightStart);
        read.setIsReverseStrand(true);
        checkClassification(classifier, read, Collections.singletonList(new BreakpointEvidence.SameStrandPair(read, readMetadata)));
        read.setIsReverseStrand(false);
        read.setMateIsReverseStrand(false);
        checkClassification(classifier, read, Collections.singletonList(new BreakpointEvidence.SameStrandPair(read, readMetadata)));
        read.setIsReverseStrand(true);
        checkClassification(classifier, read, Collections.singletonList(new BreakpointEvidence.OutiesPair(read, readMetadata)));
        read.setIsReverseStrand(false);
        read.setMatePosition(header.getSequenceDictionary().getSequence(1).getSequenceName(), rightStart);
        checkClassification(classifier, read, Collections.singletonList(new BreakpointEvidence.InterContigPair(read, readMetadata)));
    }

    void checkClassification( final ReadClassifier classifier, final GATKRead read, final List<BreakpointEvidence> expectedEvidence ) {
        final List<BreakpointEvidence> evidence = new ArrayList<>();
        classifier.apply(read).forEachRemaining(evidence::add);
        Assert.assertEquals(evidence, expectedEvidence);
    }
}
