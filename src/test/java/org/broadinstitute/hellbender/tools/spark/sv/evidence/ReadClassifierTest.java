package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

public class ReadClassifierTest extends BaseTest {
    @Test(groups = "spark")
    void restOfFragmentSizeTest() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(3, 1, 10000000, 1);
        final String groupName = header.getReadGroups().get(0).getReadGroupId();
        final int readSize = 151;
        final int fragmentLen = 400;
        final ReadMetadata.LibraryFragmentStatistics groupStats = new ReadMetadata.LibraryFragmentStatistics(fragmentLen, 175, 20);
        final Set<Integer> crossContigIgnoreSet = new HashSet<>(3);
        crossContigIgnoreSet.add(2);
        final ReadMetadata readMetadata = new ReadMetadata(crossContigIgnoreSet, header, groupStats, null, 2L, 2L, 1);
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
        read.setMatePosition(read.getContig(), read.getStart() + 2);
        checkClassification(classifier, read, Collections.emptyList());
        read.setMatePosition(read.getContig(), read.getStart() + 2 + ReadClassifier.ALLOWED_SHORT_FRAGMENT_OVERHANG);
        checkClassification(classifier, read, Collections.singletonList(new BreakpointEvidence.OutiesPair(read, readMetadata)));

        read.setIsReverseStrand(false);
        read.setMateIsReverseStrand(true);
        read.setMatePosition(read.getContig(), read.getStart() - 2);
        checkClassification(classifier, read, Collections.emptyList());
        read.setIsReverseStrand(false);

        read.setMatePosition(header.getSequenceDictionary().getSequence(1).getSequenceName(), rightStart);
        checkClassification(classifier, read, Collections.singletonList(new BreakpointEvidence.InterContigPair(read, readMetadata)));

        read.setMatePosition(header.getSequenceDictionary().getSequence(2).getSequenceName(), rightStart);
        checkClassification(classifier, read, Collections.emptyList());
    }

    private void checkClassification( final ReadClassifier classifier, final GATKRead read, final List<BreakpointEvidence> expectedEvidence ) {
        final List<BreakpointEvidence> evidence = new ArrayList<>();
        classifier.apply(read).forEachRemaining(evidence::add);
        Assert.assertEquals(evidence.size(), expectedEvidence.size());
        for ( int idx = 0; idx != evidence.size(); ++idx ) {
            Assert.assertEquals(evidence.get(idx).toString(), expectedEvidence.get(idx).toString());
        }
    }
}
