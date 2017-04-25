package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class BreakpointEvidenceTest extends BaseTest {
    @Test(groups = "spark")
    void restOfFragmentSizeTest() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(1, 1, 10000000, 1);
        final String groupName = header.getReadGroups().get(0).getReadGroupId();
        final int readSize = 151;
        final ReadMetadata.ReadGroupFragmentStatistics groupStats = new ReadMetadata.ReadGroupFragmentStatistics(401, 175, 20);
        final ReadMetadata readMetadata = new ReadMetadata(Collections.emptySet(), header, groupStats, 1, 1L, 1L, 1);
        final String templateName = "xyzzy";
        final int readStart = 1010101;
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, templateName, 0, readStart, readSize);
        read.setIsPaired(false);
        read.setIsReverseStrand(true);
        read.setReadGroup(groupName);
        final BreakpointEvidence evidence1 = new BreakpointEvidence(read, readMetadata);
        final int uncertainty = (readMetadata.getStatistics(groupName).getMedianFragmentSize()-readSize)/2;
        final int evidenceLocus = readStart - uncertainty;
        final BreakpointEvidence evidence2 = new BreakpointEvidence(read, readMetadata, evidenceLocus, uncertainty);
        Assert.assertEquals(evidence1.getLocation().getContig(), 0);
        Assert.assertEquals(evidence1.getLocation().getStart(), evidenceLocus-uncertainty);
        Assert.assertEquals(evidence1.getLocation().getEnd(), evidenceLocus+uncertainty);
        Assert.assertEquals(evidence1.getLocation().getLength(), 2*uncertainty);
        Assert.assertEquals(evidence1.getTemplateName(), templateName);
        Assert.assertEquals(evidence1.getTemplateEnd(), BreakpointEvidence.TemplateEnd.UNPAIRED);
        Assert.assertEquals(evidence1.toString(), evidence2.toString());
        read.setIsReverseStrand(false);
        final BreakpointEvidence evidence3 = new BreakpointEvidence(read, readMetadata);
        final BreakpointEvidence evidence4 = new BreakpointEvidence(read, readMetadata, readStart+readSize+uncertainty, uncertainty);
        Assert.assertEquals(evidence3.toString(), evidence4.toString());
    }

    @Test(groups = "spark")
    void serializationTest() {
        final List<BreakpointEvidence> evidenceList = new ArrayList<>(7);
        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader();
        final ReadMetadata.ReadGroupFragmentStatistics groupStatistics = new ReadMetadata.ReadGroupFragmentStatistics(400, 175, 20);
        final ReadMetadata metadata = new ReadMetadata(Collections.emptySet(), samHeader, groupStatistics, 1, 2L, 2L, 1);
        final List<GATKRead> readPair = ArtificialReadUtils.createPair(samHeader, "firstReadPair", 101, 1010, 1382, false, false);
        final GATKRead read = readPair.get(0);
        evidenceList.add(new BreakpointEvidence.SplitRead(read, metadata, true));
        evidenceList.add(new BreakpointEvidence.LargeIndel(read, metadata, read.getStart()+50));
        evidenceList.add(new BreakpointEvidence.MateUnmapped(read, metadata));
        evidenceList.add(new BreakpointEvidence.InterContigPair(read, metadata));
        evidenceList.add(new BreakpointEvidence.OutiesPair(read, metadata));
        evidenceList.add(new BreakpointEvidence.SameStrandPair(read, metadata));
        evidenceList.add(new BreakpointEvidence.WeirdTemplateSize(read, metadata));

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, evidenceList);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final List<BreakpointEvidence> evidenceList2 = (List<BreakpointEvidence>)kryo.readClassAndObject(in);
        Assert.assertEquals(evidenceList.size(), evidenceList2.size());
        for ( int idx = 0; idx != evidenceList.size(); ++idx ) {
            Assert.assertEquals(evidenceList.get(idx).toString(), evidenceList2.get(idx).toString());
        }
    }
}
