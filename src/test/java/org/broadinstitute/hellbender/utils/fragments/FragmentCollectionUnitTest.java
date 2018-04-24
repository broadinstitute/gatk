package org.broadinstitute.hellbender.utils.fragments;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class FragmentCollectionUnitTest extends GATKBaseTest {

    @Test
    public void createFromMultiSamplePileup() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, "10M");
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, "10M");
        read1.setPosition(new SimpleInterval("22", 200, 210));
        read2.setPosition(new SimpleInterval("22", 208, 218));
        read1.setMatePosition(read2);
        read2.setMatePosition(read1);
        final Locatable loc = new SimpleInterval("22", 208, 208);
        final Map<String, ReadPileup> stratified = new LinkedHashMap<>();
        stratified.put("sample1", new ReadPileup(loc, Arrays.asList(read2), 0));
        stratified.put("sample2", new ReadPileup(loc, Arrays.asList(read1), 9));
        final ReadPileup combined = new ReadPileup(loc, stratified);
        final FragmentCollection<PileupElement> elements = FragmentCollection.create(combined);
        Assert.assertTrue(elements.getSingletonReads().isEmpty());
        Assert.assertEquals(elements.getOverlappingPairs().size(), 1);
    }

}
