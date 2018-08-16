package org.broadinstitute.hellbender;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.testutils.ReadTestUtils;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public final class BwaMemIntegrationTest extends GATKBaseTest {

    private BwaMemIndex index;
    private File fastaFile;
    private File imageFile;
    private FastaSequenceIndex fastaIndex;

    private static final SAMSequenceDictionary TEST_DICTIONARY = new SAMSequenceDictionary(
            Arrays.asList(new SAMSequenceRecord("chr1", 10_000),
                          new SAMSequenceRecord("chr2", 10_000))
    );

    @BeforeClass
    public void loadIndex() throws IOException {
        final RandomDNA rdnDNA = new RandomDNA(103);
        try {
            fastaFile = rdnDNA.nextFasta(TEST_DICTIONARY, 60);
            fastaFile.deleteOnExit();
            imageFile = createTempFile(fastaFile.getAbsolutePath(), ".img");
            fastaIndex = FastaSequenceIndexCreator.buildFromFasta(fastaFile.toPath());
            BwaMemIndex.createIndexImageFromFastaFile(fastaFile.getAbsolutePath(), imageFile.getAbsolutePath());
            index = new BwaMemIndex(imageFile.getAbsolutePath());
        } catch (final RuntimeException re) {
            unloadIndex();
            throw re;
        }
    }

    @AfterClass
    public void unloadIndex() {
        if (index != null) { try { index.close(); index = null; } catch (final RuntimeException ex) {}; }
        if (imageFile != null) { try { imageFile.delete(); imageFile = null; } catch (final RuntimeException ex) {}; }
        if (fastaFile != null) { try { fastaFile.delete(); fastaFile = null; } catch (final RuntimeException ex) {}; }
    }

    @Test
    public void testPerfectUnpairedMapping() throws Exception {
        final Random rdn = new Random(13);
        final IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(fastaFile, fastaIndex);
        final List<SAMRecord> inputReads = ReadTestUtils.randomErrorFreeUnpairedReads(rdn, TEST_DICTIONARY, fasta,
                "read-", 0, 1000, 50, 200);
        final BwaMemAligner aligner = new BwaMemAligner(index);
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(inputReads, SAMRecord::getReadBases);
        Assert.assertEquals(alignments.size(), inputReads.size());
        for (int i = 0; i < alignments.size(); i++) {
            final List<BwaMemAlignment> blocks = alignments.get(i);
            final SAMRecord inputRead  = inputReads.get(i);
            Assert.assertNotNull(blocks);
            Assert.assertEquals(blocks.size(), 1);
            assertPerfectAlignmentMatch(inputRead, blocks.get(0));
        }
    }

    @Test
    public void testChimericUnpairedMapping() throws Exception {
        final Random rdn = new Random(13);
        final IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(fastaFile, fastaIndex);
        final List<SAMRecord> inputRecords = ReadTestUtils.randomErrorFreeUnpairedReads(rdn, TEST_DICTIONARY, fasta,
                "read-", 0, 1000 * 2, 100, 200);
        final List<List<SAMRecord>> chimerics = new ArrayList<>();
        for (int i = 0; i < inputRecords.size(); i += 2) {
            chimerics.add(Arrays.asList(inputRecords.get(i), inputRecords.get(i + 1)));
        }
        final BwaMemAligner aligner = new BwaMemAligner(index);
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(chimerics,
                ll -> Utils.concat(ll.get(0).getReadBases(), ll.get(1).getReadBases()));
        Assert.assertEquals(alignments.size(), chimerics.size());
        for (int i = 0; i < alignments.size(); i++) {
            final List<BwaMemAlignment> blocks = alignments.get(i);
            final List<SAMRecord> chimeric = chimerics.get(i);
            Assert.assertNotNull(blocks);
            Assert.assertTrue(blocks.size() <= 2);
            if (blocks.size() == 1) { // in occasion they record pair may be close enough and have the same orientation
                                      // resulting in a single record... we only perform simple checks in those few cases:
                Assert.assertEquals(chimeric.get(0).getReferenceIndex(), chimeric.get(1).getReferenceIndex());
                Assert.assertEquals(chimeric.get(0).getReadNegativeStrandFlag(), chimeric.get(1).getReadNegativeStrandFlag());
                Assert.assertEquals(blocks.get(0).getRefId(), chimeric.get(0).getReferenceIndex().intValue());
                final int maxDistance = chimeric.get(0).getReadLength() + chimeric.get(1).getReadLength() + 100;
                Assert.assertTrue(Math.abs(blocks.get(0).getRefStart() - chimeric.get(0).getAlignmentStart()) < maxDistance);
                Assert.assertTrue(Math.abs(blocks.get(0).getRefStart() - chimeric.get(0).getAlignmentStart()) < maxDistance);
            } else {
                BwaMemAlignment first  = blocks.get(0);
                BwaMemAlignment second = blocks.get(1);
                // We need to match the output alignments pair elements with the input records pair elements.
                // The test bellow is not bullet proof but is quite unlikely that it will fail to spot
                // a case in where the switch is needed:
                if ((first.getRefId() == second.getRefId()
                            && Math.abs(first.getRefStart() - chimeric.get(0).getStart()) > Math.abs(first.getRefStart() - chimeric.get(1).getStart()))
                        ||  (first.getRefId() != second.getRefId()
                            && second.getRefId() == TEST_DICTIONARY.getSequenceIndex(chimeric.get(0).getContig()))) {
                   first = second;
                   second = blocks.get(0);
                }
                assertChimericAlignmentMatch(chimeric.get(0), first);
                assertChimericAlignmentMatch(chimeric.get(1), second);
                // one and only one is a supplementary alignment:
                Assert.assertNotEquals(SAMFlag.SUPPLEMENTARY_ALIGNMENT.isSet(first.getSamFlag()),
                        SAMFlag.SUPPLEMENTARY_ALIGNMENT.isSet(second.getSamFlag()));
            }
        }
    }

    private void assertPerfectAlignmentMatch(final SAMRecord input, final BwaMemAlignment output) {
        assertPrimaryUnpairedOutput(input, output);
        Assert.assertEquals(output.getCigar(), "" + input.getReadBases().length + "M");
        Assert.assertEquals(output.getRefId(), input.getReferenceIndex().intValue());
        Assert.assertEquals(output.getRefStart(), input.getAlignmentStart() - 1);
        Assert.assertTrue(SAMFlag.SUPPLEMENTARY_ALIGNMENT.isUnset(output.getSamFlag()));
        Assert.assertEquals(output.getNMismatches(), 0);
    }

    private void assertChimericAlignmentMatch(final SAMRecord input, final BwaMemAlignment output) {
        assertPrimaryUnpairedOutput(input, output);
        Assert.assertEquals(output.getRefId(), input.getReferenceIndex().intValue());
        final SimpleInterval matchInterval = new SimpleInterval(
                input.getContig(),
                output.getRefStart() + 1,
                output.getRefEnd());
        Assert.assertTrue(matchInterval.contains(input));
        Assert.assertEquals(SAMFlag.READ_REVERSE_STRAND.isSet(output.getSamFlag()), input.getReadNegativeStrandFlag());
    }

    private void assertPrimaryUnpairedOutput(final SAMRecord input, final BwaMemAlignment output) {
        Assert.assertTrue(SAMFlag.SECONDARY_ALIGNMENT.isUnset(output.getSamFlag()));
        Assert.assertEquals(SAMFlag.READ_REVERSE_STRAND.isSet(output.getSamFlag()), input.getReadNegativeStrandFlag());
        Assert.assertTrue(SAMFlag.SECONDARY_ALIGNMENT.isUnset(output.getSamFlag()));
        Assert.assertTrue(SAMFlag.PROPER_PAIR.isUnset(output.getSamFlag()));
        Assert.assertTrue(SAMFlag.READ_UNMAPPED.isUnset(output.getSamFlag()));
        Assert.assertTrue(SAMFlag.READ_PAIRED.isUnset(output.getSamFlag()));
        Assert.assertTrue(output.getMapQual() >= 20);
    }

}
