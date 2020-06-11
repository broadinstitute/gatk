package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class PSBwaUtilTest extends GATKBaseTest {

    @Test
    public void testAddReferenceSequencesToHeader() {
        final SAMFileHeader header = new SAMFileHeader();
        final String dictionaryPath = hg19_chr1_1M_dict;
        PSBwaUtils.addReferenceSequencesToHeader(header, dictionaryPath);
        final SAMSequenceDictionary refDict = ReferenceUtils.loadFastaDictionary(new File(dictionaryPath));
        Assert.assertEquals(refDict.size(), header.getSequenceDictionary().size());
        for (final SAMSequenceRecord rec : refDict.getSequences()) {
            final SAMSequenceRecord recTest = header.getSequenceDictionary().getSequence(rec.getSequenceName());
            Assert.assertNotNull(recTest);
            Assert.assertEquals(rec, recTest);
        }
    }
}
