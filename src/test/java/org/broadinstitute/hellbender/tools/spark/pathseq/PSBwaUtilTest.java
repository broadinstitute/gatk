package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class PSBwaUtilTest extends GATKBaseTest {

    @Test
    public void testAddReferenceSequencesToHeader() {
        final SAMFileHeader header = new SAMFileHeader();
        final SerializableFunction<GATKRead, SimpleInterval> windowFunction = ReferenceWindowFunctions.IDENTITY_FUNCTION;
        final String referencePath = hg19MiniReference;
        PSBwaUtils.addReferenceSequencesToHeader(header, referencePath, windowFunction);
        final ReferenceMultiSparkSource ref = new ReferenceMultiSparkSource(referencePath, windowFunction);
        Assert.assertEquals(ref.getReferenceSequenceDictionary(null).size(), header.getSequenceDictionary().size());
        for (final SAMSequenceRecord rec : ref.getReferenceSequenceDictionary(null).getSequences()) {
            final SAMSequenceRecord recTest = header.getSequenceDictionary().getSequence(rec.getSequenceName());
            Assert.assertNotNull(recTest);
            Assert.assertEquals(rec, recTest);
        }
    }
}
