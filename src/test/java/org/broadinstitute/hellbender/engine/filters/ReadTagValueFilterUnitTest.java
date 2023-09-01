package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class ReadTagValueFilterUnitTest extends GATKBaseTest {

    private static final int CHR_COUNT = 1;
    private static final int CHR_START = 1;
    private static final int CHR_SIZE = 1000;
    private static final int GROUP_COUNT = 5;

    static protected GATKRead buildSAMRead(final String cigarString, Object[] attrsNameAndValue) {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(CHR_COUNT, CHR_START, CHR_SIZE, GROUP_COUNT);    final Cigar cigar = TextCigarCodec.decode(cigarString);
        final SAMRecord samRecord = ArtificialReadUtils.createArtificialSAMRecord(header, cigar);
        for ( int i = 0 ; i < attrsNameAndValue.length ; i += 2 )
            samRecord.setAttribute(attrsNameAndValue[i].toString(), attrsNameAndValue[i+1]);
        return new SAMRecordToGATKReadAdapter(samRecord);
    }

    @Test(dataProvider= "ReadTagValueFilterDataProvider")
    public void testReadTagValueFilter(final String cigarString,
                                      final Object[] attrsNameAndValue,
                                      final String tagName,
                                      final Float tagValue,
                                      final ReadTagValueFilter.Operator tagOp,
                                      final boolean expectedResult,
                                       final Class<? extends Throwable> expectedException) {

        final ReadTagValueFilter filter = new ReadTagValueFilter(tagName, tagValue, tagOp);

        final GATKRead read = buildSAMRead(cigarString, attrsNameAndValue);
        try {
            Assert.assertEquals(filter.test(read), expectedResult, cigarString);
            Assert.assertNull(expectedException);
        } catch (Throwable e) {
            if ( expectedException == null || !expectedException.isInstance(e) )
                throw e;
        }
    }

    @DataProvider(name = "ReadTagValueFilterDataProvider")
    public Iterator<Object[]> readTagValueFilterDataProvider() {
        final List<Object[]> result = new LinkedList<>();

        result.add(new Object[] {
                "100M",                             // cigar
                new Object[] {"TM", 1.0f},          // attributes
                "TM",                               // tagname
                1.0f,                               // tagvalue
                ReadTagValueFilter.Operator.EQUAL,  // tagop
                Boolean.TRUE,                       // expected
                null                                // expected exception
        });

        result.add(new Object[] {
                "100M",                             // cigar
                new Object[] {"TM", 1.0f},          // attributes
                "TM",                               // tagname
                1.0f,                               // tagvalue
                ReadTagValueFilter.Operator.LESS,   // tagop
                Boolean.FALSE,                      // expected
                null                                // expected exception
        });

        result.add(new Object[] {
                "100M",                             // cigar
                new Object[] {"TM", 1.0f},           // attributes
                "NO_SUCH",                          // tagname
                1.0f,                               // tagvalue
                ReadTagValueFilter.Operator.EQUAL,  // tagop
                Boolean.TRUE,                       // expected
                IllegalArgumentException.class      // expected exception
        });

        return result.iterator();
    }


}
