package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

public class JexlExpressionReadTagValueFilterUnitTest extends GATKBaseTest {

    @Test(dataProvider= "JexlExpressionReadTagValueFilterDataProvider")
    public void testJexlExpressionReadTagValueFilter(final String cigarString,
                                      final Object[] attrsNameAndValue,
                                      final String[] jexlExpr,
                                      final boolean expectedResult,
                                       final Class<? extends Throwable> expectedException) {

        final JexlExpressionReadTagValueFilter filter;

        // test different constructors here as well
        if ( jexlExpr.length == 1 ) {
            filter = new JexlExpressionReadTagValueFilter(jexlExpr[0]);
        } else {
            filter = new JexlExpressionReadTagValueFilter(Arrays.asList(jexlExpr));
        }

        final GATKRead read = ReadTagValueFilterUnitTest.buildSAMRead(cigarString, attrsNameAndValue);
        try {
            Assert.assertEquals(filter.test(read), expectedResult, cigarString);
            Assert.assertNull(expectedException);
        } catch (Throwable e) {
            if ( expectedException == null || !expectedException.isInstance(e) )
                throw e;
        }
    }

    @DataProvider(name = "JexlExpressionReadTagValueFilterDataProvider")
    public Iterator<Object[]> jexlExpressionReadTagValueFilterDataProvider() {
        final List<Object[]> result = new LinkedList<>();

        result.add(new Object[] {
                "100M",                             // cigar
                new Object[] {"TM", 1.0f},          // attributes
                new String[] {"TM == 1.0"},         // jexl expressions
                Boolean.TRUE,                       // expected
                null                                // expected exception
        });

        result.add(new Object[] {
                "100M",                             // cigar
                new Object[] {"TM", 1.0f},          // attributes
                new String[] {"TM < 1.0"},          // jexl expressions
                Boolean.FALSE,                      // expected
                null                                // expected exception
        });

        result.add(new Object[] {
                "100M",                             // cigar
                new Object[] {"TM", 1.0f},           // attributes
                new String[] {"NO_SUCH < 1.0"},     // jexl expressions
                Boolean.FALSE,                      // expected
                IllegalArgumentException.class      // expected exception
        });

        result.add(new Object[] {
                "100M",                             // cigar
                new Object[] {"TM", 1.0f, "TA", 2.0f},  // attributes
                new String[] {"TM >= 1.0", "TA <= 2.0"},  // jexl expressions
                Boolean.TRUE,                      // expected
                null                                // expected exception
        });

        result.add(new Object[] {
                "100M",                             // cigar
                new Object[] {"TM", 1.0f, "TA", 2.0f},  // attributes
                new String[] {"TM >= 1.0", "TA < 2.0"},  // jexl expressions
                Boolean.FALSE,                      // expected
                null                                // expected exception
        });
        return result.iterator();
    }


}
