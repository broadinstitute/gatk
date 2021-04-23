package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class AS_StrandBiasTestUnitTest extends GATKBaseTest {

    private static final Allele REF = Allele.create("T", true);
    private static final Allele ALT = Allele.create("A", false);

    /**
     * Test for issue #6766
     */
    @Test
    public void testCombineRawDataNullPointerExceptionIssue6766() {
        final List<Allele> vcAlleles = Arrays.asList(REF, ALT);
        final List<ReducibleAnnotationData<?>> combinedVCdata = new ArrayList<>();
        // The commented out code would be a normal way to create the combinedVCdata, but this did not reproduce the null pointer exception
        // we are testing for. The following one line did even though that probably would not occur if reading in from a gvcf.
//        ReducibleAnnotationData<List<Integer>> data = new ReducibleAnnotationData<>(null);
//        data.putAttribute(ALT, Arrays.asList(33640, 10));
//        data.putAttribute(REF, Arrays.asList(1,2));
//        combinedVCdata.add(data);

        // This set SB data for "." (i.e. the missing allele)
        combinedVCdata.add(new ReducibleAnnotationData<>("36000,10"));  //10 MQ60 reads

        AS_StrandBiasTest annotator = new AS_StrandOddsRatio();

        try {
            final Map<String, Object> combined = annotator.combineRawData(vcAlleles, combinedVCdata);
            final String combinedListString = (String) combined.get(annotator.getPrimaryRawKey());
            // These are both 0 because we did not set any SB data for either allele
            Assert.assertEquals(combinedListString, "0,0|0,0");
        } catch (NullPointerException npe) {
            Assert.fail("This code should not be throwing a NullPointerException");
        }

    }
}
