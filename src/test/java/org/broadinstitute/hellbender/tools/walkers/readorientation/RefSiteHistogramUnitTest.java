package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

/**
 * Created by tsato on 10/13/17.
 */
public class RefSiteHistogramUnitTest {
    @Test
    public void test() throws IOException{
        String referenceContext = "ACT";
        RefSiteHistogram histogram = new RefSiteHistogram(referenceContext);
        // depth must start at 1
        for (int i = 1; i < 400; i++){
            histogram.increment(i);
        }

        File table = File.createTempFile("histogram", ".tsv");
        RefSiteHistogram.writeRefSiteHistograms(Collections.singletonList(histogram), table);

        List<RefSiteHistogram> readHistograms = RefSiteHistogram.readRefSiteHistograms(table);
        RefSiteHistogram readHistogram = readHistograms.get(0);
        Assert.assertEquals(readHistogram.getReferenceContext(), referenceContext);
        int[] counts = readHistogram.getCounts();
        for (int j = 0; j < 199; j++){
            Assert.assertEquals(counts[j], 1);
        }

        Assert.assertEquals(counts[199], 200);
    }
}