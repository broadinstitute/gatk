package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Map;

/**
 * Created by tsato on 6/30/17.
 */
public class GCContentUnitTest extends BaseTest {
    @Test
    public void testBasic() throws Exception {
        ReferenceDataSource refSource = new ReferenceFileSource(new File(hg19_chr1_1M_Reference));
        String chr = "1";
        int start = 800000;
        SimpleInterval interval = new SimpleInterval("1", start, start);
        ReferenceContext refContext = new ReferenceContext(refSource, interval);

        // an A -> C SNP at position 50
        VariantContextBuilder variantContextBuilder = new VariantContextBuilder("Megan", "1", start / 2, start / 2,
                Arrays.asList(Allele.create("A", true), Allele.create("C")));
        VariantContext variantContext = variantContextBuilder.make();

        final GCContent gccontentAnnotation = new GCContent();
        Map<String, Object> annotation = gccontentAnnotation.annotate(refContext, variantContext, null);

        final double observed = Double.valueOf(String.valueOf(annotation.get(gccontentAnnotation.getKeyNames().get(0))));
        Assert.assertEquals(observed, 0.405, 1e-3);
    }
}