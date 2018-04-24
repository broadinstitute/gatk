package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import java.nio.file.Path;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;

/**
 * Created by davidben on 3/23/17.
 */
public class ReferenceBasesUnitTest extends GATKBaseTest {

    @Test
    public void test() {
        final Path refFasta = IOUtils.getPath(b37_reference_20_21);

        final ReferenceDataSource refDataSource = new ReferenceFileSource(refFasta);
        final ReferenceContext ref = new ReferenceContext(refDataSource, new SimpleInterval("20", 10_000_000, 10_000_200));
        final VariantContext vc = new VariantContextBuilder("source", "20", 10_000_100, 10_000_100, Collections.singleton(Allele.create((byte) 'A', true))).make();
        final String refBases = (String) new ReferenceBases().annotate(ref, vc, null)
                .get(ReferenceBases.REFERENCE_BASES_KEY);
        Assert.assertEquals(refBases, "ACTGCATCCCTTGCATTTCC");
    }

    // Asserts that the code silently failed
    @Test
    public void TestNoReferenceBehavior() {
        final VariantContext vc = new VariantContextBuilder("source", "20", 10_000_100, 10_000_100, Collections.singleton(Allele.create((byte) 'A', true))).make();
        final String refBases = (String) new ReferenceBases().annotate(null, vc, null)
                .get(ReferenceBases.REFERENCE_BASES_KEY);
        Assert.assertNull(refBases);
    }
}