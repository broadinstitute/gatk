package org.broadinstitute.hellbender.engine.spark;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;
import java.util.Collections;

public class KnownSitesCachePerfUnitTest extends GATKBaseTest {
    @Test
    @SuppressWarnings("unchecked")
    public void test() throws Exception {

        long start = System.nanoTime();
        IntervalsSkipList<GATKVariant> variants = KnownSitesCache.getVariants(Collections.singletonList("/Users/tom/workspace/gatk/src/test/resources/large/dbsnp_138.b37.1.1-65M.vcf"));
        long end = System.nanoTime();

        System.out.println("Time (millis):" + ((end - start) / 1000000));

        File file = File.createTempFile("obj", ".ser");
        serialize(variants, new BufferedOutputStream(new FileOutputStream(file)));

        long start2 = System.nanoTime();
        IntervalsSkipList<GATKVariant> variants2 = (IntervalsSkipList<GATKVariant>) deserialize(new BufferedInputStream(new FileInputStream(file)));
        long end2 = System.nanoTime();

        System.out.println("Time ser (millis):" + ((end2 - start2) / 1000000));
    }

    private static void serialize(Object obj, OutputStream out) throws IOException {
        try (ObjectOutputStream oos = new ObjectOutputStream(out)) {
            oos.writeObject(obj);
        }
    }

    private static Object deserialize(InputStream in) throws IOException, ClassNotFoundException {
        try (ObjectInputStream ois = new ObjectInputStream(in)) {
            return ois.readObject();
        }
    }
}
