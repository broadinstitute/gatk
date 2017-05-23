package org.broadinstitute.hellbender.tools.picard.vcf;

import com.beust.jcommander.internal.Lists;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

public final class GatherVcfsUnitTest {
    //tests for GatherVcfs, specifically areAllBlockCompressed()

    @DataProvider(name = "files")
    public Object[][] files(){
        final Path nullPath = null;

        return new Object[][]{
                {Lists.newArrayList(nullPath), false},
                {Lists.newArrayList(IOUtils.getPath("test.vcf")), false},
                {Lists.newArrayList(IOUtils.getPath("test.bcf")), false},
                {Lists.newArrayList(IOUtils.getPath("test.vcf.gz")), true},
                {Lists.newArrayList(IOUtils.getPath("test1.vcf"), IOUtils.getPath("test2.vcf.gz")), false}
        };
    }

    @Test(dataProvider = "files")
    public void testareAllBlockCompressed(List<Path> paths, boolean expected) throws IOException {
        Assert.assertEquals(GatherVcfs.areAllBlockCompressed(paths), expected);
    }

}
