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
        File nullFile = null;


        return new Object[][]{
                {Lists.newArrayList(nullFile), false},
                {Lists.newArrayList("test.vcf"), false},
                {Lists.newArrayList("test.bcf"), false},
                {Lists.newArrayList("test.vcf.gz"), true},
                {Lists.newArrayList("test1.vcf", "test2.vcf.gz"), false}
        };
    }

    @Test(dataProvider = "files")
    public void testareAllBlockCompressed(List<String> files, boolean expected) throws IOException {
        List<Path> paths = files.stream().map(IOUtils::getPath).collect(Collectors.toList());
        Assert.assertEquals(GatherVcfs.areAllBlockCompressed(paths), expected);
    }

}
