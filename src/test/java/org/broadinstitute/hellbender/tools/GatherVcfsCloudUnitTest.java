package org.broadinstitute.hellbender.tools;

import com.beust.jcommander.internal.Lists;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Path;
import java.util.*;

public final class GatherVcfsCloudUnitTest {
    //tests for GatherVcfsCloud, specifically areAllBlockCompressed()

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
        Assert.assertEquals(GatherVcfsCloud.areAllBlockCompressed(paths), expected);
    }

}
