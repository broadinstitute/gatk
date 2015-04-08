package org.broadinstitute.hellbender.utils.codecs.table;

import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public final class TableCodecUnitTest extends BaseTest {

    @DataProvider(name = "badNames")
    public Object[][] badNames() {
        List<Object[]> params = new ArrayList<>();
        params.add(new String[]{"a.tsv"});
        params.add(new String[]{"a.table.gz"});
        params.add(new String[]{"a.bed"});
        params.add(new String[]{"a.bcf"});
        params.add(new String[]{"a.hapmap"});
        params.add(new String[]{"a.refseq"});
        params.add(new String[]{"a.beagle"});
        return params.toArray(new Object[][]{});
    }

    @Test(dataProvider = "badNames")
    public void testBadNames(String badName){
        TableCodec tc = new TableCodec();
        Assert.assertFalse(tc.canDecode(badName), badName);
    }

    @DataProvider(name = "goodNames")
    public Object[][] goodNames() {
        List<Object[]> params = new ArrayList<>();
        params.add(new String[]{"a.table"});
        return params.toArray(new Object[][]{});
    }

    @Test(dataProvider = "goodNames")
    public void testGoodNames(String goodName){
        TableCodec tc = new TableCodec();
        Assert.assertTrue(tc.canDecode(goodName), goodName);
    }

    @Test
    public void testDecode(){
        TableCodec tc = new TableCodec();
        GenomeLocParser glp= hg19GenomeLocParser;
        final TableFeature decode = tc.decode("");

    }
}
