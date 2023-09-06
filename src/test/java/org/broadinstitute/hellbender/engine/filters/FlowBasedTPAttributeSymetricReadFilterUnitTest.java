package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.filters.flow.FlowBasedTPAttributeSymetricReadFilter;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public final class FlowBasedTPAttributeSymetricReadFilterUnitTest {

    static final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10);
    static final ReadFilter readFilter = new FlowBasedTPAttributeSymetricReadFilter();

    @DataProvider(name = "reads")
    public Object[][] readsStartEnd(){

        final GATKRead goodRead = ArtificialReadUtils.createArtificialRead(header, "AAACCC".getBytes(), "!!!!!!".getBytes(), "6M");
        goodRead.setAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME, "XYXABA".getBytes(StandardCharsets.UTF_8));

        final GATKRead badRead = ArtificialReadUtils.createArtificialRead(header, "AAACCC".getBytes(), "!!!!!!".getBytes(), "6M");
        badRead.setAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME, "XYXABB".getBytes(StandardCharsets.UTF_8));

        return new Object[][]{
                { goodRead, true },
                { badRead, false }
        };
    }


    @Test(dataProvider = "reads")
    public void testFilter(GATKRead read, boolean passesFilter) {
        Assert.assertEquals(readFilter.test(read), passesFilter);
    }
}