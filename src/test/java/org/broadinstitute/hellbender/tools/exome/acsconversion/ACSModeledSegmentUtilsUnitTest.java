package org.broadinstitute.hellbender.tools.exome.acsconversion;

import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.exome.Genome;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionSimulatedData;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by lichtens on 6/16/16.
 */
public class ACSModeledSegmentUtilsUnitTest extends BaseTest {
    static final String TEST_FILE_PATH= "src/test/resources/org/broadinstitute/hellbender/tools/exome/allelicbalancecaller/cell_line_small-sim-final.seg";

    @Test
    public void testConversion() {
        final List<ACNVModeledSegment> segs = SegmentUtils.readACNVModeledSegmentFile(new File(TEST_FILE_PATH));
        final Genome genome = new Genome(AlleleFractionSimulatedData.TRIVIAL_TARGETS, Collections.emptyList());
        final List<ACSModeledSegment> acsSegs = segs.stream().map(seg -> ACSModeledSegmentUtils.convertACNVSegmentToACSSegment(seg, 2.0, genome, true)).collect(Collectors.toList());

        for (int i = 0; i < segs.size(); i ++) {
            Assert.assertEquals(acsSegs.get(i).getTau()/2.0, segs.get(i).getSegmentMeanInCRSpace(), 1e-10);
        }
    }
}
