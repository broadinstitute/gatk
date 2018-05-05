package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.annotations.Test;

import java.util.Collections;

public class BreakpointsInferenceUnitTest {

    @Test(expectedExceptions = GATKException.class)
    public void testGetBreakpoints_ExpectException() {
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 100001, 100100), 1 ,100, TextCigarCodec.decode("100M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 100101, 100200), 101 ,200, TextCigarCodec.decode("100M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final SimpleChimera simpleChimera = new SimpleChimera(region1, region2, Collections.emptyList(), "1",
                AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME, SimpleSVDiscoveryTestDataProvider.b37_seqDict);
        simpleChimera.inferType();
    }
}
