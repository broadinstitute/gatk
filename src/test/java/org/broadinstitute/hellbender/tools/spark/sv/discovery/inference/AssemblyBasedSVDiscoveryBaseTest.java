package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import org.broadinstitute.hellbender.GATKBaseTest;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera;

/**
 * Provides a common interface to get test data and expected results.
 *
 * Intends to extended by unit tests testing methods
 * in the assembly-based SV discovery stage.
 */
abstract class AssemblyBasedSVDiscoveryBaseTest extends GATKBaseTest {

    protected final AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV forSimpleSV = new AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV();
    protected final AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints forInversionBreakpoints = new AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints();
    protected final AssemblyBasedSVDiscoveryTestDataProviderForBreakEndVariants forBreakEndVariants = new AssemblyBasedSVDiscoveryTestDataProviderForBreakEndVariants();

    AssemblyBasedSVDiscoveryBaseTest() {
    }

    protected List<AssemblyBasedSVDiscoveryTestDataForSimpleChimera> getAllTestData() {
        final List<AssemblyBasedSVDiscoveryTestDataForSimpleChimera> data = new ArrayList<>(100);

        data.addAll( forSimpleSV.getAllTestData() );
        data.addAll( forInversionBreakpoints.getAllTestData() );
        data.addAll( forBreakEndVariants.getAllTestData() );

        return data;
    }

    protected List<Tuple2<? extends AssemblyBasedSVDiscoveryTestDataForSimpleChimera, ? extends AssemblyBasedSVDiscoveryTestDataForSimpleChimera>> getAllTestDataInPairs() {
        final List<Tuple2<? extends AssemblyBasedSVDiscoveryTestDataForSimpleChimera, ? extends AssemblyBasedSVDiscoveryTestDataForSimpleChimera>> data = new ArrayList<>(50);

        data.addAll( forSimpleSV.getAllTestDataPaired() );
        data.addAll( forInversionBreakpoints.getAllTestDataPaired() );
        data.addAll( forBreakEndVariants.getAllTestDataPaired() );

        return data;
    }
}
