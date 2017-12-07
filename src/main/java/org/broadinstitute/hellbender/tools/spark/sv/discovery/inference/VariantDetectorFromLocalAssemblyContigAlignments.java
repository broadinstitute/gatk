package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryInputData;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;

import java.util.Collections;
import java.util.List;

/**
 * A base class for workflow of variant breakpoint detection from split alignments, variant type interpretation,
 * and resolving complication {@link BreakpointComplications}.
 */
interface VariantDetectorFromLocalAssemblyContigAlignments {

    @SuppressWarnings("unchecked")
    List<String> EMPTY_INSERTION_MAPPINGS = Collections.EMPTY_LIST;

    void inferSvAndWriteVCF(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                            final SvDiscoveryInputData svDiscoveryInputData);
}
