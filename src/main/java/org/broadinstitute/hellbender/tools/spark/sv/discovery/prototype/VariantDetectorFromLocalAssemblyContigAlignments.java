package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;

/**
 * A base class for workflow of variant breakpoint detection from split alignments, variant type interpretation,
 * and resolving complication {@link org.broadinstitute.hellbender.tools.spark.sv.discovery.BreakpointComplications}.
 */
interface VariantDetectorFromLocalAssemblyContigAlignments {

    void inferSvAndWriteVCF(final JavaRDD<AlignedContig> localAssemblyContigs, final String vcfOutputFileName,
                            final Broadcast<ReferenceMultiSource> broadcastReference, final String fastaReference,
                            final Logger toolLogger);
}
