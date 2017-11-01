package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.SAMSequenceDictionary;
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

    // TODO: 10/6/17 requires a ReferenceMultiSource and SAMSequenceDictionary at the same time because the 2bit reference gives a scrambled reference contig order (see #2037)
    void inferSvAndWriteVCF(final String vcfOutputFileName, final String sampleId, final JavaRDD<AlignedContig> contigs,
                            final Broadcast<ReferenceMultiSource> broadcastReference,
                            final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                            final Logger toolLogger);
}
