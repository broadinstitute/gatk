package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import org.apache.spark.api.java.JavaRDD;

/**
 * Loads various upstream assembly and alignment formats and turn into custom {@link AlignedContig} format in the discovery stage.
 * Implementations are expected to apply filter to the provided raw format contig alignments and
 * break the gapped alignments by calling into {@link ContigAlignmentsModifier#splitGappedAlignment(AlignmentInterval, int, int)}.
 */
public abstract class AlignedContigGenerator {

    public abstract JavaRDD<AlignedContig> getAlignedContigs();
}
