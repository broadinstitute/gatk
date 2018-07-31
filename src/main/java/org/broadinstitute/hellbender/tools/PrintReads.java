package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

/**
 * Write reads from SAM format file (SAM/BAM/CRAM) that pass criteria to a new file.
 *
 * <p>
 * A common use case is to subset reads by genomic interval using the -L argument.
 * Note when applying genomic intervals, the tool is literal and does not retain mates of paired-end reads outside of the interval, if any.
 * Data with missing mates will fail ValidateSamFile validation with MATE_NOT_FOUND, but certain tools may still analyze the data.
 * If needed, to rescue such mates, use either FilterSamReads or ExtractOriginalAlignmentRecordsByNameSpark.
 * </p>
 *
 * <p>
 * By default, PrintReads applies the WellformedReadFilter at the engine level.
 * What this means is that the tool does not print reads that fail the WellformedReadFilter filter.
 * You can similarly apply other engine-level filters to remove specific types of reads with the --read-filter argument.
 * See documentation category 'Read Filters' for a list of available filters.
 * To keep reads that do not pass the WellformedReadFilter, either disable the filter with --disable-read-filter or disable all default filters with --disable-tool-default-read-filters.
 * </p>
 *
 * <p>
 * The reference is strictly required when handling CRAM files.
 * </p>
 *
 * <h3> Input </h3>
 * <ul>
 *     <li> Coordinate-sorted and indexed SAM/BAM/CRAM </li>
 * </ul>
 *
 * <h3> Output </h3>
 * <ul>
 *     <li> Coordinate-sorted and indexed SAM/BAM/CRAM </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 * Print reads that pass the WellformedReadFilter and that map to chromosome 20.
 * <pre>
 * gatk PrintReads \
 *   -I input.bam \
 *   -L 20 \
 *   -O chr20.bam
 * </pre>
 * Print reads that pass the WellformedReadFilter and the NotDuplicateReadFilter.
 * <pre>
 * gatk PrintReads \
 *   -I input.bam \
 *   --read-filter NotDuplicateReadFilter \
 *   -O filtered.bam
 * </pre>
 *
 * Print reads that pass the WellformedReadFilter and that map to chromosome 20 from a BAM in a google cloud bucket to another cloud bucket.
 * <pre>
 * gatk PrintReads \
 *   -I gs://cloud-bucket/input.bam \
 *   -L 20 \
 *   -O gs://my-gcs-bucket/chr20.bam
 * </pre>
 * Print reads that pass the WellformedReadFilter and that map to chromosome 20 from a BAM in a google cloud bucket to the local system.
 * <pre>
 * gatk PrintReads \
 *   -I gs://cloud-bucket/input.bam \
 *   -L 20 \
 *   -O chr20.bam
 * </pre>
 *
 * {@GATK.walkertype ReadWalker}
 */
@CommandLineProgramProperties(
        summary = "Prints reads from the input SAM/BAM/CRAM file to the SAM/BAM/CRAM file.",
        oneLineSummary = "Print reads in the SAM/BAM/CRAM file",
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
public final class PrintReads extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    public String output;
    private SAMFileGATKReadWriter outputWriter;

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(IOUtils.getPath(output), true);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        outputWriter.addRead(read);
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
