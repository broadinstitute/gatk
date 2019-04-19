package org.broadinstitute.hellbender.tools.walkers.pacbio;

import htsjdk.samtools.*;
import ngs.ReadGroup;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.util.List;

/**
 * Restore annotations from unaligned BAM files that are discarded during conversion by samtools fastq and sorts output
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A unaligned BAM file</li>
 *     <li>A aligned BAM file</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 *     <li>An aligned BAM file with annotations restored</li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <h4>Restore annotations from unaligned BAM files that are discarded during conversion by samtools fastq and sorts output</h4>
 * <pre>
 *   gatk RepairPacBioBam \
 *     -I unaligned.bam \
 *     -A aligned.bam \
 *     -O restored.bam
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Restore annotations from unaligned BAM files that are discarded during conversion by samtools fastq and sorts output",
        oneLineSummary = "Restore annotations from unaligned BAM files that are discarded during conversion by samtools fastq and sorts output",
        programGroup = ReadDataManipulationProgramGroup.class
)
public final class RepairPacBioBam extends ReadWalker {
    @Argument(fullName = "aligned",
            shortName = "A",
            doc="aligned reads")
    public String aligned;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    public String output;

    @Argument(fullName = "sample_name",
            shortName = "S",
            doc="Sample name",
            optional = true)
    public String sample_name;

    private SAMFileWriter writer;
    private SAMRecordIterator it;
    private SAMRecord currentAlignedRead;

    @Override
    public void onTraversalStart() {
        if ( readArguments.getReadFiles().size() != 1 ) {
            throw new UserException("Specify a single unaligned BAM file to -I and a single aligned BAM file to -A");
        }

        SamReaderFactory srf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        SamReader srs = srf.open(IOUtils.getPath(aligned));

        if (getHeaderForReads().getReadGroups().size() != 1) {
            throw new UserException("Only one read group per aligned/unaligned BAM file permitted");
        }

        /*
        if (!srs.getFileHeader().getReadGroups().get(0).getId().equals(getHeaderForReads().getReadGroups().get(0).getId())) {
            throw new UserException("Read groups between unaligned and aligned BAM files are not identical ('" +
                    srs.getFileHeader().getReadGroups().get(0).getId() +
                    "' vs '" +
                    getHeaderForReads().getReadGroups().get(0).getId() +
                    "')");
        }

        if (srs.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.queryname) {
            throw new UserException("Aligned data must be unsorted");
        }
        */

        it = srs.iterator();
        if (!it.hasNext()) {
            throw new UserException("Aligned BAM file is empty");
        }

        currentAlignedRead = it.next();

        writer = createWriter(output, getHeaderForReads(), srs.getFileHeader().getSequenceDictionary());
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        SAMRecord unalignedRead = read.convertToSAMRecord(this.getHeaderForReads());

        do {
            if (currentAlignedRead != null) {
                if (!unalignedRead.getReadName().equals(currentAlignedRead.getReadName())) {
                    throw new UserException(
                            "Expected aligned read to have name '" +
                            unalignedRead.getReadName() +
                            "', but saw '" +
                            currentAlignedRead.getReadName() +
                            "'.  Are these unaligned and aligned versions of exactly the same data?"
                    );
                } else {
                    List<SAMRecord.SAMTagAndValue> attrs = unalignedRead.getAttributes();
                    attrs.addAll(currentAlignedRead.getAttributes());

                    for (SAMRecord.SAMTagAndValue tv : attrs) {
                        currentAlignedRead.setAttribute(tv.tag, tv.value);
                    }

                    writer.addAlignment(currentAlignedRead);
                }
            }

            currentAlignedRead = it.hasNext() ? it.next() : null;
        } while (currentAlignedRead != null && unalignedRead.getReadName().equals(currentAlignedRead.getReadName()));
    }

    @Override
    public void closeTool() {
        if ( writer != null ) {
            writer.close();
        }
    }

    /**
     * Creates SAMFileWriter
     * @return A SAMFileWriter.
     */
    private SAMFileWriter createWriter(String out, SAMFileHeader sfh, SAMSequenceDictionary ssd) {
        sfh.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        sfh.setSequenceDictionary(ssd);

        if (sample_name != null) {
            for (SAMReadGroupRecord rg : sfh.getReadGroups()) {
                rg.setSample(sample_name);
            }
        }

        return new SAMFileWriterFactory().makeSAMOrBAMWriter(sfh, false, new File(out));
    }
}
