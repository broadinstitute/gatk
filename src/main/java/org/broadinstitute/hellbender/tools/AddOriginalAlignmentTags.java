package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMTag;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.argparser.WorkflowProperties;
import org.broadinstitute.barclay.argparser.WorkflowOutput;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

@CommandLineProgramProperties(
        summary = "Adds Original Alignment tag and original mate contig tag",
        oneLineSummary = "Adds Original Alignment tag and original mate contig tag",
        programGroup = ReadDataManipulationProgramGroup.class
)
@ExperimentalFeature
@WorkflowProperties
public class AddOriginalAlignmentTags extends ReadWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    @WorkflowOutput(optionalCompanions = {StandardArgumentDefinitions.OUTPUT_INDEX_COMPANION})
    public GATKPath output;
    private SAMFileGATKReadWriter outputWriter;

    public final static String MATE_CONTIG_TAG_NAME = "XM";
    public final static String OA_TAG_NAME = "OA";
    public final static String OA_SEPARATOR = ",";

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(output, true);
    }
    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        addOATag(read);
        addMateContigTag(read);
        outputWriter.addRead(read);
    }
    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }

    private static void addMateContigTag(final GATKRead read) {
        read.setAttribute(MATE_CONTIG_TAG_NAME, !read.mateIsUnmapped() ? read.getMateContig() : "*");
    }

    //TODO: Once the OA tag is merged into the spec (https://github.com/samtools/hts-specs/pull/193) this should be moved to htsjdk
    private static void addOATag(final GATKRead read) {
        String OAValue;
        if(!read.isUnmapped()){
            OAValue = String.format("%s,%s,%s,%s,%s,%s;",
                    read.getContig().replace(OA_SEPARATOR, "_"),
                    read.getStart(),
                    read.isReverseStrand() ? "-" : "+",
                    read.getCigar().toString(),
                    read.getMappingQuality(),
                    read.getAttributeAsString(SAMTag.NM.name()));
        } else {
            OAValue = "*,0,*,*,0,0;";
        }

        read.setAttribute(OA_TAG_NAME, OAValue);
    }

    public static String getOAContig (final GATKRead read) {
        return(read.getAttributeAsString(OA_TAG_NAME).split(OA_SEPARATOR)[0]);
    }

}
