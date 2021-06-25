package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.util.*;

@CommandLineProgramProperties(
        summary = "Unmaps reads that have distant mates.  "+
                "This ensures that a PairWalker will always see both mates, "+
                "even when one of them is mapped far away, when given the output "+
                "of this tool along with the original inputs.",
        oneLineSummary = "Unmaps reads with distant mates.",
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
public class PrintDistantMates extends ReadWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    public GATKPath output;
    public static String DISTANT_MATE_TAG = "DM";

    private SAMFileGATKReadWriter outputWriter;

    @Override public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> readFilters = new ArrayList<>(super.getDefaultReadFilters());
        readFilters.add(ReadFilterLibrary.PAIRED);
        readFilters.add(ReadFilterLibrary.PRIMARY_LINE);
        readFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        readFilters.add(ReadFilterLibrary.MATE_DISTANT);
        return readFilters;
    }

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(output, false);
    }

    @Override
    public void apply( final GATKRead read,
                       final ReferenceContext referenceContext,
                       final FeatureContext featureContext ) {
        final GATKRead copy = doDistantMateAlterations(read);
        outputWriter.addRead(copy);
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }

    public static boolean isDistantMate( final GATKRead read ) {
        return read.hasAttribute(DISTANT_MATE_TAG);
    }

    public static GATKRead doDistantMateAlterations( final GATKRead read ) {
        final GATKRead copy = read.copy();
        final Integer nmValue = read.getAttributeAsInteger(SAMTag.NM.name());
        final String oaValue = read.getContig() + "," + read.getStart() + "," +
                (read.isReverseStrand() ? "-," : "+,") + read.getCigar() + "," +
                read.getMappingQuality() + "," + (nmValue == null ? "" : nmValue) + ";";
        copy.clearAttribute(SAMTag.NM.name());
        copy.setAttribute(SAMTag.OA.name(), oaValue);
        copy.setAttribute(DISTANT_MATE_TAG, "");
        copy.setPosition(read.getMateContig(), read.getMateStart());
        copy.setCigar(SAMRecord.NO_ALIGNMENT_CIGAR);
        copy.setMappingQuality(0);
        copy.setIsUnmapped();
        return copy;
    }

    public static GATKRead undoDistantMateAlterations( final GATKRead read ) {
        final String oaValue = read.getAttributeAsString(SAMTag.OA.name());
        if ( oaValue == null ) {
            return read;
        }
        final GATKRead copy = read.copy();
        copy.clearAttribute(DISTANT_MATE_TAG);
        copy.clearAttribute(SAMTag.OA.name());
        try {
            final String[] tokens = oaValue.split(",");
            copy.setPosition(tokens[0], Integer.parseInt(tokens[1]));
            copy.setIsReverseStrand("-".equals(tokens[2]));
            copy.setCigar(tokens[3]);
            copy.setMappingQuality(Integer.parseInt(tokens[4]));
            if ( tokens[5].length() > 1 ) {
                final String nmValue = tokens[5].substring(0, tokens[5].length() - 1);
                copy.setAttribute(SAMTag.NM.name(), Integer.parseInt(nmValue));
            }
        } catch ( final Exception exc ) {
            throw new UserException("can't recover alignment from OA tag: " + oaValue, exc);
        }
        return copy;
    }
}
