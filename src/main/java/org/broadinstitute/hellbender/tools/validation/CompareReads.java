package org.broadinstitute.hellbender.tools.validation;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.util.*;

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Compares the base qualities, cigars, alignment information, and samflags of reads between two SAM/BAM/CRAM files." +
                " The reads in the two files must have exactly the same names and appear in the same order.",
        oneLineSummary = "Compares the base qualities of two SAM/BAM/CRAM files",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
public class CompareReads extends GATKTool {
    @Argument(doc = "The first sam file against which to compare equality",
            shortName = "I1", fullName = "input1", optional = false)
    protected String input1;

    @Argument(doc = "The second sam file against which to compare equality",
            shortName = "I2", fullName = "input2", optional = false)
    protected String input2;


    @Override
    public void traverse() {
        List<String> errorMessages = new ArrayList<>();

        try(ReadsDataSource reads1 = new ReadsDataSource(IOUtils.getPath(input1));
            ReadsDataSource reads2 = new ReadsDataSource(IOUtils.getPath(input2));) {
            final Iterator<GATKRead> it1 = reads1.iterator();
            final Iterator<GATKRead> it2 = reads2.iterator();
            while (it1.hasNext() && it2.hasNext()) {
                final SAMRecord read1 = it1.next().convertToSAMRecord(reads1.getHeader());
                final SAMRecord read2 = it2.next().convertToSAMRecord(reads2.getHeader());
                final String eqMessage = readsEqualAllowAddingAttributes(read1, read2);
                if (eqMessage != null){
                    errorMessages.add(eqMessage);
                }
            }
            if (it1.hasNext() || it2.hasNext()) {
                //at least one has no more records (because the while loop is done) and at least one does have records. So we're not equal.
                throw new GATKException("Not the same number of reads");
            }
        }
        if (!errorMessages.isEmpty()) {
            System.out.println("Here are the errors generated for this file:");
            for (String e: errorMessages) {
                System.out.println("\n"+e);
            }
        }
    }

    /**
     * Compares the reads but ignores order of attributes.
     * Also allows actualRead to have a superset of attributes of expectedRead.
     */
    private static String readsEqualAllowAddingAttributes(final SAMRecord actualRead, final SAMRecord expectedRead) {
        final String actualName = actualRead.getReadName();
        final String expectedName = expectedRead.getReadName();

        String msg;

        msg = compareValues(actualName, expectedName, "name");
        if (msg != null){ return msg; }

        final String readNames = "actualName:" + actualName + " expectedName:" + expectedName;

        msg = compareValues(SAMFlag.getFlags(actualRead.getFlags()), SAMFlag.getFlags(expectedRead.getFlags()), readNames + " getFlags");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getInferredInsertSize(), expectedRead.getInferredInsertSize(), readNames + " getInferredInsertSize");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getMappingQuality(), expectedRead.getMappingQuality(), readNames + " getMappingQuality");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getMateReferenceIndex(), expectedRead.getMateReferenceIndex(), readNames + "getMateReferenceIndex");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getMateAlignmentStart(), expectedRead.getMateAlignmentStart(), readNames + "getMateAlignmentStart");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getReferenceIndex(), expectedRead.getReferenceIndex(), readNames + " getReferenceIndex");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getAlignmentStart(), expectedRead.getAlignmentStart(), readNames + " getAlignmentStart");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getCigar(), expectedRead.getCigar(), readNames + " getCigar");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getReferenceName(), expectedRead.getReferenceName(), readNames + " getReferenceName");
        if (msg != null){ return msg; }

        msg = compareValues(actualRead.getMateReferenceName(), expectedRead.getMateReferenceName(), readNames + " getMateReferenceName");
        if (msg != null){ return msg; }

        if (!Arrays.equals(actualRead.getReadBases(), expectedRead.getReadBases())){
            return "getReadBases different actualRead:" + actualName + " expectedRead:" + expectedName + " (" + Arrays.toString(actualRead.getReadBases()) + " vs " + Arrays.toString(expectedRead.getReadBases()) + ")";
        }
        if (!Arrays.equals(actualRead.getBaseQualities(), expectedRead.getBaseQualities())){
            return "getBaseQualities different actualRead:" + actualName + " expectedRead:" + expectedName + " (" + Arrays.toString(actualRead.getBaseQualities()) + " vs " + Arrays.toString(expectedRead.getBaseQualities()) + ")";
        }
        return null;
    }

    private static <T> String compareValues(final T v1, final T v2, final String label) {
        boolean eq = Objects.equals(v1, v2);
        if (eq) {
            return null;
        } else {
            final String s1 = String.valueOf(v1);
            final String s2 = String.valueOf(v2);
            return label + " differs. File 1: " + s1 + " File 2: " + s2;
        }
    }

}
