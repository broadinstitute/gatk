package org.broadinstitute.hellbender.tools.walkers.fasta;

import com.google.common.primitives.Bytes;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.FastaReferenceWriter;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;

/**
 * Create a fasta with the bases shifted by offset
 *
 * delta1 = offset - 1
 * delta2 = total - delta1
 *
 * To shift forward:
 * if you are given a position in the regular fasta (pos_r) and want the position in the shifted fasta (pos_s):
 * if pos_r > delta1 => pos_s = pos_r - delta1  ==   pos_r - offset +1
 *   otherwise          pos_s = pos_r + delta2  ==   pos_r + total - offset + 1
 *
 * To shift back:
 * if you are given a position in the shifted fasta (pos_s) and want the position in the regular fasta (pos_r):
 * if pos_s > delta2 => pos_r = pos_s - delta2  ==   pos_s - total + offset - 1
 *   otherwise          pos_r = pos_s + delta1  ==   pos_s + offset - 1
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Create a new fasta shifted by the input amount and a shift_back file based on the intervals specified",
        oneLineSummary = "Creates a shifted fasta file and shift_back file",
        programGroup = ReferenceProgramGroup.class
)
public class ShiftFasta extends GATKTool {
    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME,
            shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME,
            doc = "path to reference file, .fai should also exist in that location.")
    private String referencePath;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Path to write the output fasta to")
    protected String output;

    public static final String SHIFT_BACK_OUTPUT = "shift-back-output";
    @Argument(fullName = SHIFT_BACK_OUTPUT,
            doc = "Path to write the shift_back file to")
    protected String shiftBackOutput;

    public static final String SHIFT_OFFSET = "shift-amount";
    @Argument(fullName = SHIFT_OFFSET, doc="Number of bases to shift the reference by. For example, if 300 is specified, the new fasta will start at the 300th base")
    private int shiftOffset;

    public static final String LIFTOVER_INTERVAL = "liftover-interval";
    @Argument(fullName = LIFTOVER_INTERVAL, doc="Interval in the coordinates of the original fasta, that you plan to replace with data from the shift fasta. For example, if the reference is 1000 bases, -100:200 will make a shift back file from ")

    public static final String LINE_WIDTH_LONG_NAME = "line-width";
    @Argument(fullName= LINE_WIDTH_LONG_NAME, doc="Maximum length of sequence to write per line", optional=true)
    public int basesPerLine = FastaReferenceWriter.DEFAULT_BASES_PER_LINE;

    ReferenceDataSource refSource;
    FastaReferenceWriter refWriter;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        refSource = referenceArguments.getReferencePath() != null ? ReferenceDataSource.of(referenceArguments.getReferencePath()) : null;
        final Path path = IOUtils.getPath(output);
        try {
            refWriter = new FastaReferenceWriterBuilder()
                    .setFastaFile(path)
                    .setBasesPerLine(basesPerLine)
                    .build();
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Couldn't create " + output + ", encountered exception: " + e.getMessage(), e);
        }
    }

    public void traverse() {
        SAMSequenceDictionary refDict = refSource.getSequenceDictionary();
        long refLengthLong = refDict.getReferenceLength();
        if (refLengthLong > Integer.MAX_VALUE) {
            // TODO fix this??
            throw new UserException.BadInput("Reference length is too long");
        }
        int splitIndex = ((int) refLengthLong) - shiftOffset + 1;
        // TODO make this not mito specific
        byte[] bases = refSource.queryAndPrefetch(new SimpleInterval("chrM", 0, (int) refDict.getReferenceLength())).getBases();
        byte[] origEnd = Arrays.copyOfRange(bases, splitIndex, bases.length);
        byte[] origBegin = Arrays.copyOf(bases, splitIndex -1);
        try {
        refWriter.appendBases(origEnd).appendBases(origBegin);
        } catch (IOException e) {
            throw new UserException("Failed to write fasta due to " + e.getMessage(), e);
        }
    }

    @Override
    public Object onTraversalSuccess(){
        // TODO write shift back file
        return null;
    }

    @Override
    public void closeTool() {
        super.closeTool();
        try{
            if( refWriter != null ) {
                refWriter.close();
            }
        } catch (IllegalStateException e){
            //sink this
        } catch (IOException e) {
            throw new UserException("Failed to write fasta due to " + e.getMessage(), e);
        }
    }

    }
