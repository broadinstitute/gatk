package org.broadinstitute.hellbender.tools.walkers.fasta;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaReferenceWriter;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.ListIterator;

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
 * Example command line:
 * ShiftFasta
 * -R "<CIRCURLAR_REFERENCE.fasta>"         // the reference to shift
 * -O "<SHIFTED_REFERENCE.fasta>"           // output; the shifted fasta
 * --shift-back-output "<SHIFT_BACK.chain>" // output; the shiftback chain file to use when lifting over
 * --shift-offset-list "<SHIFT_OFFSETS>"    // optional; Specifies the offset to shift for each contig in the reference. If not specified, the offset will be half the length of the contig.
 * --interval-file-name "<SHIFT_INTERVALS>" // output; base name for output interval files (.intervals and .shifted.intervals) that should be used when calling variants against the unshifted and shifted reference.
 * --line-width 100
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        summary = "Create a new fasta starting at the shift-offset +1 position and a shift_back chain file that can be used with the Liftover tool. It will shift all contigs by default.",
        oneLineSummary = "Creates a shifted fasta file and shift_back file",
        programGroup = ReferenceProgramGroup.class
)
public class ShiftFasta extends GATKTool {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Path to write the output fasta to")
    protected String output;

    public static final String SHIFT_BACK_OUTPUT = "shift-back-output";
    @Argument(fullName = SHIFT_BACK_OUTPUT,
            doc = "Path to write the shift_back file to")
    protected String shiftBackOutput;

    public static final String SHIFT_OFFSET_LIST = "shift-offset-list";
    @Argument(fullName = SHIFT_OFFSET_LIST,
            doc="Number of bases to skip in the reference before starting the shifted reference. If the reference contains multiple contigs, a value should be specified for each contig. " +
                "For example, if 300 is specified, the new fasta will start at the 301th base (count starting at 1)." +
                "If not specified, each contig will be shifted by half the number of bases. To skip the shifting of a contig, specify 0 in the list.", optional = true)
    private List<Integer> shiftOffsets = null;

    public static final String INTERAL_FILE_NAME = "interval-file-name";
    @Argument(fullName = INTERAL_FILE_NAME,
            doc="Base name for interval files. Intervals will be midway between beginning and computed offset. If not specified, no interval files will be written.", optional = true)
    private String intervalFilename;

    public static final String LINE_WIDTH_LONG_NAME = "line-width";
    @Argument(fullName= LINE_WIDTH_LONG_NAME, doc="Maximum length of sequence to write per line", optional=true)
    public int basesPerLine = FastaReferenceWriter.DEFAULT_BASES_PER_LINE;

    private ReferenceDataSource refSource;
    private FastaReferenceWriter refWriter;
    private FileWriter chainFileWriter;
    private FileWriter intervalRegularWriter;
    private FileWriter intervalShiftedWriter;

    private int chainId = 0;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        refSource = referenceArguments.getReferencePath() != null ? ReferenceDataSource.of(referenceArguments.getReferencePath()) : null;
        final Path path = IOUtils.getPath(output);
        chainId = 1;
        try {
            refWriter = new FastaReferenceWriterBuilder()
                    .setFastaFile(path)
                    .setBasesPerLine(basesPerLine)
                    .build();
            chainFileWriter = new FileWriter(shiftBackOutput);
            if (intervalFilename != null) {
                intervalRegularWriter = new FileWriter(intervalFilename+ ".intervals");
                intervalShiftedWriter = new FileWriter(intervalFilename + ".shifted.intervals");
            }
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
        List<SAMSequenceRecord> contigs = refSource.getSequenceDictionary().getSequences();
        if (shiftOffsets != null && !shiftOffsets.isEmpty() && shiftOffsets.size() != contigs.size()) {
            throw new UserException.BadInput("Shift offset list size " + shiftOffsets.size() + " must equal number of contigs in the reference " + contigs.size());
        }
        final ListIterator<Integer> shiftOffsetsIt = shiftOffsets != null && !shiftOffsets.isEmpty()  ? shiftOffsets.listIterator() : null;
        refSource.getSequenceDictionary().getSequences().forEach(seq -> shiftContig(seq, shiftOffsetsIt));
    }

    /**
     * This method adds to a new fasta ref file that has been shifted by the amount indicated in the shiftOffsetsIt.
     * This also adds to the supporting files: chainfile, interval list for both the shifted and unshifted fasta files
     * The shift is all done in memory. This is a scaling limitation.
     * @param seq The contig or sequence within the fasta file
     * @param shiftOffsetsIt the iterator at the correct position to get the next offset or null if dividing contig by 2
     */
    protected final void shiftContig(SAMSequenceRecord seq, ListIterator<Integer> shiftOffsetsIt) {
        final int contigLength = seq.getSequenceLength();
        final String seqName = seq.getSequenceName();
        int shiftOffset = shiftOffsetsIt == null ? contigLength/2 : shiftOffsetsIt.next();
        if (shiftOffset > 0 && shiftOffset < contigLength) {
            byte[] bases = refSource.queryAndPrefetch(new SimpleInterval(seqName, 1, contigLength)).getBases();
            byte[] basesAtEnd = Arrays.copyOfRange(bases, shiftOffset, bases.length);
            byte[] basesAtStart = Arrays.copyOf(bases, shiftOffset);
            int shiftBackOffset = bases.length - shiftOffset;

            addToShiftedReference(refWriter, seqName, basesPerLine, basesAtStart, basesAtEnd);
            addToChainFile(seqName, contigLength, shiftOffset, bases, shiftBackOffset);
            if (intervalFilename != null && shiftOffsetsIt == null) {
                addToIntervalFiles(intervalRegularWriter, intervalShiftedWriter, seqName, shiftOffset, contigLength);
            }
        } else {
            logger.info("not shifting config " + seq.getContig() + " because shift offset " + shiftOffset + " is not between 1-" + contigLength );
        }
    }

    private void addToIntervalFiles(FileWriter intervalRegularWriter, FileWriter intervalShiftedWriter, String seqName, int shiftOffset, int contigLength) {
        try {
            int intervalStart = shiftOffset/2;
            int intervalEnd = intervalStart + contigLength/2 - 1;
            int shiftedIntervalStart = intervalStart;
            int shiftedIntervalEnd = intervalEnd + contigLength % 2;
            intervalRegularWriter.append(seqName + ":" + intervalStart + "-" + intervalEnd + "\n");
            intervalShiftedWriter.append(seqName + ":" + shiftedIntervalStart + "-" + shiftedIntervalEnd + "\n");
        } catch (IOException e) {
            throw new UserException("Failed to write interval files due to " + e.getMessage(), e);
        }

    }

    private void addToShiftedReference(FastaReferenceWriter refWriter, String seqName, int basesPerLine, byte[] basesAtStart, byte[] basesAtEnd) {
        try {
            refWriter.startSequence(seqName, basesPerLine);
            // swap the bases
            refWriter.appendBases(basesAtEnd).appendBases(basesAtStart);
        } catch (IOException e) {
            throw new UserException("Failed to write shifted reference due to " + e.getMessage(), e);
        }
    }

    private void addToChainFile(String seqName, int contigLength, int shiftOffset, byte[] bases, int shiftBackOffset) {
        try {
            chainFileWriter.append(createChainString(seqName, shiftBackOffset, contigLength, shiftOffset, bases.length, 0, shiftBackOffset, chainId++));
            chainFileWriter.append("\n" + shiftBackOffset + "\n\n");
            chainFileWriter.append(createChainString(seqName, shiftOffset - 1, contigLength, 0, shiftOffset, shiftBackOffset, bases.length, chainId++));
            chainFileWriter.append("\n" + shiftOffset + "\n\n");
        } catch (IOException e) {
            throw new UserException("Failed to write chainFile due to " + e.getMessage(), e);
        }
    }

    private String createChainString(String name, int score, int length, int start, int end, int shiftBackStart, int shiftBackEnd, int id) {
        String[] items = new String[] { "chain",
                Integer.toString(score),
                name,
                Integer.toString(length),
                "+",
                Integer.toString(shiftBackStart),
                Integer.toString(shiftBackEnd),
                name,
                Integer.toString(length),
                "+",
                Integer.toString(start),
                Integer.toString(end),
                Integer.toString(id)
        };
        return String.join("\t", items);
    }

    @Override
    public Object onTraversalSuccess(){
        return null;
    }

    @Override
    public void closeTool() {
        super.closeTool();
        try{
            if( refWriter != null ) {
                refWriter.close();
            }
        } catch (IOException e) {
            throw new UserException("Failed to write fasta due to " + e.getMessage(), e);
        }
        try{
            if (chainFileWriter != null) {
                chainFileWriter.close();
            }
        } catch (IOException e) {
            throw new UserException("Failed to write chain file due to " + e.getMessage(), e);
        }
        try{
            if (intervalRegularWriter != null) {
                intervalRegularWriter.close();
            }
            if (intervalShiftedWriter != null) {
                intervalShiftedWriter.close();
            }
        } catch (IOException e) {
            throw new UserException("Failed to write intervals due to " + e.getMessage(), e);
        }
    }

    }
