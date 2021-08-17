package org.broadinstitute.hellbender.tools.dragstr;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.dragstr.STRTableFile;
import org.broadinstitute.hellbender.utils.dragstr.STRTableFileBuilder;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * This tool looks for low-complexity STR sequences along the reference that are later used to estimate the Dragstr model 
 * during single sample auto calibration {@link CalibrateDragstrModel}.
 * <h3>Inputs</h3>
 * <p>
 *     This command takes as input the reference (possibly a subset of intervals) and an optional {@link STRDecimationTable decimation table} (herein referred as DT). 
 * </p>
 * <p>
 *     The DT modulates how often we sample a site for each possible period and repeat length. Since there is far more
 *     positions with short period and short repeat length sampling for those combinations should be less frequent.
 *     For further details about the format of this table and interpretation of its values please check the documentation
 *     in class {@link STRDecimationTable}.
 * </p>
 * <p>
 *    If no DT is provided, the tool uses a default one that has been tailored to work fine when run over 
 *    the entire Human genome and it should be alright with other genomes of comparable size (i.e. 1 to 10Gbps).
 *    With larger genomes, that default DT will likely result in an unecessarely large number of sampled sites 
 *    that it turn may increase the run time of tools that depend on the output. In contrast, 
 *    with smaller genomes or subsets (using targeted intervals) it might result in a number of sampled sites 
 *    too small to build accurate Dragstr model. In this case you really need to compose and provide
 *    your own DT or perhaps try out not to decimate at all ({@code --decimation NONE}). 
 * </p> 
 * <h3>Output</h3>
 * <p>
 *    The output of this command is a zip file that contains the collection of sampled sites in 
 *    binary form ({@code all.bin}), and index for that file for quick access by location interval 
 *    ({@code all.idx}). Other files in the zip provide some summary and tracking information, for example
 *    the reference sequence dictionary ({@code reference.dict}), a copy of the DT ({@code decimation.txt}) 
 *    and summarized stats ({@code summary.txt}).
 * </p>
 * <p>
 * <h3>Examples</h3>
 * <pre>

 *     # Human? just use the default.
 *     gatk ComposeSTRTableFile -R hg19.fasta -O hg19.str.zip
 *     # or ...
 *     gatk ComposeSTRTableFile -R hg19.fasta --decimation DEFAULT -O hg19.str.zip
 *
 * </pre>
 * <pre>
 *
 *     # yeast genome is roughly ~ 12Mbp long.
 *     gatk ComposeSTRTableFile -R yeast.fasta --decimation custom-yeast.dt -O yeast.str.zip
 *
 * </pre>
 * <pre>
 *
 *     #  Carsonella ruddii just about 160Kbps, prorably we don't want to decimate at all:
 *     gatk ComposeSTRTableFile -R Cruddii.fasta --decimation NONE -O yeast.str.zip
 *
 * </pre>
 * </p>
 */
@CommandLineProgramProperties(
        programGroup = ReferenceProgramGroup.class,
        oneLineSummary = "Composes a genome-wide STR location table used for DragSTR model auto-calibration",
        summary = "Composes a genome-wide STR location table used for DragSTR model auto-calibration"
)
@DocumentedFeature
public final class ComposeSTRTableFile extends GATKTool {

    public static final String REFERENCE_SEQUENCE_BUFFER_SIZE_FULL_NAME = "reference-sequence-buffer-size";
    public static final String GENERATE_SITES_TEXT_OUTPUT_FULL_NAME = "generate-sites-text-output";
    public static final int MINIMUM_REFERENCE_SEQUENCE_BUFFER_SIZE = 1024;
    public static final int MAXIMUM_REFERENCE_SEQUENCE_BUFFER_SIZE = 100_000_000;
    public static final int DEFAULT_REFERENCE_SEQUENCE_BUFFER_SIZE = 100_000;

    @Argument(fullName = "decimation", 
              doc = "decimation per period and repeat. It can be \"DEFAULT\" to use the default values (DEFAULT), " +
                    " \"NONE\" to deactivate decimation (potentially resulting in a very large output file) " + 
                    "or indicate the path to a file that contains the decimation matrix.", 
              optional = true)
    private STRDecimationTable decimationTable = STRDecimationTable.DEFAULT;

    @Argument(doc = "name of the zip file where the sites sampled will be stored",
              fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    private GATKPath outputPath = null;

    @Argument(doc = "request to generate a text formatted version of the STR table in the output zip (" + STRTableFile.SITES_TEXT_FILE_NAME + ")",
              fullName = GENERATE_SITES_TEXT_OUTPUT_FULL_NAME, optional = true)
    @Hidden
    private boolean generateSitesTextOutput = false;

    @Argument(fullName = DragstrHyperParameters.MAX_PERIOD_ARGUMENT_FULL_NAME, doc="maximum STR period sampled", optional = true, minValue = 1, maxValue = 20)
    private int maxPeriod = DragstrHyperParameters.DEFAULT_MAX_PERIOD;

    @Argument(fullName = DragstrHyperParameters.MAX_REPEATS_ARGUMENT_FULL_NAME, doc="maximum STR repeat sampled", optional = true, minValue = 1, maxValue = 100)
    private int maxRepeat = DragstrHyperParameters.DEFAULT_MAX_REPEAT_LENGTH;

    @Argument(fullName= REFERENCE_SEQUENCE_BUFFER_SIZE_FULL_NAME, doc="size of the look ahead reference sequence buffer",
            optional = true, minValue = MINIMUM_REFERENCE_SEQUENCE_BUFFER_SIZE, maxValue = MAXIMUM_REFERENCE_SEQUENCE_BUFFER_SIZE)
    @Hidden
    private int referenceSequenceBufferSize = DEFAULT_REFERENCE_SEQUENCE_BUFFER_SIZE;

    @Override
    public boolean requiresReference() {
        return true;
    }

    private static final String COMMAND_LINE_ANNOTATION_NAME = "commandLine";

    private File tempDir;
    private int[][][] nextMasks;

    public void onStartup() {
        super.onStartup();
        try {
            tempDir = File.createTempFile("gatk-sample-dragstr-sites", ".tmp");
        } catch (final IOException ex) {
            throw new GATKException("could not create temporary disk space", ex);
        }
        if (!tempDir.delete()) {
            throw new GATKException("could not create temporary disk space: could not delete tempfile");
        } else if (!tempDir.mkdir()) {
            throw new GATKException("could not create temporary disk space: could not create tempdir");
        }
    }

    public void onShutdown() {
        try {
            if (tempDir != null) {
                FileUtils.deleteDirectory(tempDir);
            }
        } catch (final IOException e) {
            throw new GATKException("issues removing temporary directory: " + tempDir, e);
        } finally {
            super.onShutdown();
        }
    }

    @Override
    public void traverse() {
        final SAMSequenceDictionary dictionary = getBestAvailableSequenceDictionary();
        initializeMasks(dictionary);
        try (final STRTableFileBuilder output = STRTableFileBuilder.newInstance(dictionary, decimationTable,
                generateSitesTextOutput, maxPeriod, maxRepeat)) {
            output.annotate(COMMAND_LINE_ANNOTATION_NAME, getCommandLine());
            final Map<String, List<SimpleInterval>> intervalsByContig = composeAndGroupTraversalIntervalsByContig(dictionary);
            for (final Map.Entry<String, List<SimpleInterval>> contigEntry : intervalsByContig.entrySet()) {
                final BufferedReferenceBases nucleotideSequence = BufferedReferenceBases.of(directlyAccessEngineReferenceDataSource(), contigEntry.getKey(), referenceSequenceBufferSize);
                final SAMSequenceRecord sequenceRecord = dictionary.getSequence(contigEntry.getKey());
                for (final SimpleInterval interval : contigEntry.getValue()) {
                    traverseInterval(sequenceRecord.getSequenceIndex(), nucleotideSequence, interval.getStart(), interval.getEnd(), decimationTable, output);
                }
            }
            output.store(outputPath);
        }
        progressMeter.stop();
    }

    private Map<String, List<SimpleInterval>> composeAndGroupTraversalIntervalsByContig(final SAMSequenceDictionary dictionary) {
        if (!intervalArgumentCollection.intervalsSpecified()) {
            return dictionary.getSequences().stream()
                    .map(s -> new SimpleInterval(s.getSequenceName(), 1, s.getSequenceLength()))
                    .collect(Collectors.groupingBy(SimpleInterval::getContig, LinkedHashMap::new, Collectors.toList()));
        } else {
            final Map<String, List<SimpleInterval>> keyUnsorted = IntervalUtils.sortAndMergeIntervals(intervalArgumentCollection.getIntervals(dictionary), dictionary, IntervalMergingRule.ALL);
            final Map<String, List<SimpleInterval>> keySorted = new LinkedHashMap<>(keyUnsorted.size());
            keyUnsorted.entrySet().stream()
                    .sorted(Comparator.comparingInt(entry -> dictionary.getSequenceIndex(entry.getKey())))
                    .forEach(entry -> keySorted.put(entry.getKey(), entry.getValue()));
            return keySorted;
        }
    }

    private void initializeMasks(final SAMSequenceDictionary dictionary) {
        nextMasks = new int[dictionary.getSequences().size()][maxPeriod + 1][maxRepeat + 1];
        for (int i = 0; i < nextMasks.length; i++) {
            for (final int[] masks : nextMasks[i]) {
                Arrays.fill(masks, i);
            }
        }
    }

    private void traverseInterval(final int seqNumber, final BufferedReferenceBases sequence, final long seqStart, final long seqEnd,
                                  final STRDecimationTable decimationTable, final STRTableFileBuilder output)
    {
        final String id = sequence.contigID();
        final long length = sequence.length();
        final byte[] unitBuffer = new byte[this.maxPeriod];
        for (long pos = seqStart; pos <= seqEnd; pos++) {
            final BestPeriodRepeat best = findBestPeriodRepeatCombination(sequence, pos, length, unitBuffer);
            if (best != null) {
               pos = best.end;
               emitOrDecimateSTR(seqNumber, id, best, decimationTable, output);
            }
        }
    }

    /**
     * Returns the STR spec starting at a particular position on the input reference base sequence.
     * @return {@code null} iff the base at the starting position is not a standard base (i.e. A, C, G, T) for example N X or ?.
     */
    private BestPeriodRepeat findBestPeriodRepeatCombination(final BufferedReferenceBases sequence, final long pos,
                                                             final long length, final byte[] unitBuffer) {
        final BestPeriodRepeat best;
        // usually maxPeriodAtPos == maxPeriod except when close to the end of the sequence
        // is the maximum period to be considered that cannot exceed min(maxPeriod, seq-length - pos + 1)
        // but is always 1 or greater.
        final int maxPeriodAtPos = sequence.copyBytesAt(pos, unitBuffer, 0, maxPeriod);

        // Efficient code for period == 1:
        final byte firstUnitBase;
        if (!Nucleotide.decode(firstUnitBase = unitBuffer[0]).isStandard()) {
            best = null;
        } else {
            long beg, end;

            // we look upstream for same bases:
            for (beg = pos - 1; beg >= 1 && Nucleotide.same(sequence.byteAt(beg), firstUnitBase); beg--) { /* nothing to be done */ }
            beg++;

            // we look downstream for same bases:
            for (end = pos + 1; end <= length && Nucleotide.same(sequence.byteAt(end), firstUnitBase); end++) { /* nothing to be done */ }
            end--;

            // Initialize best period to 1:
            best = BestPeriodRepeat.initializeToOne(beg, end);
            // Now wo do period 2 and beyond:
            int cmp; // var to hold the next position in the unit (buffer) to compare against.
            for (int period = 2; period <= maxPeriodAtPos; period++) {
                // we stop if the last base in unit is not ACGT (usually N):
                if (!Nucleotide.decode(unitBuffer[period - 1]).isStandard()) {
                    break;
                }
                // We look upstream for matching p-mers
                for (beg = pos - 1,
                             cmp = period - 1; beg >= 1 && Nucleotide.same(sequence.byteAt(beg), unitBuffer[cmp]); beg--) {
                    if (--cmp == -1) {
                        cmp = period - 1;
                    }
                }
                beg++;
                // We look downstream for matching p-mers
                for (cmp = 0, end = pos + period; end <= length && Nucleotide.same(sequence.byteAt(end), unitBuffer[cmp]); end++) {
                    if (++cmp == period) {
                        cmp = 0;
                    }
                }
                end--;
                best.updateIfBetter(period, beg, end);
            }
        }
        return best;
    }

    private void emitOrDecimateSTR(final int seqNumber, final String seqId, final BestPeriodRepeat best,
                                   final STRDecimationTable decimationTable, final STRTableFileBuilder output)
    {
        final int effectiveRepeats = Math.min(maxRepeat, best.repeats);
        final int mask = nextMasks[seqNumber][best.period][effectiveRepeats]++;
        if (!decimationTable.decimate(mask, best.period, best.repeats)) {
            final DragstrLocus locus = DragstrLocus.make(seqNumber, best.start, (byte) best.period, (short) (best.end - best.start + 1), mask);
            output.emit(locus);
            progressMeter.update(new SimpleInterval(seqId, (int) best.start, (int) best.start));
        } else {
            output.decimate(best.period, best.repeats);
        }
    }

    private static class BestPeriodRepeat {
        private int period;
        private int repeats;
        private long start;
        private long end;

        BestPeriodRepeat(final int period, final long start, final long end) {
            this.start = start;
            this.end = end;
            this.period = period;
            this.repeats = (int) (end - start + 1) / period;
        }

        private static BestPeriodRepeat initializeToOne(final long start, final long end) {
            return new BestPeriodRepeat(1, start, end);
        }

        private void updateIfBetter(final int newPeriod, final long newStart, final long newEnd) {
            final int newRepeats = (int) (newEnd - newStart + 1) / newPeriod;
            if (newRepeats > repeats || (newRepeats == repeats && newPeriod < period)) {
                start = newStart;
                end = newEnd;
                period = newPeriod;
                repeats = newRepeats;
            }
        }
    }
}
