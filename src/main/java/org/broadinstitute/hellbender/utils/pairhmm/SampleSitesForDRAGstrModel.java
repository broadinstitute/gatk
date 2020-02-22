package org.broadinstitute.hellbender.utils.pairhmm;


import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaReferenceWriter;
import org.apache.commons.compress.archivers.jar.JarArchiveOutputStream;
import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;
import org.apache.commons.compress.archivers.zip.ZipArchiveOutputStream;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.fasta.FastaReferenceMaker;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.AutoCloseableList;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.*;
import java.nio.file.Paths;
import java.util.*;

@CommandLineProgramProperties(
        programGroup = ReferenceProgramGroup.class,
        summary = "Determine the presence of STR in a reference sequence",
        oneLineSummary = "Determines the presence of STR in a reference sequence"
)
public class SampleSitesForDRAGstrModel extends GATKTool {

    private static final Logger logger = LogManager.getLogger(SampleSitesForDRAGstrModel.class);

    public static class DecimationTable {

        private static final int[][] DEFAULT_DECIMATION_MATRIX = new int[][] {
                {0}, // 0, 0, 0, 0, 0, 0, 0, 0 ...
                {0, 10, 10, 9, 8, 7, 5, 3, 1, 0},
                {0, 0, 9, 6, 3, 0}, // 0, 0, 0 ...
                {0, 0, 8, 4, 1, 0},
                {0, 0, 6, 0},
                {0, 0, 5, 0},
                {0, 0, 4, 0},
                {0}};

        public static final String NO_DECIMATION_STR = "NONE";

        public static final String DEFAULT_DECIMATION_STR = "DEFAULT";

        public static final DecimationTable DEFAULT = new DecimationTable(DEFAULT_DECIMATION_STR);

        public static final DecimationTable NONE = new DecimationTable(NO_DECIMATION_STR);

        private final long[][] decimationMask;

        private final long[][] counts;

        public DecimationTable(final String spec) {
            Utils.nonNull(spec);
            final int[][] decimation;
            if (spec.equalsIgnoreCase(NO_DECIMATION_STR)) {
                decimation = new int[][] {{0}};
            } else if (spec.equalsIgnoreCase(DEFAULT_DECIMATION_STR)) {
                decimation = DEFAULT_DECIMATION_MATRIX;
            } else {
                decimation = parseDecimationMatrixFromPath(spec);
            }
            decimationMask = calculateDecimationMask(decimation);
            counts = composeDecimationCounts(decimationMask);
        }

        public long decimationMask(final int period, final int repeats) {
            if (decimationMask.length <= period) {
                return -1;
            } else if (decimationMask[period].length <= repeats) {
                return -1;
            } else {
                return decimationMask[period][repeats];
            }
        }

        private long[][] composeDecimationCounts(final long[][] decimationMask) {
            final long[][] result = new long[decimationMask.length][];
            for (int i = 0; i < result.length; i++) {
                result[i] = new long[decimationMask[i].length];
            }
            return result;
        }

        private static int[][] parseDecimationMatrixFromPath(String spec) {
            try (final BufferedReader reader = new BufferedReader(IOUtils.makeReaderMaybeGzipped(Paths.get(spec)))) {
                final String[][] values = reader.lines()
                        .filter(str -> !str.startsWith("#") && !str.trim().isEmpty())
                        .map(str -> Arrays.stream(str.split("\\s+"))
                                   .mapToDouble(Double::parseDouble)
                                   .toArray())
                        .toArray(String[][]::new);
                return parseDecimationMatrixValues(values, spec);
            } catch (final IOException ex) {
                throw new UserException.CouldNotReadInputFile(spec, ex);
            } catch (final NumberFormatException ex){
                throw new UserException.BadInput(String.format("input decimation file %s contains non-valid values: %s", spec, ex.getMessage()));
            }
        }

        private static int[][] parseDecimationMatrixValues(final String[][] values, final String path) {
            Utils.nonNull(values);
            if (values.length == 0) {
                logger.warn("Decimation matrix path provided does not seem to contain any values, we will proceed without any decimation");
                return new int[0][];
            } else {
                int totalValues = 0;
                final int[][] result = new int[values.length][];
                for (int i = 0; i < values.length; i++) {
                    final String[] row = values[i];
                    final int[] rowValues = new int[values.length];
                    for (int j = 0; j <  row.length; j++) {
                        final String str = row[j];
                        final int value;
                        try {
                            value = Integer.parseInt(str);
                        } catch (final NumberFormatException ex) {
                            throw badDecimationValueException(str, path, i, j, "not a valid double literal");
                        }
                        if (value < 0) {
                            throw badDecimationValueException(str, path, i, j, "negatives are not allowed");
                        } else if (Double.isNaN(value)) {
                            throw badDecimationValueException(str, path, i, j, "NaN are not allowed");
                        } else if (!Double.isFinite(value)) {
                            throw badDecimationValueException(str, path, i, j, "must be finite");
                        }
                        rowValues[j] = value;
                        totalValues++;
                    }
                    result[i] = rowValues;
                }
                if (totalValues == 0) {
                    throw new UserException.BadInput("the input decimation matrix does contain any values:" + path);
                }
                return result;
            }
        }

        private static RuntimeException badDecimationValueException(final String str, final String path, final int i, final int j,
                                                                    final String details) {
            throw new UserException.BadInput(String.format("bad decimation value found in %s for period and repeats (%d, %d) with string (%s)%s",
                    path, i, j, str, details == null || details.isEmpty()? "": ": " + details));
        }

        public static long[][] calculateDecimationMask(final int[][] decimationMatrix) {
            Utils.nonNull(decimationMatrix);
            final long[][] result = new long[decimationMatrix.length][];
            for (int i = 0; i < result.length; i++) {
                final int[] row = decimationMatrix[i];
                result[i] = new long[row.length];
                for (int j = 0; j < row.length; j++) {
                    result[i][j] = (1 << row[j]) - 1;
                }
            }
            return result;
        }

        public long mask(final int period, final int repeats) {
            final int p = period >= decimationMask.length ? decimationMask.length - 1 : period;
            final long[] masks = decimationMask[p];
            if (masks.length == 0) {
                return 0;
            } else if (repeats >= masks.length) {
                return masks[masks.length - 1];
            } else {
                return masks[repeats];
            }
        }

        public boolean decimate(final int seqNumber, final int bestPeriod, final long bestPeriodRepeats) {
            if (counts.length <= bestPeriod) {
                return false;
            } else {
                final long[] periodCounts = counts[bestPeriod];
                if (periodCounts.length == 0) {
                    return false;
                } else {
                    final int effectiveRepeatCount
                            = (int) (bestPeriodRepeats < periodCounts.length ? bestPeriodRepeats : periodCounts.length - 1);
                    final long count = periodCounts[effectiveRepeatCount]++;
                    final long left = count + seqNumber;
                    final long right = decimationMask[bestPeriod][effectiveRepeatCount];
                    return ((int) left & (int) right) != 0 || ((left >> 32) & (right >> 32)) != 0;
                }
            }
        }
    }

    @Argument(fullName="decimation", doc="decimation per period and repeat. It can be \"DEFAULT\" to use the default values (default), " +
            " \"NONE\" to deactivate decimation (potentially resulting in a very large output file) or indicate the path to a file" +
            " that contains the decimation matrix.", optional = true)
    private DecimationTable decimationTable = DecimationTable.DEFAULT;

    @Argument(doc = "name of the zip file where the sites sampled will be stored",
              fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    private String outputPath = null;

    @Argument(fullName="max-period", doc="maximum STR period sampled", optional = true, minValue = 1, maxValue = 10)
    private int maxPeriod = 8;

    @Argument(fullName="max-repeats", doc="maximum STR repeat sampled", optional = true, minValue = 1, maxValue = 20)
    private int maxRepeat = 20;

    @Override
    public boolean requiresReference() {
        return true;
    }

    private File tempDir;

    private long[][] totalCounts;

    private long[][] emittedCounts;

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
        totalCounts = new long[maxPeriod + 1][maxRepeat + 1];
        emittedCounts = new long[maxPeriod + 1][maxRepeat + 1];
    }

    public void onShutdown() {
        Utils.deleteFileTree(tempDir);
        super.onShutdown();
    }

    @Override
    public void traverse() {
        final SAMSequenceDictionary dictionary = Utils.nonNull(getReferenceDictionary());
        final File tempDir = this.tempDir;
        try (final AutoCloseableList<BinaryTableWriter<DragstrLocus>> allSitesWriter
                     = AutoCloseableList.of(maxPeriod, i -> createDragstrLocusWriter(tempDir, "period_" + (i + 1) + ".bin"))) {
            if (!intervalArgumentCollection.intervalsSpecified()) {
                for (final SAMSequenceRecord sequence : dictionary.getSequences()) {
                    traverse(sequence.getSequenceIndex(), sequence, 1, sequence.getSequenceLength(), allSitesWriter, decimationTable);
                }
            } else {
                for (final SimpleInterval interval : intervalArgumentCollection.getIntervals(dictionary)) {
                    final SAMSequenceRecord sequence = dictionary.getSequence(interval.getContig());
                    traverse(sequence.getSequenceIndex(), sequence, interval.getStart(), interval.getEnd(), allSitesWriter, decimationTable);
                }
            }
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(tempDir, ex);
        }
        progressMeter.stop();
        logger.info("Finishing reference sampling. Proceeding to splitting each period case by repeat count");
        for (int i = 1; i <= maxPeriod; i++) {
            progressMeter = new ProgressMeter(secondsBetweenProgressUpdates);
            progressMeter.setRecordLabel("Splitting cases with period " + i + " by repeat count");
            progressMeter.start();
            splitPeriodLociByRepeat(i);
            progressMeter.stop();
            logger.info("Done with period " + i);
        }
        logger.info("Composing output zip");
        composeOutputZip();
    }

    private void splitPeriodLociByRepeat(final int period) {
        final File lociFile = new File(tempDir, "period_" + period + ".bin");
        final File periodDir = new File(tempDir, "" + period);
        if (!periodDir.mkdir()) {
            throw new UserException.CouldNotCreateOutputFile(periodDir, "create period loci split directory");
        }
        try (final BinaryTableReader<DragstrLocus> reader = DragstrLocus.binaryReader(new FileInputStream(lociFile));
             final AutoCloseableList<BinaryTableWriter<DragstrLocus>> writers =
                     AutoCloseableList.of(maxRepeat, i -> createDragstrLocusWriter(periodDir, "" + (i+1) + ".bin"))) {
            while (reader.hasNext()) {
                final DragstrLocus next = reader.next();
                final int effectiveRepeats = next.getRepeats() > maxRepeat ? maxRepeat : next.getRepeats();
                progressMeter.update(next.getStartInterval(getReferenceDictionary(), 0));
                final BinaryTableWriter<DragstrLocus> writer = writers.get(effectiveRepeats - 1);
                writer.write(next);

            }
        } catch (final Exception ex) {
            throw new UserException.CouldNotCreateOutputFile(lociFile, "temporary period loci splitting failed");
        }
        if (!lociFile.delete()) {
            throw new GATKException("could not delete temporary file " + lociFile);
        }
    }

    private void composeOutputZip() {
        final byte[] buffer = new byte[1 << 16];
        try (final ZipArchiveOutputStream output = new JarArchiveOutputStream(BucketUtils.createFile(outputPath))) {
            composeSummaryText(output);
            saveReferenceDictionary(output);
            composeLociFiles(buffer, output);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputPath, "error composing the zip file", ex);
        }
    }

    private void composeLociFiles(byte[] buffer, ZipArchiveOutputStream output) throws IOException {
        for (int i = 1; i <= maxPeriod; i++) {
            for (int j = 1; j <= maxRepeat; j++) {
                final File inFile = new File(new File(tempDir, "" + i), "" + j + ".bin");
                output.putArchiveEntry(new ZipArchiveEntry("" + i + "/" + j + ".bin"));
                final InputStream inStream = new FileInputStream(inFile);
                int len;
                while ((len = inStream.read(buffer)) > 0) {
                    output.write(buffer, 0, len);
                }
                inStream.close();
                output.closeArchiveEntry();
            }
        }
    }

    private void saveReferenceDictionary(final ZipArchiveOutputStream output) throws IOException {
        final SAMSequenceDictionary dictionary = getBestAvailableSequenceDictionary();
        output.putArchiveEntry(new ZipArchiveEntry("reference.dict"));
        // need to prevent closing output if the writer's close is called. We simply
        // flush instead:
        try (final Writer writer = new OutputStreamWriter(output) {
                public void close() throws IOException {
                    flush();
                }
        }) {
            final SAMSequenceDictionaryCodec codec = new SAMSequenceDictionaryCodec(writer);
            codec.encode(dictionary);
        }
        output.closeArchiveEntry();
    }

    private void composeSummaryText(ZipArchiveOutputStream output) throws IOException {
        output.putArchiveEntry(new ZipArchiveEntry("summary.txt"));
        final PrintWriter writer = new PrintWriter(new OutputStreamWriter(output));
        writer.println("period\trepeat\ttotal\temmitted\tdecimation\tactual_decimation");
        for (int i = 1; i <= maxPeriod; i++) {
            for (int j = 1; j <= maxRepeat; j++) {
                writer.println(String.format("%d\t%d\t%d\t%d\t%.2f\t%.2f",
                        i, j, totalCounts[i][j], emittedCounts[i][j],
                        Math.log(decimationTable.decimationMask(i,j) + 1) / Math.log(2),
                        (- Math.log(emittedCounts[i][j]) + Math.log(totalCounts[i][j])) / Math.log(2) ));
            }
        }
        writer.flush();
        output.closeArchiveEntry();
    }

    private static BinaryTableWriter<DragstrLocus> createDragstrLocusWriter(final File outDir, final String name) {
        try {
            final File outFile = new File(outDir, name);
            return DragstrLocus.binaryWriter(outFile);
        } catch (final IOException  e) {
            throw new UncheckedIOException(e);
        }
    }

    private void traverse(final int seqNumber, final SAMSequenceRecord sequence, final long seqStart, final long seqEnd, final AutoCloseableList<BinaryTableWriter<DragstrLocus>> output, final DecimationTable decimationTable)
        throws IOException
    {
        final String id = sequence.getSequenceName();
        final NucleotideSequence fullSequence = LazyLoadingReferenceNucleotideSequence.of(directlyAccessEngineReferenceDataSource(), sequence.getSequenceName(),10000);
        final int maxPeriod = this.maxPeriod;
        final int length = sequence.getSequenceLength();
        long pos = seqStart;
        utter: while (pos <= seqEnd) {
           final byte base = fullSequence.byteAt(pos);
           long beg, end, cmp;
           for (beg = pos - 1; beg >= 1 && fullSequence.byteAt(beg) == base; beg--);
           beg++;
           for (end = pos + 1; end <= length && fullSequence.byteAt(end) == base; end++);
           end--;
           int bestPeriod = 1;
           long bestPeriodRepeats = end - beg + 1;
           long bestEnd = end;
           long bestBeg = beg;
           for (int period = 2; period <= maxPeriod; period++) {
               for (beg = length - pos > period ? pos - 1 : length - period,
                    cmp = beg + period; beg >= 1 && fullSequence.byteAt(cmp) == fullSequence.byteAt(beg); beg--, cmp--);
               beg++;
               for (end = pos >= period ? pos + 1 : period + 1, cmp = end - period; end <= length && fullSequence.byteAt(end) == fullSequence.byteAt(cmp); end++, cmp++);
               end--;
               final long strLength = end - beg + 1;
               final long repeats =  strLength / period;
               if (repeats > bestPeriodRepeats) {
                   bestPeriod = period;
                   bestPeriodRepeats = repeats;
                   bestEnd = end;
                   bestBeg = beg;
               }
           }
           final byte[] unit = fullSequence.bytesAt(bestBeg, bestPeriod);
           pos = bestEnd + 1;

           for (int i = 0; i < bestPeriod; i++) {
               if (!Nucleotide.decode(unit[i]).isStandard()) {
                   continue utter;
               }
           }
           emitOrDecimateSTR(seqNumber, id, bestBeg, bestBeg + bestPeriod * bestPeriodRepeats - 1, bestPeriod, bestPeriodRepeats, output, decimationTable, unit);
        }
    }

    private void emitOrDecimateSTR(final int seqNumber, final String seqId, long bestBeg, long bestEnd, int bestPeriod, long bestPeriodRepeats, final AutoCloseableList<BinaryTableWriter<DragstrLocus>> output, final DecimationTable decimationTable, final byte[] unit)
        throws IOException
    {
        final int effectiveRepeats = Math.min(maxRepeat, (int) bestPeriodRepeats);
        if (!decimationTable.decimate(seqNumber, bestPeriod, bestPeriodRepeats)) {
            final DragstrLocus locus = DragstrLocus.make(seqNumber, bestBeg, unit, (int) bestPeriodRepeats);
            final BinaryTableWriter<DragstrLocus> outTable = output.get(bestPeriod - 1);
            outTable.write(locus);
            progressMeter.update(new SimpleInterval(seqId, (int) bestBeg, (int) bestBeg));
            emittedCounts[bestPeriod][effectiveRepeats]++;
        }
        totalCounts[bestPeriod][effectiveRepeats]++;
    }
}
