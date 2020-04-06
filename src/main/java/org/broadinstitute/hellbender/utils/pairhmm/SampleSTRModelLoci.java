package org.broadinstitute.hellbender.utils.pairhmm;


import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.compress.archivers.jar.JarArchiveOutputStream;
import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;
import org.apache.commons.compress.archivers.zip.ZipArchiveOutputStream;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.AutoCloseableList;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        programGroup = ReferenceProgramGroup.class,
        summary = "Determine the presence of STR in a reference sequence",
        oneLineSummary = "Determines the presence of STR in a reference sequence"
)
public class SampleSTRModelLoci extends GATKTool {

    private static final Logger logger = LogManager.getLogger(SampleSTRModelLoci.class);

    @Argument(fullName="decimation", doc="decimation per period and repeat. It can be \"DEFAULT\" to use the default values (default), " +
            " \"NONE\" to deactivate decimation (potentially resulting in a very large output file) or indicate the path to a file" +
            " that contains the decimation matrix.", optional = true)
    private STRDecimationTable decimationTable = STRDecimationTable.DEFAULT;

    @Argument(doc = "name of the zip file where the sites sampled will be stored",
              fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    private String outputPath = null;

    @Argument(doc = "name of output file to dump entries sorted by coordinates",
              fullName = "sorted-output", optional = true)
    @Hidden
    private String sortedOutputPath = null;

    @Argument(fullName="max-period", doc="maximum STR period sampled", optional = true, minValue = 1, maxValue = 10)
    private int maxPeriod = 8;

    @Argument(fullName="max-repeats", doc="maximum STR repeat sampled", optional = true, minValue = 1, maxValue = 20)
    private int maxRepeat = 20;

    @Argument(fullName="down-sample", doc="target number of STR sites per period, repeat-count combination. The default of 0 indicates no down-sampling", optional = true, minValue = 0)
    private int downSample = 0;


    @Override
    public boolean requiresReference() {
        return true;
    }

    private File tempDir;

    private long[][] totalCounts;
    private int[][][] nextMasks;

    private long[][] emittedCounts;

    private PrintWriter sortedOutputWriter;

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

        sortedOutputWriter = sortedOutputPath == null ? null : new PrintWriter(new OutputStreamWriter(BucketUtils.createFile(sortedOutputPath)));
        totalCounts = new long[maxPeriod + 1][maxRepeat + 1];
        emittedCounts = new long[maxPeriod + 1][maxRepeat + 1];
    }

    public void onShutdown() {
        try {
            if (tempDir != null) {
                Utils.deleteFileTree(tempDir);
            }
            if (sortedOutputWriter != null) {
                sortedOutputWriter.close();
            }
        } finally {
            super.onShutdown();
        }
    }

    @Override
    public void traverse() {
        final SAMSequenceDictionary dictionary = Utils.nonNull(getReferenceDictionary());
        final File tempDir = this.tempDir;
        initializeMasks(dictionary);
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
            if (sortedOutputWriter != null) {
                sortedOutputWriter.flush();
            }
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(tempDir, ex);
        }
        progressMeter.stop();
        logger.info("Finished reference sampling. Proceeding to splitting each period case by repeat count");
        for (int i = 1; i <= maxPeriod; i++) {
            progressMeter = new ProgressMeter(secondsBetweenProgressUpdates);
            progressMeter.setRecordLabel("Splitting cases with period " + i + " by repeat count");
            progressMeter.start();
            splitPeriodLociByRepeat(i);
            progressMeter.stop();
            logger.info("Done with period " + i);
        }
        if (downSample > 0) {
            logger.info("Downsampling period, repeat-count case to a maximum of " + downSample + " each.");
            for (int i = 1; i <= maxPeriod; i++) {
                for (int r = 1; r <= maxRepeat; r++) {
                    downSampleLociByRepeat(i, r);
                }
            }
        }
        logger.info("Composing output zip");
        composeOutputZip();
    }

    private void downSampleLociByRepeat(int period, int repeatCount) {
        final File periodDir = new File(tempDir, "" + period);
        final File inFile = new File(periodDir, repeatCount + ".bin");
        final List<DragstrLocus> all;
        try (final BinaryTableReader<DragstrLocus> reader = DragstrLocus.binaryReader(new FileInputStream(inFile))) {
            all = reader.stream().collect(Collectors.toList());
        } catch (final Exception ex) {
            throw new UserException.CouldNotCreateOutputFile(inFile, "temporary period loci downsample failed");
        }
        if (all.size() > downSample) {
            final long bit = decimationTable.decimationBit(period, repeatCount);
            long mask = 1 << bit;
            List<DragstrLocus> result = all;
            while (result.size() > downSample) {
                final long filterMask = mask;
                result = result.stream()
                        .filter(locus -> (locus.getMask() & filterMask) == 0L)
                        .collect(Collectors.toList());
                mask <<= 1;
            }
            logger.info("Down-sampling period " + period + " repeat count " + repeatCount + " from " + all.size() + " (> " + downSample + ") to " + result.size());
            try (final BinaryTableWriter<DragstrLocus> writer = DragstrLocus.binaryWriter(inFile)) {
                writer.writeAll(result);
            } catch (final Exception ex) {
                throw new UserException.CouldNotCreateOutputFile(inFile, "temporary period loci downsample failed");
            }
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
            printDecimationMatrix(output);
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

    private void printDecimationMatrix(final ZipArchiveOutputStream output) throws IOException {
        output.putArchiveEntry(new ZipArchiveEntry("decimation.txt"));
        final PrintStream writer = new PrintStream(output);
        decimationTable.printDecimationMatrix(writer);
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

    private void traverse(final int seqNumber, final SAMSequenceRecord sequence, final long seqStart, final long seqEnd, final AutoCloseableList<BinaryTableWriter<DragstrLocus>> output, final STRDecimationTable decimationTable)
        throws IOException
    {
        final String id = sequence.getSequenceName();
        final NucleotideSequence fullSequence = LazyLoadingReferenceNucleotideSequence.of(directlyAccessEngineReferenceDataSource(), sequence.getSequenceName(),10000);
        final int length = sequence.getSequenceLength();
        final byte[] unitBuffer = new byte[this.maxPeriod];
        long pos = seqStart;
        utter: while (pos <= seqEnd) {
           // usually maxPeriodAtPos == maxPeriod except when close to the end of the sequence
           // is the maximum period to be considered that cannot exceed min(maxPeriod, seq-length - pos + 1)
           // but is always 1 or greater.
           final int maxPeriodAtPos = fullSequence.copyBytesAt(unitBuffer, 0, pos, maxPeriod);

           // Efficient code for period == 1:
           final byte firstUnitBase;
           if (!Nucleotide.decode(firstUnitBase = unitBuffer[0]).isStandard()) {
               pos++; continue;
           }
           long beg, end;
           // we look upstream for same bases:
           for (beg = pos - 1; beg >= 1 && Nucleotide.same(fullSequence.byteAt(beg), firstUnitBase); beg--);
           beg++;
           // we look downstream for same bases:
           for (end = pos + 1; end <= length && Nucleotide.same(fullSequence.byteAt(end), firstUnitBase); end++);
           end--;
           // Initialize best period to 1:
           final BestPeriod best = BestPeriod.initialize(1, beg, end);
           // Now wo do period 2 and beyond:
           int cmp; // var to hold the next position in the unit (buffer) to compare against.
           for (int period = 2; period <= maxPeriodAtPos; period++) {
               // we stop if the last base in unit is not ACGT (usually N):
               if (!Nucleotide.decode(unitBuffer[period - 1]).isStandard()) {
                   break;
               }
               // We look upstream for matching p-mers
               for (beg = pos - 1,
                   cmp = period - 1; beg >= 1 && Nucleotide.same(fullSequence.byteAt(beg), unitBuffer[cmp]); beg--) {
                   if (--cmp == -1) {
                       cmp = period - 1;
                   }
               }
               beg++;
               // We look downstream for matching p-mers
               for (cmp = 0, end = pos + period; end <= length && Nucleotide.same(fullSequence.byteAt(end), unitBuffer[cmp]); end++) {
                   if (++cmp == period) {
                       cmp = 0;
                   }
               }
               end--;
               best.updateIfBetter(period, beg, end);
           }
           final byte[] unit = Arrays.copyOfRange(unitBuffer, 0, best.period);
           pos = best.end + 1;
           emitOrDecimateSTR(seqNumber, id, best, output, decimationTable, unit);
        }
    }

    private void emitOrDecimateSTR(final int seqNumber, final String seqId, final BestPeriod best, final AutoCloseableList<BinaryTableWriter<DragstrLocus>> output, final STRDecimationTable decimationTable, final byte[] unit)
        throws IOException
    {
        final int effectiveRepeats = Math.min(maxRepeat, best.repeats);
        final int mask = nextMasks[seqNumber][best.period][effectiveRepeats]++;
        if (!decimationTable.decimate(mask, best.period, best.repeats)) {
            final DragstrLocus locus = DragstrLocus.make(seqNumber, best.start, (byte) best.period, (short) (best.end - best.start + 1), mask);
            final BinaryTableWriter<DragstrLocus> outTable = output.get(best.period - 1);
            outTable.write(locus);
            if (sortedOutputPath != null) {
                sortedOutputWriter.println("refId " + seqNumber + ", position " + (best.start - 1));
                sortedOutputWriter.println("  mask: " + mask);
                sortedOutputWriter.println("  length: " + (best.end - best.start + 1));
                sortedOutputWriter.println("  period: " + best.period + ", length " + effectiveRepeats);
            }
            progressMeter.update(new SimpleInterval(seqId, (int) best.start, (int) best.start));
            emittedCounts[best.period][effectiveRepeats]++;
        }
        totalCounts[best.period][effectiveRepeats]++;
    }

    private static class BestPeriod {
        private int period;
        private int repeats;
        private long start;
        private long end;

        public BestPeriod(final int period, final long start, final long end) {
            this.start = start;
            this.end = end;
            this.period = period;
            this.repeats = (int) (end - start + 1) / period;
        }

        private static BestPeriod initialize(final int period, final long start, final long end) {
            return new BestPeriod(period, start, end);
        }

        private boolean updateIfBetter(final int newPeriod, final long newStart, final long newEnd) {
            final int newRpeats = (int) (newEnd - newStart + 1) / newPeriod;
            if (newRpeats > repeats || (newRpeats == repeats && newPeriod < period)) {
                start = newStart;
                end = newEnd;
                period = newPeriod;
                repeats = newRpeats;
                return true;
            } else {
                return false;
            }
        }
    }
}
