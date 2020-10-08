package org.broadinstitute.hellbender.utils.dragstr;

import com.google.common.io.Files;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.dragstr.DragstrLocus;
import org.broadinstitute.hellbender.tools.dragstr.DragstrLocusUtils;
import org.broadinstitute.hellbender.utils.BinaryTableWriter;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.ZipUtils;
import org.broadinstitute.hellbender.tools.dragstr.STRDecimationTable;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.*;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Utility class to compose the contents of the STR Table file.
 * <p>
 *     The builder uses a temporary folder to compose the content of the final str-table-file.
 * </p>
 * <p>
 *     In order to create an actual str-table-file, once all the information has been added to the builder,
 *     you need to invoke {@link #store}.
 * </p>
 * <p>
 *     In order to clear resources after the work is done you need to invoke {@link #close}.
 * </p>
 */
public final class STRTableFileBuilder implements AutoCloseable {
    private boolean closed;
    private final File dir;
    private final Map<String, String> annotations;
    private int maxPeriod;
    private int maxRepeatLength;
    private final BinaryTableWriter<DragstrLocus> sitesWriter;
    private final TableWriter<DragstrLocus>     textSitesWriter;
    private final SAMSequenceDictionary dictionary;
    private final STRDecimationTable decimationTable;
    private final long[][] emittedCounts;
    private final long[][] totalCounts;

    private STRTableFileBuilder(final File dir, final boolean generateTextSitesFile, final SAMSequenceDictionary dictionary,
                        final STRDecimationTable decimationTable, final int maxPeriod, final int maxRepeatLength) {
        this.maxPeriod = maxPeriod;
        this.maxRepeatLength = maxRepeatLength;
        this.emittedCounts = new long[maxPeriod + 1][maxRepeatLength + 1];
        this.totalCounts = new long[maxPeriod + 1][maxRepeatLength + 1];
        this.dir = dir;
        this.dictionary = dictionary;
        this.decimationTable = decimationTable;
        this.annotations = new LinkedHashMap<>();
        try {
            sitesWriter = DragstrLocusUtils.binaryWriter(new File(dir, STRTableFile.SITES_FILE_NAME), new File(dir, STRTableFile.SITES_INDEX_FILE_NAME));
        } catch (final FileNotFoundException ex) {
            throw new GATKException("possible bug, the parent directory " + dir + " must exists at this point", ex);
        }
        try {
            textSitesWriter = generateTextSitesFile ? DragstrLocusUtils.textWriter(new FileOutputStream(new File(dir, STRTableFile.SITES_TEXT_FILE_NAME)), dictionary) : null;
        } catch (IOException e) {
            throw new GATKException("possible bug", e);
        }
        writeReferenceDictionary(dir, dictionary);
        writeDecimationTable(dir, decimationTable);
    }

    /**
     * Creates a new builder.
     * @return never {@code null}.
     */
    public static STRTableFileBuilder newInstance(final SAMSequenceDictionary dictionary, final STRDecimationTable decimationTable, final boolean generateTextSitesFile, final int maxPeriod, final int maxRepeatLength) {
        Utils.validateArg(maxPeriod >= 1, "max period must be positive");
        Utils.validateArg(maxRepeatLength >= 1, "max repeat length must be positive");
        Utils.nonNull(decimationTable, "decimation table must not be negative");
        Utils.nonNull(dictionary, "dictionary must not be negative");
        final File tempDir = Files.createTempDir();
        return new STRTableFileBuilder(tempDir, generateTextSitesFile, dictionary, decimationTable, maxPeriod, maxRepeatLength);
    }

    /**
     * Add an annotation that would go to the header of the summary file.
     * @param name the annotation name.
     * @param value the annotation value.
     */
    public void annotate(final String name, final String value) {
        Utils.nonNull(name);
        Utils.nonNull(value);
        annotations.put(name, value);
    }

    private static void writeReferenceDictionary(final File dir, final SAMSequenceDictionary dictionary) {
        final File dictionaryFile = new File(dir, STRTableFile.REF_DICTIONARY_FILE_NAME);
        try (final Writer dictWriter = new PrintWriter(new FileWriter(dictionaryFile))) {
            final SAMSequenceDictionaryCodec codec = new SAMSequenceDictionaryCodec(dictWriter);
            codec.encode(dictionary);
        } catch (final IOException e) {
            throw new GATKException("issues writing dictionary file in stage directory " + dir, e);
        }

    }

    private static void writeDecimationTable(final File dir, final STRDecimationTable decimationTable) {
        final File decimationTableFile = new File(dir, STRTableFile.DECIMATION_TABLE_FILE_NAME);
        try (final PrintWriter deciWriter = new PrintWriter(new FileWriter(decimationTableFile))) {
            decimationTable.print(deciWriter);
        } catch (final IOException e) {
            throw new GATKException("issues writing dictionary file in stage directory " + dir, e);
        }
    }

    /**
     * Note that a site with a period and repeat-length has been decimated.
     * <p>
     *     The str table file (writer) would updated counts and summary accordingly.
     * </p>
     */
    public void decimate(final int period, final int repeatLength) {
        checkIsNotClosed();
        final int effectiveRepeatLength = Math.min(maxRepeatLength, repeatLength);
        final int effectivePeriod = Math.min(maxPeriod, period);
        totalCounts[effectivePeriod][effectiveRepeatLength]++;
    }

    /**
     * Emits a locus in the str table.
     * @param locus the locus to emit.
     * @throws GATKException if any low-level issue  occurred while emitting the locus.
     */
    public void emit(final DragstrLocus locus) throws GATKException {
        checkIsNotClosed();
        checkLocusIsValid(locus);
        final int effectiveRepeatLength = Math.min(maxRepeatLength, locus.getRepeats());
        final int effectivePeriod = Math.min(maxPeriod, locus.getPeriod());
        totalCounts[effectivePeriod][effectiveRepeatLength]++;
        emittedCounts[effectivePeriod][effectiveRepeatLength]++;
        try {
            sitesWriter.write(locus);
            if (textSitesWriter != null) {
                textSitesWriter.writeRecord(locus);
            }
        } catch (final IOException ex) {
            throw new GATKException("issues writing loci to the staging files in " + dir, ex);
        }
    }

    private void checkLocusIsValid(final DragstrLocus locus) {
        Utils.nonNull(locus, "the locus cannot be null");
        final SAMSequenceRecord seq = dictionary.getSequence(locus.getChromosomeIndex());
        Utils.nonNull(seq, "the locus chr idx is out of range");
        Utils.validateArg(locus.getStart() >= 1, "the start coordinate must be positive");
        Utils.validateArg(locus.getEnd() <= seq.getSequenceLength(), "the end position is beyond the seq's end");
    }

    private void writeSummary() {
        final File summaryFile = new File(dir, STRTableFile.SUMMARY_FILE_NAME);
        try (final PrintWriter writer = new PrintWriter(new FileWriter(summaryFile))) {
            writer.println("##########################################################################################");
            writer.println("# STRTableSummary");
            writer.println("# ---------------------------------------");
            writer.println("# maxPeriod = " + maxPeriod);
            writer.println("# maxRepeatLength = " + maxRepeatLength);
            for (final String name : annotations.keySet()) {
                writer.println("# " + name + " = " + annotations.get(name));
            }
            writer.println("##########################################################################################");
            writer.println(String.join("\t", "period", "repeatLength", "totalCounts", "emittedCounts", "intendedDecimation", "actualDecimation"));
            for (int period = 1; period <= maxPeriod; period++) {
                for (int repeatLength = period == 1 ? 1 : 2; repeatLength <= maxRepeatLength; repeatLength++) {
                    final long total = totalCounts[period][repeatLength];
                    final long emitted = emittedCounts[period][repeatLength];
                    final int decimation = decimationTable.decimationBit(period, repeatLength);
                    final double actualDecimation = total > 0 ? (MathUtils.INV_LOG_2 * (Math.log(total) - Math.log(emitted))): 0;
                    writer.println(Utils.join("\t", period, repeatLength, total, emitted, decimation, Math.round(actualDecimation * 100) / 100.0));
                }
            }
        } catch (final IOException e) {
            throw new GATKException("unexpected issues writing summary file in " + dir);
        }
    }

    public void store(final GATKPath path) {
        checkIsNotClosed();
        try {
            sitesWriter.flush();
            if (textSitesWriter != null) textSitesWriter.flush();
        } catch (final IOException ex) {
            throw new GATKException("problems flushing the str-table-file content to " + dir, ex);
        }
        writeSummary();
        try {
            ZipUtils.zip(dir, path);
        } catch (final GATKException ex) {
            throw new GATKException("problems flushing the str-table-file content from " + dir + " to " + path, ex.getCause());
        }
    }

    private void checkIsNotClosed() {
        if (closed) {
            throw new IllegalStateException("the writer is already closed");
        }
    }

    public void close() {
        if (!closed) {
            closed = true;
            try {
                if (sitesWriter != null) sitesWriter.close();
                if (textSitesWriter != null) textSitesWriter.close();
                if (dir.exists()) FileUtils.deleteDirectory(dir);
            } catch (final IOException ex) {
                throw new GATKException("issues finishing writing the sites files in the stage directory " + dir);
            }
        }
    }
}
