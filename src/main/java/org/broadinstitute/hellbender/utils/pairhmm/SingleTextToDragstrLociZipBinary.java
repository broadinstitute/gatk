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
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.*;
import java.nio.file.Paths;
import java.util.Arrays;

@CommandLineProgramProperties(
        programGroup = ReferenceProgramGroup.class,
        summary = "Determine the presence of STR in a reference sequence",
        oneLineSummary = "Determines the presence of STR in a reference sequence"
)
public class SingleTextToDragstrLociZipBinary extends GATKTool {

    private static final Logger logger = LogManager.getLogger(SingleTextToDragstrLociZipBinary.class);

    @Argument(doc = "text", fullName = "text")
    private String inputPath = null;

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
        try (final TableReader<DragstrLocus> inputReader = DragstrLocus.textReader(BucketUtils.openFile(inputPath), dictionary);
                final AutoCloseableList<BinaryTableWriter<DragstrLocus>> allSitesWriter
                     = AutoCloseableList.of(maxPeriod, i -> createDragstrLocusWriter(tempDir, "period_" + (i + 1) + ".bin"))) {
            inputReader.stream().forEach(
                    locus -> {try {
                allSitesWriter.get(locus.getPeriod() - 1).write(locus);
            } catch (final IOException ex) {
                        throw new UncheckedIOException(ex);
            }});
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
        writer.println("period\trepeat\temmitted");
        for (int i = 1; i <= maxPeriod; i++) {
            for (int j = 1; j <= maxRepeat; j++) {
                writer.println(String.format("%d\t%d\t%d",
                        i, j, totalCounts[i][j]));
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
}
