package org.broadinstitute.hellbender.tools.filediagnostics;

import htsjdk.samtools.CRAMCRAIIndexer;
import htsjdk.samtools.cram.CRAIIndex;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.io.*;
import java.nio.file.Files;

/**
 * Analyzer for CRAM (.crai) index files.
 */
public class CRAIAnalyzer extends HTSAnalyzer {

    final OutputStream fos;

    public CRAIAnalyzer(final GATKPath inputPath, final GATKPath outputPath) {
        super(inputPath, outputPath);
        try {
            fos = Files.newOutputStream(outputPath.toPath());
        } catch (final IOException e) {

            throw new RuntimeIOException(e);
        }
    }

    protected void emitln(final String s) {
        try {
            fos.write(s.getBytes());
            fos.write('\n');
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Run the analyzer for the file.
     */
    protected void doAnalysis() {
        try (final InputStream is = inputPath.getInputStream()) {
            final CRAIIndex craiIndex = CRAMCRAIIndexer.readIndex(is);
            emitln("\nSeqId AlignmentStart AlignmentSpan ContainerOffset SliceOffset SliceSize\n");
            craiIndex.getCRAIEntries().stream().forEach(e -> emitln(e.toString()));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void close() throws IOException {
        if (fos != null) {
            fos.close();
        }
    }

}

