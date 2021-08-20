package org.broadinstitute.hellbender.tools.gvs.ingest;

import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class RefRangesTsvWriter implements RefRangesWriter {
    private BufferedWriter writer;
    private final static char SEPARATOR = IngestConstants.SEPARATOR;

    public RefRangesTsvWriter(String outputFile) throws IOException{
        writer = Files.newBufferedWriter(Paths.get(outputFile));

        writer.append("location");
        writer.append(SEPARATOR);
        writer.append("sample_id");
        writer.append(SEPARATOR);
        writer.append("length");
        writer.append(SEPARATOR);
        writer.append("state");
        writer.append("\n");
    }

    public void write(long location, long sampleId, int length, String state) throws IOException {
        writer.append(String.valueOf(location));
        writer.append(SEPARATOR);
        writer.append(String.valueOf(sampleId));
        writer.append(SEPARATOR);
        writer.append(String.valueOf(length));
        writer.append(SEPARATOR);
        writer.append(state);
        writer.append("\n");
    }

    public void close() throws IOException {
        writer.flush();
        writer.close();
    }
}
