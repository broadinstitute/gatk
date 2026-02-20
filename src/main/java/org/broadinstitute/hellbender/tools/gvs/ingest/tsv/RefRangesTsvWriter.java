package org.broadinstitute.hellbender.tools.gvs.ingest.tsv;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.ingest.RefRangesWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class RefRangesTsvWriter implements RefRangesWriter {
    private final BufferedWriter writer;
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

    @Override
    public void writeCompressed(long packedData, long sampleId) throws IOException {
        throw new GATKException("TSVWriter doesn't support compressed reference blocks");
    }

    public void close() throws IOException {
        writer.flush();
        writer.close();
    }
}
