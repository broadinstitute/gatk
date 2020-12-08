package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.broadinstitute.hellbender.tools.variantdb.IngestConstants;

public class PetTsvWriter implements Closeable {
    private BufferedWriter writer;
    private final static char SEPARATOR = IngestConstants.SEPARATOR;

    public PetTsvWriter(String outputFile) throws IOException{
        writer = Files.newBufferedWriter(Paths.get(outputFile));

        writer.append("location");
        writer.append(SEPARATOR);
        writer.append("sample");
        writer.append(SEPARATOR);
        writer.append("state");
    }

    public void addRow(long location, long sampleId, String state) throws IOException {              
        writer.append(String.valueOf(location));
        writer.append(SEPARATOR);
        writer.append(String.valueOf(sampleId));
        writer.append(SEPARATOR);
        writer.append(state);
    }

    public void close() throws IOException {
        writer.flush();
        writer.close();
    }
}
