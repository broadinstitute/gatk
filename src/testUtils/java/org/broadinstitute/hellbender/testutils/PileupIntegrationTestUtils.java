package org.broadinstitute.hellbender.testutils;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

public final class PileupIntegrationTestUtils {

    public static File createTempFile() throws IOException {
        final File out = File.createTempFile("out", ".txt");
        out.delete();
        out.deleteOnExit();
        return out;
    }

    public static File createFeaturesPileupTestFile(String sourcePath, String templatePath) throws IOException { ;
        final File filledFile = createTempFile();
        String template = new String(Files.readAllBytes(Paths.get(templatePath)), StandardCharsets.UTF_8);
        template = template.replaceAll("&source", sourcePath);
        Files.write(Paths.get(filledFile.getPath()), template.getBytes(StandardCharsets.UTF_8));
        return filledFile;
    }

}
