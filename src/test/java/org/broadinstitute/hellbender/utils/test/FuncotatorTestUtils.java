package org.broadinstitute.hellbender.utils.test;

import org.apache.commons.compress.archivers.ArchiveEntry;
import org.apache.commons.compress.archivers.ArchiveException;
import org.apache.commons.compress.archivers.ArchiveInputStream;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import org.apache.commons.compress.utils.IOUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;

public class FuncotatorTestUtils {

    private FuncotatorTestUtils() {}

    /**
     * TODO: Finish docs (LTL)
     *
     * Much of this code taken from: http://commons.apache.org/proper/commons-compress/examples.html
     *
     * @param fastaTarGz Should include the dict file and index, though this is not enforced.  Never {@code null}
     * @param destDir
     * @return
     */
    public static String extractFastaTarGzToTemp(final File fastaTarGz, final Path destDir) throws IOException, ArchiveException {
        Utils.nonNull(fastaTarGz);
        String result = null;
        try (final ArchiveInputStream i = new TarArchiveInputStream(new GzipCompressorInputStream(new FileInputStream(fastaTarGz)))) {
            ArchiveEntry entry = null;
            while ((entry = i.getNextEntry()) != null) {
                if (!i.canReadEntryData(entry)) {
                    // log something?
                    continue;
                }
                final String name = destDir + "/" + entry.getName();
                if (entry.getName().endsWith(".fasta")) {
                    result = name;
                }
                final File f = new File(name);
                if (entry.isDirectory()) {
                    if (!f.isDirectory() && !f.mkdirs()) {
                        throw new IOException("failed to create directory " + f);
                    }
                } else {
                    File parent = f.getParentFile();
                    if (!parent.isDirectory() && !parent.mkdirs()) {
                        throw new IOException("failed to create directory " + parent);
                    }
                    try (OutputStream o = Files.newOutputStream(f.toPath())) {
                        IOUtils.copy(i, o);
                    }
                }
            }
        }
        return result;
    }
}
