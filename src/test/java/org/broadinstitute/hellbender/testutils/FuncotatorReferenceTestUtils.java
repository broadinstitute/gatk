package org.broadinstitute.hellbender.testutils;

import htsjdk.samtools.util.Lazy;
import org.apache.commons.compress.archivers.ArchiveEntry;
import org.apache.commons.compress.archivers.ArchiveInputStream;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import org.apache.commons.compress.utils.IOUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.function.Supplier;

public class FuncotatorReferenceTestUtils {
    private static final Logger logger = LogManager.getLogger(FuncotatorReferenceTestUtils.class);
    private FuncotatorReferenceTestUtils() {} // Cannot instantiate this class.
    private static Lazy<String> hg19Chr3Ref;
    private static Lazy<String> hg19Chr19Ref;
    private static Lazy<String> b37Chr3Ref;
    private static Lazy<String> b37Chr2Ref;
    private static Lazy<String> hg38Chr3Ref;

    static {
        hg19Chr3Ref = new Lazy<>(createInitializer(new File(FuncotatorTestConstants.HG19_CHR3_REFERENCE_FILE_NAME)));
        hg19Chr19Ref = new Lazy<>(createInitializer(new File(FuncotatorTestConstants.HG19_CHR19_REFERENCE_FILE_NAME)));
        b37Chr3Ref = new Lazy<>(createInitializer(new File(FuncotatorTestConstants.HG19_3_REFERENCE_FILE_NAME)));
        b37Chr2Ref = new Lazy<>(createInitializer(new File(FuncotatorTestConstants.HG19_2_REFERENCE_FILE_NAME)));
        hg38Chr3Ref = new Lazy<>(createInitializer(new File(FuncotatorTestConstants.HG38_3_REFERENCE_FILE_NAME)));
    }

    private static Supplier<String> createInitializer(final File fastaTarGz) {
        return () -> FuncotatorReferenceTestUtils.extractFastaTarGzToTemp(fastaTarGz,
                org.broadinstitute.hellbender.utils.io.IOUtils.createTempDir("funcotatorTmpFolder").toPath());
    }
    /**
     * Extract a fasta tar.gz (should include only one fasta, the index, and a dict).  Though this is not enforced.
     *
     * Much of this code taken from: http://commons.apache.org/proper/commons-compress/examples.html
     *
     * @param fastaTarGz Should include one fasta file, the associated dict file and the associated index, though this is not enforced.  Never {@code null}
     * @param destDir Where to store te untarred files.  Never {@code null}
     * @return a path (as string) to the fasta file.  If multiple fasta files were in the tar gz, it will return that last one seen.
     */
    private static String extractFastaTarGzToTemp(final File fastaTarGz, final Path destDir) {
        Utils.nonNull(fastaTarGz);
        Utils.nonNull(destDir);
        String result = null;
        try (final ArchiveInputStream i = new TarArchiveInputStream(new GzipCompressorInputStream(new FileInputStream(fastaTarGz)))) {
            ArchiveEntry entry = null;
            while ((entry = i.getNextEntry()) != null) {
                if (!i.canReadEntryData(entry)) {
                    logger.warn("Could not read an entry in the archive: " + entry.getName() + " from " + fastaTarGz.getAbsolutePath());
                }
                final String name = destDir + "/" + entry.getName();
                if (entry.getName().endsWith(".fasta")) {
                    result = name;
                }
                final File f = new File(name);
                if (entry.isDirectory()) {
                    if (!f.isDirectory() && !f.mkdirs()) {
                        throw new IOException("Failed to create directory " + f);
                    }
                } else {
                    File parent = f.getParentFile();
                    if (!parent.isDirectory() && !parent.mkdirs()) {
                        throw new IOException("Failed to create directory " + parent);
                    }
                    try (OutputStream o = Files.newOutputStream(f.toPath())) {
                        IOUtils.copy(i, o);
                    }
                }
            }
        } catch (final IOException ioe) {
            throw new GATKException("Could not untar test reference: " + fastaTarGz, ioe);
        }
        return result;
    }

    /**
     * @return a path (as String) to the hg19 chromosome 3 only reference.  ("chr3")
     */
    public static String retrieveHg19Chr3Ref() {
        return hg19Chr3Ref.get();
    }

   /**
     * @return a path (as String) to the hg19 chromosome 19 only reference.  ("chr19")
     */
    public static String retrieveHg19Chr19Ref() {
        return hg19Chr19Ref.get();
    }

    /**
     * @return a path (as String) to the b37 chromosome 3 only reference.  ("3")
     */
    public static String retrieveB37Chr3Ref() { return b37Chr3Ref.get(); }

    /**
     * @return a path (as String) to the b37 chromosome 2 only reference.  ("2")
     */
    public static String retrieveB37Chr2Ref() { return b37Chr2Ref.get(); }

    /**
     * @return a path (as String) to the hg38 chromosome 3 only reference.  ("chr3" in hg38)
     */
    public static String retrieveHg38Chr3Ref() { return hg38Chr3Ref.get(); }
}
