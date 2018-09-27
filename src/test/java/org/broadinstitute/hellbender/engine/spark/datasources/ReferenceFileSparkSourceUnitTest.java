package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import java.io.*;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceFileSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.Test;

public class ReferenceFileSparkSourceUnitTest extends GATKBaseTest {

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void testMissingReferenceFile() throws IOException {
        new ReferenceFileSparkSource(
                GATKBaseTest.getSafeNonExistentFile("NonExistentReference.fasta")
                        .getAbsolutePath());
    }

    @Test
    public void testDatasourcesReferenceInPath() throws IOException {
        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path refPath = jimfs.getPath("reference.fasta");
            Files.createFile(refPath);
            final Path faiFile = jimfs.getPath("reference.fasta.fai");
            Files.createFile(faiFile);
            final Path dictPath = jimfs.getPath("reference.dict");
            Files.createFile(dictPath);

            new ReferenceFileSparkSource(refPath);
        }
    }

    @Test
    public void testDatasourcesReferenceSerializes() throws IOException, ClassNotFoundException {
        Path refPath = GATKBaseTest.createTempFile("reference", ".fasta").toPath();
        GATKBaseTest.createTempFile("reference", ".fasta.fai");
        GATKBaseTest.createTempFile("reference", ".dict");

        final ReferenceFileSparkSource referenceFileSource = new ReferenceFileSparkSource(refPath);

        // Can we serialize it?
        final ByteArrayOutputStream baos = new ByteArrayOutputStream();
        ObjectOutputStream os = new ObjectOutputStream(baos);
        os.writeObject(referenceFileSource);

        ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(baos.toByteArray()));
        ReferenceFileSparkSource otherSide = (ReferenceFileSparkSource)is.readObject();
        // After deserialization, will it crash?
        otherSide.getReferenceSequenceDictionary(null);

    }

}
