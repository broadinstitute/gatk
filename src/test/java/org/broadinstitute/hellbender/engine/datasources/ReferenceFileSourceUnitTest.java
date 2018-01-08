package org.broadinstitute.hellbender.engine.datasources;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import java.io.*;
import java.nio.ByteBuffer;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.Test;

public class ReferenceFileSourceUnitTest extends GATKBaseTest {

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void testMissingReferenceFile() throws IOException {
        new ReferenceFileSource(
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

            new ReferenceFileSource(refPath);
        }
    }

    @Test
    public void testDatasourcesReferenceSerializes() throws IOException, ClassNotFoundException {
        Path refPath = GATKBaseTest.createTempFile("reference", ".fasta").toPath();
        GATKBaseTest.createTempFile("reference", ".fasta.fai");
        GATKBaseTest.createTempFile("reference", ".dict");

        final ReferenceFileSource referenceFileSource = new ReferenceFileSource(refPath);

        // Can we serialize it?
        final ByteArrayOutputStream baos = new ByteArrayOutputStream();
        ObjectOutputStream os = new ObjectOutputStream(baos);
        os.writeObject(referenceFileSource);

        ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(baos.toByteArray()));
        ReferenceFileSource otherSide = (ReferenceFileSource)is.readObject();
        // After deserialization, will it crash?
        otherSide.getReferenceSequenceDictionary(null);

    }

}
