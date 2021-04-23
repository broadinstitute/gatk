package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.SparkTestUtils;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;

public class ReferenceFileSparkSourceUnitTest extends GATKBaseTest {

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void testMissingReferenceFile() {
        new ReferenceFileSparkSource(
                new GATKPath(GATKBaseTest.getSafeNonExistentFile("NonExistentReference.fasta").toString()));
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

            new ReferenceFileSparkSource(new GATKPath(refPath.toUri().toString()));
        }
    }

    @Test
    public void testDatasourcesReferenceSerializes() throws IOException, ClassNotFoundException {
        Path refPath = GATKBaseTest.createTempFile("reference", ".fasta").toPath();
        GATKBaseTest.createTempFile("reference", ".fasta.fai");
        GATKBaseTest.createTempFile("reference", ".dict");

        final ReferenceFileSparkSource referenceFileSource = new ReferenceFileSparkSource(new GATKPath(refPath.toUri().toString()));

        //can we serialize it?
        ReferenceFileSparkSource otherSide = SparkTestUtils.roundTripThroughJavaSerialization(referenceFileSource);

        // After deserialization, will it crash?
        otherSide.getReferenceSequenceDictionary(null);

    }

}
