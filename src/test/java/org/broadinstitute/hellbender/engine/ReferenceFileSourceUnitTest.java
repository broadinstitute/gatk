package org.broadinstitute.hellbender.engine;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.GATKBaseTest;

import org.testng.annotations.Test;

import java.io.IOException;

public class ReferenceFileSourceUnitTest extends GATKBaseTest {

    @Test(expectedExceptions = UserException.MissingReferenceFaiFile.class)
    public void testReferenceInPathWithoutFai() throws IOException {
        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path refPath = jimfs.getPath("reference.fasta");
            Files.createFile(refPath);
            final Path dictPath = jimfs.getPath("reference.dict");
            Files.createFile(dictPath);

            new ReferenceFileSource(refPath);
        }
    }

    @Test(expectedExceptions = UserException.MissingReferenceDictFile.class)
    public void testReferenceInPathWithoutDict() throws IOException {
        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path refPath = jimfs.getPath("reference.fasta");
            Files.createFile(refPath);
            final Path faiFile = jimfs.getPath("reference.fasta.fai");
            Files.createFile(faiFile);

            new ReferenceFileSource(refPath);
        }
    }

    @Test
    public void testReferenceInPath() throws IOException {
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

}
