package org.broadinstitute.hellbender.utils.nio;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import java.io.BufferedWriter;
import java.io.IOException;
import java.net.URI;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class PathLineIteratorUnitTest extends GATKBaseTest {

    private static final String[] opus = {"Hello world", "What's new?"};

    @Test
    public void testLineIteratorInForEach() throws IOException {
        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path path = jimfs.getPath("test.txt");
            try (BufferedWriter bufferedWriter = Files.newBufferedWriter(path)) {
                bufferedWriter.write(String.join("\n", opus));
            }
            ArrayList<String> got = new ArrayList<>();
            try (PathLineIterator lines = new PathLineIterator(path)) {
                for (String s : lines) {
                    got.add(s);
                }
            }
            Assert.assertEquals(got.toArray(), opus);
        }
    }

    @Test
    public void testLineIteratorAsResource() throws IOException {
        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path path = jimfs.getPath("test.txt");
            try (BufferedWriter bufferedWriter = Files.newBufferedWriter(path)) {
                bufferedWriter.write(String.join("\n", opus));
            }
            ArrayList<String> got = new ArrayList<>();
            try (PathLineIterator lines = new PathLineIterator(path)) {
                lines.forEach(s -> got.add(s));
            }
            Assert.assertEquals(got.toArray(), opus);
        }
    }

    @Test(groups = {"bucket"})
    public void testPathLineIteratorWithGCS() throws IOException {
        // this file, potentially unlike the others in the set, is not marked as "Public link".
        final String privateFilePath = "org/broadinstitute/hellbender/utils/nio/private_file.txt";
        final String privateFile = getGCPTestInputPath() + privateFilePath;

        try {
            Path path = Paths.get(URI.create((privateFile)));
            try (PathLineIterator lines = new PathLineIterator(path)) {
                int count = 0;
                for (String s : lines) {
                    count++;
                }
                Assert.assertEquals(count, 1, "Unexpected line count in " + privateFile);
            }
        } catch (Exception x) {
            System.err.println("Unable to read " + privateFile);
            helpDebugAuthError();
            throw x;
        }
    }
}