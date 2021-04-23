package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

public class OptionalTextOutputArgumentCollectionUnitTest extends GATKBaseTest {

    @Test
    public void testPrintToOutput() throws IOException {
        final OptionalTextOutputArgumentCollection out = new OptionalTextOutputArgumentCollection();
        final File tempOutput = createTempFile("testPrintToOutput", ".txt");
        out.output = new GATKPath(tempOutput.getAbsolutePath());

        out.print("test");

        final byte[] result = Files.readAllBytes(out.output.toPath());
        Assert.assertEquals(result.length, 4);
        Assert.assertEquals(result, "test".getBytes());
    }

    @Test
    public void testPrintlnToOutput() throws IOException {
        final OptionalTextOutputArgumentCollection out = new OptionalTextOutputArgumentCollection();
        final File tempOutput = createTempFile("testPrintlnToOutput", ".txt");
        out.output = new GATKPath(tempOutput.getAbsolutePath());

        out.println("test");

        final byte[] result = Files.readAllBytes(out.output.toPath());
        Assert.assertEquals(result.length, 5);
        Assert.assertEquals(result, "test\n".getBytes());
    }
}
