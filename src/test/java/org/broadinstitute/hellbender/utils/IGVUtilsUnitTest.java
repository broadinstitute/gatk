package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;

public class IGVUtilsUnitTest extends GATKBaseTest {

    @Test
    public void testPrintIGVFormatHeader() throws IOException {
        try ( final ByteArrayOutputStream baos = new ByteArrayOutputStream();
              final PrintStream out = new PrintStream(baos) ) {

            IGVUtils.printIGVFormatHeader(out, "line", "Track1", "Track2");

            Assert.assertEquals(baos.toString(), "#track graphType=line\nChromosome\tStart\tEnd\tFeature\tTrack1\tTrack2\n");
        }
    }

    @Test
    public void testPrintIGVFormatRow() throws IOException {
        try ( final ByteArrayOutputStream baos = new ByteArrayOutputStream();
              final PrintStream out = new PrintStream(baos) ) {

            IGVUtils.printIGVFormatRow(out, new SimpleInterval("1", 5, 10), "myFeature", 4.0, 7.0);

            // Our interval should get converted from 1-based closed form to 0-based closed-open form,
            // and our values should get formatted to 5 digits after the decimal point.
            Assert.assertEquals(baos.toString(), "1\t4\t10\tmyFeature\t4.00000\t7.00000\n");
        }
    }
}
