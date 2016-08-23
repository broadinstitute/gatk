package org.broadinstitute.hellbender.engine.spark.datasources;

import org.testng.annotations.Test;

import java.io.IOException;
import java.net.URI;
import java.nio.file.Paths;

public class NioProviderExceptionUnitTest {
    @Test()
    public void test() throws IOException {
        Paths.get(URI.create("hdfs://nn/path"));
    }
}
