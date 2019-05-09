package org.broadinstitute.hellbender.utils.tsv;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;

import static org.testng.Assert.*;

public class SimpleCSVWriterWrapperWithHeaderUnitTest extends GATKBaseTest {


    @Test (groups = "bucket")
    public void testWriteToBucketPath() throws IOException {
        Path bucketPath = IOUtils.getPath(BucketUtils.getTempFilePath(
                getGCPTestStaging() +"testWriteToBucketPath", ".tsv"));
        Path localPath = IOUtils.getPath(createTempFile("testWriteToBucketPath",  ".tsv").getPath() );

        SimpleCSVWriterWrapperWithHeader bucketWriter = new SimpleCSVWriterWrapperWithHeader(bucketPath, '\t');
        SimpleCSVWriterWrapperWithHeader localWriter = new SimpleCSVWriterWrapperWithHeader(localPath, '\t');

        String[] header = new String[]{"a","b","c"};

        for (int i = 0; i < 100; i++) {
            SimpleCSVWriterWrapperWithHeader.SimpleCSVWriterLineBuilder bucketLine = bucketWriter.getNewLineBuilder();
            SimpleCSVWriterWrapperWithHeader.SimpleCSVWriterLineBuilder localLine = localWriter.getNewLineBuilder();
            Arrays.stream(header).forEach(column -> {
                double rand = Math.random();
                bucketLine.setColumn(column, Double.toString(rand));
                localLine.setColumn(column, Double.toString(rand));
            });
        }

        bucketWriter.close();
        localWriter.close();

        IntegrationTestSpec.assertEqualTextFiles(bucketPath, localPath, null);
    }

}