package org.broadinstitute.hellbender.tools.funcotator;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Class to test the {@link FuncotatorDataSourceBundler}.
 * Created by Hailey on 8/3/21.
 */

public class FuncotatorDataSourceBundlerUnitTest extends CommandLineProgramTest{

    @Test
    public void testPath() {
        String orgName = "bacteria";
        String baseURL = "http://ftp.ensemblgenomes.org/pub/current/bacteria/gtf/bacteria_collection_128/";
        String speciesName = "absiella_dolichum_dsm_3991_gca_000154285";
        Path testPath = FuncotatorDataSourceBundler.getPath(baseURL, orgName, speciesName);
        Path outputDestination = getOutputLocation(testPath);
    }

    @Test
    private Path getOutputLocation(final Path dataSourcesPath) {
        final File tmpDir = createTempDir("FuncotatorDataSourceDownloaderIntegrationTest_testDownloadDummySmallDataSources");
        final Path tmpDirPath = tmpDir.toPath();
        final Path outputDataSourcesPath = tmpDirPath.resolve(IOUtils.getPath(FuncotatorTestConstants.DUMMY_DATA_SOURCES_TAR_GZ).getFileName());
        final String outputFile = outputDataSourcesPath.toFile().getAbsolutePath();
        if ( outputFile == null ) {
            return IOUtils.getPath(dataSourcesPath.getFileName().toString());
        }
        else {
            return outputDataSourcesPath;
        }
    }
}
