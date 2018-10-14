package org.broadinstitute.hellbender.utils.nio;

import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.nio.file.Path;

/**
 * Class for test utilities shared with the NioFileCopier* classes.
 * Created by jonn on 9/4/18.
 */
public class NioFileCopierWithProgessMeterTestUtils {

    //==================================================================================================================
    // Public Static Members:

    static Path getSourcePathFromPseudoUrl(final String source) {
        final Path sourcePath;
        if (source.startsWith("GS:")) {
            sourcePath = IOUtils.getPath(BaseTest.getGCPTestInputPath() + source.substring(3));
        }
        else {
            final String sourceBaseName = FilenameUtils.getBaseName(source);
            final String sourceExtension = FilenameUtils.getExtension(source);
            sourcePath = IOUtils.createTempFile(sourceBaseName, sourceExtension).toPath();
        }
        return sourcePath;
    }

    static Path getDestPathFromPseudoUrl(final String dest) {
        final Path destPath;
        if (dest.startsWith("GS:")) {
            final String destBaseName = FilenameUtils.getBaseName(dest.substring(3));
            final String destExtension = FilenameUtils.getExtension(dest.substring(3));
            destPath = BucketUtils.getPathOnGcs(BucketUtils.getTempFilePath(BaseTest.getGCPTestStaging() + destBaseName, destExtension));
        }
        else {
            destPath = BaseTest.getSafeNonExistentPath(dest);
        }
        return destPath;
    }
}
