package org.broadinstitute.hellbender.utils.nio;

import com.google.cloud.AuthCredentials;
import com.google.cloud.RetryParams;
import com.google.cloud.storage.StorageException;
import com.google.cloud.storage.StorageOptions;
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration;
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.nio.file.*;
import java.util.Arrays;

/**
 * Test GCS access via the NIO APIs.
 */
public final class GcsNioIntegrationTest extends BaseTest {

    final String privateFilePath = "org/broadinstitute/hellbender/utils/nio/private_file.txt";
    final String privateFilePath2 = "org/broadinstitute/hellbender/utils/nio/private_file_2.txt";

    @Test
    public void testGcsEnabled() {
        FileSystem fs = FileSystems.getFileSystem(URI.create("gs://domain-registry-alpha"));
    }

    @Test
    public void openPublicFile() throws IOException {
        try (InputStream inputStream = Files.newInputStream(Paths.get(URI.create(
                "gs://pgp-harvard-data-public/hu011C57/GS000018120-DID/GS000015172-ASM/manifest.all.sig")))) {
            int firstByte = inputStream.read();
        }
    }


    /**
     * When we give no explicit credentials, then the system will still work if either
     * - the GOOGLE_APPLICATION_CREDENTIALS environment variable is set
     * - the user ran "gcloud auth login" first
     * - the code is running on a Google Cloud Compute machine.
     * (see http://gcloud-python.readthedocs.org/en/latest/gcloud-auth.html)
     */
    @Test(groups = {"cloud"})
    public void openPrivateFileUsingDefaultCredentials() throws IOException {
        // this file, potentially unlike the others in the set, is not marked as "Public link".
        final String privateFile = getGCPTestInputPath() + privateFilePath;

        Path path = Paths.get(URI.create((privateFile)));
        int firstByte = Files.newInputStream(path).read();
    }

    /**
     * Opening the private file even when the user is not logged in on gcloud should work
     * when we provide explicit credentials.
     */
    @Test(groups = {"cloud"})
    public void openPrivateFileWithExplicitCredentials() throws IOException {
        // this file, potentially unlike the others in the set, is not marked as "Public link".
        final String privateFile = getGCPTestInputPath() + privateFilePath;
        final String BUCKET = BucketUtils.getBucket(privateFile);
        final String pathWithoutBucket = BucketUtils.getPathWithoutBucket(privateFile);

        FileSystem fs = getAuthenticatedGcs(BUCKET);
        Path path = fs.getPath(pathWithoutBucket);
        int firstByte = Files.newInputStream(path).read();
    }

    /**
     * Using explicit credentials only works on that access, they are not kept.
     * This test will fail if default credentials are available.
     * That means that you must NOT set $GOOGLE_APPLICATION_CREDENTIALS
     * Yet you must set getGoogleServiceAccountKeyPath() (you may have to switch it to a different
     * environment variable).
     */
    @Test(enabled = false, groups = {"cloud"}, expectedExceptions = {StorageException.class})
    public void explicitCredentialsAreNotKept() throws IOException {
        // this file, potentially unlike the others in the set, is not marked as "Public link".
        final String privateFile = getGCPTestInputPath() + privateFilePath;
        final String privateFile2 = getGCPTestInputPath() + privateFilePath2;
        final String BUCKET = BucketUtils.getBucket(privateFile);
        final String pathWithoutBucket = BucketUtils.getPathWithoutBucket(privateFile);

        FileSystem fs = getAuthenticatedGcs(BUCKET);
        Path path = fs.getPath(pathWithoutBucket);
        int firstByte = Files.newInputStream(path).read();

        // now let's open another private file. It shouldn't work.
        path = Paths.get(URI.create((privateFile2)));
        firstByte = Files.newInputStream(path).read();
    }

    private FileSystem getAuthenticatedGcs(String bucket) throws IOException {
        byte[] creds = Files.readAllBytes(Paths.get(getGoogleServiceAccountKeyPath()));
        return BucketUtils.getAuthenticatedGcs(getGCPTestProject(), bucket, creds);
    }


}
