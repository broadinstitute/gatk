package org.broadinstitute.hellbender.utils.nio;

import com.google.cloud.storage.StorageException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.nio.ByteBuffer;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.FileSystem;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Random;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;

/**
 * Test GCS access via the NIO APIs.
 */
public final class GcsNioIntegrationTest extends GATKBaseTest {

    final String privateFilePath = "org/broadinstitute/hellbender/utils/nio/private_file.txt";
    final String privateFilePath2 = "org/broadinstitute/hellbender/utils/nio/private_file_2.txt";
    final String largeFilePath = "large/human_g1k_v37.20.21.fasta";

    @Test(groups={"bucket"})
    public void testGcsEnabled() {
        FileSystem fs = FileSystems.getFileSystem(URI.create("gs://domain-registry-alpha"));
    }

    @Test(groups={"bucket"})
    public void openPublicFile() throws IOException {
        try (InputStream inputStream = Files.newInputStream(Paths.get(URI.create(
                "gs://pgp-harvard-data-public/hu011C57/GS000018120-DID/GS000015172-ASM/manifest.all.sig")))) {
            int firstByte = inputStream.read();
        }
    }

    /**
     * Ensure we can write to the staging folder. If this fails, it may indicate a misconfiguration.
     *
     * @throws IOException
     */
    @Test(groups={"bucket"})
    public void writePrivateFile() throws IOException {
        try {
            final String dest = BucketUtils.getTempFilePath(
                    getGCPTestStaging() +"GcsNioIntegrationTest-writePrivateFile-test", ".txt");
            final Path outputPath = BucketUtils.getPathOnGcs(dest);
            System.out.println("Writing to " + dest);
            try (OutputStream os = Files.newOutputStream(outputPath)) {
                os.write(42);
            }
        } catch (shaded.cloud_nio.com.google.api.client.http.HttpResponseException forbidden) {
            helpDebugAuthError();
            throw forbidden;
        }
    }


    /**
     * When we give no explicit credentials, then the system will still work if either
     * - the GOOGLE_APPLICATION_CREDENTIALS environment variable is set
     * - the user ran "gcloud auth login" first
     * - the code is running on a Google Cloud Compute machine.
     * (see http://gcloud-python.readthedocs.org/en/latest/gcloud-auth.html)
     */
    @Test(groups = {"bucket"})
    public void openPrivateFileUsingDefaultCredentials() throws IOException {
        // this file, potentially unlike the others in the set, is not marked as "Public link".
        final String privateFile = getGCPTestInputPath() + privateFilePath;

        try {
            Path path = Paths.get(URI.create((privateFile)));
            int firstByte = Files.newInputStream(path).read();
        } catch (Exception x) {
            System.err.println("Unable to open " + privateFile);
            helpDebugAuthError();
            throw x;
        }
    }

    // TODO(jpmartin):uncomment once getAuthenticatedGcs is back
    /**
     * Opening the private file even when the user is not logged in on gcloud should work
     * when we provide explicit credentials.
     *
    @Test(groups = {"bucket"})
    public void openPrivateFileWithExplicitCredentials() throws IOException {
        // this file, potentially unlike the others in the set, is not marked as "Public link".
        final String privateFile = getGCPTestInputPath() + privateFilePath;
        final String BUCKET = BucketUtils.getBucket(privateFile);
        final String pathWithoutBucket = BucketUtils.getPathWithoutBucket(privateFile);

        try {
            FileSystem fs = getAuthenticatedGcs(BUCKET);
            Path path = fs.getPath(pathWithoutBucket);
            int firstByte = Files.newInputStream(path).read();
        } catch (Exception x) {
            System.err.println("Unable to open " + privateFile);
            helpDebugAuthError();
            throw x;
        }
    }

    /**
     * Using explicit credentials only works on that access, they are not kept.
     * This test will fail if default credentials are available.
     * That means that you must NOT set $GOOGLE_APPLICATION_CREDENTIALS
     * Yet you must set getGoogleServiceAccountKeyPath() (you may have to switch it to a different
     * environment variable).
     **/
    @Test(enabled = false, groups = {"bucket"}, expectedExceptions = {StorageException.class})
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


    @Test(groups = {"bucket"})
    public void testCloseWhilePrefetching() throws Exception {
        final String large = getGCPTestInputPath() + largeFilePath;
        SeekableByteChannel chan = new SeekableByteChannelPrefetcher(
            Files.newByteChannel(Paths.get(URI.create(large))), 10*1024*1024);
        // read just 1 byte, get the prefetching going
        ByteBuffer one = ByteBuffer.allocate(1);
        chan.read(one);
        // closing must not throw an exception, even if the prefetching
        // thread is active.
        chan.close();
    }

}
