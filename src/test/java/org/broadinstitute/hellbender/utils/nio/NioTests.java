package org.broadinstitute.hellbender.utils.nio;

import com.google.gcloud.AuthCredentials;
import com.google.gcloud.RetryParams;
import com.google.gcloud.storage.StorageException;
import com.google.gcloud.storage.StorageOptions;
import com.google.gcloud.storage.contrib.nio.CloudStorageConfiguration;
import com.google.gcloud.storage.contrib.nio.CloudStorageFileSystem;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.nio.file.FileSystem;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

/**
 * Test GCS access via the NIO APIs.
 */
public final class NioTests extends BaseTest {

  @Test
  public void testGcsEnabled() {
    FileSystem fs = FileSystems.getFileSystem(URI.create("gs://domain-registry-alpha"));
  }

  @Test
  public void openPublicFile() throws IOException {
    InputStream inputStream = Files.newInputStream(Paths.get(URI.create("gs://pgp-harvard-data-public/hu011C57/GS000018120-DID/GS000015172-ASM/manifest.all.sig")));
    int firstByte = inputStream.read();
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
    final String privateFile = getGCPTestInputPath() + "org/broadinstitute/hellbender/utils/private_file.txt";

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
    final String privateFile = getGCPTestInputPath() + "org/broadinstitute/hellbender/utils/private_file.txt";
    final String[] split = privateFile.split("/");
    final String BUCKET = split[2];
    final String pathWithoutBucket = String.join("/", Arrays.copyOfRange(split,3, split.length));

    FileSystem fs = getAuthenticatedGcs(BUCKET);
    Path path = fs.getPath(pathWithoutBucket);
    int firstByte = Files.newInputStream(path).read();
  }

  /**
   * Using explicit credentials only works on that access, they are not kept.
   * This test will fail if default credentials are available.
   */
  @Test(enabled = false, groups = {"cloud"}, expectedExceptions = { StorageException.class })
  public void explicitCredentialsAreNotKept() throws IOException {
    // this file, potentially unlike the others in the set, is not marked as "Public link".
    final String privateFile = getGCPTestInputPath() + "org/broadinstitute/hellbender/utils/private_file.txt";
    final String privateFile2 = getGCPTestInputPath() + "org/broadinstitute/hellbender/utils/private_file_2.txt";
    final String[] split = privateFile.split("/");
    final String BUCKET = split[2];
    final String pathWithoutBucket = String.join("/", Arrays.copyOfRange(split,3, split.length));

    FileSystem fs = getAuthenticatedGcs(BUCKET);
    Path path = fs.getPath(pathWithoutBucket);
    int firstByte = Files.newInputStream(path).read();

    // now let's open another private file. It shouldn't work.
    path = Paths.get(URI.create((privateFile2)));
    firstByte = Files.newInputStream(path).read();
  }

  private FileSystem getAuthenticatedGcs(String bucket) throws IOException {
    // 1. Read credentials from disk, as per
    //    https://github.com/GoogleCloudPlatform/gcloud-java#authentication
    StorageOptions storageOptions;
    // try-with-resources to close the key file as soon as we're done
    try (InputStream keyStream = Files.newInputStream(Paths.get(getServiceAccountKeyPath()))) {
      storageOptions = StorageOptions.builder()
          .projectId(getGCPTestProject())
          .authCredentials((AuthCredentials.createForJson(keyStream)))
          // generous timeouts, to avoid tests failing when not warranted.
          .connectTimeout(60000)
          .readTimeout(60000)
          .retryParams(RetryParams.builder()
              .retryMaxAttempts(10)
              .retryMinAttempts(6)
              .maxRetryDelayMillis(30000)
              .totalRetryPeriodMillis(120000)
              .initialRetryDelayMillis(250)
              .build())
          .build();
    }

    // 2. Create GCS filesystem object with those credentials
    return CloudStorageFileSystem.forBucket(bucket, CloudStorageConfiguration.DEFAULT, storageOptions);
  }

}
