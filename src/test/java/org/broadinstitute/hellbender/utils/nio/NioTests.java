package org.broadinstitute.hellbender.utils.nio;

import com.google.gcloud.AuthCredentials;
import com.google.gcloud.RetryParams;
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

  @Test
  public void OpenPrivateFile() throws IOException {
    // this file, potentially unlike the others in the set, is not marked as "Public link".
    final String privateFile = getGCPTestInputPath() + "org/broadinstitute/hellbender/utils/private_file.txt";
    final String BUCKET = privateFile.split("/")[2];
    FileSystem fs = getAuthenticatedGcs(BUCKET);
    Path path = fs.getPath(privateFile);
    int firstByte = Files.newInputStream(path).read();
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
