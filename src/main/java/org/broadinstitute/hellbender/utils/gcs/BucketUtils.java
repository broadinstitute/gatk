package org.broadinstitute.hellbender.utils.gcs;

import com.google.cloud.http.HttpTransportOptions;
import com.google.cloud.storage.StorageOptions;
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration;
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration.Builder;
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystemProvider;
import com.google.common.base.Strings;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.util.TabixUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import shaded.cloud_nio.com.google.api.gax.retrying.RetrySettings;
import shaded.cloud_nio.com.google.auth.oauth2.GoogleCredentials;
import shaded.cloud_nio.org.threeten.bp.Duration;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.UUID;

/**
 * Utilities for dealing with google buckets.
 */
public final class BucketUtils {
    public static final String GCS_SCHEME = "gs";
    public static final String GCS_PREFIX = GCS_SCHEME + "://";
    public static final String HDFS_PREFIX = "hdfs://";

    // slashes omitted since hdfs paths seem to only have 1 slash which would be weirder to include than no slashes
    public static final String FILE_PREFIX = "file:";

    public static final Logger logger = LogManager.getLogger("org.broadinstitute.hellbender.utils.gcs");

    private BucketUtils(){} //private so that no one will instantiate this class

    public static boolean isCloudStorageUrl(final String path) {
        Utils.nonNull(path);
        return path.startsWith(GCS_PREFIX);
    }

    /**
     * Return true if this {@code GATKPathSpecifier} represents a gcs URI.
     * @param pathSpec specifier to inspect
     * @return true if this {@code GATKPathSpecifier} represents a gcs URI.
     */
    public static boolean isCloudStorageUrl(final GATKPathSpecifier pathSpec) {
        Utils.nonNull(pathSpec);
        return pathSpec.getScheme().equals(GCS_SCHEME);
    }

    public static boolean isCloudStorageUrl(final Path path) {
        // the initial "" protects us against a null scheme
        return ("" + path.toUri().getScheme() + "://").equals(GCS_PREFIX);
    }

    /**
     * Returns true if the given path is a HDFS (Hadoop filesystem) URL.
     */
    public static boolean isHadoopUrl(String path) {
        return path.startsWith(HDFS_PREFIX);
    }

    /**
     * Returns true if the given path is a GCS or HDFS (Hadoop filesystem) URL.
     */
    public static boolean isRemoteStorageUrl(String path) {
        return isCloudStorageUrl(path) || isHadoopUrl(path);
    }

    /**
     * Changes relative local file paths to be absolute file paths. Paths with a scheme are left unchanged.
     * @param path the path
     * @return an absolute file path if the original path was a relative file path, otherwise the original path
     */
    public static String makeFilePathAbsolute(String path){
        if (isCloudStorageUrl(path) || isHadoopUrl(path) || isFileUrl(path)){
            return path;
        } else {
            return new File(path).getAbsolutePath();
        }
    }

    /**
     * Get a temporary file path based on the prefix and extension provided.
     * This file (and possible indexes associated with it) will be scheduled for deletion on shutdown
     *
     * @param prefix a prefix for the file name
     *               for remote paths this should be a valid URI to root the temporary file in (ie. gs://hellbender/staging/)
     *               there is no guarantee that this will be used as the root of the tmp file name, a local prefix may be placed in the tmp folder for example
     * @param extension and extension for the temporary file path, the resulting path will end in this
     * @return a path to use as a temporary file, on remote file systems which don't support an atomic tmp file reservation a path is chosen with a long randomized name
     *
     */
    public static String getTempFilePath(String prefix, String extension){
        if (isCloudStorageUrl(prefix) || (isHadoopUrl(prefix))){
            final String path = randomRemotePath(prefix, "", extension);
            IOUtils.deleteOnExit(IOUtils.getPath(path));
            IOUtils.deleteOnExit(IOUtils.getPath(path + Tribble.STANDARD_INDEX_EXTENSION));
            IOUtils.deleteOnExit(IOUtils.getPath(path + TabixUtils.STANDARD_INDEX_EXTENSION));
            IOUtils.deleteOnExit(IOUtils.getPath(path + ".bai"));
            IOUtils.deleteOnExit(IOUtils.getPath(path + ".md5"));
            IOUtils.deleteOnExit(IOUtils.getPath(path.replaceAll(extension + "$", ".bai"))); //if path ends with extension, replace it with .bai
            return path;
        } else {
            return IOUtils.createTempFile(prefix, extension).getAbsolutePath();
        }
    }

    /**
     * Picks a random name, by putting some random letters between "prefix" and "suffix".
     *
     * @param stagingLocation The folder where you want the file to be. Must start with "gs://" or "hdfs://"
     * @param prefix The beginning of the file name
     * @param suffix The end of the file name, e.g. ".tmp"
     */
    public static String randomRemotePath(String stagingLocation, String prefix, String suffix) {
        if (isCloudStorageUrl(stagingLocation)) {
            // Go through URI because Path.toString isn't guaranteed to include the "gs://" prefix.
            return getPathOnGcs(stagingLocation).resolve(prefix + UUID.randomUUID().toString() + suffix).toUri().toString();
        } else if (isHadoopUrl(stagingLocation)) {
            return new org.apache.hadoop.fs.Path(stagingLocation, prefix + UUID.randomUUID().toString() + suffix).toString();
        } else {
            throw new IllegalArgumentException("Staging location is not remote: " + stagingLocation);
        }
    }

    public static boolean isFileUrl(String path) {
        return path.startsWith(FILE_PREFIX);
    }

    /**
     * Given a path of the form "gs://bucket/folder/folder/file", returns "bucket".
     */
    public static String getBucket(String path) {
        return path.split("/")[2];
    }

    /**
     * Given a path of the form "gs://bucket/folder/folder/file", returns "folder/folder/file".
     */
    public static String getPathWithoutBucket(String path) {
        final String[] split = path.split("/");
        return String.join("/", Arrays.copyOfRange(split, 3, split.length));
    }

    /**
     * Sets max_reopens, requester_pays, and generous timeouts as the global default.
     * These will apply even to library code that creates its own paths to access with NIO.
     *
     * @param maxReopens If the GCS bucket channel errors out, how many times it will attempt to
     *                   re-initiate the connection.
     * @param requesterProject Project to bill when accessing "requester pays" buckets. If unset,
     *                         these buckets cannot be accessed.
     */
    public static void setGlobalNIODefaultOptions(int maxReopens, String requesterProject) {
        CloudStorageFileSystemProvider.setDefaultCloudStorageConfiguration(getCloudStorageConfiguration(maxReopens, requesterProject));
        CloudStorageFileSystemProvider.setStorageOptions(setGenerousTimeouts(StorageOptions.newBuilder()).build());
    }

    /**
     * String -> Path. This *should* not be necessary (use Paths.get(URI.create(...)) instead) , but it currently is
     * on Spark because using the fat, shaded jar breaks the registration of the GCS FilesystemProvider.
     * To transform other types of string URLs into Paths, use IOUtils.getPath instead.
     */
    public static Path getPathOnGcs(String gcsUrl) {
        // use a split limit of -1 to preserve empty split tokens, especially trailing slashes on directory names
        final String[] split = gcsUrl.split("/", -1);
        final String BUCKET = split[2];
        final String pathWithoutBucket = String.join("/", Arrays.copyOfRange(split, 3, split.length));
        return CloudStorageFileSystem.forBucket(BUCKET).getPath(pathWithoutBucket);
    }

    /**
     * The config we want to use.
     *
     * @param maxReopens If the GCS bucket channel errors out, how many times it will attempt to
     *                   re-initiate the connection.
     * @param requesterProject Project to bill when accessing "requester pays" buckets. If unset,
     *                         these buckets cannot be accessed.
     *
     **/
    public static CloudStorageConfiguration getCloudStorageConfiguration(int maxReopens, String requesterProject) {
        Builder builder = CloudStorageConfiguration.builder()
            // if the channel errors out, re-open up to this many times
            .maxChannelReopens(maxReopens);
        if (!Strings.isNullOrEmpty(requesterProject)) {
            // enable requester pays and indicate who pays
            builder = builder.autoDetectRequesterPays(true).userProject(requesterProject);
        }

        //this causes the gcs filesystem to treat files that end in a / as a directory
        //true is the default but this protects against future changes in behavior
        builder.usePseudoDirectories(true);
        return builder.build();
    }

    private static StorageOptions.Builder setGenerousTimeouts(StorageOptions.Builder builder) {
        return builder
            .setTransportOptions(HttpTransportOptions.newBuilder()
                .setConnectTimeout(120000)
                .setReadTimeout(120000)
                .build())
            .setRetrySettings(RetrySettings.newBuilder()
                .setMaxAttempts(15)
                .setMaxRetryDelay(Duration.ofMillis(256_000L))
                .setTotalTimeout(Duration.ofMillis(4000_000L))
                .setInitialRetryDelay(Duration.ofMillis(1000L))
                .setRetryDelayMultiplier(2.0)
                .setInitialRpcTimeout(Duration.ofMillis(180_000L))
                .setRpcTimeoutMultiplier(1.0)
                .setMaxRpcTimeout(Duration.ofMillis(180_000L))
                .build());
    }

    /**
     * Get an authenticated GCS-backed NIO FileSystem object representing the selected projected and bucket.
     * Credentials are found automatically when running on Compute/App engine, logged into gcloud, or
     * if the GOOGLE_APPLICATION_CREDENTIALS env. variable is set. In that case leave credentials null.
     * Otherwise, you must pass the contents of the service account credentials file.
     * See https://github.com/GoogleCloudPlatform/gcloud-java#authentication
     *
     * Note that most of the time it's enough to just open a file via
     * Files.newInputStream(Paths.get(URI.create( path ))).
     **/
    public static FileSystem getAuthenticatedGcs(String projectId, String bucket, byte[] credentials) throws IOException {
        StorageOptions.Builder builder = StorageOptions.newBuilder()
                .setProjectId(projectId);
        if (null != credentials) {
            builder = builder.setCredentials(GoogleCredentials.fromStream(new ByteArrayInputStream(credentials)));
        }
        // generous timeouts, to avoid tests failing when not warranted.
        StorageOptions storageOptions = setGenerousTimeouts(builder).build();

        // 2. Create GCS filesystem object with those credentials
        return CloudStorageFileSystem.forBucket(bucket, CloudStorageConfiguration.DEFAULT, storageOptions);
    }
}
