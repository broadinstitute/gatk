package org.broadinstitute.hellbender.utils.gcs;

import com.google.cloud.storage.StorageOptions;
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration;
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration.Builder;
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystemProvider;
import com.google.cloud.storage.contrib.nio.SeekableByteChannelPrefetcher;
import com.google.common.base.Strings;
import com.google.common.io.ByteStreams;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import shaded.cloud_nio.com.google.api.gax.retrying.RetrySettings;
import shaded.cloud_nio.com.google.auth.oauth2.GoogleCredentials;
import shaded.cloud_nio.com.google.cloud.http.HttpTransportOptions;
import shaded.cloud_nio.org.threeten.bp.Duration;

import java.io.*;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.UUID;
import java.util.function.Function;

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

    public static boolean isCloudStorageUrl(final java.nio.file.Path path) {
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
    //NOTE: this should no longer be needed once everything uses GATKPathSpecifier, since the URI underlying
    // a GATKPathSpecifier always has a protocol; for a local file the protocol is "file", and the file protocol
    // requires absolute paths
    public static String makeFilePathAbsolute(String path){
        if (isCloudStorageUrl(path) || isHadoopUrl(path) || isFileUrl(path)){
            return path;
        } else {
            return new File(path).getAbsolutePath();
        }
    }

    /**
     * Open a file for reading regardless of whether it's on GCS, HDFS or local disk.
     *
     * If the file ends with .gz will attempt to wrap it in an appropriate unzipping stream
     *
     * @param path the GCS, HDFS or local path to read from. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     * @return an InputStream that reads from the specified file.
     */
    public static InputStream openFile(String path) {
        try {
            Utils.nonNull(path);
            InputStream inputStream;
            if (isCloudStorageUrl(path)) {
                java.nio.file.Path p = getPathOnGcs(path);
                inputStream = Files.newInputStream(p);
            } else if (isHadoopUrl(path)) {
                Path file = new org.apache.hadoop.fs.Path(path);
                FileSystem fs = file.getFileSystem(new Configuration());
                inputStream = fs.open(file);
            } else {
                //TODO: move this to a higher level (or make GATKPathSpecifier.inputStream do the block
                // compressed wrapping, or else make this
                 inputStream = new GATKPathSpecifier(path).getInputStream();
            }

            if(IOUtil.hasBlockCompressedExtension(path)){
                return IOUtils.makeZippedInputStream(new BufferedInputStream(inputStream));
            } else {
                return inputStream;
            }
        } catch (IOException x) {
            throw new UserException.CouldNotReadInputFile(path, x);
        }
    }

    /**
     * Open a binary file for writing regardless of whether it's on GCS, HDFS or local disk.
     * For writing to GCS it'll use the application/octet-stream MIME type.
     *
     * @param path the GCS or local path to write to. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     * @return an OutputStream that writes to the specified file.
     */
    public static OutputStream createFile(String path) {
        Utils.nonNull(path);
        try {
            if (isCloudStorageUrl(path)) {
                java.nio.file.Path p = getPathOnGcs(path);
                return Files.newOutputStream(p);
            } else if (isHadoopUrl(path)) {
                Path file = new Path(path);
                FileSystem fs = file.getFileSystem(new Configuration());
                return fs.create(file);
            } else {
                return new FileOutputStream(path);
            }
        } catch (IOException x) {
            throw new UserException.CouldNotCreateOutputFile("Could not create file at path:" + path + " due to " + x.getMessage(), x);
        }
    }

    /**
     * Copies a file. Can be used to copy e.g. from GCS to local.
     *
     * @param sourcePath the path to read from. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     * @param destPath the path to copy to. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     * @throws IOException
     */
    public static void copyFile(String sourcePath, String destPath) throws IOException {
        try (
            InputStream in = openFile(sourcePath);
            OutputStream fout = createFile(destPath)) {
            ByteStreams.copy(in, fout);
        }
    }

    /**
     * Deletes a file: local, GCS or HDFS.
     *  @param pathToDelete the path to delete. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     *
     */
    public static void deleteFile(String pathToDelete) throws IOException {
        if (isCloudStorageUrl(pathToDelete)) {
            java.nio.file.Path p = getPathOnGcs(pathToDelete);
            Files.delete(p);
        } else if (isHadoopUrl(pathToDelete)) {
            Path file = new Path(pathToDelete);
            FileSystem fs = file.getFileSystem(new Configuration());
            fs.delete(file, false);
        } else {
            boolean ok = new File(pathToDelete).delete();
            if (!ok) throw new IOException("Unable to delete '"+pathToDelete+"'");
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
            IOUtils.deleteOnExit(IOUtils.getPath(path + FileExtensions.TRIBBLE_INDEX));
            IOUtils.deleteOnExit(IOUtils.getPath(path + FileExtensions.TABIX_INDEX));
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
            return new Path(stagingLocation, prefix + UUID.randomUUID().toString() + suffix).toString();
        } else {
            throw new IllegalArgumentException("Staging location is not remote: " + stagingLocation);
        }
    }

    /**
     * Returns true if we can read the first byte of the file.
     *  @param path The folder where you want the file to be (local, GCS or HDFS).
     *
     */
    public static boolean fileExists(String path) {
        final boolean MAYBE = false;
        try (InputStream inputStream = openFile(path)) {
            int ignored = inputStream.read();
        } catch (UserException.CouldNotReadInputFile notthere) {
            // file isn't there
            return false;
        } catch (FileNotFoundException x) {
            // file isn't there
            return false;
        } catch (IOException x) {
            // unexpected problem while reading the file. The file may exist, but it's not accessible.
            return MAYBE;
        }
        return true;
    }

    /**
     * Returns the file size of a file pointed to by a GCS/HDFS/local path
     *
     * @param path The URL to the file whose size to return
     * @return the file size in bytes
     * @throws IOException
     */
    public static long fileSize(String path) throws IOException {
        if (isCloudStorageUrl(path)) {
            java.nio.file.Path p = getPathOnGcs(path);
            return Files.size(p);
        } else if (isHadoopUrl(path)) {
            Path hadoopPath = new Path(path);
            FileSystem fs = hadoopPath.getFileSystem(new Configuration());
            return fs.getFileStatus(hadoopPath).getLen();
        } else {
            return new File(path).length();
        }
    }

    /**
     * Returns the total file size of all files in a directory, or the file size if the path specifies a file.
     * Note that sub-directories are ignored - they are not recursed into.
     * Only supports HDFS and local paths.
     *
     * @param path The URL to the file or directory whose size to return
     * @return the total size of all files in bytes
     */
    public static long dirSize(String path) {
        try {
            // GCS case (would work with local too)
            if (isCloudStorageUrl(path)) {
                java.nio.file.Path p = getPathOnGcs(path);
                if (Files.isRegularFile(p)) {
                    return Files.size(p);
                }
                return Files.list(p).mapToLong(
                    q -> {
                        try {
                            return (Files.isRegularFile(q) ? Files.size(q) : 0);
                        } catch (IOException e) {
                            throw new RuntimeIOException(e);
                        }
                    }
                ).sum();
            }
            // local file or HDFS case
            Path hadoopPath = new Path(path);
            FileSystem fs = new Path(path).getFileSystem(new Configuration());
            FileStatus status = fs.getFileStatus(hadoopPath);
            if (status == null) {
                throw new UserException.CouldNotReadInputFile(path, "File not found.");
            }
            long size = 0;
            if (status.isDirectory()) {
                for (FileStatus st : fs.listStatus(status.getPath())) {
                    if (st.isFile()) {
                        size += st.getLen();
                    }
                }
            } else {
                size += status.getLen();
            }
            return size;
        } catch (RuntimeIOException | IOException e) {
            throw new UserException("Failed to determine total input size of " + path + "\n Caused by:" + e.getMessage(), e);
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
    public static java.nio.file.Path getPathOnGcs(String gcsUrl) {
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
    public static java.nio.file.FileSystem getAuthenticatedGcs(String projectId, String bucket, byte[] credentials) throws IOException {
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

    /**
     * Wrap a SeekableByteChannel with a prefetcher.
     * @param bufferSizeMB buffer size in mb which the prefetcher should fetch ahead.
     * @param channel a channel that needs prefetching
     */
    public static SeekableByteChannel addPrefetcher(final int bufferSizeMB, final SeekableByteChannel channel) {
        try {
            return SeekableByteChannelPrefetcher.addPrefetcher(bufferSizeMB, channel);
        } catch (final IOException ex) {
            throw new GATKException("Unable to initialize the prefetcher: " + ex);
        }
    }

    /**
     * Creates a wrapping function which adds a prefetcher if the buffer size is > 0 if it's <= 0 then this wrapper returns the
     * original channel.
     * @param cloudPrefetchBuffer the prefetcher buffer size in MB
     */
    public static Function<SeekableByteChannel, SeekableByteChannel> getPrefetchingWrapper(final int cloudPrefetchBuffer) {
        return cloudPrefetchBuffer > 0 ? rawChannel -> addPrefetcher(cloudPrefetchBuffer, rawChannel) : Utils.identityFunction();
    }
}
