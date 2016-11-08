package org.broadinstitute.hellbender.utils.gcs;

import com.google.api.services.storage.Storage;
import com.google.cloud.AuthCredentials;
import com.google.cloud.RetryParams;
import com.google.cloud.dataflow.sdk.options.GcsOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.util.GcsUtil;
import com.google.cloud.dataflow.sdk.util.Transport;
import com.google.cloud.dataflow.sdk.util.gcsfs.GcsPath;
import com.google.cloud.storage.StorageOptions;
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration;
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import com.google.common.io.ByteStreams;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.util.TabixUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.*;
import java.nio.channels.Channels;
import java.nio.file.NoSuchFileException;
import java.util.Arrays;
import java.util.List;
import java.util.UUID;

/**
 * Utilities for dealing with google buckets
 */
public final class BucketUtils {
    public static final String GCS_PREFIX = "gs://";
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
     * Open a file for reading regardless of whether it's on GCS, HDFS or local disk.
     *
     * If the file ends with .gz will attempt to wrap it in an appropriate unzipping stream
     *
     * @param path the GCS, HDFS or local path to read from. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     * @param popts the pipeline's options, with authentication information.
     * @return an InputStream that reads from the specified file.
     */
    public static InputStream openFile(String path, PipelineOptions popts) {
        try {
            Utils.nonNull(path);
            InputStream inputStream;
            if (BucketUtils.isCloudStorageUrl(path)) {
                Utils.nonNull(popts, "Cannot load from a GCS path without authentication.");
                inputStream = Channels.newInputStream(new GcsUtil.GcsUtilFactory().create(popts).open(GcsPath.fromUri(path)));
            } else if (isHadoopUrl(path)) {
                Path file = new org.apache.hadoop.fs.Path(path);
                FileSystem fs = file.getFileSystem(new Configuration());
                inputStream = fs.open(file);
            } else {
                inputStream = new FileInputStream(path);
            }

            if(AbstractFeatureReader.hasBlockCompressedExtension(path)){
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
     * @param popts the pipeline's options, with authentication information.
     * @return an OutputStream that writes to the specified file.
     */
    public static OutputStream createFile(String path, PipelineOptions popts) {
        Utils.nonNull(path);
        try {
            if (isCloudStorageUrl(path)) {
                return Channels.newOutputStream(new GcsUtil.GcsUtilFactory().create(popts).create(GcsPath.fromUri(path), "application/octet-stream"));
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
     * Open a binary file for writing regardless of whether it's on GCS, HDFS or local disk.
     * For writing to GCS it'll use the application/octet-stream MIME type.
     *
     * @param path the GCS , HDFS, or local path to write to. If HDFS, it must start with "hdfs://".
     *             If GCS, it must start with "gs://" and you must be using API Key authentication.
     * @param auth authentication information.
     * @return an OutputStream that writes to the specified file.
     */
    public static OutputStream createFile(String path, AuthHolder auth) {
        PipelineOptions popts = auth.asPipelineOptionsDeprecated();
        return createFile(path, popts);
    }

    /**
     * Open a binary file for writing regardless of whether it's on HDFS or local disk.
     *
     * @param path the local path to write to.
     * @return an OutputStream that writes to the specified file.
     */
    public static OutputStream createNonGCSFile(String path) {
        return createFile(path, (PipelineOptions)null);
    }

    /**
     * Copies a file. Can be used to copy e.g. from GCS to local.
     *
     * @param sourcePath the path to read from. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     * @param popts the pipeline's options, with authentication information.
     * @param destPath the path to copy to. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     * @throws IOException
     */
    public static void copyFile(String sourcePath, PipelineOptions popts, String destPath) throws IOException {
        try (
            InputStream in = openFile(sourcePath, popts);
            OutputStream fout = createFile(destPath, popts)) {
            ByteStreams.copy(in, fout);
        }
    }

    /**
     * Deletes a file: local, GCS or HDFS.
     *
     * @param pathToDelete the path to delete. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     * @param popts the pipeline's options, with authentication information.
     */
    public static void deleteFile(String pathToDelete, PipelineOptions popts) throws IOException {
        if (BucketUtils.isCloudStorageUrl(pathToDelete)) {
            GcsPath path = GcsPath.fromUri(pathToDelete);
            GcsOptions gcsOptions = popts.as(GcsOptions.class);
            Storage storage = Transport.newStorageClient(gcsOptions).build();
            storage.objects().delete(path.getBucket(), path.getObject()).execute();
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
     * This file (and possible indexes associated with it will be scheduled for deleton on shutdown
     *
     * @param prefix a prefix for the file name
     *               for remote paths this should be a valid URI to root the temporary file in (ie. gcs://hellbender/staging/)
     *               there is no guarantee that this will be used as the root of the tmp file name, a local prefix may be placed in the tmp folder for exapmle
     * @param extension and extension for the temporary file path, the resulting path will end in this
     * @param authHolder authentication for remote file paths, may be null
     * @return a path to use as a temporary file, on remote file systems which don't support an atomic tmp file reservation a path is chosen with a long randomized name
     *
     */
    public static String getTempFilePath(String prefix, String extension, AuthHolder authHolder){
        if (BucketUtils.isCloudStorageUrl(prefix) || (BucketUtils.isHadoopUrl(prefix))){
            final String path = randomRemotePath(prefix, "", extension);
            deleteOnExit(path, authHolder);
            deleteOnExit(path + Tribble.STANDARD_INDEX_EXTENSION, authHolder);
            deleteOnExit(path + TabixUtils.STANDARD_INDEX_EXTENSION, authHolder);
            deleteOnExit(path + ".bai", authHolder);
            deleteOnExit(path.replaceAll(extension + "$", ".bai"), authHolder); //if path ends with extension, replace it with .bai
            return path;
        } else {
            return IOUtils.createTempFile(prefix, extension).getAbsolutePath();
        }
    }

    /**
     * Schedule a file to be deleted on JVM shutdown.
     * @param fileToDelete the path to the file to be deleted
     * @param authHolder authentication for remote files
     */
    public static void deleteOnExit(String fileToDelete, AuthHolder authHolder){
        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                try {
                    deleteFile(fileToDelete, authHolder.asPipelineOptionsDeprecated());
                } catch (IOException e) {
                    logger.warn("Failed to delete file: " + fileToDelete+ ".", e);
                }
            }
        });
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
            return GcsPath.fromUri(stagingLocation).resolve(prefix + UUID.randomUUID().toString() + suffix).toString();
        } else if (isHadoopUrl(stagingLocation)) {
            return new Path(stagingLocation, prefix + UUID.randomUUID().toString() + suffix).toString();
        } else {
            throw new IllegalArgumentException("Staging location is not remote: " + stagingLocation);
        }
    }

    /**
     * Returns true if we can read the first byte of the file.
     *
     * @param path The folder where you want the file to be (local, GCS or HDFS).
     * @param popts the pipeline's options, with authentication information.
     */
    public static boolean fileExists(String path, PipelineOptions popts) {
        final boolean MAYBE = false;
        try {
            InputStream inputStream = BucketUtils.openFile(path, popts);
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
     * @param popts PipelineOptions for GCS (if relevant; otherwise pass null)
     * @return the file size in bytes
     * @throws IOException
     */
    public static long fileSize(String path, PipelineOptions popts) throws IOException {
        if (isCloudStorageUrl(path)) {
            return new GcsUtil.GcsUtilFactory().create(popts).fileSize(GcsPath.fromUri(path));
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
     * @param popts PipelineOptions for GCS (if relevant; otherwise pass null)
     * @return the total size of all files in bytes
     */
    public static long dirSize(String path, PipelineOptions popts) {
        try {
            // GCS case
            if (isCloudStorageUrl(path)) {
                final GcsUtil gcsUtil = new GcsUtil.GcsUtilFactory().create(popts);
                final List<GcsPath> pathsInDir = gcsUtil.expand(GcsPath.fromUri(path));
                return pathsInDir.stream().mapToLong((path1) -> {
                    try {
                        return gcsUtil.fileSize(path1);
                    } catch( final NoSuchFileException e) {
                        return 0;
                    } catch( final IOException e) {
                        throw new RuntimeIOException(e);
                    }
                }).sum();

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
        final String BUCKET = split[2];
        return String.join("/", Arrays.copyOfRange(split, 3, split.length));

    }

    /**
     * String -> Path. This *should* not be necessary (use Paths.get(URI.create(...)) instead) , but it currently is
     * on Spark because using the fat, shaded jar breaks the registration of the GCS FilesystemProvider.
     */
    public static java.nio.file.Path getPathOnGcs(String gcsUrl) {
        final String[] split = gcsUrl.split("/");
        final String BUCKET = split[2];
        final String pathWithoutBucket = String.join("/", Arrays.copyOfRange(split, 3, split.length));
        return CloudStorageFileSystem.forBucket(BUCKET).getPath(pathWithoutBucket);
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
     */
    public static java.nio.file.FileSystem getAuthenticatedGcs(String projectId, String bucket, byte[] credentials) throws IOException {
        StorageOptions.Builder builder = StorageOptions.newBuilder()
                .setProjectId(projectId);
        if (null != credentials) {
            builder = builder.setAuthCredentials((AuthCredentials.createForJson(new ByteArrayInputStream(credentials))));
        }
        // generous timeouts, to avoid tests failing when not warranted.
        StorageOptions storageOptions = builder.setConnectTimeout(60000)
            .setReadTimeout(60000)
            .setRetryParams(RetryParams.newBuilder()
                    .setRetryMaxAttempts(10)
                    .setRetryMinAttempts(6)
                    .setMaxRetryDelayMillis(30000)
                    .setTotalRetryPeriodMillis(120000)
                    .setInitialRetryDelayMillis(250)
                    .build())
            .build();

        // 2. Create GCS filesystem object with those credentials
        return CloudStorageFileSystem.forBucket(bucket, CloudStorageConfiguration.DEFAULT, storageOptions);
    }
}
