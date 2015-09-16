package org.broadinstitute.hellbender.utils.dataflow;

import com.google.api.client.util.ByteStreams;
import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.options.GcsOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.util.GcsUtil;
import com.google.cloud.dataflow.sdk.util.Transport;
import com.google.cloud.dataflow.sdk.util.gcsfs.GcsPath;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.channels.Channels;
import java.util.UUID;

/**
 * Utilities for dealing with google buckets
 */
public final class BucketUtils {
    public static final String GCS_PREFIX = "gs://";
    public static final String HDFS_PREFIX = "hdfs://";

    private BucketUtils(){} //private so that no one will instantiate this class

    public static boolean isCloudStorageUrl(String path) {
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

            if(path.endsWith(".gz")){
                return IOUtils.makeZippedInputStream(new BufferedInputStream(inputStream));
            } else {
                return inputStream;
            }
        } catch (IOException x) {
            throw new UserException.CouldNotReadInputFile(path, x);
        }
    }

    /**
     * Open a file for reading regardless of whether it's on GCS, HDFS or local disk.
     *
     * @param path the GCS, HDFS or local path to read from. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     * @param apiKey a GCS api key for accessing gcs paths, may be null if you will never open a gcs file
     * @return an InputStream that reads from the specified file.
     */
    public static InputStream openFile(String path, String apiKey) {
        return openFile(path, apiKey != null ? getAuthenticatedGCSOptions(apiKey) : null);
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
            throw new UserException.CouldNotCreateOutputFile(path, x);
        }
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
     * Create a new GCS options with a given api key.
     * @param apiKey a valid GCS api key
     */
    public static GCSOptions getAuthenticatedGCSOptions(String apiKey) {
        Utils.nonNull(apiKey);
        GCSOptions options = PipelineOptionsFactory.as(GCSOptions.class);
        options.setApiKey(apiKey);
        return options;
    }
}
