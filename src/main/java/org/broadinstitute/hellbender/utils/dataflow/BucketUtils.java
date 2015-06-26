package org.broadinstitute.hellbender.utils.dataflow;

import com.google.api.client.util.ByteStreams;
import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.options.GcsOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.util.GcsUtil;
import com.google.cloud.dataflow.sdk.util.Transport;
import com.google.cloud.dataflow.sdk.util.gcsfs.GcsPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;

import java.io.*;
import java.nio.channels.Channels;
import java.nio.channels.Pipe;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.GeneralSecurityException;
import java.util.Random;
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
     * Open a file for reading regardless of whether it's on GCS or local disk.
     *
     * @param path the GCS or local path to read from. If GCS, it must start with "gs://".
     * @param popts the pipeline's options, with authentication information.
     * @return an InputStream that reads from the specified file.

     */
    public static InputStream openFile(String path, PipelineOptions popts) {
        try {
            Utils.nonNull(path);
            if (BucketUtils.isCloudStorageUrl(path)) {
                Utils.nonNull(popts);
                return Channels.newInputStream(new GcsUtil.GcsUtilFactory().create(popts).open(GcsPath.fromUri(path)));
            } else {
                return new FileInputStream(path);
            }
        } catch (IOException x) {
            throw new UserException.CouldNotReadInputFile(path, x);
        }
    }

    /**
     * Open a binary file for writing regardless of whether it's on GCS or local disk.
     * For writing to GCS it'll use the application/octet-stream MIME type.
     *
     * @param path the GCS or local path to write to. If GCS, it must start with "gs://".
     * @param popts the pipeline's options, with authentication information.
     * @return an OutputStream that writes to the specified file.
     */
    public static OutputStream createFile(String path, PipelineOptions popts) {
        try {
            if (isCloudStorageUrl(path)) {
                return Channels.newOutputStream(new GcsUtil.GcsUtilFactory().create(popts).create(GcsPath.fromUri(path), "application/octet-stream"));
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
     * @param sourcePath the path to read from. If GCS, it must start with "gs://".
     * @param popts the pipeline's options, with authentication information.
     * @param destPath the path to copy to. If GCS, it must start with "gs://".
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
     * Deletes a file, local or on GCS.
     *
     * @param pathToDelete the path to delete. If GCS, it must start with "gs://".
     * @param popts the pipeline's options, with authentication information.
     */
    public static void deleteFile(String pathToDelete, PipelineOptions popts) throws IOException, GeneralSecurityException {
        if (BucketUtils.isCloudStorageUrl(pathToDelete)) {
            GcsPath path = GcsPath.fromUri(pathToDelete);
            GcsOptions gcsOptions = popts.as(GcsOptions.class);
            Storage storage = Transport.newStorageClient(gcsOptions).build();
            storage.objects().delete(path.getBucket(), path.getObject()).execute();
        } else {
            boolean ok = new File(pathToDelete).delete();
            if (!ok) throw new IOException("Unable to delete '"+pathToDelete+"'");
        }
    }

    /**
     * Picks a random name, by putting some random letters between "prefix" and "suffix".
     *
     * @param stagingLocation The folder where you want the file to be. Must start with "gs://"
     * @param prefix The beginning of the file name
     * @param suffix The end of the file name, e.g. ".tmp"
     */
    public static String randomGcsPath(String stagingLocation, String prefix, String suffix) {
        return GcsPath.fromUri(stagingLocation).resolve(prefix + UUID.randomUUID().toString() + suffix).toString();
    }

    /**
     * Returns true if we can read the first byte of the file.
     *
     * @param path The folder where you want the file to be (local or GCS).
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
}
