package org.broadinstitute.hellbender.utils.dataflow;

import com.google.appengine.tools.cloudstorage.GcsFilename;
import com.google.appengine.tools.cloudstorage.GcsInputChannel;
import com.google.appengine.tools.cloudstorage.GcsOutputChannel;
import com.google.appengine.tools.cloudstorage.GcsService;
import com.google.appengine.tools.cloudstorage.GcsServiceFactory;
import com.google.appengine.tools.cloudstorage.RetryParams;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.util.GcsUtil;
import com.google.cloud.dataflow.sdk.util.gcsfs.GcsPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.*;
import java.nio.channels.Channels;

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
            if (BucketUtils.isCloudStorageUrl(path)) {
                return Channels.newInputStream(new GcsUtil.GcsUtilFactory().create(popts).open(GcsPath.fromUri(path)));
            } else {
                return new FileInputStream(path);
            }
        } catch (Exception x) {
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
        } catch (Exception x) {
            throw new UserException.CouldNotCreateOutputFile(path, x);
        }
    }
}
