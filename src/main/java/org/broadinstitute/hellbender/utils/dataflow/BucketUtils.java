package org.broadinstitute.hellbender.utils.dataflow;

import com.google.appengine.tools.cloudstorage.GcsFilename;
import com.google.appengine.tools.cloudstorage.GcsInputChannel;
import com.google.appengine.tools.cloudstorage.GcsOutputChannel;
import com.google.appengine.tools.cloudstorage.GcsService;
import com.google.appengine.tools.cloudstorage.GcsServiceFactory;
import com.google.appengine.tools.cloudstorage.RetryParams;

import java.io.*;
import java.nio.channels.Channels;

/**
 * Utilities for dealing with google buckets
 */
public final class BucketUtils {
    public static final String GCS_PREFIX = "gs://";

    private BucketUtils(){} //private so that no one will instantiate this class

    public static boolean isCloudStorageUrl(String path) {
        return path.startsWith(GCS_PREFIX);
    }
}
