package org.broadinstitute.hellbender.utils;

import com.google.appengine.tools.cloudstorage.GcsFilename;
import com.google.appengine.tools.cloudstorage.GcsInputChannel;
import com.google.appengine.tools.cloudstorage.GcsOutputChannel;
import com.google.appengine.tools.cloudstorage.GcsService;
import com.google.appengine.tools.cloudstorage.GcsServiceFactory;
import com.google.appengine.tools.cloudstorage.RetryParams;

import java.io.*;
import java.nio.channels.Channels;

public class BucketUtils {
    //no instances of this class
    private BucketUtils(){};

    private static InputStream openStreamFromGcs(GcsFilename fileName)
            throws IOException, ClassNotFoundException {
        final GcsService gcsService = GcsServiceFactory.createGcsService(RetryParams.getDefaultInstance());
        GcsInputChannel readChannel = gcsService.openPrefetchingReadChannel(fileName, 0, 1024 * 1024);
        return Channels.newInputStream(readChannel);
    }
}
