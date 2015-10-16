package org.broadinstitute.hellbender.utils.gcs;

import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.options.DataflowPipelineOptions;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.utils.GenomicsFactory;

import java.io.*;
import java.security.GeneralSecurityException;

/**
 * We need authentication options from GCSOptions, and Dataflow options from DataflowPipelineOptions.
 * Neither inherits from the other, so we have to put them together like this.
 *
 * This also includes code to save the OfflineAuth onto the pipelineoptions, so we can get to them later.
 * You see, GCSOptions will store the API Key and the path to client-secrets.json, but the latter
 * doesn't help us if we call createCredentials on the worker (because it won't be able to find
 * the client-secrets file). So here instead we create the OfflineAuth first (which has the
 * contents of the file, if one was specified) and stash that in the options.
 * This guarantees that anyone can call GATKGCSOptions.Methods.getOfflineAuth and
 * succeed.
 *
 * NOTE: We're deprecating those, use AuthHolder instead. It's serializable which makes it easier
 * to move around, and it can also create a storage client (plus also a Genomics client!)
 */
public interface GATKGCSOptions extends GCSOptions, DataflowPipelineOptions {

    public static class Methods {
        public static void setOfflineAuth(GATKGCSOptions opts, GenomicsFactory.OfflineAuth auth) throws IOException {
            ByteArrayOutputStream os = new ByteArrayOutputStream();
            try (ObjectOutputStream oos = new ObjectOutputStream(os)) {
                oos.writeObject(auth);
                oos.flush();
            }
            opts.setSerializedOfflineAuth(os.toByteArray());
        }
        public static GenomicsFactory.OfflineAuth getOfflineAuth(GATKGCSOptions opts) throws IOException, ClassNotFoundException, GeneralSecurityException {
            byte[] serialized = opts.getSerializedOfflineAuth();
            if (null==serialized && opts.getApiKey()!=null) {
                // fall back to using the API key only (even if a secrets file was also specified).
                GenomicsFactory.Builder builder =
                        GenomicsFactory.builder(opts.getAppName()).setNumberOfRetries(opts.getNumberOfRetries());
                return builder.build().getOfflineAuthFromApiKey(opts.getApiKey());
            }
            try (ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(serialized))) {
                return (GenomicsFactory.OfflineAuth)(is.readObject());
            }
        }

        public static Storage.Objects createStorageClient(GATKGCSOptions opts) throws GeneralSecurityException, IOException, ClassNotFoundException {
            GenomicsFactory.OfflineAuth auth = getOfflineAuth(opts);
            return GCSOptions.Methods.createStorageClient(opts.as(GCSOptions.class), auth);
        }

    }

    // you don't need to call those directly, use the helper methods in Methods instead.
    void setSerializedOfflineAuth(byte[] auth);
    byte[] getSerializedOfflineAuth();
}