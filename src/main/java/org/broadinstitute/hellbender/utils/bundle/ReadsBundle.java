package org.broadinstitute.hellbender.utils.bundle;

import com.google.common.base.Charsets;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collections;

public class ReadsBundle {
    public static final String BUNDLE_EXTENSION = ".json";
    private static final Logger LOG = LogManager.getLogger(ReadsBundle.class);
    private static final String SCHEMA_VERSION = "0.1.0";
    private static final Gson GSON = new GsonBuilder().setPrettyPrinting().create();

    public static final String FILE_TYPE_KEY = "FILE_TYPE";
    public static final String READS_BUNDLE_SCHEMA_TYPE = "ReadsBundle";
    private final BundleFile reads;
    private final BundleFile index;
    private final String schemaVersion;
    //TODO make this part of the json without being a field
    private final String schemaType;

    public ReadsBundle(final GATKPath reads, ReadsType readsInputType, final GATKPath index, ReadsIndexType indexType){
        this(new BundleFile(reads.getRawInputString(), readsInputType.toString()), new BundleFile(index.getRawInputString(), indexType.toString()));
    }

    public ReadsBundle(final BundleFile reads, final BundleFile index) {
        this.schemaVersion = SCHEMA_VERSION;
        this.schemaType = READS_BUNDLE_SCHEMA_TYPE;
        this.reads = reads;
        this.index = index;
    }

    public static ReadsBundle fromPath(GATKPath path){
        try(final InputStreamReader inputStreamReader = new InputStreamReader(new BufferedInputStream(path.getInputStream()), Charsets.UTF_8)) {
            final ReadsBundle readsBundle = GSON.fromJson(inputStreamReader, ReadsBundle.class);
            checkSchema(readsBundle, path.toString());
            return readsBundle;
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile("Failed to load reads bundle json from: " +path.getRawInputString(), e);
        }
    }



    private static void checkSchema(final ReadsBundle readsBundle, String source) {
        if(!READS_BUNDLE_SCHEMA_TYPE.equals(readsBundle.schemaType)){
            throw new UserException("Failed while loading ReadsBundle from " + source +
                    ".\n Json does not include a correct schemaType key.  Expected " + READS_BUNDLE_SCHEMA_TYPE + " but found " +readsBundle.schemaType + " ." );
        }
        if (!SCHEMA_VERSION.equals(readsBundle.schemaVersion)) {
            LOG.warn("Loaded a reads bundle with an unexpected schema version.  Expected:"+ SCHEMA_VERSION +" but found: " + readsBundle.schemaVersion);
        }
    }

    public static ReadsBundle fromJson(String json){
        final ReadsBundle readsBundle = GSON.fromJson(json, ReadsBundle.class);
        checkSchema(readsBundle, "a string");
        return readsBundle;
    }

    public static boolean looksLikeAReadsBundle(final GATKPath rawReadPath) {
        return rawReadPath.getURI().getPath().endsWith(BUNDLE_EXTENSION);
    }

    public String toJson(){
        return GSON.toJson(this);
    }

    public GATKPath getReads() {
        return getPathSpecifierWithTypeTag(reads);
    }

    public GATKPath getIndex() {
        return getPathSpecifierWithTypeTag(index);
    }

    private static GATKPath getPathSpecifierWithTypeTag(final BundleFile index) {
        final GATKPath path = new GATKPath(index.getPath());
        path.setTagAttributes(Collections.singletonMap(FILE_TYPE_KEY, index.getFileType()));
        return path;
    }
}

