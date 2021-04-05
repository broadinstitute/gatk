package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.io.IOPath;
import htsjdk.beta.plugin.bundle.BundleResourceType;
import htsjdk.beta.plugin.bundle.Bundle;
import htsjdk.beta.plugin.bundle.IOPathResource;
import htsjdk.beta.plugin.bundle.BundleResource;
import htsjdk.samtools.SamFiles;
import htsjdk.utils.ValidationUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;
import java.io.StringWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

//TODO:
// propagate GATKPath tag and attributes on JSON deserialization ?
// subContentType is always inferred, never explicitly provided...
// only IOPath resources can be serialized
// Barclay won’t be able to print out bundle contents on the command line...
// What if there is a clash between the attributes specified ON an bundle input, and those specified IN the bundle ?

/**
 * A reads bundle may optionally have an index, but its not required.
 */
public class ReadsBundle extends Bundle implements Serializable {
    private static final long serialVersionUID = 1L;
    private static final Logger logger = LogManager.getLogger(ReadsBundle.class);

    public static final String BUNDLE_EXTENSION = ".json";

    private final GATKPath cachedReadsPath;

    public ReadsBundle(final GATKPath reads) {
        super(Collections.singletonList(toInputResource(Utils.nonNull(reads), BundleResourceType.READS)));
        cachedReadsPath = getReadsPathFromBundle();
    }

    public ReadsBundle(final GATKPath reads, final GATKPath index) {
        super(Arrays.asList(
                toInputResource(Utils.nonNull(reads), BundleResourceType.READS),
                toInputResource(Utils.nonNull(index), BundleResourceType.INDEX)));
        cachedReadsPath = getReadsPathFromBundle();
    }

    public ReadsBundle(final String jsonString){
        super(jsonString, GATKPath::new);
        cachedReadsPath = getReadsPathFromBundle();
        //TODO: bridge/propagate the InputResource tag/attributes back to the GATKPath in this bundle...

    }

    public static ReadsBundle getReadsBundleFromJSON(final GATKPath jsonPath) {
        return new ReadsBundle(getJSONStringFromPath(jsonPath));
    }

    public GATKPath getReads() { return cachedReadsPath; }

    public Optional<GATKPath> getIndex() {
        final Supplier<GATKException> ex = () -> new GATKException("Index resource is present with a null path");
        final Optional<BundleResource> inputResource = get(BundleResourceType.INDEX);
        // its OK for there to be no index resource, but if there *is* an index resource, it must contain
        // a non-null path...
        return inputResource.isPresent() ?
                Optional.of((GATKPath) inputResource.get().getIOPath().orElseThrow(ex)) :
                Optional.empty();
    }

    public static List<ReadsBundle> fromPathLists(List<Path> reads, List<Path> indexes){
        Utils.nonNull(reads);
        Utils.nonNull(indexes);
        return fromLists(Utils.map(reads, IOUtils::toGATKPath), Utils.map(indexes, IOUtils::toGATKPath));
    }

    public static List<ReadsBundle> fromLists(List<GATKPath> reads, List<GATKPath> indexes){
        Utils.nonNull(reads);
        Utils.nonNull(indexes);
        Utils.validate(reads.size() == indexes.size(),
                String.format("reads (%d) and indexes (%d) must be the same length", reads.size(), indexes.size()));

        return IntStream.range(0, reads.size())
                .mapToObj(i -> new ReadsBundle(reads.get(i), indexes.get(i)))
                .collect(Collectors.toList());
    }

    public static ReadsBundle resolveIndex(GATKPath reads){
        final Path index = SamFiles.findIndex(reads.toPath());
        if (index == null) {
            return new ReadsBundle(reads);
        }
        return new ReadsBundle(reads, IOUtils.toGATKPath(index));
    }

    public static boolean looksLikeAReadsBundle(final GATKPath rawReadPath) {
        return rawReadPath.getURI().getPath().endsWith(BUNDLE_EXTENSION);
    }

    // try to get the reads for the side effect of validating that the bundle actually contains
    // a reads resource
    private GATKPath getReadsPathFromBundle() {
        final Supplier<GATKException> ex = () -> new GATKException("No reads in reads bundle");
        final Optional<IOPath> ioPath = get(BundleResourceType.READS).orElseThrow(ex).getIOPath();
        return (GATKPath) ioPath.orElseThrow(ex);
    }

    private static String getJSONStringFromPath(final GATKPath path) {
        try (final InputStream is = new BufferedInputStream(path.getInputStream())) {
            final StringWriter jsonStringWriter = new StringWriter();
            //TODO: the UTF-8 encoding ofcthese should be codified somewhere else...
            org.apache.commons.io.IOUtils.copy(is, jsonStringWriter, "UTF-8");
            return jsonStringWriter.toString();
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile("Failed to load reads bundle json from: " +path.getRawInputString(), e);
        }
    }

    //TODO: move to htsjdk?
    private static IOPathResource toInputResource(final GATKPath gatkPath, final String providedContentType) {
        ValidationUtils.nonNull(gatkPath, "gatkPath");
        final Optional<Pair<String, String>> typePair = getContentTypes(gatkPath);
        if (typePair.isPresent()) {
            if (providedContentType != null && !typePair.get().getLeft().equals(providedContentType)) {
                logger.warn(String.format(
                        "Provided content type \"%s\" for \"%s\" doesn't match derived content type \"%s\"",
                        providedContentType,
                        gatkPath.getRawInputString(),
                        typePair.get().getLeft()));
            }
            return new IOPathResource(
                    gatkPath,
                    providedContentType,  // prefer the provided content type
                    typePair.get().getRight());
        } else {
            return new IOPathResource(
                    gatkPath,
                    providedContentType);
        }
    }

    private static Optional<Pair<String, String>> getContentTypes(final GATKPath gatkPath) {
        ValidationUtils.nonNull(gatkPath, "gatkPath");
        final Optional<String> extension = gatkPath.getExtension();
        if (extension.isPresent()) {
        }
        return Optional.empty();//TODO: finish this...
    }

}
