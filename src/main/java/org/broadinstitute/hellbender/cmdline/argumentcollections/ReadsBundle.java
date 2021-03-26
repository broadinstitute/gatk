package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.io.IOPath;
import htsjdk.plugin.bundle.BundleResourceType;
import htsjdk.plugin.bundle.InputBundle;
import htsjdk.plugin.bundle.InputIOPathResource;
import htsjdk.plugin.bundle.InputResource;
import htsjdk.samtools.SamFiles;
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
import java.util.List;
import java.util.Optional;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

//TODO:
// reference input triplet
// propagate GATKPath tag and attributes
// Note: This class only handles *input* pairs. Output pairs would require a separate class since the
// interface types for output bundle resources differ from input resources.

/**
 * A reads bundle may optionally have an index, but its not required.
 */
public class ReadsBundle extends InputBundle implements Serializable {
    private static final long serialVersionUID = 1L;
    public static final String BUNDLE_EXTENSION = ".json";

    final private transient GATKPath cachedReadsPath;

    //TODO: bridge/propagate tag/attributes...
    public ReadsBundle(final GATKPath reads) {
        super(Utils.nonNull(reads), BundleResourceType.READS);
        cachedReadsPath = getReadsPathFromBundle();
    }

    public ReadsBundle(final GATKPath reads, final GATKPath index) {
        super(Arrays.asList(
                new InputIOPathResource(Utils.nonNull(reads), BundleResourceType.READS),
                new InputIOPathResource(Utils.nonNull(index), BundleResourceType.INDEX)));
        cachedReadsPath = getReadsPathFromBundle();
    }

    public ReadsBundle(final String jsonString){
        // GATKPath::new ensures that the IOPath resources in the resulting bundle
        // contain GATKPath objects
        super(jsonString, GATKPath::new);
        cachedReadsPath = getReadsPathFromBundle();
    }

    public static ReadsBundle getReadsBundleFromJSON(final GATKPath jsonPath) {
        return new ReadsBundle(getJSONStringFromPath(jsonPath));
    }

    public GATKPath getReads() { return cachedReadsPath; }

    public Optional<GATKPath> getIndex() {
        final Supplier<GATKException> ex = () -> new GATKException("Index resource is present with a null path");
        // its OK for there to be no index resource...
        final Optional<InputResource> inputResource = get(BundleResourceType.INDEX);
        // ...but if there *is* an index resources, it must contain a non-null path...
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
            org.apache.commons.io.IOUtils.copy(is, jsonStringWriter, "UTF-8");
            return jsonStringWriter.toString();
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile("Failed to load reads bundle json from: " +path.getRawInputString(), e);
        }
    }

//    public GATKPath getReads() {
//        return getPathSpecifierWithTypeTag(reads);
//    }
//
//    public GATKPath getIndex() {
//        return getPathSpecifierWithTypeTag(index);
//    }

//    private static GATKPath getPathSpecifierWithTypeTag(final BundleFile index) {
//        final GATKPath path = new GATKPath(index.getPath());
//        path.setTagAttributes(Collections.singletonMap(FILE_TYPE_KEY, index.getFileType()));
//        return path;
//    }

}
