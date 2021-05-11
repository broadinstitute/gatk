package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.beta.plugin.bundle.BundleJSON;
import htsjdk.beta.plugin.bundle.BundleResource;
import htsjdk.beta.plugin.reads.ReadsBundle;
import htsjdk.samtools.SamFiles;
import htsjdk.utils.ValidationUtils;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.Serializable;
import java.nio.file.Path;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

//TODO:
// propagate GATKPath tag and attributes on JSON deserialization to the GATKPath in the bundle ? or leave it on
// the tool arg ?
// subContentType is always inferred, never explicitly provided...
// Barclay won’t be able to print out bundle contents on the command line...

/**
 * A reads bundle may optionally have an index, but its not required.
 */
public class GATKReadsBundle extends ReadsBundle<GATKPath> implements Serializable {
    private static final long serialVersionUID = 1L;

    public GATKReadsBundle(final GATKPath reads) { super(reads); }

    public GATKReadsBundle(final GATKPath reads, final GATKPath index) {
        super(reads, index);
    }

    /**
     * @param resources collection of resources. must have
     */
    public GATKReadsBundle(final Collection<BundleResource> resources) {
        super(resources);
    }

    public static GATKReadsBundle resolveIndex(GATKPath reads){
        final Path index = SamFiles.findIndex(reads.toPath());
        if (index == null) {
            return new GATKReadsBundle(reads);
        }
        return new GATKReadsBundle(reads, IOUtils.toGATKPath(index));
    }

    public static GATKReadsBundle getGATKReadsBundleFromString(final String jsonString) {
        return new GATKReadsBundle(BundleJSON.toBundle(jsonString, GATKPath::new).getResources());
    }

    public static GATKReadsBundle getGATKReadsBundleFromPath(final GATKPath jsonPath) {
        return getGATKReadsBundleFromString(htsjdk.beta.plugin.IOUtils.getStringFromPath(jsonPath));
    }

    public static List<GATKReadsBundle> fromPathLists(List<Path> reads, List<Path> indexes) {
        ValidationUtils.nonNull(reads, "reads");
        ValidationUtils.nonNull(indexes, "indexes");
        return fromLists(
                reads.stream().map(IOUtils::toGATKPath).collect(Collectors.toList()),
                indexes.stream().map(IOUtils::toGATKPath).collect(Collectors.toList()));
    }

    public static List<GATKReadsBundle> fromLists(List<GATKPath> reads, List<GATKPath> indexes) {
        ValidationUtils.nonNull(reads, "reads");
        ValidationUtils.nonNull(indexes, "indexes");
        ValidationUtils.validateArg(reads.size() == indexes.size(),
                String.format("reads (%d) and indexes (%d) must be the same length", reads.size(), indexes.size()));

        return IntStream.range(0, reads.size())
                .mapToObj(i -> new GATKReadsBundle(reads.get(i), indexes.get(i)))
                .collect(Collectors.toList());
    }

}
