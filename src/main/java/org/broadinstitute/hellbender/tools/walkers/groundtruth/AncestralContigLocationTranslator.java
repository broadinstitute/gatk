package org.broadinstitute.hellbender.tools.walkers.groundtruth;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.Tuple;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

 class AncestralContigLocationTranslator {

    // locals
    final private GATKPath                                    basePath;
    final private Map<String, SingleFileLocationTranslator>    translators = new LinkedHashMap<>();

    AncestralContigLocationTranslator(GATKPath basePath) {
        this.basePath = basePath;
    }

    protected Tuple<SimpleInterval, SimpleInterval> translate(final Locatable loc) throws IOException {
        return new Tuple<>(translate(GroundTruthReadsBuilder.C_MATERNAL, loc),
                translate(GroundTruthReadsBuilder.C_PATERNAL, loc));
    }

    private SimpleInterval translate(final String ancestor, final Locatable loc) throws IOException {

        int         start = translate(ancestor, loc.getContig(), loc.getStart());
        int         end = translate(ancestor, loc.getContig(), loc.getEnd());

        if ( end > start ) {
            return new SimpleInterval(loc.getContig() + "_" + ancestor, start, end);
        } else {
            throw new LocationTranslationException("location " + loc + " failed to translate for " + ancestor + ", start:" + start + " ,end:" + end);
        }
    }

    private int translate(final String ancestor, final String contig, final int from) throws IOException {

        // check-for/create translator
        final String                          key = ancestor + "." + contig + ".csv";
        if ( !translators.containsKey(key) ) {
            final GATKPath        path = new GATKPath(basePath.getURIString() + key);
            translators.put(key, new SingleFileLocationTranslator(path));
        }

        // translate
        return translators.get(key).translate(from);
    }
}
