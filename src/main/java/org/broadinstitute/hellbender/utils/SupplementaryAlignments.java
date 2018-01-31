package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMTag;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.AbstractList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by valentin on 9/22/17.
 */
public final class SupplementaryAlignments extends AbstractList<AlignmentInterval> {

    private static final SupplementaryAlignments EMPTY = new SupplementaryAlignments(Collections.emptyList());

    private final List<AlignmentInterval> intervals;

    private SupplementaryAlignments(final List<AlignmentInterval> intervals) {
        this.intervals = intervals;
    }

    public static SupplementaryAlignments of(final GATKRead read) {
        return of(read, SAMTag.SA.name());
    }

    public static SupplementaryAlignments of(final GATKRead read, final String tag) {
        Utils.nonNull(read);
        Utils.nonNull(tag);
        final String str = read.getAttributeAsString(tag);
        if (str == null) {
            return EMPTY;
        } else {
            final String[] parts = str.split("\\s*;\\s*");
            return new SupplementaryAlignments(Stream.of(parts)
                    .filter(s -> !s.isEmpty())
                    .map(AlignmentInterval::new)
                    .collect(Collectors.toList()));
        }
    }

    @Override
    public AlignmentInterval get(final int index) {
        return intervals.get(index);
    }

    @Override
    public int size() {
        return intervals.size();
    }
}
