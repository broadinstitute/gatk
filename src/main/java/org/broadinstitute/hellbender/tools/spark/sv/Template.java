package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.TemplateFragmentOrdinal;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalLocator;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.ArrayUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Represent a template sequence as reconstructed from the corresponding
 * fragments such as reads pairs.
 */
public final class Template implements Serializable {

    private static final long serialVersionUID = 1L;

    private List<Fragment> fragments;

    private String name;

    public static class Fragment implements Serializable {

        private static final long serialVersionUID = 1L;

        private final byte[] bases;
        private final int[] qualities;
        private final int length;
        private int mappingQuality;
        private final String name;
        private final TemplateFragmentOrdinal number;
        private final List<AlignmentInterval> mapping;

        public static Fragment of(final SVFastqUtils.FastqRead read) {
            return new Fragment(read.getName(), read.getFragmentOrdinal(),
                    Collections.unmodifiableList(new ArrayList<>(read.getMapping().getAllIntervals())),
                    read.getBases(), ArrayUtils.toInts(read.getQuals(), false));
        }

        private Fragment(final String name, final TemplateFragmentOrdinal number, final List<AlignmentInterval> mapping, final byte[] bases, final int[] qualities) {
            this.name = name;
            this.mapping = mapping;
            this.bases = Utils.nonNull(bases);
            this.qualities = Utils.nonNull(qualities);
            this.mappingQuality = SAMRecord.UNKNOWN_MAPPING_QUALITY;
            this.length = bases.length;
            if (this.length != qualities.length) {
                throw new IllegalArgumentException("the input bases and qualities must have the same length");
            }
            this.number = number;
        }

        public byte[] bases() {
            return bases.clone();
        }

        public void copyBases(final byte[] dest, final int offset, final int from, final int to) {
            Utils.nonNull(dest);
            Utils.validate(from <= to, "from cannot be larger than to");
            final int copyLength = to - from;
            Utils.validate(copyLength + offset < dest.length, "the copy will go beyond the end of the destination array");
            Utils.validIndex(offset, dest.length);
            ParamUtils.isPositiveOrZero(offset, "the destination offset cannot be negative");
            System.arraycopy(bases, from, dest, offset, copyLength);
        }

        public int[] qualities() {
            return qualities.clone();
        }

        public int length() {
            return length;
        }

        public List<AlignmentInterval> alignmentIntervals() { return this.mapping; }

        public int getMappingQuality() {
            return mappingQuality;
        }

        public void setMappingQuality(final int mq) {
            mappingQuality = mq;
        }
    }

    private Template(final String name, final List<Fragment> fragments) {
        this.name = name;
        this.fragments = Collections.unmodifiableList(fragments);
    }

    public static <T> Template create(final String name, final Iterable<T> fragmentPrecursors, final Function<T, Fragment> fragmentTranslator) {
        Utils.nonNull(name);
        Utils.nonNull(fragmentPrecursors);
        Utils.nonNull(fragmentTranslator);
        final List<Fragment> fragments = StreamSupport.stream(Spliterators.spliteratorUnknownSize(fragmentPrecursors.iterator(), Spliterator.NONNULL), false)
                .map(fragmentTranslator)
                .collect(Collectors.toList());
        return new Template(name, fragments);
    }

    public List<Fragment> fragments() {
        return fragments;
    }

    public String name() {
        return name;
    }

    public boolean equals(final Template other) {
        return name.equals(other.name);
    }

    public int hashCode() {
        return name.hashCode();
    }

    public GATKRead asUnmappedRead(final SAMFileHeader header) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name());
        record.setReadUnmappedFlag(true);
        return new SAMRecordToGATKReadAdapter(record);
    }

    public void calculateMaximumMappingQualities(final SVIntervalTree<?> targetIntervals, SVIntervalLocator locator,
                                                 final InsertSizeDistribution insertSizeDistribution) {
        if (fragments.isEmpty()) {
            return;
        }
        int definedMappingQualities = 0;
        for (int i = 0; i < fragments.size(); i++) {
            int maxInZoneQual = -1;
            final Fragment fragment = fragments.get(i);
            for (final AlignmentInterval mappingInterval : fragment.alignmentIntervals()) {
                if (mappingInterval.mapQual == SAMRecord.UNKNOWN_MAPPING_QUALITY) continue;
                //final int gap = targetIntervals.smallestGapLength(locator.toSVInterval(mappingInterval.referenceSpan));
                //final boolean inZone = gap < insertSizeDistribution.maximum() + fragment.length;
                final boolean inZone = targetIntervals.hasOverlapper(locator.toSVInterval(mappingInterval.referenceSpan));
                if (inZone) {
                    maxInZoneQual = Math.max(maxInZoneQual, mappingInterval.mapQual);
                } else {
                    maxInZoneQual = 0; break;
                }
            }
            if (maxInZoneQual == -1) {
                fragment.setMappingQuality(maxInZoneQual);
            } else {
                definedMappingQualities++;
                fragment.setMappingQuality(Math.max(0, maxInZoneQual));
            }
        }
        if (definedMappingQualities == 0) {
            for (final Fragment fragment : fragments) {
                fragment.setMappingQuality(60); // default qual for unmapped fragments.
            }
        } else if (definedMappingQualities < fragments.size()) {
            final int maxQual = fragments.stream().mapToInt(Fragment::getMappingQuality)
                .map(qual -> qual == SAMRecord.UNKNOWN_MAPPING_QUALITY ? -1 : qual).max().orElse(0);
            for (final Fragment fragment: fragments) {
                if (fragment.getMappingQuality() == SAMRecord.UNKNOWN_MAPPING_QUALITY) {
                    fragment.setMappingQuality(maxQual);
                }
            }
        }
    }
}
