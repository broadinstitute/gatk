package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.TemplateFragmentOrdinal;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalLocator;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.ArrayUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.report.GATKReportColumnFormat;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Represent a template sequence as reconstructed from the corresponding
 * fragments such as reads pairs.
 */
public class Template implements Serializable {

    private static final long serialVersionUID = 1L;

    private List<Fragment> fragments;

    private String name;

    public static class Fragment implements Serializable {

        private static final long serialVersionUID = 1L;

        private final byte[] bases;
        private final int[] qualities;
        private final int length;
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
            this.length = bases.length;
            if (this.length != qualities.length) {
                throw new IllegalArgumentException("the input bases and qualities must have the same length");
            }
            this.number = number;
        }

        public byte[] bases() {
            return bases.clone();
        }

        public int[] qualities() {
            return qualities.clone();
        }

        public int length() {
            return length;
        }

        public List<AlignmentInterval> alignmentIntervals() { return this.mapping; }

        public GATKRead toUnmappedRead(final SAMFileHeader header, final boolean paired) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadName(name);
            record.setBaseQualities(encodeQuals());
            record.setReadBases(bases);
            record.setReadUnmappedFlag(true);
            if (paired) {
                record.setReadPairedFlag(true);
                record.setMateUnmappedFlag(true);
                record.setFirstOfPairFlag(number == TemplateFragmentOrdinal.PAIRED_FIRST);
                record.setSecondOfPairFlag(number == TemplateFragmentOrdinal.PAIRED_SECOND);
            }
            record.setValidationStringency(ValidationStringency.STRICT);
            return SAMRecordToGATKReadAdapter.headerlessReadAdapter(record);
        }

        private byte[] encodeQuals() {
            final byte[] result = new byte[qualities.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = (byte) qualities[i];
            }
            return result;
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


    public int[] maximumMappingQuality(final SVIntervalTree<?> targetIntervals, SVIntervalLocator locator,
                                       final InsertSizeDistribution insertSizeDistribution) {
        final int[] result = new int[fragments.size()];
        int definedMappingQualities  = 0;
        int maxQual = 0;
        for (int i = 0; i < result.length; i++) {
            final Fragment fragment = fragments.get(i);
            result[i] = -1;
            for (final AlignmentInterval mappingInterval : fragment.alignmentIntervals()) {
                if (!targetIntervals.hasOverlapper(locator.toSVInterval(mappingInterval.referenceSpan))) {
                    result[i] = 0;
                    definedMappingQualities++;
                    continue;
                } else {
                    if (result[i] == -1) {
                        definedMappingQualities++;
                    }
                    result[i] = Math.max(result[i], mappingInterval.mapQual);
                    if (result[i] > maxQual) {
                        maxQual = result[i];
                    }
                }
            }
        }
        if (definedMappingQualities == result.length) {
            return result;
        } else if (definedMappingQualities == 0) {
            Arrays.fill(result, 10);
            return result;
        } else if (maxQual == 0) {
            Arrays.fill(result, 0);
            return result;
        } else {
            int minGap = Integer.MAX_VALUE;
            final int readLength = fragments.stream().mapToInt(Fragment::length).min().getAsInt();
            for (int i = 0; i < result.length; i++) {
                if (result[i] <= 0) continue;
                minGap = Math.min(minGap, fragments.get(i).alignmentIntervals().stream()
                                            .mapToInt(ai -> targetIntervals.smallestGapLength(locator.toSVInterval(ai.referenceSpan)))
                                            .min().getAsInt());
            }
            final double logDistanceProbability = insertSizeDistribution.logCumulativeProbability(minGap + 2 * readLength);
            final double distancePenalty = -10.0 * MathUtils.log1mexp(logDistanceProbability) / Math.log(10);
            final int unmappedMappingQuality = (int) Math.max(0, maxQual - distancePenalty);
            for (int i = 0; i < result.length; i++) {
                if (result[i] == -1) { result[i] = unmappedMappingQuality; }
            }
            return result;
        }

    }


}
