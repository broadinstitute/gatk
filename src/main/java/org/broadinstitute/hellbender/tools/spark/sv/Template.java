package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.TemplateFragmentOrdinal;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.io.Serializable;
import java.util.Collections;
import java.util.List;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.Function;
import java.util.stream.Collectors;
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
        public final String name;
        private final TemplateFragmentOrdinal number;
        public final String mapping;

        public Fragment(final String name, final TemplateFragmentOrdinal number, final String mapping, final byte[] bases, final int[] qualities) {
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
            return bases;
        }

        public int[] qualities() {
            return qualities;
        }

        public int length() {
            return length;
        }

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


}
