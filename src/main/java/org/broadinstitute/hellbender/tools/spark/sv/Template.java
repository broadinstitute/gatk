package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

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
public class Template {

    private List<Fragment> fragments;

    private String name;

    public static class Fragment {

        private final byte[] bases;
        private final int[] qualities;
        private final int length;

        public Fragment(final byte[] bases, final int[] qualities) {
            this.bases = Utils.nonNull(bases);
            this.qualities = Utils.nonNull(qualities);
            this.length = bases.length;
            if (this.length != qualities.length) {
                throw new IllegalArgumentException("the input bases and qualities must have the same length");
            }
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
    }

    private Template(final String name, final List<Fragment> fragments) {
        this.name = name;
        this.fragments = fragments;
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


}
