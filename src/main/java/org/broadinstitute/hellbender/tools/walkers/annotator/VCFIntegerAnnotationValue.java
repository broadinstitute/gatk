package org.broadinstitute.hellbender.tools.walkers.annotator;

import org.broadinstitute.hellbender.utils.Nucleotide;

import java.util.List;
import java.util.NoSuchElementException;

public interface VCFIntegerAnnotationValue {

    int length();

    default long getFirst() {
        if (length() < 1) {
            throw new NoSuchElementException("");
        } else {
            return get(0);
        }
    }

    default long getOnly() {
        if (length() != 1) {
            throw new IllegalStateException("");
        } else {
            return get(0);
        }
    }

    long get(int index);

    default void copyRange(final long[] dest, final int offset, final int firstIndex, final int afterLastIndex) {
        final int length = length();
        for (int i = offset, j = firstIndex; j < afterLastIndex; j++, i++) {
            dest[i] = get(j);
        }
    }

    default void copyRange(final int[] dest, final int offset, final int firstIndex, final int afterLastIndex) {
        final int length = length();
        for (int i = offset, j = firstIndex; j < afterLastIndex; j++, i++) {
            dest[i] = (int) get(j);
        }
    }

    default long[] toArray() {
        final int length = length();
        final long[] result = new long[length];
        for (int i = 0; i < length; i++) {
            result[i] = get(i);
        }
        return result;
    }

    default int[] toIntArray() {
        final int length = length();
        final int[] result = new int[length];
        for (int i = 0; i < length; i++) {
            result[i] = (int) get(i);
        }
        return result;
    }


    static VCFIntegerAnnotationValue empty() {
        return new VCFIntegerAnnotationValue() {
            @Override
            public int length() { return 1; }

            @Override
            public long get(int index) {
                throw new NoSuchElementException();
            }
        };
    }

    static VCFIntegerAnnotationValue of(final long value) {
        return new VCFIntegerAnnotationValue() {
            @Override
            public int length() {
                return 1;
            }

            @Override
            public long get(int index) {
                return value;
            }
        };
    }

    static VCFIntegerAnnotationValue of(final long[] values) {
        return new VCFIntegerAnnotationValue() {
            @Override
            public int length() {
                return values.length;
            }

            @Override
            public long get(int index) {
                return values[index];
            }
        };
    }

    static VCFIntegerAnnotationValue of(final int[] values) {
        return new VCFIntegerAnnotationValue() {
            @Override
            public int length() {
                return values.length;
            }

            @Override
            public long get(int index) {
                return values[index];
            }
        };
    }

    static VCFIntegerAnnotationValue of(final List<?> values) {
        return new VCFIntegerAnnotationValue() {
            @Override
            public int length() {
                return values.size();
            }

            @Override
            public long get(int index) {
                final Object obj = values.get(index);
                if (obj instanceof Number) {
                    return ((Number)obj).longValue();
                } else {
                    return Long.parseLong(obj.toString());
                }
            }
        };
    }

    static VCFIntegerAnnotationValue of(final String values) {
        final String[] array = values.split(",");
        return new VCFIntegerAnnotationValue() {
            @Override
            public int length() {
                return array.length;
            }

            @Override
            public long get(int index) {
                return Long.parseLong(array[index]);
            }
        };
    }

    static VCFIntegerAnnotationValue of(final Object obj) {
        final Class<?> clazz = obj.getClass();
        if (clazz.isArray()) {

        } else {
            
        }
    }

}
