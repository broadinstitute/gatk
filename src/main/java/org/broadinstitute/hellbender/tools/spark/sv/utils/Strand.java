package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.NoSuchElementException;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

// TODO htsjdk has it own Strand annotation enum. These classes could be merged if that Strand would me updated
// TODO so that one can get the enum constant char encoding; currently one can only do the transformation the other way.
@DefaultSerializer(Strand.Serializer.class)
public enum Strand {
    POSITIVE('+'),
    NEGATIVE('-');

    public static final Pattern PATTERN = Pattern.compile("\\+|\\-");

    private final char charEncoding;

    Strand(final char ce) {
        charEncoding = ce;
    }

    public static Strand decode(final String ce) {
        Utils.nonNull(ce);
        if (ce.length() == 1)
            return decode(ce.charAt(0));
        else
            throw new NoSuchElementException(String.format("there is no strand designation for encoding %s valid encodings are: %s.",
                    ce, Stream.of(values()).map(Strand::encodeAsString).collect(Collectors.joining(", "))));
    }

    public static Strand decode(final char ce) {
        if (ce == POSITIVE.charEncoding)
            return POSITIVE;
        else if (ce == NEGATIVE.charEncoding)
            return NEGATIVE;
        else
            throw new NoSuchElementException("there is no strand designation for encoding " + ce + " valid encodings are: " +
                    Stream.of(values()).map(s -> "" + s.charEncoding).collect(Collectors.joining(", ")) + ".");
    }

    @Override
    public String toString() { return encodeAsString(); };

    String encodeAsString() { return "" + charEncoding; }

    public void serialize(final Kryo kryo, final Output output) {
        output.write(this.ordinal());
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<Strand> {

        @Override
        public void write(final Kryo kryo, final Output output, final Strand strand) {
            strand.serialize(kryo, output);
        }

        @Override
        public Strand read(final Kryo kryo, final Input input, final Class<Strand> kclass) {
            return Strand.values()[input.readInt()];
        }
    }
}
