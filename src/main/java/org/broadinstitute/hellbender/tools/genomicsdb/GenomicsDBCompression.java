package org.broadinstitute.hellbender.tools.genomicsdb;

import org.broadinstitute.barclay.argparser.CommandLineException;

import java.util.Arrays;
import java.util.Locale;
import java.util.stream.Collectors;

/**
 * The codec and level that GenomicsDB uses to compress the tiles of a new workspace, parsed from a specification of
 * the form {@code <codec>} or {@code <codec>:<level>}. A codec given without a level uses level 1: the fastest gzip
 * and zstd level, and lz4's default mode.
 *
 * <p>The codec is recorded in the workspace's array schema, so it only needs to be chosen when a workspace is
 * created; readers need no option. Zstandard is loaded at runtime from the system's {@code libzstd}, so it must be
 * installed on every machine that imports into or reads a workspace compressed with it.</p>
 *
 * @param codec the codec
 * @param level the level passed to the codec (0 for {@link Codec#NONE})
 */
record GenomicsDBCompression(Codec codec, int level) {

    /** The level used when a specification names only a codec. */
    static final int DEFAULT_LEVEL = 1;

    /**
     * The codecs GenomicsDB's TileDB can use for tiles, with TileDB's numeric code and accepted levels for each.
     * TileDB stores a level in the array schema as one signed byte, so no level may exceed 127.
     */
    enum Codec {
        NONE(0, 0, 0),
        GZIP(1, 1, 9),
        ZSTD(2, 1, 22),
        /** For lz4 the level is the acceleration of LZ4's fast mode: 1 is LZ4's default mode, and higher is faster. */
        LZ4(3, 1, 127);

        private final int tiledbCode;
        private final int minLevel;
        private final int maxLevel;

        Codec(final int tiledbCode, final int minLevel, final int maxLevel) {
            this.tiledbCode = tiledbCode;
            this.minLevel = minLevel;
            this.maxLevel = maxLevel;
        }

        /** The numeric codec identifier that GenomicsDB's import configuration expects. */
        int tiledbCode() {
            return tiledbCode;
        }
    }

    /**
     * Parses a compression specification.
     *
     * @param argumentName the argument the specification came from, for error messages
     * @param spec {@code <codec>} or {@code <codec>:<level>}; the codec name is case-insensitive
     * @return the codec and level, with level {@link #DEFAULT_LEVEL} when none is given (0 for {@code none})
     * @throws CommandLineException.BadArgumentValue if the codec is unknown, the level is not an integer in the
     *         codec's range, or a level is given for {@code none}
     */
    static GenomicsDBCompression parse(final String argumentName, final String spec) {
        final String[] parts = spec.split(":", -1);
        if (parts.length > 2) {
            throw new CommandLineException.BadArgumentValue(argumentName, spec, "expected <codec> or <codec>:<level>");
        }
        final Codec codec;
        try {
            codec = Codec.valueOf(parts[0].trim().toUpperCase(Locale.ROOT));
        } catch (final IllegalArgumentException e) {
            throw new CommandLineException.BadArgumentValue(argumentName, spec, "the codec must be one of " + codecNames());
        }
        if (parts.length == 1) {
            return new GenomicsDBCompression(codec, codec == Codec.NONE ? 0 : DEFAULT_LEVEL);
        }
        if (codec == Codec.NONE) {
            throw new CommandLineException.BadArgumentValue(argumentName, spec, "none takes no level");
        }
        final String levelOutOfRange = String.format("the level for %s must be an integer between %d and %d",
                codecName(codec), codec.minLevel, codec.maxLevel);
        final int level;
        try {
            level = Integer.parseInt(parts[1].trim());
        } catch (final NumberFormatException e) {
            throw new CommandLineException.BadArgumentValue(argumentName, spec, levelOutOfRange);
        }
        if (level < codec.minLevel || level > codec.maxLevel) {
            throw new CommandLineException.BadArgumentValue(argumentName, spec, levelOutOfRange);
        }
        return new GenomicsDBCompression(codec, level);
    }

    private static String codecName(final Codec codec) {
        return codec.name().toLowerCase(Locale.ROOT);
    }

    private static String codecNames() {
        return Arrays.stream(Codec.values()).map(GenomicsDBCompression::codecName).collect(Collectors.joining(", "));
    }
}
