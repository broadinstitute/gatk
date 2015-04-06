package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.IOUtil;

import org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers.FileFaker;

import java.io.File;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public abstract class ParameterizedFileUtil {
    public static final String PER_TILE_PATTERN_STRING = "s_(\\d+)_(\\d{1,5})";
    /**
     * The file extension for this class, file extension does not have the standard meaning
     * in this instance.  It means, all the characters that come after the identifying portion of
     * the file (after lane, tile, and end that is).  So _qseq.txt and .filter are both file extensions
     */
    protected final String extension;

    /**
     * A pattern that will match files of this type for this lane
     */
    protected Pattern matchPattern;

    protected final int lane;
    protected List<Integer> tiles;
    /**
     * If you think of the file system as a tree, this is the deepest directory(node) on the tree that
     * still contains all of the files for this given type (e.g. If we're talking about BCLs the directory
     * structure is:
     * <p/>
     * BaseCall Dir
     * |
     * L001
     * |     |        |
     * C1.1 C2.1 ... Cn.1
     * |     |        |
     * bcl Files ... bclFiles
     * <p/>
     * L001 is the base because it contains every BCL file in the run (though those files are nested in
     * other folders).
     */
    protected final File base;
    protected final FileFaker faker;

    public ParameterizedFileUtil(final boolean laneTileRegex, final String extension, final File base,
                                 final FileFaker faker, final int lane) {
        this(extension, base, faker, lane);
        if (laneTileRegex) {
            matchPattern = Pattern.compile(escapePeriods(makeLaneTileRegex(processTxtExtension(extension), lane)));
        } else {
            matchPattern = Pattern.compile(escapePeriods(makeLaneRegex(extension, lane)));
        }
    }

    public ParameterizedFileUtil(final String pattern, final String extension, final File base, final FileFaker faker,
                                 final int lane) {
        this(extension, base, faker, lane);
        this.matchPattern = Pattern.compile(pattern);
    }

    private ParameterizedFileUtil(final String extension, final File base, final FileFaker faker,
                                  final int lane) {
        this.faker = faker;
        this.extension = extension;
        this.base = base;
        this.lane = lane;
    }

    /**
     * Determine whether or not files are available
     *
     * @return return true if files are found matching this types pattern, false otherwise
     */
    public abstract boolean filesAvailable();

    /**
     * Return a list of all tiles available for this file format and run
     *
     * @return A List of tile integers
     */
    public List<Integer> getTiles() {
        return tiles;
    }

    /**
     * Given the expected tiles/expected cycles for this file type, return a list of error messages describing any
     * missing/or malformed files
     *
     * @param expectedTiles  An ordered list of tile numbers
     * @param expectedCycles An ordered list of cycle numbers that may contain gaps
     * @return A list of error messages for this format
     */
    public abstract List<String> verify(List<Integer> expectedTiles, int[] expectedCycles);

    /**
     * Given the expected tiles/expected cycles for this file type create a set of fake files such that the
     * verification criteria are met.
     *
     * @param expectedTiles An ordered list of tile numbers
     * @param cycles        An ordered list of cycle numbers that may contain gaps
     * @param format        The format of the files that are to be faked
     * @return A list of error messages for this format
     */
    public abstract List<String> fakeFiles(List<Integer> expectedTiles, int[] cycles,
                                           IlluminaFileUtil.SupportedIlluminaFormat format);

    /**
     * Returns only lane and tile information as PerTileFt's do not have End information.
     *
     * @param fileName Filename to analyze for data
     * @return A LaneTile object with the discovered Lane and Tile information and a null end field.
     */
    protected Integer fileToTile(final String fileName) {
        final Matcher matcher = matchPattern.matcher(fileName);
        if (!matcher.matches()) {
            return null;
        }
        return Integer.parseInt(matcher.group(1));
    }

    /**
     * Return a regex string for finding Lane and Tile given a file extension pattern
     */
    public static String makeLaneTileRegex(final String fileNameEndPattern, final int lane) {
        if (lane < 0) {
            throw new IlluminaParserException("Lane (" + lane + ") cannot be negative");
        }
        return "^" + "s_" + lane + "_(\\d{1,5})" + fileNameEndPattern + "$";
    }

    private String makeLaneRegex(final String fileNameEndPattern, final int lane) {
        return "^s_" + lane + fileNameEndPattern + "$";
    }

    /**
     * The period separator is expected in the file extension, since some do not start with it
     */
    private String escapePeriods(final String preEscaped) {
        return preEscaped
                .replaceAll("\\.", "\\."); //In the first one the \\ is inside a regex in the second it's NOT
    }

    /**
     * For filename patterns that end with .txt tack on the option .gz extension
     */
    private String processTxtExtension(final String fileNameEndPattern) {
        if (fileNameEndPattern.endsWith(".txt")) {
            return fileNameEndPattern + "(\\.gz|\\.bz2)?";
        } else {
            return fileNameEndPattern;
        }
    }

    /**
     * Return all files that match pattern of the given file type in the given base directory
     */
    protected IlluminaFileMap getTiledFiles(final File baseDirectory, final Pattern pattern) {
        final IlluminaFileMap fileMap = new IlluminaFileMap();
        if (baseDirectory.exists()) {
            IOUtil.assertDirectoryIsReadable(baseDirectory);
            final File[] files = IOUtil.getFilesMatchingRegexp(baseDirectory, pattern);
            for (final File file : files) {
                if (file.length() > 0) {
                    fileMap.put(fileToTile(file.getName()), file);
                }
            }
        }

        return fileMap;
    }

}
