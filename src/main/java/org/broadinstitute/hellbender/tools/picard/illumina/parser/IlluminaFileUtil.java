package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers.*;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.TileMetricsOutReader;

import java.io.File;
import java.util.*;
import java.util.regex.Pattern;

import static htsjdk.samtools.util.CloserUtil.close;
import static htsjdk.samtools.util.IOUtil.assertFileIsReadable;
import static java.lang.String.valueOf;
import static java.util.regex.Pattern.compile;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaFileUtil.SupportedIlluminaFormat.*;

/**
 * General utils for dealing with IlluminaFiles as well as utils for specific, support formats.
 * This class contains utils that span across multiple Illumina files but it's primary intent
 * was to provide support for basic file types.  Each supported file type can be accessed
 * via a factory method (make<filetype>Ft).  When IlluminaFileUtil is created it is parameterized
 * by basecallDir and lane and all IlluminaFileTypes created by IlluminaFileUtil will also be
 * parameterized in this fashion.
 *
 * @author jburke@broadinstitute.org
 */
public class IlluminaFileUtil {
    public static final Pattern CYCLE_SUBDIRECTORY_PATTERN = compile("^C(\\d+)\\.1$");

    public enum SupportedIlluminaFormat {
        Bcl,
        Locs,
        Clocs,
        Pos,
        Filter,
        Barcode,
        MultiTileFilter,
        MultiTileLocs,
        MultiTileBcl
    }

    private final File basecallLaneDir;
    private final File intensityLaneDir;
    private final File basecallDir;
    private final File barcodeDir;
    private final File intensityDir;
    private final int lane;

    private final File tileMetricsOut;
    private final Map<SupportedIlluminaFormat, ParameterizedFileUtil> utils = new HashMap<>();

    public IlluminaFileUtil(final File basecallDir, final int lane) {
        this(basecallDir, null, lane);
    }


    public IlluminaFileUtil(final File basecallDir, File barcodeDir, final int lane) {
        this.lane = lane;
        this.basecallDir = basecallDir;
        this.barcodeDir = barcodeDir;
        this.intensityDir = basecallDir.getParentFile();
        final File dataDir = intensityDir.getParentFile();
        this.basecallLaneDir = new File(basecallDir, longLaneStr(lane));
        this.intensityLaneDir = new File(intensityDir, longLaneStr(lane));
        final File interopDir = new File(dataDir.getParentFile(), "InterOp");
        tileMetricsOut = new File(interopDir, "TileMetricsOut.bin");
    }


    /**
     * Return the lane we're inspecting
     */
    public int getLane() {
        return lane;
    }

    /**
     * Given a file type, get the Parameterized File Util object associated with it
     */
    public ParameterizedFileUtil getUtil(final SupportedIlluminaFormat format) {
        ParameterizedFileUtil parameterizedFileUtil = utils.get(format);
        if (parameterizedFileUtil == null) {
            switch (format) {
                case Bcl:
                    final ParameterizedFileUtil bclFileUtil = new PerTilePerCycleFileUtil(".bcl", basecallLaneDir, new BclFileFaker(), lane);
                    final ParameterizedFileUtil gzBclFileUtil = new PerTilePerCycleFileUtil(".bcl.gz", basecallLaneDir, new BclFileFaker(), lane);
                    if (bclFileUtil.filesAvailable() && !gzBclFileUtil.filesAvailable()) {
                        parameterizedFileUtil = bclFileUtil;
                    } else if (!bclFileUtil.filesAvailable() && gzBclFileUtil.filesAvailable()) {
                        parameterizedFileUtil = gzBclFileUtil;
                    } else if (!bclFileUtil.filesAvailable() && !gzBclFileUtil.filesAvailable()) {
                        parameterizedFileUtil = bclFileUtil;
                    } else {
                        throw new IlluminaParserException("Not all BCL files in " + basecallLaneDir.getAbsolutePath() + " have the same extension!");
                    }
                    utils.put(Bcl, parameterizedFileUtil);
                    break;
                case Locs:
                    parameterizedFileUtil = new PerTileFileUtil(".locs", intensityLaneDir, new LocsFileFaker(), lane);
                    utils.put(Locs, parameterizedFileUtil);
                    break;
                case Clocs:
                    parameterizedFileUtil = new PerTileFileUtil(".clocs", intensityLaneDir, new ClocsFileFaker(), lane);
                    utils.put(Clocs, parameterizedFileUtil);
                    break;
                case Pos:
                    parameterizedFileUtil = new PerTileFileUtil("_pos.txt", intensityDir, new PosFileFaker(), lane);
                    utils.put(Pos, parameterizedFileUtil);
                    break;
                case Filter:
                    parameterizedFileUtil = new PerTileFileUtil(".filter", basecallLaneDir, new FilterFileFaker(), lane);
                    utils.put(Filter, parameterizedFileUtil);
                    break;
                case Barcode:
                    parameterizedFileUtil = new PerTileFileUtil("_barcode.txt", barcodeDir != null ? barcodeDir : basecallDir, new BarcodeFileFaker(), lane);
                    utils.put(Barcode, parameterizedFileUtil);
                    break;
                case MultiTileFilter:
                    parameterizedFileUtil = new MultiTileFilterFileUtil(basecallLaneDir, lane);
                    utils.put(MultiTileFilter, parameterizedFileUtil);
                    break;
                case MultiTileLocs:
                    parameterizedFileUtil = new MultiTileLocsFileUtil(new File(intensityDir, basecallLaneDir.getName()), basecallLaneDir, lane);
                    utils.put(MultiTileLocs, parameterizedFileUtil);
                    break;
                case MultiTileBcl:
                    parameterizedFileUtil = new MultiTileBclFileUtil(basecallLaneDir, lane);
                    utils.put(MultiTileBcl, parameterizedFileUtil);
                    break;
            }
        }
        return parameterizedFileUtil;
    }

    /**
     * Return the list of tiles we would expect for this lane based on the metrics found in InterOp/TileMetricsOut.bin
     */
    public List<Integer> getExpectedTiles() {
        assertFileIsReadable(tileMetricsOut);
        //Used just to ensure predictable ordering
        final TreeSet<Integer> expectedTiles = new TreeSet<>();

        final Iterator<TileMetricsOutReader.IlluminaTileMetrics> tileMetrics = new TileMetricsOutReader(tileMetricsOut);
        while (tileMetrics.hasNext()) {
            final TileMetricsOutReader.IlluminaTileMetrics tileMetric = tileMetrics.next();

            if (tileMetric.getLaneNumber() == lane) {
                if (!expectedTiles.contains(tileMetric.getTileNumber())) {
                    expectedTiles.add(tileMetric.getTileNumber());
                }
            }
        }

        close(tileMetrics);
        return new ArrayList<>(expectedTiles);
    }

    /**
     * Get the available tiles for the given formats, if the formats have tile lists that differ then
     * throw an exception, if any of the format
     */
    public List<Integer> getActualTiles(final List<SupportedIlluminaFormat> formats) {
        if (formats == null) {
            throw new IlluminaParserException("Format list provided to getTiles was null!");
        }

        if (formats.size() == 0) {
            throw new IlluminaParserException("0 Formats were specified.  You need to specify at least SupportedIlluminaFormat to use getTiles");
        }

        final List<Integer> tiles = getUtil(formats.get(0)).getTiles();
        for (int i = 1; i < formats.size(); i++) {
            final List<Integer> fmTiles = getUtil(formats.get(i)).getTiles();
            if (tiles.size() != fmTiles.size() || !tiles.containsAll(fmTiles)) {
                throw new IlluminaParserException("Formats do not have the same number of tiles! " + summarizeTileCounts(formats));
            }
        }

        return tiles;
    }

    /*
     * Return a string representing the Lane in the format "L00<lane>"
     *
     * @param lane The lane to transform
     * @return A long string representation of the name
     */
    public static String longLaneStr(final int lane) {
        String lstr = valueOf(lane);
        final int zerosToAdd = 3 - lstr.length();

        for (int i = 0; i < zerosToAdd; i++) {
            lstr = "0" + lstr;
        }
        return "L" + lstr;
    }


    private String liToStr(final List<Integer> intList) {
        if (intList.size() == 0) {
            return "";
        }

        String summary = valueOf(intList.get(0));
        for (int i = 1; i < intList.size(); i++) {
            summary += ", " + valueOf(intList.get(i));
        }

        return summary;
    }

    private String summarizeTileCounts(final List<SupportedIlluminaFormat> formats) {
        String summary;
        ParameterizedFileUtil pfu = getUtil(formats.get(0));
        List<Integer> tiles = pfu.getTiles();
        summary = pfu.extension + "(" + liToStr(tiles) + ")";

        for (final SupportedIlluminaFormat format : formats) {
            pfu = getUtil(format);
            tiles = pfu.getTiles();

            summary += ", " + pfu.extension + "(" + liToStr(tiles) + ")";
        }

        return summary;
    }
}
