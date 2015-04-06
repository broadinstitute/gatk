package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.IOUtil;

import org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers.BciFileFaker;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers.FileFaker;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers.FilterFileFaker;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers.MultiTileLocsFileFaker;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * For file types for which there is one file per lane, with fixed record size, and all the tiles in it,
 * so the s_<lane>.bci file can be used to figure out where each tile starts and ends.
 */
public abstract class MultiTileFileUtil<OUTPUT_RECORD extends IlluminaData> extends ParameterizedFileUtil {
    protected final File bci;
    protected TileIndex tileIndex;
    protected File dataFile;

    MultiTileFileUtil(final String extension, final File base, final File bciDir, final FileFaker fileFaker,
                      final int lane) {
        super(false, extension, base, fileFaker, lane);
        bci = new File(bciDir, "s_" + lane + ".bci");
        if (bci.exists()) {
            tileIndex = new TileIndex(bci);
        } else {
            tileIndex = null;
        }
        final File[] filesMatchingRegexp = IOUtil.getFilesMatchingRegexp(base, matchPattern);
        if (filesMatchingRegexp == null || filesMatchingRegexp.length == 0) {
            dataFile = null;
        } else if (filesMatchingRegexp.length == 1) {
            dataFile = filesMatchingRegexp[0];
        } else {
            throw new IlluminaParserException("More than one filter file found in " + base.getAbsolutePath());
        }
    }

    @Override
    public boolean filesAvailable() {
        return tileIndex != null && dataFile != null && dataFile.exists();
    }

    @Override
    public List<Integer> getTiles() {
        if (tileIndex == null) {
            return Collections.emptyList();
        }
        return tileIndex.getTiles();
    }

    /**
     * expectedCycles are not checked in this implementation.
     */
    @Override
    public List<String> verify(final List<Integer> expectedTiles, final int[] expectedCycles) {
        if (tileIndex == null) {
            return Collections.singletonList("Tile index(" + bci.getAbsolutePath() + ") does not exist!");
        }
        return tileIndex.verify(expectedTiles);
    }

    @Override
    public List<String> fakeFiles(final List<Integer> expectedTiles, final int[] expectedCycles,
                                  final IlluminaFileUtil.SupportedIlluminaFormat format) {
        //we need to fake a bci file for the tile index
        final BciFileFaker bciFileFaker = new BciFileFaker();
        try {
            bciFileFaker.fakeBciFile(bci, expectedTiles);
            tileIndex = new TileIndex(bci);
            faker.fakeFile(base, expectedTiles, lane, extension);
            final File[] filesMatchingRegexp = IOUtil.getFilesMatchingRegexp(base, matchPattern);
            if (filesMatchingRegexp == null || filesMatchingRegexp.length == 0) {
                dataFile = null;
            } else if (filesMatchingRegexp.length == 1) {
                dataFile = filesMatchingRegexp[0];
            } else {
                throw new IlluminaParserException("More than one filter file found in " + base.getAbsolutePath());
            }
        } catch (final IOException e) {
            return Collections.singletonList("Could not create tile index file: " + bci.getAbsolutePath());
        }
        return tileIndex.verify(expectedTiles);
    }

    abstract IlluminaParser<OUTPUT_RECORD> makeParser(List<Integer> requestedTiles);

}

class MultiTileFilterFileUtil extends MultiTileFileUtil<PfData> {

    /**
     * @param basecallLaneDir location of .filter file and also .bci file
     */
    MultiTileFilterFileUtil(final File basecallLaneDir, final int lane) {
        super(".filter", basecallLaneDir, basecallLaneDir, new FilterFileFaker(), lane);
    }

    @Override
    IlluminaParser<PfData> makeParser(final List<Integer> requestedTiles) {
        return new MultiTileFilterParser(tileIndex, requestedTiles, dataFile);
    }
}

class MultiTileLocsFileUtil extends MultiTileFileUtil<PositionalData> {

    MultiTileLocsFileUtil(final File basecallLaneDir, final File bciDir, final int lane) {
        super(".locs", basecallLaneDir, bciDir, new MultiTileLocsFileFaker(), lane);
    }

    @Override
    IlluminaParser<PositionalData> makeParser(final List<Integer> requestedTiles) {
        return new MultiTileLocsParser(tileIndex, requestedTiles, dataFile, lane);
    }
}

