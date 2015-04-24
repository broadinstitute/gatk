package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.IOUtil;

import org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers.FileFaker;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclReader;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;

public class PerTilePerCycleFileUtil extends ParameterizedFileUtil {

    private final CycleIlluminaFileMap cycleFileMap;
    private final Set<Integer> detectedCycles = new TreeSet<>();

    public PerTilePerCycleFileUtil(final String extension,
                                   final File base, final FileFaker faker, final int lane) {
        super(true, extension, base, faker, lane);
        //sideEffect, assigned to numCycles
        this.cycleFileMap = getPerTilePerCycleFiles();
    }

    /**
     * For the given tiles, populate a CycleIlluminaFileMap that contains all these tiles and will iterate through
     * all the files for these tiles in expectedBase
     * Side Effect: Assigns numCycles
     *
     * @return A CycleIlluminaFileMap with the listed (or all) tiles for at least expectedCycles number of cycles(or total available
     * cycles if expectedCycles is null)
     */
    protected CycleIlluminaFileMap getPerTilePerCycleFiles() {
        final CycleIlluminaFileMap cycledMap = new CycleIlluminaFileMap();

        final File laneDir = base;
        final File[] tempCycleDirs;
        tempCycleDirs = IOUtil.getFilesMatchingRegexp(laneDir, IlluminaFileUtil.CYCLE_SUBDIRECTORY_PATTERN);
        if (tempCycleDirs == null || tempCycleDirs.length == 0) {
            return cycledMap;
        }

        for (final File tempCycleDir : tempCycleDirs) {
            detectedCycles.add(getCycleFromDir(tempCycleDir));
        }

        final Set<Integer> uniqueTiles = new HashSet<>();

        for (final File cycleDir : tempCycleDirs) {
            final IlluminaFileMap fileMap = getTiledFiles(cycleDir, matchPattern);
            uniqueTiles.addAll(fileMap.keySet());
            cycledMap.put(getCycleFromDir(cycleDir), fileMap);
        }

        this.tiles = Collections.unmodifiableList(new ArrayList<>(uniqueTiles));
        return cycledMap;
    }

    public CycleIlluminaFileMap getFiles() {
        return cycleFileMap;
    }

    public CycleIlluminaFileMap getFiles(final List<Integer> tiles) {
        return cycleFileMap.keep(tiles, detectedCycles);
    }

    /**
     * Returns a cycleIlluminaFileMap with all available tiles but limited to the cycles passed in.  Any cycles that are missing
     * cycle files or directories will be removed from the cycle list that is kept.
     *
     * @param cycles Cycles that should be present in the output CycleIlluminaFileMap
     * @return A CycleIlluminaFileMap with all available tiles but at most the cycles passed in by the cycles parameter
     */
    public CycleIlluminaFileMap getFiles(final int[] cycles) {
        //Remove any cycles that were discovered to be NON-EXISTENT when this util was instantiated
        final Set<Integer> filteredCycles = removeNonExistentCycles(cycles);
        return cycleFileMap.keep(tiles, filteredCycles);
    }

    /**
     * Returns a cycleIlluminaFileMap that contains only the tiles and cycles specified (and fewer if the original CycleIlluminaFileMap, created
     * on util instantiation, doesn't contain any of these tiles/cycles).
     *
     * @param cycles Cycles that should be present in the output CycleIlluminaFileMap
     * @return A CycleIlluminaFileMap with at most the tiles/cycles listed in the parameters
     */
    public CycleIlluminaFileMap getFiles(final List<Integer> tiles, final int[] cycles) {
        //Remove any cycles that were discovered to be NON-EXISTENT when this util was instantiated
        final Set<Integer> filteredCycles = removeNonExistentCycles(cycles);
        return cycleFileMap.keep(tiles, filteredCycles);
    }

    private Set<Integer> removeNonExistentCycles(final int[] cycles) {

        final TreeSet<Integer> inputCyclesSet = new TreeSet<>();
        for (final Integer inputCycle : cycles) {
            inputCyclesSet.add(inputCycle);
        }

        inputCyclesSet.retainAll(detectedCycles);

        return inputCyclesSet;
    }

    public Set<Integer> getDetectedCycles() {
        return detectedCycles;
    }

    /**
     * Discover all files of this type in expectedBase that match pattern and construct a list of tiles
     * available based on these files.  The same number of tiles is expected in each cycle dir.
     *
     * @return A list of tile integers for all tiles available
     */
    public List<Integer> getTiles() {
        return tiles;
    }

    public boolean filesAvailable() {
        boolean filesAvailable = false;
        for (final IlluminaFileMap fileMap : cycleFileMap.values()) {
            if (!fileMap.isEmpty()) {
                filesAvailable = true;
                break;
            }
        }
        return filesAvailable;
    }

    @Override
    public List<String> verify(final List<Integer> expectedTiles, final int[] expectedCycles) {
        final List<String> failures = new LinkedList<>();
        final Map<Integer, Long> tileToFileLengthMap = new HashMap<>();

        if (!base.exists()) {
            failures.add("Base directory(" + base.getAbsolutePath() + ") does not exist!");
        } else {
            final CycleIlluminaFileMap cfm = getFiles(expectedTiles, expectedCycles);
            for (final int currentCycle : expectedCycles) {
                final IlluminaFileMap fileMap = cfm.get(currentCycle);
                if (fileMap != null) {
                    for (final int tile : expectedTiles) {
                        final File cycleFile = fileMap.get(tile);
                        if (cycleFile != null) {
                            if (tileToFileLengthMap.get(tile) == null) {
                                tileToFileLengthMap.put(tile, cycleFile.length());
                            } else if (!extension.equals(".bcl.gz") && tileToFileLengthMap.get(tile) != cycleFile.length()) {

                                // TODO: The gzip bcl files might not be the same length despite having the same content,
                                // for now we're punting on this but this should be looked into at some point
                                failures.add("File type " + extension
                                        + " has cycles files of different length.  Current cycle ("
                                        + currentCycle + ") " +
                                        "Length of first non-empty file (" + tileToFileLengthMap.get(tile)
                                        + ") length of current cycle (" + cycleFile.length() + ")"
                                        + " File(" + cycleFile.getAbsolutePath() + ")");
                            }
                        } else {
                            failures.add("File type " + extension + " is missing a file for cycle " + currentCycle + " and tile " + tile);
                        }
                    }
                } else {
                    failures.add("Missing file for cycle " + currentCycle + " in directory " + base.getAbsolutePath()
                            + " for file type " + extension);
                }
            }

        }


        return failures;
    }

    @Override
    public List<String> fakeFiles(final List<Integer> expectedTiles, final int[] expectedCycles,
                                  final IlluminaFileUtil.SupportedIlluminaFormat format) {
        final List<String> failures = new LinkedList<>();

        if (!base.exists()) {
            base.mkdirs();
        }

        final Set<Integer> missingCycleSet = new TreeSet<>();
        for (final Integer cycle : expectedCycles) {
            missingCycleSet.add(cycle);
        }

        missingCycleSet.removeAll(detectedCycles);

        for (final Integer cycle : missingCycleSet) {
            final File cycleDirectory = new File(base, "C" + cycle + ".1");
            if (cycleDirectory.mkdirs()) {
                detectedCycles.add(cycle);
            }
        }

        final CycleIlluminaFileMap cfm = getPerTilePerCycleFiles();
        final Map<Integer, Integer> tileToSizeMap = new HashMap<>();
        for (final int currentCycle : expectedCycles) {
            final IlluminaFileMap fileMap = cfm.get(currentCycle);

            if (fileMap == null) {
                for (final Integer tile : expectedTiles) {
                    final File fileToFake = new File(base + File.separator + getFileForCycle(currentCycle, tile));
                    try {
                        if (tileToSizeMap.containsKey(tile)) {
                            faker.fakeFile(fileToFake, tileToSizeMap.get(tile));
                        }
                        else{
                            faker.fakeFile(fileToFake, 1);
                        }
                    } catch (final IOException e) {
                        failures.add("Could not create fake file: " + e.getMessage());
                    }
                }
            } else {
                for (final int tile : expectedTiles) {
                    final File cycleFile = fileMap.get(tile);
                    if (cycleFile != null && !tileToSizeMap.containsKey(tile)) {
                        tileToSizeMap.put(tile, (int) BclReader.getNumberOfClusters(cycleFile));
                    }
                    try {
                        if (cycleFile == null) {
                            final File fileToFake = new File(base + File.separator + getFileForCycle(currentCycle, tile));
                            if (tileToSizeMap.containsKey(tile)) {
                                faker.fakeFile(fileToFake, tileToSizeMap.get(tile));
                            } else {
                                faker.fakeFile(fileToFake, 1);
                            }
                        }
                    } catch (final IOException e) {
                        failures.add("Could not create fake file: " + e.getMessage());
                    }
                }
            }

        }

        for (final Integer cycle : missingCycleSet) {
            failures.add("Missing cycle directory " + cycle + " in directory " + base.getAbsolutePath()
                    + " for file type " + extension);
        }
        return failures;
    }

    private String getFileForCycle(final int currentCycle, final int tile) {
        return "C" + currentCycle + ".1" + File.separator + "s_" + lane + "_" + tile + extension;
    }

    private static int getCycleFromDir(final File tempCycleDir) {
        final String fileName = tempCycleDir.getName();

        final Matcher matcher = IlluminaFileUtil.CYCLE_SUBDIRECTORY_PATTERN.matcher(fileName);
        if (!matcher.matches()) {
            throw new IlluminaParserException("Invalid cycle directory name " + tempCycleDir.getName());
        }

        return Integer.parseInt(matcher.group(1));
    }
}
