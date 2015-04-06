package org.broadinstitute.hellbender.tools.picard.illumina;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProcessExecutor;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.IlluminaProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataProviderFactory;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaFileUtil;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.OutputMapping;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ParameterizedFileUtil;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadStructure;

import java.io.File;
import java.util.*;

/**
 * Program to check a lane of an Illumina output directory.  This program checks that files exist, are non-zero in length, for every tile/cycle and
 * specified data type.  If NO data type is specified then the default data types used by IlluminaBasecallsToSam are used.
 */
@CommandLineProgramProperties(
        usage = "Check that the files to provide the data specified by DATA_TYPES are available, exist, and are reasonably sized for every tile/cycle.  " +
                "Reasonably sized means non-zero sized for files that exist per tile and equal size for binary files that exist per cycle/per tile. " +
                "CheckIlluminaDirectory DOES NOT check that the individual records in a file are well-formed.",
        usageShort = "Asserts the validity of the data in the specified Illumina basecalling data",
        programGroup = IlluminaProgramGroup.class
)
public class CheckIlluminaDirectory extends PicardCommandLineProgram {
    private static final Log log = Log.getInstance(CheckIlluminaDirectory.class);

    @Argument(doc = "The basecalls output directory. ", shortName = "B")
    public File BASECALLS_DIR;

    @Argument(doc = "The data types that should be checked for each tile/cycle.  If no values are provided then the data types checked are those " +
            "required by IlluminaBaseCallsToSam (which is a superset of those used in ExtractIlluminaBarcodes).  These data types vary slightly depending on " +
            "whether or not the run is barcoded so READ_STRUCTURE should be the same as that which will be passed to IlluminaBasecallsToSam.  If this option " +
            "is left unspecified then both ExtractIlluminaBarcodes and IlluminaBaseCallsToSam should complete successfully UNLESS the " +
            "individual records of the files themselves are spurious.",
            shortName = "DT", optional = true)
    public final Set<IlluminaDataType> DATA_TYPES = new TreeSet<IlluminaDataType>();

    @Argument(doc = ReadStructure.PARAMETER_DOC + " Note:  If you want to check whether or not a future IlluminaBasecallsToSam or ExtractIlluminaBarcodes " +
            "run will fail then be sure to use the exact same READ_STRUCTURE that you would pass to these programs for this run.",
            shortName = "RS")
    public String READ_STRUCTURE;

    @Argument(doc = "The number of the lane(s) to check. ", shortName = StandardArgumentDefinitions.LANE_SHORT_NAME)
    public List<Integer> LANES;

    @Argument(doc = "The number(s) of the tile(s) to check. ", shortName = "T", optional = true)
    public List<Integer> TILE_NUMBERS;

    @Argument(doc = "A flag to determine whether or not to create fake versions of the missing files.", shortName = "F",
            optional = true)
    public Boolean FAKE_FILES = false;

    @Argument(doc = "A flag to create symlinks to the loc file for the X Ten for each tile.", shortName = "X",
            optional = true)
    public Boolean LINK_LOCS = false;

    /**
     * @return a status code indicating the number of errors detected
     */
    @Override
    protected Object doWork() {
        final ReadStructure readStructure = new ReadStructure(READ_STRUCTURE);
        if (DATA_TYPES.isEmpty()) {
            DATA_TYPES.addAll(Arrays.asList(IlluminaBasecallsConverter.DATA_TYPES_NO_BARCODE));
        }

        final List<Integer> failingLanes = new ArrayList<Integer>();
        int totalFailures = 0;

        final int[] expectedCycles = new OutputMapping(readStructure).getOutputCycles();
        log.info("Checking lanes(" + StringUtil.join(",", LANES) + " in basecalls directory (" + BASECALLS_DIR
                .getAbsolutePath() + ")\n");
        log.info("Expected cycles: " + StringUtil.intValuesToString(expectedCycles));

        for (final Integer lane : LANES) {
            IlluminaFileUtil fileUtil = new IlluminaFileUtil(BASECALLS_DIR, lane);
            final List<Integer> expectedTiles = fileUtil.getExpectedTiles();
            if (!TILE_NUMBERS.isEmpty()) {
                expectedTiles.retainAll(TILE_NUMBERS);
            }

            if (LINK_LOCS) {
                createLocFileSymlinks(fileUtil, lane);
                //we need to create a new file util because it stores a cache to the files it found on
                //construction and this doesn't inclue the recently created symlinks
                fileUtil = new IlluminaFileUtil(BASECALLS_DIR, lane);
            }

            log.info("Checking lane " + lane);
            log.info("Expected tiles: " + StringUtil.join(", ", expectedTiles));

            final int numFailures = verifyLane(fileUtil, expectedTiles, expectedCycles, DATA_TYPES, FAKE_FILES);

            if (numFailures > 0) {
                log.info("Lane " + lane + " FAILED " + " Total Errors: " + numFailures);
                failingLanes.add(lane);
                totalFailures += numFailures;
            } else {
                log.info("Lane " + lane + " SUCCEEDED ");
            }
        }

        int status = 0;
        if (totalFailures == 0) {
            log.info("SUCCEEDED!  All required files are present and non-empty.");
        } else {
            status = totalFailures;
            log.info("FAILED! There were " + totalFailures + " in the following lanes: " + StringUtil
                    .join(", ", failingLanes));
        }

        return status;
    }

    private void createLocFileSymlinks(final IlluminaFileUtil fileUtil, final int lane) {
        final File baseFile = new File(BASECALLS_DIR.getParentFile().getAbsolutePath() + File.separator + "s.locs");
        final File newFileBase = new File(baseFile.getParent() + File.separator + IlluminaFileUtil
                .longLaneStr(lane) + File.separator);
        if (baseFile.exists()) {
            boolean success = true;
            if (!newFileBase.exists()) {
                success = newFileBase.mkdirs();
            }
            if (success) {
                for (final Integer tile : fileUtil.getExpectedTiles()) {
                    final String newName =
                            newFileBase + File.separator + String.format("s_%d_%d.locs", lane, tile);
                    final ProcessExecutor.ExitStatusAndOutput output =
                            ProcessExecutor.executeAndReturnInterleavedOutput(new String[]{"ln", "-fs", baseFile.getAbsolutePath(), newName});
                    if (output.exitStatus != 0) {
                        throw new RuntimeIOException("Could not create symlink: " + output.stdout);
                    }
                }
            } else {
                throw new RuntimeIOException(String.format("Could not create lane directory: %s.", newFileBase.getAbsolutePath()));
            }
        } else {
            throw new RuntimeIOException(String.format("Locations file %s does not exist.", baseFile.getAbsolutePath()));
        }

    }

    /**
     * Use fileUtil to find the data types that would be used by IlluminaDataProvider.  Verify that for the expected
     * tiles/cycles/data types that all the files needed to provide their data is present.  This method logs every
     * error that is found (excluding file faking errors) and returns the number of errors found
     *
     * @param fileUtil      A file util paramterized with the directory/lane to check
     * @param expectedTiles The tiles we expect to be available/well-formed
     * @param cycles        The cycles we expect to be available/well-formed
     * @param dataTypes     The data types we expect to be available/well-formed
     * @return The number of errors found/logged for this directory/lane
     */
    private static final int verifyLane(final IlluminaFileUtil fileUtil, final List<Integer> expectedTiles,
                                        final int[] cycles,
                                        final Set<IlluminaDataType> dataTypes, final boolean fakeFiles) {
        if (expectedTiles.isEmpty()) {
            throw new UserException("0 input tiles were specified!  Check to make sure this lane is in the InterOp file!");
        }

        if (cycles.length == 0) {
            throw new UserException("0 output cycles were specified!");
        }

        int numFailures = 0;

        //find what request IlluminaDataTypes we have files for and select the most preferred file format available for that type
        final Map<IlluminaFileUtil.SupportedIlluminaFormat, Set<IlluminaDataType>> formatToDataTypes =
                IlluminaDataProviderFactory.determineFormats(dataTypes, fileUtil);

        //find if we have any IlluminaDataType with NO available file formats and, if any exist, increase the error count
        final Set<IlluminaDataType> unmatchedDataTypes =
                IlluminaDataProviderFactory.findUnmatchedTypes(dataTypes, formatToDataTypes);
        if (!unmatchedDataTypes.isEmpty()) {
            if (fakeFiles) {
                for (final IlluminaDataType dataType : unmatchedDataTypes) {
                    final IlluminaFileUtil.SupportedIlluminaFormat format =
                            IlluminaDataProviderFactory.findPreferredFormat(dataType, fileUtil);
                    fileUtil.getUtil(format).fakeFiles(expectedTiles, cycles, format);

                }
            }
            log.info("Could not find a format with available files for the following data types: " + StringUtil
                    .join(", ", new ArrayList<IlluminaDataType>(unmatchedDataTypes)));
            numFailures += unmatchedDataTypes.size();
        }

        for (final IlluminaFileUtil.SupportedIlluminaFormat format : formatToDataTypes.keySet()) {
            final ParameterizedFileUtil util = fileUtil.getUtil(format);
            final List<String> failures = util.verify(expectedTiles, cycles);
            //if we have failures and we want to fake files then fake them now.
            if (!failures.isEmpty() && fakeFiles) {
                //fake files
                util.fakeFiles(expectedTiles, cycles, format);

            }
            numFailures += failures.size();
            for (final String failure : failures) {
                log.info(failure);
            }
        }

        return numFailures;
    }

    @Override
    protected String[] customCommandLineValidation() {
        IOUtil.assertDirectoryIsReadable(BASECALLS_DIR);
        final List<String> errors = new ArrayList<String>();

        for (final Integer lane : LANES) {
            if (lane < 1) {
                errors.add(
                        "LANES must be greater than or equal to 1.  LANES passed in " + StringUtil.join(", ", LANES));
                break;
            }
        }

        if (errors.isEmpty()) {
            return null;
        } else {
            return errors.toArray(new String[errors.size()]);
        }
    }
}
