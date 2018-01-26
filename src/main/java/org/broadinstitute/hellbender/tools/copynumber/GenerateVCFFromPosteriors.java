package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.ChunkedCopyNumberPosteriorCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LocatableCopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVPostProcessingEngine;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVNamingConstants;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberStateCollection;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Create a VCF given the output of {@link GermlineCNVCaller}.
 *
 * <p> This tool takes a list of {@link GermlineCNVCaller} call directories containing posterior records for
 * different copy number states and outputs a VCF for a specified sample. Note that the list should be sorted
 * according order of the intervals in corresponding chunks. The VCF will contain an ALT allele for every
 * non-reference copy number state specified in the posterior files. The reference allele corresponds to copy number 2.
 * </p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>A list of paths to {@link GermlineCNVCaller} chunked calls directories</li>
 *     <li>Name of the sample directory (it must be the same across all chunks)</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A VCF file with CNV calls</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *   gatk GenerateVCFFromPosteriors \
 *     --called-chunk-directory path/to/chunk_1
 *     --called-chunk-directory path/to/chunk_2
 *     --sample-name-directory SAMPLE_1
 *     --output output.vcf
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Create a VCF given the output of GermlineCNVCaller.",
        oneLineSummary = "Create a VCF given the output of GermlineCNVCaller.",
        programGroup = CopyNumberProgramGroup.class
)
public final class GenerateVCFFromPosteriors extends GATKTool {
    private static final Logger logger = LogManager.getLogger(GenerateVCFFromPosteriors.class);

    public static final String POSTERIOR_CALL_DIRECTORY_FULL_NAME = "chunk-path";
    public static final String SAMPLE_DIRECTORY_NAME_FULL_NAME = "sample-directory";

    private static final String COMMENT_PREFIX = "@";

    @Argument(
            doc = "List of paths to GermlineCNVCaller call directories, sorted in order of the intervals contained.",
            fullName = POSTERIOR_CALL_DIRECTORY_FULL_NAME,
            minElements = 1
    )
    private List<File> orderedChunkDirectoryList = null;

    @Argument(
            doc = "Name of the sample directory (must be contained in all call directories).",
            fullName = SAMPLE_DIRECTORY_NAME_FULL_NAME
    )
    private String sampleDirectoryName = null;

    @Argument(
            doc = "Output VCF file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputFile = null;

    @Override
    public void traverse() {
        final SimpleIntervalCollection firstChunkSimpleIntervalCollection =
                new SimpleIntervalCollection(new File(getIntervalFileFromChunkDirectory(
                        orderedChunkDirectoryList.get(0).getAbsolutePath()).getAbsolutePath()));
        final SAMSequenceDictionary samSequenceDictionary = firstChunkSimpleIntervalCollection
                .getMetadata().getSequenceDictionary();
        Triple<List<LocatableCopyNumberPosteriorDistribution>, IntegerCopyNumberStateCollection, String> firstCopyNumberPosteriorFileInfo =
                readChunkedPosteriorFileFromDirectory(orderedChunkDirectoryList.get(0), samSequenceDictionary);
        final IntegerCopyNumberStateCollection copyNumberStateCollection = firstCopyNumberPosteriorFileInfo.getMiddle();
        final String sampleName = firstCopyNumberPosteriorFileInfo.getRight();

        if (firstCopyNumberPosteriorFileInfo.getMiddle().size() < 3) {
            throw new UserException.BadInput("There should be at least 3 copy number states.");
        }

        final VariantContextWriter outputWriter = createVCFWriter(outputFile);
        final GermlineCNVPostProcessingEngine germlineCNVPostProcessingEngine = new GermlineCNVPostProcessingEngine(
                outputWriter, copyNumberStateCollection, sampleName, samSequenceDictionary);

        //TODO pass the command line invocation string to the header composer method
        germlineCNVPostProcessingEngine.composeVariantContextHeader(getDefaultToolVCFHeaderLines());
        int currentChunk = 0;

        for(File chunkRootDirectory: orderedChunkDirectoryList) {
            currentChunk++;
            logger.info(String.format("Analyzing copy number posterior chunk %d", currentChunk));

            final Triple<List<LocatableCopyNumberPosteriorDistribution>, IntegerCopyNumberStateCollection, String> chunkedCopyNumberPosteriorFileInfo =
            readChunkedPosteriorFileFromDirectory(chunkRootDirectory, samSequenceDictionary);
            //check that the copy number state collection is equivalent across all chunks
            if (!copyNumberStateCollection.equals(chunkedCopyNumberPosteriorFileInfo.getMiddle())) {
                throw new UserException.BadInput("Copy number collection differs across chunked posterior outputs.");
            }
            //check that the sample name is consistent across all chunks
            if (!chunkedCopyNumberPosteriorFileInfo.getRight().equals(sampleName)) {
                throw new UserException.BadInput("Sample name differs across chunked posterior outputs.");
            }

            final List<LocatableCopyNumberPosteriorDistribution> locatableRecordsList =
                    chunkedCopyNumberPosteriorFileInfo.getLeft();

            germlineCNVPostProcessingEngine.writeChunkedVariantContext(locatableRecordsList);
        }
        outputWriter.close();
    }

    /**
     * Reads and processes a {@link GermlineCNVCaller} analysis chunk for a single sample.
     *
     * @param chunkRootDirectory path to chunk directory where posteriors are stores
     * @param samSequenceDictionary SAM sequence dictionary
     * @return a {@link Triple} object containing list of posterior locatable records, copy number state collection
     * and sample name
     */
    private Triple<List<LocatableCopyNumberPosteriorDistribution>, IntegerCopyNumberStateCollection, String> readChunkedPosteriorFileFromDirectory(
            final File chunkRootDirectory,
            final SAMSequenceDictionary samSequenceDictionary) {
        //get a list of intervals associated with chunk currently being processed
        final SimpleIntervalCollection chunkSimpleIntervalCollection =
                new SimpleIntervalCollection(new File(getIntervalFileFromChunkDirectory(
                        chunkRootDirectory.getAbsolutePath()).getAbsolutePath()));

        //check that sequence dictionaries are equivalent among all chunks
        if (!chunkSimpleIntervalCollection.getMetadata().getSequenceDictionary().equals(samSequenceDictionary)) {
            throw new UserException.BadInput("Not all sequence dictionaries are equivalent among interval" +
                    " lists in chunk directories.");
        }
        final File chunkPosteriorFile = getPosteriorFileFromChunkDirectory(chunkRootDirectory.getAbsolutePath(),
                sampleDirectoryName);

        final List<String> copyNumberStateColumns = extractPosteriorFileColumns(chunkPosteriorFile);
        final IntegerCopyNumberStateCollection chunkCopyNumberStateCollection =
                new IntegerCopyNumberStateCollection(copyNumberStateColumns);

        final ChunkedCopyNumberPosteriorCollection chunkedCopyNumberPosteriorCollection =
                new ChunkedCopyNumberPosteriorCollection(chunkPosteriorFile, chunkCopyNumberStateCollection);

        //Combine together chunked interval list and chunked posterior collection
        // into list of LocatableCopyNumberPosteriorDistribution
        final List<LocatableCopyNumberPosteriorDistribution> locatableCopyNumberPosteriorDistributionList =
                IntStream.range(0, chunkedCopyNumberPosteriorCollection.size())
                        .mapToObj(i -> new LocatableCopyNumberPosteriorDistribution(
                                chunkSimpleIntervalCollection.getIntervals().get(i),
                                chunkedCopyNumberPosteriorCollection.getRecords().get(i)))
                        .collect(Collectors.toList());

        return new ImmutableTriple<>(locatableCopyNumberPosteriorDistributionList, chunkCopyNumberStateCollection,
                chunkedCopyNumberPosteriorCollection.getMetadata().getSampleName());
    }

    /**
     * Extracts column names from a copy-number log posterior TSV file
     */
    private static List<String> extractPosteriorFileColumns(final File copyNumberPosteriorFile) {
        List<String> columns = null;
        try (final XReadLines reader = new XReadLines(copyNumberPosteriorFile)) {
            while (reader.hasNext()) {
                String nextLine = reader.next();
                if (!nextLine.startsWith(COMMENT_PREFIX)) {
                    columns = Arrays.asList(nextLine.split(TableUtils.COLUMN_SEPARATOR_STRING));
                    break;
                }
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(copyNumberPosteriorFile);
        }
        if (columns == null) {
            throw new UserException.BadInput("Copy number posterior file does not have a header.");
        }
        return columns;
    }

    /**
     * Get the posterior file from chunk directory
     */
    private static File getPosteriorFileFromChunkDirectory(final String chunkDirectoryPath, final String sampleIndex) {
        final String posteriorCallsDirectory = chunkDirectoryPath + File.separator + sampleIndex;
        return new File(posteriorCallsDirectory, GermlineCNVNamingConstants.COPY_NUMBER_POSTERIOR_FILE_NAME);
    }

    /**
     * Get the intervals file from chunk directory
     */
    private static File getIntervalFileFromChunkDirectory(final String chunkDirectoryPath) {
        return new File(chunkDirectoryPath, GermlineCNVNamingConstants.INTERVAL_LIST_FILE_NAME);
    }
}
