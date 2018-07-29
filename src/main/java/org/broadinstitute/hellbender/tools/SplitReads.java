package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.readersplitters.LibraryNameSplitter;
import org.broadinstitute.hellbender.tools.readersplitters.ReadGroupIdSplitter;
import org.broadinstitute.hellbender.tools.readersplitters.ReaderSplitter;
import org.broadinstitute.hellbender.tools.readersplitters.SampleNameSplitter;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.stream.Collectors;

/**
 * Outputs reads from a SAM/BAM/CRAM by read group, sample and library name
 *
 * <p>
 *     Note: Not to be confused with a tool that splits reads for RNA-Seq analysis
 * </p>
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A single BAM file</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 *     <li>A collection of BAM files each corresponding to a group, sample and/or library of the original BAM file</li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <h4>Split reads in BAM file by sample name, read group and library name</h4>
 * <pre>
 *   gatk SplitReads \
 *     -I input.bam \
 *     -O outputDirectory \
 *     --split-sample \
 *     --split-read-group \
 *     --split-library-name
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Outputs reads from a SAM/BAM/CRAM by read group, sample and library name",
        oneLineSummary = "Outputs reads from a SAM/BAM/CRAM by read group, sample and library name",
        programGroup = ReadDataManipulationProgramGroup.class
)
public final class SplitReads extends ReadWalker {

    public static final String SAMPLE_SHORT_NAME = "SM";
    public static final String READ_GROUP_SHORT_NAME = "RG";
    public static final String LIBRARY_NAME_SHORT_NAME = "LB";
    public static final String SAMPLE_LONG_NAME = "split-sample";
    public static final String READ_GROUP_LONG_NAME = "split-read-group";
    public static final String LIBRARY_NAME_LONG_NAME = "split-library-name";
    public static final String UNKNOWN_OUT_PREFIX = "unknown";


    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "The directory to output SAM/BAM/CRAM files."
    )
    public File OUTPUT_DIRECTORY = new File("");

    @Argument(
            fullName = SAMPLE_LONG_NAME,
            shortName = SAMPLE_SHORT_NAME,
            doc = "Split file by sample."
    )
    public boolean SAMPLE;

    @Argument(
            fullName = READ_GROUP_LONG_NAME,
            shortName = READ_GROUP_SHORT_NAME,
            doc = "Split file by read group."
    )
    public boolean READ_GROUP;

    @Argument(
            fullName = LIBRARY_NAME_LONG_NAME,
            shortName = LIBRARY_NAME_SHORT_NAME,
            doc = "Split file by library."
    )
    public boolean LIBRARY_NAME;

    private final List<ReaderSplitter<?>> splitters = new ArrayList<>();
    private Map<String, SAMFileGATKReadWriter> outs = null;

    @Override
    public void onTraversalStart() {
        IOUtil.assertDirectoryIsWritable(OUTPUT_DIRECTORY);
        if ( readArguments.getReadFiles().size() != 1 ) {
            throw new UserException("This tool only accepts a single SAM/BAM/CRAM as input");
        }

        if (SAMPLE) {
            splitters.add(new SampleNameSplitter());
        }
        if (READ_GROUP) {
            splitters.add(new ReadGroupIdSplitter());
        }
        if (LIBRARY_NAME) {
            splitters.add(new LibraryNameSplitter());
        }
        outs = createWriters(splitters);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        outs.computeIfAbsent(getKey(splitters, read), this::createUnknownOutOnDemand).addRead(read);
    }

    @Override
    public void closeTool() {
        if ( outs != null ) {
            outs.values().forEach(writer -> writer.close());
        }
    }

    // Create an output stream on demand for holding any reads that do not have a value for one or more of the
    // attributes we're grouping by
    private SAMFileGATKReadWriter createUnknownOutOnDemand(String attributeValue) {
        if (!attributeValue.equals("."+UNKNOWN_OUT_PREFIX)) {
            // the only attribute value we should ever discover at runtime is the string ".unknown" which is
            // synthesized by "getkey" below when a splitter returns null because we're splitting on some
            // attribute for which a given read/group has no value; anything else indicates a coding error
            throw new GATKException.ShouldNeverReachHereException("Unrecognized attribute value found: " + attributeValue);
        }
        final SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();
        final SAMFileHeader samFileHeaderIn = getHeaderForReads();

        return prepareSAMFileWriter(samFileWriterFactory, samFileHeaderIn, attributeValue);
    }

    //  Create a new output file and prepare and return the corresponding SAMFileGATKReadWriter.
    private SAMFileGATKReadWriter prepareSAMFileWriter(
            SAMFileWriterFactory samFileWriterFactory,
            SAMFileHeader samFileHeaderIn,
            final String keyName) {
        final String base = FilenameUtils.getBaseName(readArguments.getReadFiles().get(0).getName());
        final String extension = "." + FilenameUtils.getExtension(readArguments.getReadFiles().get(0).getName());
        final File outFile = new File(OUTPUT_DIRECTORY, base + keyName + extension);
        return createSAMWriter(outFile, true);
    }

    /**
     * Creates SAMFileWriter instances for the reader splitters based on the input file.
     * @param splitters Reader splitters.
     * @return A map of file name keys to SAMFileWriter.
     */
    private Map<String, SAMFileGATKReadWriter> createWriters(final List<ReaderSplitter<?>> splitters) {
        final Map<String, SAMFileGATKReadWriter> outs = new LinkedHashMap<>();

        final SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();
        final SAMFileHeader samFileHeaderIn = getHeaderForReads();

        // Build up a list of key options at each level.
        final List<List<?>> splitKeys = splitters.stream()
                .map(splitter -> splitter.getSplitsBy(samFileHeaderIn))
                .collect(Collectors.toList());

        // For every combination of keys, add a SAMFileWriter.
        addKey(splitKeys, 0, "", key -> {
            outs.put(key, prepareSAMFileWriter(samFileWriterFactory, samFileHeaderIn, key));
        });

        return outs;
    }

    /**
     * Recursively builds up a key, then when it reaches the bottom of the list, calls the adder on the generated key.
     * @param listKeys A outer list, where each inner list contains the output options for that level.
     * @param listIndex The current recursive index within the listKeys.
     * @param key The built up key recursively
     * @param adder Function to run on the recursively generated key once the bottom of the outer list is reached.
     */
    private void addKey(final List<List<?>> listKeys, final int listIndex,
                        final String key, final Consumer<String> adder) {
        if (listIndex < listKeys.size()) {
            for (final Object newKey : listKeys.get(listIndex)) {
                addKey(listKeys, listIndex + 1, key + "." + newKey, adder);
            }
        } else {
            adder.accept(key);
        }
    }

    /**
     * Traverses the splitters generating a key for this particular record.
     * @param splitters The list of splitters.
     * @param record The record to analyze.
     * @return The generated key that may then be used to find the appropriate SAMFileWriter.
     */
    private String getKey(final List<ReaderSplitter<?>> splitters, final GATKRead record) {
        // if a read is missing the value for the target split, return the constant "unknown" which will
        // result in a new output stream being created on demand to hold uncategorized reads

        return splitters.stream()
                .map(s -> {
                    final Object key = s.getSplitBy(record, getHeaderForReads());
                    return key == null ? UNKNOWN_OUT_PREFIX : key.toString();
                })
                .reduce("", (acc, item) -> acc + "." + item);
    }
}
