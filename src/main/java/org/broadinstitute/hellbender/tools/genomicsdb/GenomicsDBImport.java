package org.broadinstitute.hellbender.tools.genomicsdb;

import com.google.common.annotations.VisibleForTesting;
import com.intel.genomicsdb.ChromosomeInterval;
import com.intel.genomicsdb.GenomicsDBCallsetsMapProto;
import com.intel.genomicsdb.GenomicsDBImportConfiguration;
import com.intel.genomicsdb.GenomicsDBImporter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.nio.SeekableByteChannelPrefetcher;

import java.io.File;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;


/**
 * This tool imports GVCFs to GenomicsDB. To run this tool,
 * 1. A single interval must be provided
 * 2. The tool accepts multiple GVCFs each of which must contain data
 *    for one sample
 * 3. The path to the GenomicsDB workspace must be specified
 * 4. User may optionally specify paths to which to write JSON files
 *
 * To read data from GenomicsDB, use the query interface GenomicsDBFeatureReader
 */
@CommandLineProgramProperties(
    summary = "Import VCFs to GenomicsDB",
    oneLineSummary = "Import VCFs to GenomicsDB",
    programGroup = VariantProgramGroup.class
)
public final class GenomicsDBImport extends GATKTool {

    private static final long DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE = 16*1024L;
    private static final long DEFAULT_SEGMENT_SIZE = 1048576L;
    private static final int DEFAULT_ZERO_BATCH_SIZE = 0;

    public static final String WORKSPACE_ARG_NAME = "genomicsDBWorkspace";
    public static final String SEGMENT_SIZE_ARG_NAME = "genomicsDBSegmentSize";
    public static final String OVERWRITE_WORKSPACE_NAME = "overwriteExistingGenomicsDBWorkspace";

    public static final String VCF_BUFFER_SIZE_ARG_NAME = "genomicsDBVCFBufferSize";

    public static final String BATCHSIZE_ARG_NAME = "batchSize";
    public static final String CONSOLIDATE_ARG_NAME = "consolidate";
    public static final String SAMPLE_NAME_MAP_LONG_NAME = "sampleNameMap";

    @Argument(fullName = WORKSPACE_ARG_NAME,
              shortName = WORKSPACE_ARG_NAME,
              doc = "Workspace for GenomicsDB. Has to be a POSIX file system path")
    private String workspace;

    @Argument(fullName = SEGMENT_SIZE_ARG_NAME,
              shortName = SEGMENT_SIZE_ARG_NAME,
              doc = "Buffer size in bytes allocated for GenomicsDB attributes during " +
                    "import. Should be large enough to hold data from one site. " +
                    " Defaults to " + DEFAULT_SEGMENT_SIZE,
              optional = true)
    private long segmentSize = DEFAULT_SEGMENT_SIZE;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
              shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
              doc = "GVCF files to be imported to GenomicsDB. Each file must contain" +
                    "data for only a single sample. Either this or " + SAMPLE_NAME_MAP_LONG_NAME +
                    " must be specified.",
              optional = true,
              mutex = {SAMPLE_NAME_MAP_LONG_NAME})
    private List<String> variantPaths;

    @Argument(fullName = VCF_BUFFER_SIZE_ARG_NAME,
              shortName = VCF_BUFFER_SIZE_ARG_NAME,
              doc = "Buffer size in bytes to store variant contexts." +
                    " Larger values are better as smaller values cause frequent disk writes." +
                    " Defaults to " + DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE,
              optional = true,
              minValue = 1024L)
    private long vcfBufferSizePerSample = DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE;

    @Argument(fullName = OVERWRITE_WORKSPACE_NAME,
              shortName = OVERWRITE_WORKSPACE_NAME,
              doc = "Will overwrite given workspace if it exists. " +
                    "Otherwise a new workspace is created. " +
                    "Defaults to false",
              optional = true)
    private Boolean overwriteExistingWorkspace = false;

    @Argument(fullName = BATCHSIZE_ARG_NAME,
              shortName = BATCHSIZE_ARG_NAME,
              doc = "Batch size controls the number of samples for which readers are open at once " +
                    "and therefore provides a way to minimize memory consumption. However, it can take longer to complete. " +
                    "Use the consolidate flag if more than a hundred batches were used. This will improve feature read time. " +
                    "batchSize=0 means no batching (i.e. readers for all samples will be opened at once) " +
                    "Defaults to " + DEFAULT_ZERO_BATCH_SIZE,
              optional = true)
    private int batchSize = DEFAULT_ZERO_BATCH_SIZE;

    @Argument(fullName = CONSOLIDATE_ARG_NAME,
              shortName = CONSOLIDATE_ARG_NAME,
              doc = "Boolean flag to enable consolidation. If importing data in batches, a new fragment is created for " +
                    "each batch. In case thousands of fragments are created, GenomicsDB feature readers will try " +
                    "to open ~20x as many files. Also, internally GenomicsDB would consume more memory to maintain " +
                    "bookkeeping data from all fragments. Use this flag to merge all fragments into one. " +
                    "Merging can potentially improve read performance, however overall benefit might not be noticeable " +
                    "as the top Java layers have significantly higher overheads. This flag has no effect if only one " +
                    "batch is used. Defaults to false",
              optional = true)
    private Boolean doConsolidation = false;

    @Advanced
    @Argument(fullName = SAMPLE_NAME_MAP_LONG_NAME,
            shortName = SAMPLE_NAME_MAP_LONG_NAME,
            doc = "Path to file containing a mapping of sample name to file name in tab delimited format.  If this is " +
                    "specified then the header from the first sample will be used as the header template rather than " +
                    "merging the headers, and the sample names will taken from this file.  This is a performance optimization " +
                    "that relaxes the normal checks for consistent headers.  Using vcfs with incompatible headers may result " +
                    "in silent data corruption.",
            optional = true,
            mutex = {StandardArgumentDefinitions.VARIANT_LONG_NAME})
    private String sampleNameToPathFile;

    @Override
    public boolean requiresIntervals() { return true; }

    @Override
    public boolean requiresReads() { return false; }

    @Override
    public boolean requiresReference() {
      return false;
    }

    @Override
    public int getDefaultCloudPrefetchBufferSize() {
        // Since this tool is typically run with many inputs, we want a smaller NIO buffer size
        // than normal:
        return 2;
    }

    @Override
    public int getDefaultCloudIndexPrefetchBufferSize() {
        // Since this tool is typically run with many inputs, we want a smaller NIO buffer size for the index
        // than normal:
        return 2;
    }

    // Intervals from command line (singleton for now)
    private List<ChromosomeInterval> intervals;

    // List of all sample names
    private final Set<String> sampleNames = new LinkedHashSet<>();

    // Linked hash map between sample names and corresponding GVCF file name
    private final Map<String, String> sampleNameToVcfUri = new LinkedHashMap<>();

    // Needed as smartMergeHeaders() returns a set of VCF header lines
    private Set<VCFHeaderLine> mergedHeaderLines = null;

    // VCFHeader created from the merged header
    private SAMSequenceDictionary mergedHeaderSequenceDictionary;

    // Path to vidmap file to be written by GenomicsDBImporter
    private File vidMapJSONFile;

    // Path to callsetmap file to be written by GenomicsDBImporter
    private File callsetMapJSONFile;

    // GenomicsDB callset map protobuf structure containing all callset names
    // used to write the callset json file on traversal success
    private GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB;

    /**
     * Before traversal starts, create the feature readers
     * for all the input GVCFs, create the merged header and
     * initialize the interval
     */
    @Override
    public void onStartup() {
        assertVariantPathsOrSampleNameFileWasSpecified();
        initializeHeaderAndSampleMappings();
        initializeIntervals();
        super.onStartup();
    }

    private void assertVariantPathsOrSampleNameFileWasSpecified(){
        if ( (variantPaths == null || variantPaths.isEmpty()) && sampleNameToPathFile == null) {
            throw new CommandLineException.MissingArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME,
                                                       "One of " + StandardArgumentDefinitions.VARIANT_LONG_NAME + " or " + SAMPLE_NAME_MAP_LONG_NAME + " must be specified" );
        }
    }

    private void initializeHeaderAndSampleMappings() {
        // Only one of -V and -sampleNameToPathFile may be specified
        if (sampleNameToPathFile == null) {
            // -V was specified
            final List<VCFHeader> headers = new ArrayList<>(variantPaths.size());

            for (final String variantPath : variantPaths) {
                final  VCFHeader header = loadHeaderFromVCFUri(variantPath);

                // A GVCF file must contain only one sample, throw an exception otherwise
                assertGVCFHasOnlyOneSample(variantPath, header);
                headers.add(header);

                final String sampleName = header.getGenotypeSamples().get(0);
                sampleNames.add(sampleName);
                final String previousPath = sampleNameToVcfUri.put(sampleName, variantPath);
                if (previousPath != null) {
                    throw new UserException("Duplicate sample: " + sampleName + ". Sample was found in both "
                                                    + variantPath + " and " + previousPath + ".");
                }
            }
            mergedHeaderLines = VCFUtils.smartMergeHeaders(headers, true);
            mergedHeaderSequenceDictionary = new VCFHeader(mergedHeaderLines).getSequenceDictionary();
        } else {
            // -sampleNameMap was specified
            sampleNameToVcfUri.putAll(loadSampleNameMapFile(IOUtils.getPath(sampleNameToPathFile)));
            sampleNames.addAll(sampleNameToVcfUri.keySet());

            final String firstHeaderPath = sampleNameToVcfUri.entrySet().iterator().next().getValue();
            final VCFHeader header = loadHeaderFromVCFUri(firstHeaderPath);
            mergedHeaderLines = header.getMetaDataInInputOrder();
            mergedHeaderSequenceDictionary = header.getSequenceDictionary();
        }
    }

    private VCFHeader loadHeaderFromVCFUri(final String variantPath) {
        try(final AbstractFeatureReader<VariantContext, LineIterator> reader = getReaderFromVCFUri(variantPath)) {
            return (VCFHeader) reader.getHeader();
        } catch (final IOException e) {
            throw new UserException("Error while reading vcf header from " + variantPath, e);
        }
    }

    private static void assertGVCFHasOnlyOneSample(String variantPath, VCFHeader header) {
        // A GVCF file must contain only one sample, throw an exception otherwise
        final int numberOfSamples = header.getNGenotypeSamples();
        if (numberOfSamples != 1) {
            throw new UserException("Input GVCF: " + variantPath + " was expected to contain a single sample but actually contained " + numberOfSamples + ".");
        }
    }

    @VisibleForTesting
    static Map<String, String> loadSampleNameMapFile(final Path sampleToFileMapPath) {
        try {
            final Map<String, String> sampleToFilename = Files.lines(sampleToFileMapPath)
                    .map(line -> {
                        final String[] split = line.split("\\s");
                        if (split.length != 2) {
                            throw new UserException.BadInput(
                                    "Expected a file of format\nSample\tFile\n but found line: " + line);
                        }
                        return split;
                    })
                    .collect(Collectors.toMap(split -> split[0],
                                             split -> split[1],
                                              (a, b) -> { throw new UserException.BadInput("Found two mappings for the same sample:\n" + a + "\n" + b); },
                                              HashMap::new));

            if (sampleToFilename.entrySet().isEmpty()) {
                throw new UserException.BadInput(
                        "At least 1 sample is required but none were found in the sample mapping file");
            }
            return sampleToFilename;
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(sampleToFileMapPath, "exception while reading sample->filename mapping file",  e);
        }
    }


    /**
     * Before traversal, fix configuration parameters and initialize
     * GenomicsDB. Hard-coded to handle only VCF files and headers
     */
    @Override
    public void onTraversalStart() {

        File workspaceDir = overwriteOrCreateWorkspace();

        vidMapJSONFile = new File(workspaceDir + "/" + GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME);
        callsetMapJSONFile = new File (workspaceDir + "/" + GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME);

        logger.info("Vid Map JSON file will be written to " + vidMapJSONFile);
        logger.info("Callset Map JSON file will be written to " + callsetMapJSONFile);
        logger.info("Importing to array - " + workspace + "/" + GenomicsDBConstants.DEFAULT_ARRAY_NAME);

        // Passing in false here so that sample names will be sorted.
        // This is needed for consistent ordering across partitions/machines
        callsetMappingPB = GenomicsDBImporter.generateSortedCallSetMap(new ArrayList<>(sampleNames), false);
    }

    /**
     * A complete traversal from start to finish. This method will import all samples
     * specified in the input GVCF files.
     */
    @Override
    public void traverse() {

        final int updatedBatchSize = (batchSize == DEFAULT_ZERO_BATCH_SIZE) ? sampleNames.size() : batchSize;
        final int sampleCount = sampleNames.size();
        final int totalBatchCount = (sampleCount/updatedBatchSize) + (sampleCount%updatedBatchSize==0 ? 0 : 1);

        GenomicsDBImporter importer;

        for (int i = 0, batchCount = 1; i < sampleCount; i += updatedBatchSize, ++batchCount) {

            final Map<String, FeatureReader<VariantContext>> sampleToReaderMap =
                getFeatureReaders(new ArrayList<>(sampleNames), sampleNameToVcfUri, updatedBatchSize, sampleCount, i);

            logger.info("Importing batch " + batchCount + " with " + sampleToReaderMap.size() + " samples");
            final long variantContextBufferSize = vcfBufferSizePerSample * sampleToReaderMap.size();

            final GenomicsDBImportConfiguration.ImportConfiguration importConfiguration =
                    createImportConfiguration(workspace, GenomicsDBConstants.DEFAULT_ARRAY_NAME,
                            variantContextBufferSize, segmentSize,
                            i, (i+updatedBatchSize-1));

            try {
                importer = new GenomicsDBImporter(sampleToReaderMap, mergedHeaderLines, intervals.get(0), importConfiguration);
            } catch (final IOException e) {
                throw new UserException("Error initializing GenomicsDBImporter in batch " + batchCount, e);
            }

            try {
                importer.importBatch();
            } catch (final IOException e) {
                throw new UserException("GenomicsDB import failed in batch " + batchCount, e);
            }

            closeReaders(sampleToReaderMap);
            logger.info("Done importing batch " + batchCount + "/" + totalBatchCount);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        if (batchSize==DEFAULT_ZERO_BATCH_SIZE) {
            logger.info("Import completed!");
        } else {
            logger.info("Import of all batches to GenomicsDB completed!");
        }

        // Write the vid and callset map JSON files
        try {
            GenomicsDBImporter.writeVidMapJSONFile(vidMapJSONFile.getAbsolutePath(), mergedHeaderLines);
        } catch (final FileNotFoundException fe) {
            throw new UserException("Unable to write vid map JSON file " + vidMapJSONFile.getAbsolutePath(), fe);
        }
        try {
            GenomicsDBImporter.writeCallsetMapJSONFile(callsetMapJSONFile.getAbsolutePath(), callsetMappingPB);
        } catch (final FileNotFoundException fe) {
            throw new UserException("Unable to write callset map JSON file " + callsetMapJSONFile.getAbsolutePath(), fe);
        }

        if (doConsolidation) {
            logger.info("GenomicsDB consolidation started");
            GenomicsDBImporter.consolidateTileDBArray(workspace, GenomicsDBConstants.DEFAULT_ARRAY_NAME);
            logger.info("GenomicsDB consolidation completed");
        }

        return true;
    }

    /**
     * Method to create feature readers for input files or GCS URLs
     * in the current batch
     *
     * @param sampleNames  List of all sample names
     * @param sampleNameToFileName  Sample name to file name mapping
     * @param batchSize  Current batch size
     * @param sampleCount  total samples in this import
     * @param lowerSampleIndex  0-based Lower bound of sample index -- inclusive
     * @return  Feature readers to be imported in the current batch
     */
    private Map<String, FeatureReader<VariantContext>> getFeatureReaders(final List<String> sampleNames, final Map<String, String> sampleNameToFileName, final int batchSize, int sampleCount, int lowerSampleIndex) {
        final Map<String, FeatureReader<VariantContext>> sampleToReaderMap = new LinkedHashMap<>();

        for(int i = lowerSampleIndex; i < sampleCount && i < lowerSampleIndex+batchSize; ++i) {
            final String sampleName = sampleNames.get(i);
            assert sampleNameToFileName.containsKey(sampleName);

            final AbstractFeatureReader<VariantContext, LineIterator> reader = getReaderFromVCFUri(sampleNameToFileName.get(sampleName));

            assert sampleName.equals(((VCFHeader) reader.getHeader()).getGenotypeSamples().get(0));
            sampleToReaderMap.put(sampleName, reader);
        }
        return sampleToReaderMap;
    }

    /**
     * Creates a feature reader object from a given VCF URI (can also be
     * a local file path) and returns it
     * @param variantPath  URI or file path
     * @return  Feature reader
     */
    private AbstractFeatureReader<VariantContext, LineIterator> getReaderFromVCFUri(String variantPath) {
        final String variantURI = IOUtils.getPath(variantPath).toAbsolutePath().toUri().toString();
        final Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper = (cloudPrefetchBuffer > 0 ? is -> SeekableByteChannelPrefetcher.addPrefetcher(cloudPrefetchBuffer, is) : Function.identity());
        final Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper = (cloudIndexPrefetchBuffer > 0 ? is -> SeekableByteChannelPrefetcher.addPrefetcher(cloudIndexPrefetchBuffer, is) : Function.identity());
        return AbstractFeatureReader.getFeatureReader(variantURI, null, new VCFCodec(), true, cloudWrapper, cloudIndexWrapper);
    }

    /**
     * Creates a GenomicsDB configuration data structure
     * instead of sending a long list of parameters to the constructor call
     *
     * @param workspace  GenomicsDB workspace
     * @param arrayName  GenomicsDB array
     * @param variantContextBufferSize  Buffer size to store VCF records for all samples
     * @param segmentSize  Buffer size to store columnar data to be serialized to disk
     * @param lbSampleIndex  Lower bound of sample index -- inclusive (0-based)
     * @param ubSampleIndex  Upper bound of sample index -- inclusive (0-based)
     * @return  GenomicsDB import configuration object
     */
    private static GenomicsDBImportConfiguration.ImportConfiguration createImportConfiguration(
        final String workspace,
        final String arrayName,
        final long variantContextBufferSize,
        final long segmentSize,
        final  long lbSampleIndex,
        final long ubSampleIndex) {

        GenomicsDBImportConfiguration.Partition.Builder pBuilder =
            GenomicsDBImportConfiguration.Partition.newBuilder();

        // Since, there is one partition for this import, the
        // begin column partition index is 0
        GenomicsDBImportConfiguration.Partition partition =
            pBuilder
                .setWorkspace(workspace)
                .setArray(arrayName)
                .setBegin(0)
                .build();

        GenomicsDBImportConfiguration.GATK4Integration.Builder gBuilder =
            GenomicsDBImportConfiguration.GATK4Integration.newBuilder();

        GenomicsDBImportConfiguration.GATK4Integration gatk4Parameters =
            gBuilder
                .setLowerSampleIndex(lbSampleIndex)
                .setUpperSampleIndex(ubSampleIndex)
                .build();

        GenomicsDBImportConfiguration.ImportConfiguration.Builder cBuilder =
            GenomicsDBImportConfiguration.ImportConfiguration.newBuilder();

        return cBuilder
                .addColumnPartitions(0, partition)
                .setGatk4IntegrationParameters(gatk4Parameters)
                .setSizePerColumnPartition(variantContextBufferSize)
                .setSegmentSize(segmentSize)
                .build();
    }

    /**
     * Close all readers in the current batch
     *
     * @param sampleToReaderMap  Map of sample names to readers
     */
    private void closeReaders(Map<String, FeatureReader<VariantContext>> sampleToReaderMap) {
        for (Map.Entry<String, FeatureReader<VariantContext>> reader : sampleToReaderMap.entrySet()) {
            try {
                reader.getValue().close();
            } catch (IOException e) {
                throw new GATKException("FeatureReader close() failed for " + reader.getKey(), e);
            }
        }
    }

    /**
     * Input argument "overwriteExistingWorkspace" defaults to false.
     * The tool creates a new workspace if it doesn't exist. Deletes
     * an existing workspace if argument is true
     *
     * @return  The workspace directory
     */
    private File overwriteOrCreateWorkspace() {
        File workspaceDir = new File(workspace);

        if (overwriteExistingWorkspace) {
            IOUtils.tryDelete(workspaceDir);
        }

        if (!workspaceDir.exists()) {
            int ret = GenomicsDBImporter.createTileDBWorkspace(workspaceDir.getAbsolutePath());

            if (ret > 0) {
                checkIfValidWorkspace(workspaceDir);
                logger.info("Importing data to GenomicsDB workspace: " + workspaceDir);
            } else if (ret < 0) {
                throw new UserException("Error creating GenomicsDB workspace: " + workspaceDir);
            }
        } else {
            // Check whether its a valid workspace
            checkIfValidWorkspace(workspaceDir);
        }
        return workspaceDir;
    }

    private void checkIfValidWorkspace(File workspaceDir) {
        File tempFile = new File(workspaceDir.getAbsolutePath() + "/__tiledb_workspace.tdb");
        if (!tempFile.exists()) {
            throw new UserException(workspaceDir.getAbsolutePath() + " is not a valid GenomicsDB workspace");
        }
    }

    /**
     * Loads our intervals using the best available sequence
     * dictionary (as returned by {@link #getBestAvailableSequenceDictionary})
     * to parse/verify them. Does nothing if no intervals were specified.
     */
    private void initializeIntervals() {
        if (intervalArgumentCollection.intervalsSpecified()) {

            final SAMSequenceDictionary intervalDictionary = getBestAvailableSequenceDictionary();
            if (intervalDictionary == null) {
                throw new UserException("We require at least one input source that " +
                    "has a sequence dictionary (reference or reads) when intervals are specified");
            }

            intervals = new ArrayList<>();

            List<SimpleInterval> simpleIntervalList =
                intervalArgumentCollection.getIntervals(intervalDictionary);

            if (simpleIntervalList.size() > 1) {
                throw new UserException("More than one interval specified. The tool takes only one");
            }

            for (SimpleInterval simpleInterval : simpleIntervalList) {
                intervals.add(new ChromosomeInterval(simpleInterval.getContig(),
                  simpleInterval.getStart(), simpleInterval.getEnd()));
            }
        } else {
            throw new UserException("No intervals specified");
        }
    }

    /**
     * Overriding getBestAvailableSequenceDictionary() to prefer the mergedVCFHeader's
     * sequence directory, if present, over any other dictionaries
     *
     * @return  Sequence directory from merged header, or super.getBestAvailableSequenceDictionary()
     *          if none
     */
    @Override
    public SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        final SAMSequenceDictionary sequenceDictionary = mergedHeaderSequenceDictionary;
        if (sequenceDictionary == null) {
            return super.getBestAvailableSequenceDictionary();
        } else {
            return sequenceDictionary;
        }
    }
}

