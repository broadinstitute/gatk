package org.broadinstitute.hellbender.tools.genomicsdb;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.intel.genomicsdb.importer.GenomicsDBImporter;
import com.intel.genomicsdb.importer.model.ChromosomeInterval;
import com.intel.genomicsdb.model.Coordinates;
import com.intel.genomicsdb.model.GenomicsDBCallsetsMapProto;
import com.intel.genomicsdb.model.GenomicsDBImportConfiguration;
import com.intel.genomicsdb.model.ImportConfig;
import com.intel.genomicsdb.model.BatchCompletionCallbackFunctionArgument;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.TribbleException;
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
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.nio.SeekableByteChannelPrefetcher;

import java.io.File;
import java.io.IOException;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;
import java.util.function.Function;
import java.util.stream.Collectors;


/**
 * Import single-sample GVCFs into GenomicsDB before joint genotyping.
 *
 * <p>The GATK4 Best Practice Workflow for SNP and Indel calling uses GenomicsDBImport to merge GVCFs from multiple samples.
 * GenomicsDBImport offers the same functionality as CombineGVCFs and comes from the <i>Intel-Broad Center for Genomics</i>.
 * The datastore transposes sample-centric variant information across genomic loci to make data more accessible to tools.
 * </p>
 *
 * <p>To query the contents of the GenomicsDB datastore, use
 * <a href='https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php'>SelectVariants</a>.
 * See <a href='https://software.broadinstitute.org/gatk/documentation/article?id=10061'>Tutorial#10061</a> to get started. </p>
 *
 * <p>Details on GenomicsDB are at
 * <a href='https://github.com/Intel-HLS/GenomicsDB/wiki'>https://github.com/Intel-HLS/GenomicsDB/wiki</a>.
 * In brief, GenomicsDB is a utility built on top of TileDB. TileDB is a format for efficiently representing sparse data.
 * Genomics data is typically sparse in that each sample has few variants with respect to the entire reference genome.
 * GenomicsDB contains code to specialize TileDB for genomics applications, such as VCF parsing and INFO field annotation
 * calculation.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * One or more GVCFs produced by in HaplotypeCaller with the `-ERC GVCF` or `-ERC BP_RESOLUTION` settings, containing
 * the samples to joint-genotype.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A GenomicsDB workspace
 * </p>
 *
 *  <h3>Usage examples</h3>
 *
 *  Provide each sample GVCF separately.
 *  <pre>
 *    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
 *      -V data/gvcfs/mother.g.vcf.gz \
 *      -V data/gvcfs/father.g.vcf.gz \
 *      -V data/gvcfs/son.g.vcf.gz \
 *      --genomicsdb-workspace-path my_database \
 *      --tmp-dir=/path/to/large/tmp \
 *      -L 20
 *  </pre>
 *
 *  Provide sample GVCFs in a map file.
 *
 *  <pre>
 *    gatk --java-options "-Xmx4g -Xms4g" \
 *       GenomicsDBImport \
 *       --genomicsdb-workspace-path my_database \
 *       --batch-size 50 \
 *       -L chr1:1000-10000 \
 *       --sample-name-map cohort.sample_map \
 *       --tmp-dir=/path/to/large/tmp \
 *       --reader-threads 5
 *  </pre>
 *
 *  The sample map is a tab-delimited text file with sample_name--tab--path_to_sample_vcf per line. Using a sample map
 *  saves the tool from having to download the GVCF headers in order to determine the sample names. Sample names in
 *  the sample name map file may have non-tab whitespace, but may not begin or end with whitespace.
 *
 *  <pre>
 *  sample1      sample1.vcf.gz
 *  sample2      sample2.vcf.gz
 *  sample3      sample3.vcf.gz
 *  </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 *     <li>IMPORTANT: The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB, as the native TileDB library requires additional memory on top of the Java memory. Failure to leave enough memory for the native code can result in confusing error messages!</li>
 *     <li>At least one interval must be provided</li>
 *     <li>Input GVCFs cannot contain multiple entries for a single genomic position</li>
 *     <li>The --genomicsdb-workspace-path must point to a non-existent or empty directory.</li>
 *     <li>GenomicsDBImport uses temporary disk storage during import. The amount of temporary disk storage required can exceed the space available, especially when specifying a large number of intervals. The command line argument `--tmp-dir` can be used to specify an alternate temporary storage location with sufficient space..</li>
 * </ul>
 *
 * <h3>Developer Note</h3>
 * To read data from GenomicsDB, use the query interface {@link com.intel.genomicsdb.reader.GenomicsDBFeatureReader}
 */
@DocumentedFeature
@CommandLineProgramProperties(
    summary = "Import VCFs to GenomicsDB",
    oneLineSummary = "Import VCFs to GenomicsDB",
    programGroup = ShortVariantDiscoveryProgramGroup.class
)
public final class GenomicsDBImport extends GATKTool {

    private static final long DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE = 16*1024L;
    private static final long DEFAULT_SEGMENT_SIZE = 1048576L;
    private static final int DEFAULT_ZERO_BATCH_SIZE = 0;

    public static final String WORKSPACE_ARG_LONG_NAME = "genomicsdb-workspace-path";
    public static final String SEGMENT_SIZE_ARG_LONG_NAME = "genomicsdb-segment-size";
    public static final String OVERWRITE_WORKSPACE_LONG_NAME = "overwrite-existing-genomicsdb-workspace";

    public static final String VCF_BUFFER_SIZE_ARG_NAME = "genomicsdb-vcf-buffer-size";

    public static final String BATCHSIZE_ARG_LONG_NAME = "batch-size";
    public static final String CONSOLIDATE_ARG_NAME = "consolidate";
    public static final String SAMPLE_NAME_MAP_LONG_NAME = "sample-name-map";
    public static final String VALIDATE_SAMPLE_MAP_LONG_NAME = "validate-sample-name-map";
    public static final String VCF_INITIALIZER_THREADS_LONG_NAME = "reader-threads";
    public static final String MAX_NUM_INTERVALS_TO_IMPORT_IN_PARALLEL = "max-num-intervals-to-import-in-parallel";
    public static final int INTERVAL_LIST_SIZE_WARNING_THRESHOLD = 100;

    @Argument(fullName = WORKSPACE_ARG_LONG_NAME,
              doc = "Workspace for GenomicsDB. Must be a POSIX file system path, but can be a relative path." +
                      " Must be an empty or non-existent directory.")
    private String workspace;

    @Argument(fullName = SEGMENT_SIZE_ARG_LONG_NAME,
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
                    " Defaults to " + DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE + " which was empirically determined to work" +
                    " well for many inputs.",
              optional = true,
              minValue = 1024L,
              minRecommendedValue = 10 * 1024)
    private long vcfBufferSizePerSample = DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE;

    @Argument(fullName = OVERWRITE_WORKSPACE_LONG_NAME,
              doc = "Will overwrite given workspace if it exists. " +
                    "Otherwise a new workspace is created. " +
                    "Defaults to false",
              optional = true)
    private Boolean overwriteExistingWorkspace = false;

    @Argument(fullName = BATCHSIZE_ARG_LONG_NAME,
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
            doc = "Path to file containing a mapping of sample name to file uri in tab delimited format.  If this is " +
                    "specified then the header from the first sample will be treated as the merged header rather than " +
                    "merging the headers, and the sample names will be taken from this file.  This may be used to rename " +
                    "input samples. This is a performance optimization that relaxes the normal checks for consistent " +
                    "headers.  Using vcfs with incompatible headers may result in silent data corruption.",
            optional = true,
            mutex = {StandardArgumentDefinitions.VARIANT_LONG_NAME})
    private String sampleNameMapFile;

    @Argument(fullName = VALIDATE_SAMPLE_MAP_LONG_NAME,
            shortName = VALIDATE_SAMPLE_MAP_LONG_NAME,
            doc = "Boolean flag to enable checks on the sampleNameMap file. If true, tool checks whether" +
                "feature readers are valid and shows a warning if sample names do not match with the headers. " +
                "Defaults to false",
            optional = true)
    private Boolean validateSampleToReaderMap = false;

    @Advanced
    @Argument(fullName = VCF_INITIALIZER_THREADS_LONG_NAME,
            shortName = VCF_INITIALIZER_THREADS_LONG_NAME,
            doc = "How many simultaneous threads to use when opening VCFs in batches; higher values may improve performance " +
                    "when network latency is an issue",
            optional = true,
            minValue = 1)
    private int vcfInitializerThreads = 1;

    @Advanced
    @Argument(fullName = MAX_NUM_INTERVALS_TO_IMPORT_IN_PARALLEL,
            shortName = MAX_NUM_INTERVALS_TO_IMPORT_IN_PARALLEL,
            doc = "Max number of intervals to import in parallel; higher values may improve performance, but require more" +
                  " memory and a higher number of file descriptors open at the same time",
            optional = true,
            minValue = 1)
    private int maxNumIntervalsToImportInParallel = 1;

    //executor service used when vcfInitializerThreads > 1
    private ExecutorService inputPreloadExecutorService;

    @Override
    public boolean requiresIntervals() { return true; }

    @Override
    public int getDefaultCloudPrefetchBufferSize() {
        // Empirical testing has shown that this tool performs best at scale with cloud buffering
        // disabled. With cloud buffering on and thousands of concurrent GenomicsDBImport tasks,
        // we do too many simultaneous GCS accesses (since the prefetcher spawns a new thread for each
        // reader upon a query) and start seeing intermittent failures, even with aggressive retries.
        return 0;
    }

    @Override
    public int getDefaultCloudIndexPrefetchBufferSize() {
        // Empirical testing has shown that this tool performs best at scale with cloud buffering
        // disabled. With cloud buffering on and thousands of concurrent GenomicsDBImport tasks,
        // we do too many simultaneous GCS accesses (since the prefetcher spawns a new thread for each
        // reader upon a query) and start seeing intermittent failures, even with aggressive retries.
        return 0;
    }

    @Override
    public String getProgressMeterRecordLabel() { return "batches"; }

    // Intervals from command line (singleton for now)
    private List<ChromosomeInterval> intervals;

    // Sorted mapping between sample names and corresponding GVCF file name
    //
    // IMPORTANT: This must be sorted or it will result in sample name swaps in the output database.
    // This happens because the callset json is generated independently from the import process
    // each imported batch is then sorted, so if we have an unsorted list we'll end up with different global vs batch
    // sorting.
    // We preemptively sort here so we will have consistent sorting.
    private SortedMap<String, Path> sampleNameToVcfPath = new TreeMap<>();

    // Needed as smartMergeHeaders() returns a set of VCF header lines
    private Set<VCFHeaderLine> mergedHeaderLines = null;

    // sequence dictionary created from the merged header
    private SAMSequenceDictionary mergedHeaderSequenceDictionary;

    // Path to vidmap file to be written by GenomicsDBImporter
    private File vidMapJSONFile;

    // Path to callsetmap file to be written by GenomicsDBImporter
    private File callsetMapJSONFile;

    // Path to combined VCF header file to be written by GenomicsDBImporter
    private File vcfHeaderFile;

    // GenomicsDB callset map protobuf structure containing all callset names
    // used to write the callset json file on traversal success
    private GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB;

    //in-progress batchCount
    private int batchCount = 1;

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
        if ( (variantPaths == null || variantPaths.isEmpty()) && sampleNameMapFile == null) {
            throw new CommandLineException.MissingArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME,
                                                       "One of --" + StandardArgumentDefinitions.VARIANT_LONG_NAME + " or --" + SAMPLE_NAME_MAP_LONG_NAME + " must be specified" );
        }
    }

    /**
     * sets the values of mergedHeaderLines, mergedHeaderSequenceDictionary, and sampleNameToVcfPath
     */
    private void initializeHeaderAndSampleMappings() {
        // Only one of -V and --sampleNameMapFile may be specified
        if (sampleNameMapFile == null) {
            // -V was specified
            final List<VCFHeader> headers = new ArrayList<>(variantPaths.size());
            for (final String variantPathString : variantPaths) {
                final Path variantPath = IOUtils.getPath(variantPathString);
                final  VCFHeader header = getHeaderFromPath(variantPath);
                Utils.validate(header != null, "Null header was found in " + variantPath + ".");
                assertGVCFHasOnlyOneSample(variantPathString, header);
                headers.add(header);

                final String sampleName = header.getGenotypeSamples().get(0);
                final Path previousPath = sampleNameToVcfPath.put(sampleName, variantPath);
                if (previousPath != null) {
                    throw new UserException("Duplicate sample: " + sampleName + ". Sample was found in both "
                                                    + variantPath.toUri() + " and " + previousPath.toUri() + ".");
                }
            }
            mergedHeaderLines = VCFUtils.smartMergeHeaders(headers, true);
            mergedHeaderSequenceDictionary = new VCFHeader(mergedHeaderLines).getSequenceDictionary();

        } else {
            // --sampleNameMap was specified

            //it's VERY IMPORTANT that this map is Sorted according to String's natural ordering, if it is not
            //the resulting database will have incorrect sample names
            //see https://github.com/broadinstitute/gatk/issues/3682 for more information
            sampleNameToVcfPath = loadSampleNameMapFileInSortedOrder(IOUtils.getPath(sampleNameMapFile));
            final Path firstHeaderPath = sampleNameToVcfPath.entrySet().iterator().next().getValue();
            final VCFHeader header = getHeaderFromPath(firstHeaderPath);
            //getMetaDataInInputOrder() returns an ImmutableSet - LinkedHashSet is mutable and preserves ordering
            mergedHeaderLines = new LinkedHashSet<VCFHeaderLine>(header.getMetaDataInInputOrder());
            mergedHeaderSequenceDictionary = header.getSequenceDictionary();
        }

        mergedHeaderLines.addAll(getDefaultToolVCFHeaderLines());

        if ( mergedHeaderSequenceDictionary == null) {
            throw new UserException("The merged vcf header has no sequence dictionary. Please provide a header that contains a sequence dictionary.");
        }
    }

    private VCFHeader getHeaderFromPath(final Path variantPath) {
        try(final AbstractFeatureReader<VariantContext, LineIterator> reader = getReaderFromPath(variantPath)) {
            return (VCFHeader) reader.getHeader();
        } catch (final IOException e) {
            throw new UserException("Error while reading vcf header from " + variantPath.toUri(), e);
        }
    }

    private static void assertGVCFHasOnlyOneSample(final String variantPath, final VCFHeader header) {
        // A GVCF file must contain only one sample, throw an exception otherwise
        final int numberOfSamples = header.getNGenotypeSamples();
        if (numberOfSamples != 1) {
            throw new UserException("Input GVCF: " + variantPath + " was expected to contain a single sample but actually contained " + numberOfSamples + " samples.");
        }
    }

    /**
     * load a tab delimited new line separated file of sample name to URI mapping:
     * this maintains the keys in the same order that they appeared in the file
     *
     * this tool should only call {@link #loadSampleNameMapFileInSortedOrder(Path)},
     * this version is exposed for the benefit of {@link org.broadinstitute.hellbender.tools.FixCallSetSampleOrdering}
     *
     * ex:
     *
     * Sample1\tpathToSample1.vcf\n
     * Sample2\tpathTosample2.vcf\n
     * ...
     *
     * The sample names must be unique.
     * @param sampleToFileMapPath path to the mapping file
     * @return map of sample name to corresponding file, the map will be ordered according to the order in the input file
     */
    public static LinkedHashMap<String, Path> loadSampleNameMapFile(final Path sampleToFileMapPath) {
        try {
            final List<String> lines = Files.readAllLines(sampleToFileMapPath);
            if (lines.isEmpty()) {
                throw new UserException.BadInput( "At least 1 sample is required but none were found in the sample mapping file");
            }

            final LinkedHashMap<String, Path> sampleToFilename = new LinkedHashMap<>();
            for ( final String line : lines) {
                final String[] split = line.split("\\t",-1);
                if (split.length != 2) {
                    throw new UserException.BadInput("Expected a file with 2 fields per line in the format\nSample\tFile\n but found line: \""
                            + line +"\" with "+split.length+" fields");
                }
                if ( !split[0].trim().equals(split[0]) || split[0].trim().isEmpty()
                        || split[1].trim().isEmpty()) {
                    throw new UserException.BadInput("Expected a file of format\nSample\tFile\n but found line: '" + line + "'\nValid sample names must be non-empty strings that cannot begin or end with whitespace and valid file names must be non-empty and not all whitespace");
                }
                final String sample = split[0];
                final String path = split[1].trim();
                final Path oldPath = sampleToFilename.put(sample, IOUtils.getPath(path));
                if (oldPath != null){
                    throw new UserException.BadInput("Found two mappings for the same sample: " + sample + "\n" + path + "\n" + oldPath.toUri() );
                }
            }
            return sampleToFilename;
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(sampleToFileMapPath, "exception while reading sample->filename mapping file",  e);
        }
    }

    /**
     * load a tab delimited new line separated file of sample name to URI mapping:
     *
     * ex:
     * Sample1\tpathToSample1.vcf\n
     * Sample2\tpathTosample2.vcf\n
     * ...
     *
     * The sample names must be unique.
     * @param sampleToFileMapPath path to the mapping file
     * @return map of sample name to corresponding file, sorted by sample name
     */
    public static SortedMap<String, Path> loadSampleNameMapFileInSortedOrder(final Path sampleToFileMapPath){
        return new TreeMap<>(loadSampleNameMapFile(sampleToFileMapPath));
    }

    /**
     * Before traversal, fix configuration parameters and initialize
     * GenomicsDB. Hard-coded to handle only VCF files and headers
     */
    @Override
    public void onTraversalStart() {

        final File workspaceDir = overwriteOrCreateWorkspace();

        vidMapJSONFile = new File(workspaceDir + "/" + GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME);
        callsetMapJSONFile = new File(workspaceDir + "/" + GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME);
        vcfHeaderFile = new File(workspaceDir + "/" + GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME);

        logger.info("Vid Map JSON file will be written to " + vidMapJSONFile);
        logger.info("Callset Map JSON file will be written to " + callsetMapJSONFile);
        logger.info("Complete VCF Header will be written to " + vcfHeaderFile);
        logger.info("Importing to array - " + workspace + "/" + GenomicsDBConstants.DEFAULT_ARRAY_NAME);

        initializeInputPreloadExecutorService();
    }

    private void initializeInputPreloadExecutorService() {
        if( vcfInitializerThreads > 1) {
            final ThreadFactory threadFactory = new ThreadFactoryBuilder()
                    .setNameFormat("readerInitializer-thread-%d")
                    .setDaemon(true)
                    .build();
            this.inputPreloadExecutorService = Executors.newFixedThreadPool(vcfInitializerThreads, threadFactory);
        } else {
            inputPreloadExecutorService = null;
        }
    }

    private Map<String, FeatureReader<VariantContext>> createSampleToReaderMap(
            final Map<String, Path> sampleNameToVcfPath, final int batchSize, final int index) {
        // TODO: fix casting since it's really ugly
        return inputPreloadExecutorService != null ?
                getFeatureReadersInParallel((SortedMap<String, Path>) sampleNameToVcfPath, batchSize, index)
                : getFeatureReadersSerially(sampleNameToVcfPath, batchSize, index);
    }

    private Void logMessageOnBatchCompletion(final BatchCompletionCallbackFunctionArgument arg) {
        progressMeter.update(intervals.get(0));
        logger.info("Done importing batch " + arg.batchCount + "/" + arg.totalBatchCount);
        this.batchCount = arg.batchCount + 1;
        return null;
    }

    private List<GenomicsDBImportConfiguration.Partition> generatePartitionListFromIntervals(List<ChromosomeInterval> chromosomeIntervals) {
        return chromosomeIntervals.stream().map(interval -> {
            GenomicsDBImportConfiguration.Partition.Builder partitionBuilder = GenomicsDBImportConfiguration.Partition.newBuilder();
            Coordinates.ContigPosition.Builder contigPositionBuilder = Coordinates.ContigPosition.newBuilder();
            Coordinates.GenomicsDBColumn.Builder columnBuilder = Coordinates.GenomicsDBColumn.newBuilder();
            //begin
            contigPositionBuilder.setContig(interval.getContig()).setPosition(interval.getStart());
            columnBuilder.setContigPosition(contigPositionBuilder.build());
            partitionBuilder.setBegin(columnBuilder.build());
            //end
            contigPositionBuilder.setPosition(interval.getEnd());
            columnBuilder.setContigPosition(contigPositionBuilder.build());
            partitionBuilder.setEnd(columnBuilder.build());
            partitionBuilder.setWorkspace(workspace);
            partitionBuilder.setGenerateArrayNameFromPartitionBounds(true);
            return partitionBuilder.build();
        }).collect(Collectors.toList());
    }

    private ImportConfig createImportConfig(final int batchSize) {
        final List<GenomicsDBImportConfiguration.Partition> partitions = generatePartitionListFromIntervals(intervals);
        GenomicsDBImportConfiguration.ImportConfiguration.Builder importConfigurationBuilder =
                GenomicsDBImportConfiguration.ImportConfiguration.newBuilder();
        importConfigurationBuilder.addAllColumnPartitions(partitions);
        importConfigurationBuilder.setSizePerColumnPartition(vcfBufferSizePerSample);
        importConfigurationBuilder.setFailIfUpdating(true);
        importConfigurationBuilder.setSegmentSize(segmentSize);
        importConfigurationBuilder.setConsolidateTiledbArrayAfterLoad(doConsolidation);
        ImportConfig importConfig = new ImportConfig(importConfigurationBuilder.build(), validateSampleToReaderMap, true,
                batchSize, mergedHeaderLines, sampleNameToVcfPath, this::createSampleToReaderMap);
        importConfig.setOutputCallsetmapJsonFile(callsetMapJSONFile.getAbsolutePath());
        importConfig.setOutputVidmapJsonFile(vidMapJSONFile.getAbsolutePath());
        importConfig.setOutputVcfHeaderFile(vcfHeaderFile.getAbsolutePath());
        importConfig.setUseSamplesInOrder(true);
        importConfig.setFunctionToCallOnBatchCompletion(this::logMessageOnBatchCompletion);
        return importConfig;
    }

    /**
     * A complete traversal from start to finish. This method will import all samples
     * specified in the input GVCF files.
     */
    @Override
    public void traverse() {
        // Force the progress meter to update after every batch
        progressMeter.setRecordsBetweenTimeChecks(1L);

        final int sampleCount = sampleNameToVcfPath.size();
        final int updatedBatchSize = (batchSize == DEFAULT_ZERO_BATCH_SIZE) ? sampleCount : batchSize;
        final ImportConfig importConfig = createImportConfig(updatedBatchSize);

        GenomicsDBImporter importer;
        try {
            importer = new GenomicsDBImporter(importConfig);
            importer.executeImport(maxNumIntervalsToImportInParallel);
        } catch (final IOException e) {
            throw new UserException("Error initializing GenomicsDBImporter", e);
        } catch (final IllegalArgumentException iae) {
            throw new GATKException("Null feature reader found in sampleNameMap file: " + sampleNameMapFile, iae);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        if (batchSize == DEFAULT_ZERO_BATCH_SIZE) {
            logger.info("Import completed!");
        } else {
            logger.info("Import of all batches to GenomicsDB completed!");
        }
        return true;
    }

    /**
     * Method to create feature readers for input files or GCS URLs
     * in the current batch
     *
     * @param sampleNametoPath  Sample name to file name mapping
     * @param batchSize  Current batch size
     * @param lowerSampleIndex  0-based Lower bound of sample index -- inclusive
     * @return  Feature readers to be imported in the current batch, sorted by sample name
     */
    private SortedMap<String, FeatureReader<VariantContext>> getFeatureReadersInParallel(
            final SortedMap<String, Path> sampleNametoPath, final int batchSize, final int lowerSampleIndex) {
        final SortedMap<String, FeatureReader<VariantContext>> sampleToReaderMap = new TreeMap<>();
        logger.info("Starting batch input file preload");
        final Map<String, Future<FeatureReader<VariantContext>>> futures = new LinkedHashMap<>();
        final List<String> sampleNames = new ArrayList<>(sampleNametoPath.keySet());
        for(int i = lowerSampleIndex; i < sampleNametoPath.size() && i < lowerSampleIndex+batchSize; ++i) {
            final String sampleName = sampleNames.get(i);
            futures.put(sampleName, inputPreloadExecutorService.submit(() -> {
                final Path variantPath = sampleNametoPath.get(sampleName);
                try {
                    return new InitializedQueryWrapper(getReaderFromPath(variantPath), intervals.get(0));
                } catch (final IOException e) {
                    throw new UserException.CouldNotReadInputFile("Couldn't read file: " + variantPath.toUri(), e);
                }
            }));
        }

        futures.forEach((sampleName, future) -> {
            try {
                final FeatureReader<VariantContext> reader = future.get();
                sampleToReaderMap.put(sampleName, reader);
            } catch (InterruptedException | ExecutionException e) {
                throw new UserException.CouldNotReadInputFile("Failure while waiting for FeatureReader to initialize ",
                                                              e);
            }
        });
        logger.info("Finished batch preload");
        logger.info("Importing batch " + this.batchCount + " with " + sampleToReaderMap.size() + " samples");
        return sampleToReaderMap;
    }

    private SortedMap<String, FeatureReader<VariantContext>> getFeatureReadersSerially(final Map<String, Path> sampleNameToPath,
                                                                                 final int batchSize, final int lowerSampleIndex){
        final SortedMap<String, FeatureReader<VariantContext>> sampleToReaderMap = new TreeMap<>();
        final List<String> sampleNames = new ArrayList<>(sampleNameToPath.keySet());
        for(int i = lowerSampleIndex; i < sampleNameToPath.size() && i < lowerSampleIndex+batchSize; ++i) {
            final String sampleName = sampleNames.get(i);
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getReaderFromPath(sampleNameToPath.get(sampleName));
            sampleToReaderMap.put(sampleName, reader);
        }
        logger.info("Importing batch " + this.batchCount + " with " + sampleToReaderMap.size() + " samples");
        return sampleToReaderMap;
    }

    /**
     * Creates a feature reader object from a given VCF URI (can also be
     * a local file path) and returns it
     * @return  Feature reader
     * @param variantPath
     */
    private AbstractFeatureReader<VariantContext, LineIterator> getReaderFromPath(final Path variantPath) {
        final String variantURI = variantPath.toAbsolutePath().toUri().toString();
        final Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper = (cloudPrefetchBuffer > 0 ? is -> SeekableByteChannelPrefetcher.addPrefetcher(cloudPrefetchBuffer, is) : Function.identity());
        final Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper = (cloudIndexPrefetchBuffer > 0 ? is -> SeekableByteChannelPrefetcher.addPrefetcher(cloudIndexPrefetchBuffer, is) : Function.identity());
        try {
            return AbstractFeatureReader.getFeatureReader(variantURI, null, new VCFCodec(), true, cloudWrapper, cloudIndexWrapper);
        } catch (final TribbleException e){
            throw new UserException("Failed to create reader from " + variantURI, e);
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
        final File workspaceDir = new File(workspace);

        if (overwriteExistingWorkspace) {
            IOUtils.tryDelete(workspaceDir);
        }

        if (!workspaceDir.exists()) {
            final int ret = GenomicsDBImporter.createTileDBWorkspace(workspaceDir.getAbsolutePath());
            if (ret > 0) {
                checkIfValidWorkspace(workspaceDir);
                logger.info("Importing data to GenomicsDB workspace: " + workspaceDir);
            } else if (ret < 0) {
                throw new UnableToCreateGenomicsDBWorkspace("Error creating GenomicsDB workspace: " + workspaceDir);
            }
            return workspaceDir;
        } else {
            throw new UnableToCreateGenomicsDBWorkspace("The workspace you're trying to create already exists. ( " + workspaceDir.getAbsolutePath() + " ) " +
                                                  "Writing into an existing workspace can cause data corruption. " +
                                                  "Please choose an output path that doesn't already exist. ");
        }
    }

    static class UnableToCreateGenomicsDBWorkspace extends UserException {
        private static final long serialVersionUID = 1L;

        UnableToCreateGenomicsDBWorkspace(final String message){
            super(message);
        }
    }

    private static void checkIfValidWorkspace(final File workspaceDir) {
        final File tempFile = new File(workspaceDir.getAbsolutePath() + "/__tiledb_workspace.tdb");
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
            final List<SimpleInterval> simpleIntervalList = intervalArgumentCollection.getIntervals(intervalDictionary);
            if (simpleIntervalList.size() > INTERVAL_LIST_SIZE_WARNING_THRESHOLD) {
                logger.warn(String.format(
                        "A large number of intervals were specified. " +
                        "Using more than %d intervals in a single import is not recommended and can cause performance to suffer. " +
                        "It is recommended that intervals be aggregated together.",
                        INTERVAL_LIST_SIZE_WARNING_THRESHOLD)
                );
            }
            simpleIntervalList.forEach(interval -> intervals.add(new ChromosomeInterval(interval.getContig(),
                    interval.getStart(), interval.getEnd())));
        } else {
            throw new UserException("No intervals specified");
        }
    }

    @Override
    public void onShutdown(){
        if(inputPreloadExecutorService != null) {
            inputPreloadExecutorService.shutdownNow();
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

    /**
     * This class is a hack to force parallel loading of the headers and indexes of remote gvcf files.
     * It initializes a feature reader and starts a query.  This causes the header and index to be read, and also causes any
     * pre-fetching to begin if enabled.  It is very narrowly crafted and should not be used for other purposes.
     */
    private static final class InitializedQueryWrapper implements FeatureReader<VariantContext> {
        private final FeatureReader<VariantContext> reader;
        private final SimpleInterval interval;
        private CloseableTribbleIterator<VariantContext> query;

        private InitializedQueryWrapper(final FeatureReader<VariantContext> reader, final Locatable interval) throws IOException {
            this.reader = reader;
            this.interval = new SimpleInterval(interval);
            this.query = reader.query(interval.getContig(), interval.getStart(), interval.getEnd());
        }

        @Override
        public CloseableTribbleIterator<VariantContext> query(final String chr, final int start, final int end) {
            final SimpleInterval queryInterval = new SimpleInterval(chr, start, end);
            if( !interval.equals(queryInterval)){
                throw new GATKException("Cannot call query with different interval, expected:" + this.interval + " queried with: " + queryInterval);
            }
            if( query != null ){
                final CloseableTribbleIterator<VariantContext> tmp = query;
                query = null;
                return tmp;
            } else {
                throw new GATKException("Cannot call query twice on this wrapper.");
            }
        }

        @Override
        public CloseableTribbleIterator<VariantContext> iterator() {
            throw new UnsupportedOperationException("iterator() not supported, this should not have been called and indicates an issue with GenomicsDB integration");
        }

        @Override
        public void close() throws IOException {
            reader.close();
        }

        @Override
        public List<String> getSequenceNames() {
            throw new UnsupportedOperationException("getSequenceNames() not supported, this should not have been called and indicates an issue with GenomicsDB integration");
        }

        @Override
        public Object getHeader() {
            return reader.getHeader();
        }
    }
}

