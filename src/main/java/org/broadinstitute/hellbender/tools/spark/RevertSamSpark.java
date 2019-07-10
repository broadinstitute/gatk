package org.broadinstitute.hellbender.tools.spark;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.LinkedListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.table.TableCodec;
import org.broadinstitute.hellbender.utils.codecs.table.TableFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import scala.Tuple2;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Reverts a SAM file by optionally restoring original quality scores and by removing
 * all alignment information.
 * <p>
 * <p>
 * This tool removes or restores certain properties of the SAM records, including alignment information.
 * It can be used to produce an unmapped BAM (uBAM) from a previously aligned BAM. It is also capable of
 * restoring the original quality scores of a BAM file that has already undergone base quality score recalibration
 * (BQSR) if the original qualities were retained during the calibration (in the OQ tag).
 * <p>
 * <h3>Usage Examples</h3>
 * <h4>Output to a single file</h4>
 * <pre>
 * gatk RevertSamSpark  \\
 *      -I input.bam \\
 *      -O reverted.bam
 * </pre>
 * <p>
 * <h4>Output by read group into multiple files with sample map</h4>
 * <pre>
 * gatk RevertSamSpark \\
 *      -I input.bam \\
 *      --output-by-readgroup \\
 *      --output-map reverted_bam_paths.tsv
 * </pre>
 * <p>
 * <h4>Output by read group with no output map</h4>
 * <pre>
 * gatk RevertSamSpark \\
 *      -I input.bam \\
 *      --output-by-readgroup \\
 *      -O /write/reverted/read/group/bams/in/this/dir
 * </pre>
 * This will output a BAM (Can be overridden with outputByReadgroupFileFormat option.)
 * <br/>
 * Note: If the program fails due to a SAM validation error, consider setting the VALIDATION_STRINGENCY option to
 * LENIENT or SILENT if the failures are expected to be obviated by the reversion process
 * (e.g. invalid alignment information will be obviated when the keepAlignmentInformation option is used).
 */

@DocumentedFeature
@CommandLineProgramProperties(
        summary = RevertSamSpark.USAGE_DETAILS,
        oneLineSummary = RevertSamSpark.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@BetaFeature
public class RevertSamSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    static final String USAGE_SUMMARY = "Reverts SAM or BAM files to a previous state.";
    static final String USAGE_DETAILS = "This tool removes or restores certain properties of the SAM records, including alignment " +
            "information, which can be used to produce an unmapped BAM (uBAM) from a previously aligned BAM. It is also capable of " +
            "restoring the original quality scores of a BAM file that has already undergone base quality score recalibration (BQSR) if the" +
            "original qualities were retained.\n" +
            "<h3>Examples</h3>\n" +
            "<h4>Example with single output:</h4>\n" +
            "gatk RevertSamSpark \\\n" +
            "     -I input.bam \\\n" +
            "     -O reverted.bam\n" +
            "\n" +
            "<h4>Example outputting by read group with output map:</h4>\n" +
            "gatk RevertSamSpark \\\n" +
            "     -I input.bam \\\n" +
            "     --output-by-readgroup \\\n" +
            "     --output-map reverted_bam_paths.tsv\n" +
            "\n" +
            "Will output a BAM/SAM file per read group.\n" +
            "<h4>Example outputting by read group without output map:</h4>\n" +
            "gatk RevertSamSpark \\\n" +
            "     -I input.bam \\\n" +
            "     --output-by-readgroup \\\n" +
            "     -O /write/reverted/read/group/bams/in/this/dir\n" +
            "\n" +
            "Will output a BAM file per read group." +
            " Output format can be overridden with the outputByReadgroupFileFormat option.\n" +
            "Note: If the program fails due to a SAM validation error, consider setting the VALIDATION_STRINGENCY option to " +
            "LENIENT or SILENT if the failures are expected to be obviated by the reversion process " +
            "(e.g. invalid alignment information will be obviated when the keep-alignment-information option is used).\n" +
            "";
    public static final String OUTPUT_MAP_READ_GROUP_FIELD_NAME = "READ_GROUP_ID";
    public static final String OUTPUT_MAP_OUTPUT_FILE_FIELD_NAME = "OUTPUT";

    @Override
    public boolean requiresReads() { return true; }

    @Argument(mutex = {OUTPUT_MAP_LONG_NAME}, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "The output SAM/BAM file to create, or an output directory if '--output-by-readgroup' is set.")
    public String output;

    public static final String OUTPUT_MAP_LONG_NAME = "output-map";
    @Argument(mutex = {StandardArgumentDefinitions.OUTPUT_LONG_NAME},
            fullName = OUTPUT_MAP_LONG_NAME,
            doc = "Tab separated file with two columns, OUTPUT_MAP_READ_GROUP_FIELD_NAME and OUTPUT_MAP_OUTPUT_FILE_FIELD_NAME, providing file mapping only used if '--output-by-readgroup' is set.")
    public String outputMap;

    public static final String OUTPUT_BY_READGROUP_LONG_NAME = "output-by-readgroup";
    @Argument(fullName = OUTPUT_BY_READGROUP_LONG_NAME, doc = "When true, outputs each read group in a separate file.")
    public boolean outputByReadGroup = false;

    @Argument(doc = "WARNING: This option is potentially destructive. If enabled will discard reads in order to produce " +
            "a consistent output BAM. Reads discarded include (but are not limited to) paired reads with missing " +
            "mates, duplicated records, records with mismatches in length of bases and qualities. This option should " +
            "only be enabled if the output sort order is queryname and will always cause sorting to occur.")
    public boolean sanitize = false;

    public static final String KEEP_FIRST_DUPLICATE_LONG_NAME = "keep-first-duplicate";
    @Argument(doc = "If 'sanitize' only one record when we find more than one record with the same name for R1/R2/unpaired reads respectively. " +
            "For paired end reads, keeps only the first R1 and R2 found respectively, and discards all unpaired reads. " +
            "Duplicates do not refer to the duplicate flag in the FLAG field, but instead reads with the same name.",
            fullName = KEEP_FIRST_DUPLICATE_LONG_NAME)
    public boolean keepFirstDuplicate = false;

    public static final String OUTPUT_BY_READGROUP_FILE_FORMAT_LONG_NAME = "output-by-readgroup-file-format";
    @Argument(fullName = OUTPUT_BY_READGROUP_FILE_FORMAT_LONG_NAME, doc = "When using --output-by-readgroup, the output file format can be set to a certain format.")
    public FileType outputByReadgroupFileFormat = FileType.dynamic;

    @Argument(shortName = StandardArgumentDefinitions.SORT_ORDER_SHORT_NAME, fullName = StandardArgumentDefinitions.SORT_ORDER_LONG_NAME, doc = "The sort order to create the reverted output file with, defaults to whatever is specified in the current file", optional = true)
    public SAMFileHeader.SortOrder sortOrder = SAMFileHeader.SortOrder.queryname;

    public static final String DONT_RESTORE_ORIGINAL_QUALITIES_LONG_NAME = "dont-restore-original-qualities";
    @Argument( fullName = DONT_RESTORE_ORIGINAL_QUALITIES_LONG_NAME, doc = "Set to prevent the tool from setting the OQ field to the QUAL where available.", optional = true)
    public boolean dontRestoreOriginalQualities = false;

    public static final String DONT_REMOVE_DUPLICATE_INFORMATION_LONG_NAME = "remove-duplicate-information";
    @Argument(fullName = DONT_REMOVE_DUPLICATE_INFORMATION_LONG_NAME, doc = "By default we remove duplicate read flags from all reads.  Note that if this is true, " +
            " the output may have the unusual but sometimes desirable trait of having unmapped reads that are marked as duplicates.")
    public boolean keepDuplicateInformation = false;

    public static final String KEEP_ALIGNMENT_INFORMATION = "keep-alignment-information";
    @Argument(fullName = KEEP_ALIGNMENT_INFORMATION, doc = "Don't remove any of the alignment information from the file.")
    public boolean keepAlignmentInformation = false;

    public static final String ATTRIBUTE_TO_CLEAR_LONG_NAME = "attributes-to-clear";
    @Argument(fullName = ATTRIBUTE_TO_CLEAR_LONG_NAME, doc = "When removing alignment information, the set of optional tags to remove.", optional = true)
    public Set<String> attributesToClear = new HashSet<String>();

    public static final String REMOVE_DEFAULT_ATTRIBUTE_TO_CLEAR_LONG_NAME = "remove-default-attributes-to-clear";
    @Argument(fullName = REMOVE_DEFAULT_ATTRIBUTE_TO_CLEAR_LONG_NAME, doc = "When removing alignment information, the set of optional tags to remove.", optional = true)
    public boolean removeDefaults = false;

    public static final String SAMPLE_ALIAS_ARG = "sample-alias";
    @Argument(fullName = SAMPLE_ALIAS_ARG, doc = "The sample alias to use in the reverted output file.  This will override the existing " +
            "sample alias in the file and is used only if all the read groups in the input file have the " +
            "same sample alias.", shortName = StandardArgumentDefinitions.SAMPLE_ALIAS_SHORT_NAME, optional = true)
    public String sampleAlias;

    public static final String LIBRARY_NAME_ARG = "library-name";
    @Argument(fullName = LIBRARY_NAME_ARG, doc = "The library name to use in the reverted output file.  This will override the existing " +
            "sample alias in the file and is used only if all the read groups in the input file have the " +
            "same library name.", optional = true)
    public String libraryName;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    final public static List<String> DEFAULT_ATTRIBUTES_TO_CLEAR = Collections.unmodifiableList(new ArrayList<String>(){
        private static final long serialVersionUID = 1L;{
        add(SAMTag.NM.name());
        add(SAMTag.UQ.name());
        add(SAMTag.PG.name());
        add(SAMTag.MD.name());
        add(SAMTag.MQ.name());
        add(SAMTag.SA.name()); // Supplementary alignment metadata
        add(SAMTag.MC.name()); // Mate Cigar
        add(SAMTag.AS.name());
    }});

    public enum FileType implements CommandLineParser.ClpEnum {
        sam("Generate SAM files."),
        bam("Generate BAM files."),
        cram("Generate CRAM files."),
        dynamic("Generate files based on the extention of input.");

        final String description;

        FileType(String description) {
            this.description = description;
        }

        @Override
        public String getHelpDoc() {
            return description;
        }
    }

    /**
     * Enforce that output ordering is queryname when sanitization is turned on since it requires a queryname sort.
     * Also checks to ensure that the user has chosen a valid subset of arguments pertaining to output and sanitization.
     */
    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<>();
        validateOutputParams(outputByReadGroup, output, outputMap);

        if (!sanitize && keepFirstDuplicate) {
            errors.add("'keepFirstDuplicate' cannot be used without 'sanitize'");
        }

        if (!errors.isEmpty()) {
            return errors.toArray(new String[errors.size()]);
        }
        return null;
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        Broadcast<SAMFileHeader> headerBroadcast = ctx.broadcast(getHeaderForReads());
        JavaRDD<GATKRead> reads = getReads();

        ////////////////////////////////////////////////////////////////////////////
        // Grab the input header and remap values where appropriate
        ////////////////////////////////////////////////////////////////////////////
        SAMFileHeader localHeader = headerBroadcast.getValue();
        validateHeaderOverrides(localHeader, sampleAlias, libraryName);
        if (sampleAlias != null) {
            localHeader.getReadGroups().forEach(rg -> rg.setSample(sampleAlias));
        }
        if (libraryName != null) {
            localHeader.getReadGroups().forEach(rg -> rg.setLibrary(libraryName));
        }

        ////////////////////////////////////////////////////////////////////////////
        // Map the readgroups in the header to appropriate
        ////////////////////////////////////////////////////////////////////////////
        Map<String, Path> writerMap = getOutputMap(outputMap,
                                                  output,
                                                  getDefaultExtension(readArguments.getReadFiles().get(0).toString(), outputByReadgroupFileFormat),
                                                  localHeader.getReadGroups(),
                                                  outputByReadGroup);

        ////////////////////////////////////////////////////////////////////////////
        // Construct appropriate headers for the output files
        ////////////////////////////////////////////////////////////////////////////
        final Map<String, SAMFileHeader> headerMap = getReadGroupHeaderMap(localHeader, writerMap);

        // Revert the reads based on the given attributes
        List<String> attributesToRevert = removeDefaults ? DEFAULT_ATTRIBUTES_TO_CLEAR : new ArrayList<>();
        attributesToRevert.addAll(attributesToClear);
        JavaRDD<GATKRead> readsReverted = revertReads(reads, attributesToRevert);

        ////////////////////////////////////////////////////////////////////////////
        // Sanitize the reads, sorting them into appropriate order if necessary
        ////////////////////////////////////////////////////////////////////////////
        if (sanitize) {
            Map<String, FastqQualityFormat> readGroupFormatMap = createReadGroupFormatMap(readsReverted, headerBroadcast, !dontRestoreOriginalQualities);

            readsReverted = sanitize(readGroupFormatMap, readsReverted, localHeader, keepFirstDuplicate);
        }

        // Write the one or many read output files
        for (Map.Entry<String, Path> rmap: writerMap.entrySet()) {
            //TODO what to do if the readgroup isn't present
            final String key = rmap.getKey();
            JavaRDD<GATKRead> filteredreads = rmap.getKey()==null? readsReverted :
                                                                    readsReverted.filter(r -> r.getReadGroup().equals(key));
            writeReads(ctx, rmap.getValue().toString(), filteredreads, headerMap.get(rmap.getKey()), false); //TODO proper header map
        }
    }

    /**
     * Runs the QualityEncodingDetector over a sampling of each readGroup present in the file to detect what the encoding format
     * for base quality those reads are.
     *
     * @param reads Reads RDD over which to iterate and detect readgroups
     * @param inHeader Header describing the readgroups present in the bam
     * @param restoreOriginalQualities Whether to use the OQ tag for determining the map
     * @return the best guess at the quality encoding format present for each readgroup based on the first {@link QualityEncodingDetector#DEFAULT_MAX_RECORDS_TO_ITERATE} reads in each readgroup.
     */
    private Map<String, FastqQualityFormat> createReadGroupFormatMap( final JavaRDD<GATKRead> reads,
                                                                      final Broadcast<SAMFileHeader> inHeader,
                                                                      final boolean restoreOriginalQualities) {
        final Map<String, FastqQualityFormat> output = new HashMap<>();

        inHeader.getValue().getReadGroups().stream().forEach(rg -> {
            // For each readgroup filter down to just the reads in that group
            final String key = rg.getId();
            JavaRDD<GATKRead> filtered = reads.filter(r -> r.getReadGroup().equals(key));

            // NOTE: this method has the potential to be expensive as it may end up pulling on the first partition many times, and potentially
            //       end up iterating over the entire genome in the case where there are readgroups missing from the bam
            if (!filtered.isEmpty()) {

                // take the number of reads required by QualityEncodingDetector to determine quality score map
                CloseableIterator<SAMRecord> iterator = new CloseableIterator<SAMRecord>() {
                    Iterator<GATKRead> delegateIterator = filtered.take((int) QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE).iterator();

                    @Override
                    public void close() {
                        delegateIterator = null;
                    }

                    @Override
                    public boolean hasNext() {
                        return delegateIterator != null && delegateIterator.hasNext();
                    }

                    @Override
                    public SAMRecord next() {
                        if (!hasNext()) {
                            throw new NoSuchElementException("hasNext should be called before next");
                        }
                        return delegateIterator.next().convertToSAMRecord(inHeader.getValue());
                    }
                };

                // Save what the quality format for each readgroup was.
                output.put(rg.getId(), QualityEncodingDetector.detect(
                        QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE,
                        iterator,
                        restoreOriginalQualities));
            }
        });
        return output;
    }

    /**
     * If this is run, we want to be careful to remove copied reads from the bam.
     *
     * In order to do this we group each read by its readname and randomly select one read labeled as first in pair
     * and one read labeled as second in pair to treat as the representative reads, throwing away the rest.
     */
    private JavaRDD<GATKRead> sanitize(final Map<String, FastqQualityFormat> readGroupToFormat, final JavaRDD<GATKRead> reads, final SAMFileHeader header, final boolean keepFirstDuplicate) {
        JavaRDD<GATKRead> sortedReads = SparkUtils.querynameSortReadsIfNecessary(
                reads.filter(r -> r.getLength() == r.getBaseQualityCount()),
                getRecommendedNumReducers(), header);
        JavaPairRDD<String, Iterable<GATKRead>> readsByGroup = spanReadsByKey(sortedReads);

        return readsByGroup.flatMap(group -> {
            final List<GATKRead> out = Lists.newArrayList();

            List<GATKRead> primaryReads = Utils.stream(group._2()).collect(Collectors.toList());

            // Get the number of R1s, R2s, and unpaired reads respectively.
            int firsts = 0, seconds = 0, unpaired = 0;
            GATKRead firstRecord = null, secondRecord = null, unpairedRecord = null;
            for (final GATKRead rec : primaryReads) {
                if (!rec.isPaired()) {
                    if (unpairedRecord == null) {
                        unpairedRecord = rec;
                    }
                    ++unpaired;
                } else {
                    if (rec.isFirstOfPair()) {
                        if (firstRecord == null) {
                            firstRecord = rec;
                        }
                        ++firsts;
                    }
                    if (rec.isSecondOfPair()) {
                        if (secondRecord == null) {
                            secondRecord = rec;
                        }
                        ++seconds;
                    }
                }
            }

            // If we have paired reads, then check that there is exactly one first of pair and one second of pair.
            // Otherwise, check that we have only one unpaired read.
            if (firsts > 0 || seconds > 0) { // if we have any paired reads
                if (firsts != 1 || seconds != 1) { // if we do not have exactly one R1 and one R2
                    if (keepFirstDuplicate && firsts >= 1 && seconds >= 1) { // if we have at least one R1 and one R2, we can discard all but an arbitrary one
                        primaryReads = Arrays.asList(firstRecord, secondRecord);
                    }
                    // Otherwise don't admit anything from this group
                    else {
                        return out.iterator();
                    }
                }
            } else if (unpaired > 1) { // only unpaired reads, and we have too many
                if (keepFirstDuplicate) {
                    primaryReads = Collections.singletonList(unpairedRecord);
                }
                // Otherwise remove these reads entirely
                else {
                    return out.iterator();
                }
            }

            // If we've made it this far spit the records into the output!
            for (final GATKRead rec : primaryReads) {
                // The only valid quality score encoding scheme is standard; if it's not standard, change it.
                final FastqQualityFormat recordFormat = readGroupToFormat.get(rec.getReadGroup());
                if (recordFormat != null && !recordFormat.equals(FastqQualityFormat.Standard)) {
                    final byte[] quals = rec.getBaseQualities();
                    for (int i = 0; i < quals.length; i++) {
                        quals[i] -= SolexaQualityConverter.ILLUMINA_TO_PHRED_SUBTRAHEND;
                    }
                    rec.setBaseQualities(quals);
                }
                out.add(rec);
            }
            return out.iterator();
        });

    }

    /**
     * Method which takes an RDD of reads that is guaranteed to have every readname group placed together on the same
     * partition and maps those so a JavaPairRDD with the readname as the key.
     */
    private static JavaPairRDD<String, Iterable<GATKRead>> spanReadsByKey(final JavaRDD<GATKRead> reads) {
        JavaPairRDD<String, GATKRead> nameReadPairs = reads.mapToPair(read -> new Tuple2<>(read.getName(), read));
        return SparkUtils.spanByKey(nameReadPairs).flatMapToPair(namedRead -> {
            // for each name, separate reads by key (group name)
            List<Tuple2<String, Iterable<GATKRead>>> out = Lists.newArrayList();
            ListMultimap<String, GATKRead> multi = LinkedListMultimap.create();
            for (GATKRead read : namedRead._2()) {
                multi.put(read.getName(), read);
            }
            for (String key : multi.keySet()) {
                // list from Multimap is not serializable by Kryo, so put in a new array list
                out.add(new Tuple2<>(key, Lists.newArrayList(multi.get(key))));
            }
            return out.iterator();
        });
    }


    private Map<String, SAMFileHeader> getReadGroupHeaderMap(SAMFileHeader inHeader, Map<String, Path> writerMap) {
        final Map<String, SAMFileHeader> headerMap;
        if (outputByReadGroup) {
            if (inHeader.getReadGroups().isEmpty()) {
                throw new UserException("The header is missing its read group map");
            }

            assertAllReadGroupsMapped(writerMap, inHeader.getReadGroups());
            headerMap = new HashMap<>();
            for (final SAMReadGroupRecord readGroup : inHeader.getReadGroups()) {
                final SAMFileHeader header = createOutHeader(inHeader, sortOrder, !keepAlignmentInformation);
                header.addReadGroup(readGroup);
                headerMap.put(readGroup.getId(), header);
            }
        } else {
            final SAMFileHeader singleOutHeader = createOutHeader(inHeader, sortOrder, !keepAlignmentInformation);
            inHeader.getReadGroups().forEach(singleOutHeader::addReadGroup);
            headerMap = Collections.singletonMap(null, singleOutHeader);
        }
        return headerMap;
    }

    private static SAMFileHeader createOutHeader(
            final SAMFileHeader inHeader,
            final SAMFileHeader.SortOrder sortOrder,
            final boolean removeAlignmentInformation) {
        final SAMFileHeader outHeader = new SAMFileHeader();
        outHeader.setSortOrder(sortOrder);
        if (!removeAlignmentInformation) {
            outHeader.setSequenceDictionary(inHeader.getSequenceDictionary());
            outHeader.setProgramRecords(inHeader.getProgramRecords());
        }
        return outHeader;
    }

    @VisibleForTesting
    static String getDefaultExtension(final String input, final FileType setting) {
        if (setting == FileType.dynamic) {
            if (input.endsWith(FileExtensions.SAM)) {
                return FileExtensions.SAM;
            }
            if (input.endsWith(FileExtensions.CRAM)) {
                throw new UserException.UnimplementedFeature("Input file is a cram. This is currently unsupported for this tool");
            }
            return FileExtensions.BAM;
        } else {
            return "." + setting.toString();
        }
    }

    /**
     * Takes an individual SAMRecord and applies the set of changes/reversions to it that
     * have been requested by program level options.
     */
    public JavaRDD<GATKRead> revertReads(JavaRDD<GATKRead> reads, List<String> attributesToClear) {
        Broadcast<List<String>> attrBroadcast = JavaSparkContext.fromSparkContext(reads.context()).broadcast(attributesToClear);

        if (!dontRestoreOriginalQualities) {
            reads = reads.map(r -> {
                final byte[] oq = ReadUtils.getOriginalBaseQualities(r);
                if (oq != null) {
                    r.setBaseQualities(oq);
                    r.setAttribute("OQ", (String)null);
                }
                return r;
            });
        }

        if (!keepDuplicateInformation) {
            reads = reads.map(r -> {r.setIsDuplicate(false); return r;});
        }

        if (!keepAlignmentInformation) {
            reads = reads.map(rec -> {
                if (rec.isReverseStrand()) {
                    rec.reverseComplement();
                    rec.setIsReverseStrand(false);
                }

                // Remove all alignment based information about the read itself
                rec.setIsUnplaced();
                rec.setCigar(SAMRecord.NO_ALIGNMENT_CIGAR);

                rec.setFragmentLength(0);
                rec.setIsSecondaryAlignment(false);
                rec.setIsProperlyPaired(false);

                // Then remove any mate flags and info related to alignment
                rec.setMateIsUnplaced();

                // And then remove any tags that are calculated from the alignment
                attrBroadcast.getValue().forEach(tag -> rec.setAttribute(tag, (String) null));
                return rec;
            });
        }

        return reads;
    }

    @VisibleForTesting
    static Map<String, Path> getOutputMap(
            final String outputMapFile,
            final String outputDir,
            final String defaultExtension,
            final List<SAMReadGroupRecord> readGroups,
            final boolean outputByReadgroup) {

        if (outputByReadgroup) {
            final Map<String, Path> outputMap;
            if (outputMapFile != null) {
                try {
                    outputMap = createOutputMapFromFile(outputMapFile);
                } catch (IOException e) {
                    throw new UserException("Encountered an error reading output map file", e);
                }
            } else {
                outputMap = createOutputMapFromHeader(readGroups, outputDir, defaultExtension);
            }
            return outputMap;
        } else {
            return Collections.singletonMap(null, IOUtils.getPath(outputDir));
        }
    }

    // Names the files based on the locations laid out in the readgroup map
    private static Map<String, Path> createOutputMapFromFile(final String outputMapFile) throws IOException {
        final Map<String, Path> outputMap = new HashMap<>();
        try (final FeatureReader<TableFeature>  parser = AbstractFeatureReader.getFeatureReader(outputMapFile, new TableCodec(null), false);) {
            for (final TableFeature row : parser.iterator()) {
                final String id = row.get(OUTPUT_MAP_READ_GROUP_FIELD_NAME);
                final String output = row.get(OUTPUT_MAP_OUTPUT_FILE_FIELD_NAME);
                final Path outputPath = IOUtils.getPath(output);
                outputMap.put(id, outputPath);
            }
        }
        return outputMap;
    }

    // Names the files based on the readgroups individually presented in the header
    private static Map<String, Path> createOutputMapFromHeader(final List<SAMReadGroupRecord> readGroups, final String outputDir, final String extension) {
        final Map<String, Path> outputMap = new HashMap<>();
        for (final SAMReadGroupRecord readGroup : readGroups) {
            final String id = readGroup.getId();
            final String fileName = id + extension;
            final Path outputPath = Paths.get(outputDir, fileName);
            outputMap.put(id, outputPath);
        }
        return outputMap;
    }

    /**
     * Methods used for validating parameters to RevertSam.
     */
    static List<String> validateOutputParams(final boolean outputByReadGroup, final String output, final String outputMap) {
        final List<String> errors = new ArrayList<>();
        try {
            if (outputByReadGroup) {
                errors.addAll(validateOutputParamsByReadGroup(output, outputMap));
            } else {
                errors.addAll(validateOutputParamsNotByReadGroup(output, outputMap));
            }
        } catch (IOException e) {
            throw new UserException.BadInput("Error while validating input file", e);
        }
        return errors;
    }

    @SuppressWarnings("unchecked")
    static List<String>  validateOutputParamsByReadGroup(final String output, final String outputMap) throws IOException {
        final List<String> errors = new ArrayList<>();
        if (output != null) {
            if (!Files.isDirectory(IOUtil.getPath(output))) {
                errors.add("When '--output-by-readgroup' is set and output is provided, it must be a directory: " + output);
            }
            return errors;
        }
        // output is null if we reached here
        if (outputMap == null) {
            errors.add("Must provide either output or outputMap when '--output-by-readgroup' is set.");
            return errors;
        }
        if (!Files.isReadable(IOUtil.getPath(outputMap))) {
            errors.add("Cannot read outputMap " + outputMap);
            return errors;
        }
        final FeatureReader<TableFeature>  parser = AbstractFeatureReader.getFeatureReader(outputMap, new TableCodec(null),false);
        if (!isOutputMapHeaderValid((List<String>)parser.getHeader())) {
            errors.add("Invalid header: " + outputMap + ". Must be a tab-separated file with +"+OUTPUT_MAP_READ_GROUP_FIELD_NAME+"+ as first column and output as second column.");
        }
        return errors;
    }

    static List<String> validateOutputParamsNotByReadGroup(final String output, final String outputMap) throws IOException {
        final List<String> errors = new ArrayList<>();
        if (outputMap != null) {
            errors.add("Cannot provide outputMap when '--output-by-read' isn't set. Provide output instead.");
        }
        if (output == null) {
            errors.add("output is required when '--output-by-read'");
            return errors;
        }
        if (Files.isDirectory(IOUtil.getPath(output))) {
            errors.add("output " + output + " should not be a directory when '--output-by-read'");
        }
        return errors;
    }

    /**
     * If we are going to override sampleAlias or libraryName, make sure all the read
     * groups have the same values.
     */
    static void validateHeaderOverrides(
            final SAMFileHeader inHeader,
            final String sampleAlias,
            final String libraryName) {

        final List<SAMReadGroupRecord> rgs = inHeader.getReadGroups();
        if (sampleAlias != null || libraryName != null) {
            boolean allSampleAliasesIdentical = true;
            boolean allLibraryNamesIdentical = true;
            for (int i = 1; i < rgs.size(); i++) {
                if (!rgs.get(0).getSample().equals(rgs.get(i).getSample())) {
                    allSampleAliasesIdentical = false;
                }
                if (!rgs.get(0).getLibrary().equals(rgs.get(i).getLibrary())) {
                    allLibraryNamesIdentical = false;
                }
            }
            if (sampleAlias != null && !allSampleAliasesIdentical) {
                throw new UserException("Read groups have multiple values for sample.  " +
                        "A value for sampleAlias cannot be supplied.");
            }
            if (libraryName != null && !allLibraryNamesIdentical) {
                throw new UserException("Read groups have multiple values for library name.  " +
                        "A value for library name cannot be supplied.");
            }
        }
    }

    static void assertAllReadGroupsMapped(final Map<String, Path> outputMap, final List<SAMReadGroupRecord> readGroups) {
        for (final SAMReadGroupRecord readGroup : readGroups) {
            final String id = readGroup.getId();
            final Path output = outputMap.get(id);
            if (output == null) {
                throw new GATKException("Read group id " + id + " not found in outputMap " + outputMap);
            }
        }
    }

    static boolean isOutputMapHeaderValid(final List<String> columnLabels) {
        return columnLabels.size() >= 2 &&
                OUTPUT_MAP_READ_GROUP_FIELD_NAME.equals(columnLabels.get(0)) &&
                OUTPUT_MAP_OUTPUT_FILE_FIELD_NAME.equals(columnLabels.get(1));
    }

}
