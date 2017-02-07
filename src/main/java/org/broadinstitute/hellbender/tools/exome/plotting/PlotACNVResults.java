package org.broadinstitute.hellbender.tools.exome.plotting;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Create plots of denoised and segmented copy-ratio and minor-allele-fraction estimates.",
        oneLineSummary = "Create plots of denoised and segmented copy-ratio and minor-allele-fraction estimates.",
        programGroup = CopyNumberProgramGroup.class
)
public final class PlotACNVResults extends CommandLineProgram {
    private static final String CNV_PLOTTING_R_LIBRARY = "CNVPlottingLibrary.R";
    private static final String ACNV_PLOTTING_R_SCRIPT = "ACNVResultsPlotting.R";

    private static final String CONTIG_DELIMITER = "CONTIG_DELIMITER";  //used to delimit contig names and lengths passed to the R script
    private static final int DEFAULT_MINIMUM_CONTIG_LENGTH = 1000000;   //can be used to filter out mitochondrial contigs, unlocalized contigs, etc.

    //CLI arguments
    protected static final String OUTPUT_PREFIX_LONG_NAME = "outputPrefix";
    protected static final String OUTPUT_PREFIX_SHORT_NAME = "pre";

    protected static final String MINIMUM_CONTIG_LENGTH_LONG_NAME = "minimumContigLength";
    protected static final String MINIMUM_CONTIG_LENGTH_SHORT_NAME = "minContigLength";

    @Argument(
            doc = "File containing het SNP positions, ref counts, and alt counts, produced by GetHetCoverage/GetBayesianHetCoverage.",
            fullName =  ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpCountsFile;

    @Argument(
            doc = "File containing tangent-normalized coverage of targets, produced by NormalizeSomaticReadCounts.",
            fullName =  ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File tangentFile;

    @Argument(
            doc = "File containing segmented regions of the genome, produced by AllelicCNV.",
            fullName =  ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File segmentsFile;

    @Argument(
            doc = "File containing the reference sequence dictionary (used to determine relative contig lengths). " +
                    "Contigs will be plotted in the order given. " +
                    "Contig names should not include \"" + CONTIG_DELIMITER + "\"." +
                    "Data on contigs not in the dictionary will not be plotted and a warning will be thrown.",
            shortName = StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME,
            optional = false
    )
    protected File sequenceDictionaryFile;

    @Argument(
            doc = "Threshold length (in bp) for contigs to be plotted. " +
                    "Contigs with lengths less than this threshold will not be plotted. " +
                    "This can be used to filter out mitochondrial contigs, unlocalized contigs, etc.",
            fullName =  MINIMUM_CONTIG_LENGTH_LONG_NAME,
            shortName = MINIMUM_CONTIG_LENGTH_SHORT_NAME,
            optional = false
    )
    protected int minContigLength = DEFAULT_MINIMUM_CONTIG_LENGTH;

    @Argument(
            doc = "Prefix for output filenames.",
            fullName =  OUTPUT_PREFIX_LONG_NAME,
            shortName = OUTPUT_PREFIX_SHORT_NAME,
            optional = false
    )
    protected String outputPrefix;

    @Argument(
            doc = "Output directory.",
            fullName =  StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false
    )
    protected String outputDir;

    @Override
    protected Object doWork() {
        checkRegularReadableUserFiles();

        //get sample name from input files (consistency check is performed)
        final String sampleName = getSampleName();

        //load contig names and lengths from the sequence dictionary into a LinkedHashMap
        final SAMSequenceDictionary sequenceDictionary = ReferenceUtils.loadFastaDictionary(sequenceDictionaryFile);
        Utils.validateArg(sequenceDictionary.getSequences().stream().map(SAMSequenceRecord::getSequenceName).noneMatch(n -> n.contains(CONTIG_DELIMITER)),
                String.format("Contig names cannot contain \"%s\".", CONTIG_DELIMITER));
        final Map<String, Integer> contigLengthMap = sequenceDictionary.getSequences().stream()
                .filter(s -> s.getSequenceLength() >= minContigLength)
                .collect(Collectors.toMap(SAMSequenceRecord::getSequenceName, SAMSequenceRecord::getSequenceLength,
                        (c, l) -> {
                            throw new IllegalArgumentException(String.format("Duplicate contig in sequence dictionary: %s", c));
                        },
                        LinkedHashMap::new));
        Utils.validateArg(contigLengthMap.size() > 0,
                "There must be at least one contig above the threshold length in the sequence dictionary.");
        logger.info("Contigs above length threshold: " + contigLengthMap.toString());

        //check that contigs in input files are present in sequence dictionary and that data points are valid given lengths
        validateContigs(contigLengthMap);

        //generate the plots
        final List<String> contigNames = new ArrayList<>(contigLengthMap.keySet());
        final List<Integer> contigLengths = new ArrayList<>(contigLengthMap.values());
        writeSegmentedAlleleFractionPlot(sampleName, contigNames, contigLengths);
        return "SUCCESS";
    }

    private void checkRegularReadableUserFiles() {
        Utils.regularReadableUserFile(snpCountsFile);
        Utils.regularReadableUserFile(tangentFile);
        Utils.regularReadableUserFile(segmentsFile);
        Utils.regularReadableUserFile(sequenceDictionaryFile);
        if (!new File(outputDir).exists()) {
            throw new UserException(String.format("Output directory %s does not exist.", outputDir));
        }
    }

    private String getSampleName() {
        //TODO add check of snpCountsFile once sample name has been moved to headers (and appropriate test)
        final String tangentSampleName = ReadCountCollectionUtils.getSampleNameForCLIsFromReadCountsFile(tangentFile);
        final String segmentsSampleName = SegmentUtils.getSampleNameForCLIsFromSegmentFile(segmentsFile);
        Utils.validateArg(new HashSet<>(Arrays.asList(tangentSampleName, segmentsSampleName)).size() == 1,
                "Sample names in input files do not all match.");
        return tangentSampleName;
    }

    private void validateContigs(final Map<String, Integer> contigLengthMap) {
        final Set<String> contigNames = contigLengthMap.keySet();

        //validate contig names and lengths in SNP counts file
        final AllelicCountCollection snpCounts = new AllelicCountCollection(snpCountsFile);
        final Set<String> snpCountsContigNames = snpCounts.getCounts().stream().map(AllelicCount::getContig).collect(Collectors.toSet());
        if (!contigNames.containsAll(snpCountsContigNames)) {
            logger.warn("Contigs present in the SNP counts file are missing from the sequence dictionary and will not be plotted.");
        }
        final Map<String, Integer> snpCountsContigMaxPositionMap = snpCounts.getCounts().stream().filter(c -> contigNames.contains(c.getContig()))
                .collect(Collectors.toMap(AllelicCount::getContig, AllelicCount::getEnd, Integer::max));
        snpCountsContigMaxPositionMap.keySet().forEach(c -> Utils.validateArg(snpCountsContigMaxPositionMap.get(c) <= contigLengthMap.get(c),
                "Position present in the SNP-counts file exceeds contig length in the sequence dictionary."));

        //validate contig names and lengths in tangent file
        final ReadCountCollection tangent;
        try {
            tangent = ReadCountCollectionUtils.parse(tangentFile);
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(tangentFile, e);
        }
        final Set<String> tangentContigNames = tangent.targets().stream().map(Target::getContig).collect(Collectors.toSet());
        if (!contigNames.containsAll(tangentContigNames)) {
            logger.warn("Contigs present in the tangent-normalized coverage file are missing from the sequence dictionary and will not be plotted.");
        }
        final Map<String, Integer> tangentContigMaxPositionMap = tangent.targets().stream().filter(t -> contigNames.contains(t.getContig()))
                .collect(Collectors.toMap(Target::getContig, Target::getEnd, Integer::max));
        tangentContigMaxPositionMap.keySet().forEach(c -> Utils.validateArg(tangentContigMaxPositionMap.get(c) <= contigLengthMap.get(c),
                "Position present in the tangent-normalized coverage file exceeds contig length in the sequence dictionary."));

        //validate contig names and lengths in segments file
        final List<ACNVModeledSegment> segments = SegmentUtils.readACNVModeledSegmentFile(segmentsFile);
        final Set<String> segmentsContigNames = segments.stream().map(ACNVModeledSegment::getContig).collect(Collectors.toSet());
        if (!contigNames.containsAll(segmentsContigNames)) {
            logger.warn("Contigs present in the segments file are missing from the sequence dictionary and will not be plotted.");
        }
        final Map<String, Integer> segmentsContigMaxPositionMap = segments.stream().filter(s -> contigNames.contains(s.getContig()))
                .collect(Collectors.toMap(ACNVModeledSegment::getContig, ACNVModeledSegment::getEnd, Integer::max));
        segmentsContigMaxPositionMap.keySet().forEach(c -> Utils.validateArg(segmentsContigMaxPositionMap.get(c) <= contigLengthMap.get(c),
                "Position present in the segments file exceeds contig length in the sequence dictionary."));
    }

    /**
     * @param sampleName Sample name derived from input files
     * @param contigNames List containing contig names
     * @param contigLengths List containing contig lengths (same order as contigNames)
     */
    private void writeSegmentedAlleleFractionPlot(final String sampleName,
                                                  final List<String> contigNames,
                                                  final List<Integer> contigLengths) {
        final String contigNamesArg = contigNames.stream().collect(Collectors.joining(CONTIG_DELIMITER));                            //names separated by delimiter
        final String contigLengthsArg = contigLengths.stream().map(Object::toString).collect(Collectors.joining(CONTIG_DELIMITER));  //lengths separated by delimiter
        final String outputDirArg = addTrailingSlashIfNecessary(outputDir);
        final RScriptExecutor executor = new RScriptExecutor();

        //this runs the R statement "source("CNVPlottingLibrary.R")" before the main script runs
        executor.addScript(new Resource(CNV_PLOTTING_R_LIBRARY, PlotSegmentedCopyRatio.class));
        executor.addScript(new Resource(ACNV_PLOTTING_R_SCRIPT, PlotACNVResults.class));
        //--args is needed for Rscript to recognize other arguments properly
        executor.addArgs("--args",
                "--sample_name=" + sampleName,
                "--snp_counts_file=" + snpCountsFile,
                "--tangent_file=" + tangentFile,
                "--segments_file=" + segmentsFile,
                "--contig_names=" + contigNamesArg,
                "--contig_lengths=" + contigLengthsArg,
                "--output_dir=" + outputDirArg,
                "--output_prefix=" + outputPrefix);
        executor.exec();
    }

    private static String addTrailingSlashIfNecessary(final String outputDir) {
        return outputDir.endsWith(File.separator) ? outputDir : outputDir + File.separator;
    }
}
