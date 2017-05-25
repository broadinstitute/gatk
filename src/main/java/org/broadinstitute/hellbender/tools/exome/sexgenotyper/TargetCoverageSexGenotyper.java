package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableReader;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import javax.annotation.Nonnull;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Infers sex genotypes from raw target coverage data and germline contig ploidy annotations.
 *
 * <p>
 *     In addition to the raw target coverage file, the user must provide a tab-separated contig
 *     germline ploidy annotation file. The format is described in {@link ContigGermlinePloidyAnnotationTableReader}.
 *     Note: Every contig that appears in the read counts table must be annotated.
 * </p>
 *
 * <p>
 *     The provided target coverage file must contain both AUTOSOMAL and ALLOSOMAL targets. The tool
 *     infers the read depth density from AUTOSOMAL targets and uses this information to calculate the
 *     likelihood of sex genotypes.
 * </p>
 *
 * <p>
 *     Note: due to the uncertainty in the alignment of short reads, a small number of reads may be
 *     erroneously mapped to contigs with 0 actual ploidy. For example, female Homo sapiens samples may
 *     have a small number of reads aligned to the Y contig. The user must specify the typical mapping
 *     error probability for the tool to properly account for these errors.
 * </p>
 *
 * <p>
 *     Note: setting {@link TargetCoverageSexGenotyper#baselineMappingErrorProbability} to 0 (or to
 *     an unreasonably small number) will bias genotyping toward classes that cover more
 *     targets, e.g. XY over XX.
 * </p>
 *
 * <h2>Output File Format</h2>
 * <p>
 *     The inferred genotypes will be written to a tab-separated file such as:
 *
 *     <pre>
 *         SAMPLE_NAME      SEX_GENOTYPE      SEX_XX             SEX_XY
 *         sample_1       SEX_XX            [log likelihood]   [log likelihood]
 *         sample_2       SEX_XY            [log likelihood]   [log likelihood]
 *                                             ...
 *     </pre>
 * </p>
 *
 * <p>
 *     To sex genotype, remove targets that overlap pseudoautosomal regions (PAR) of allosteric contigs, e.g. human X and Y.
 *     For example, if the read counts include PAR targets, then omit these from the analysis by excluding them from the --targets file.
 *     If these regions are included, then sex genotypes can be unreliable.
 * </p>
 *
 * <p>
 *     Allosteric contig mapping errors, e.g. human X contig reads that map to Y,
 *     as well as partial blacklisting of PAR regions, can bias inferred sexes towards genotypes that cover more targets such as XY.
 *     Inspect results for such systematic biases, e.g. with a histogram of log likelihoods.
 *     This recommendation applies especially to whole genome sequencing (WGS) samples.
 * </p>
 * <h3>Example</h3>
 *
 * <pre>
 * java -Xmx4g -jar $gatk_jar TargetCoverageSexGenotyper \
 *   --input combined_read_counts.tsv \
 *   --contigAnnotations grch37_contig_annotations.tsv \
 *   --targets targets.tsv \
 *   --output inferred_sex_genotypes.tsv
 * </pre>
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "This tool infers sample sex genotypes from raw target read counts by calculating the likelihoods of all provided " +
                "genotypes and choosing the most likely one. The required inputs are (1) a table of raw target read counts " +
                "from one or more samples, and (2) a table of annotated contigs that includes a CONTIG column, " +
                "a CLASS column (AUTOSOMAL, or ALLOSOMAL), and one additional column for each sex genotype " +
                "that lists the expected germline ploidy of the contig. Sex genotypes may have arbitrary names and are " +
                "identified with their given column name. All contigs that appear in the read count table " +
                "must be annotated in this file. The output is a tab-separated file that includes sample names, their " +
                "inferred sex genotypes, and the log likelihood of each sex genotype.",
        oneLineSummary = "Infers sample sex genotypes from raw read counts.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public class TargetCoverageSexGenotyper extends CommandLineProgram {

    private final Logger logger = LogManager.getLogger(TargetCoverageSexGenotyper.class);

    public static final String INPUT_CONTIG_ANNOTATIONS_LONG_NAME = "contigAnnotations";
    public static final String INPUT_CONTIG_ANNOTATIONS_SHORT_NAME = "annots";

    public static final String BASELINE_MAPPING_ERROR_PROBABILITY_LONG_NAME = "baselineMappingError";
    public static final String BASELINE_MAPPING_ERROR_PROBABILITY_SHORT_NAME = "mapErr";

    @Argument(
            doc = "Input raw read count collection tab-separated file for either cohort or single sample.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            optional = false
    )
    protected File inputRawReadCountsFile;

    @Argument(
            doc = "Output sample sex genotype tab-separated file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false
    )
    protected File outputSampleGenotypesFile;

    @Argument(
            doc = "Input germline contig ploidy annotations file. For an example file, see the GATK Resource Bundle.",
            fullName = INPUT_CONTIG_ANNOTATIONS_LONG_NAME,
            shortName = INPUT_CONTIG_ANNOTATIONS_SHORT_NAME,
            optional = false
    )
    protected File inputContigAnnotsFile;

    @Argument(
            doc = "Baseline mapping error probability or an estimate of the typical mapping error.",
            fullName = BASELINE_MAPPING_ERROR_PROBABILITY_LONG_NAME,
            shortName = BASELINE_MAPPING_ERROR_PROBABILITY_SHORT_NAME,
            optional = true
    )
    protected double baselineMappingErrorProbability = 1e-4;

    @Argument(
            doc = "Input targets file. A tab-separated (.tsv) list of targets (not a BED file) with at least" +
                    " the following header columns: contig, start, stop, name. It must be a subset or can be all of the" +
                    " targets included in the read counts table.",
            fullName = ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME,
            optional = true
    )
    protected File inputTargetListFile = null;

    @ArgumentCollection
    protected OptionalIntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @Override
    protected Object doWork() {
        /* check args */
        IOUtils.canReadFile(inputRawReadCountsFile);
        IOUtils.canReadFile(inputContigAnnotsFile);
        Utils.validateArg(baselineMappingErrorProbability > 0 && baselineMappingErrorProbability < 1,
                "Must have 0 < baseline mapping error probability < 1.");

        /* read input read counts */
        final ReadCountCollection rawReadCounts;
        try {
            logger.info("Parsing raw read count collection file...");
            rawReadCounts = ReadCountCollectionUtils.parse(inputRawReadCountsFile);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the raw read count collection file");
        }

        /* parse contig genotype ploidy annotations */
        final List<ContigGermlinePloidyAnnotation> contigGermlinePloidyAnnotationList;
        logger.info("Parsing contig genotype ploidy annotations file...");
        contigGermlinePloidyAnnotationList = ContigGermlinePloidyAnnotationTableReader.readContigGermlinePloidyAnnotationsFromFile(inputContigAnnotsFile);

        /* parse target list */
        final List<Target> inputTargetList;
        if (inputTargetListFile != null) {
            inputTargetList = TargetTableReader.readTargetFile(inputTargetListFile);
        } else {
            logger.info("No target list was provided; using all targets present in the read count collection");
            inputTargetList = rawReadCounts.targets();
        }
        final List<SimpleInterval> genotypingIntervals;
        if (intervalArgumentCollection.intervalsSpecified()) {
            final SAMSequenceDictionary sequenceDictionary = generateMaximallyCoveringSAMSequenceDictionary(inputTargetList);
            try {
                genotypingIntervals = intervalArgumentCollection.getIntervals(sequenceDictionary);
            } catch (final UserException.MalformedGenomeLoc ex) {
                throw new UserException.BadInput("Genotyping intervals could not be resolved properly; make sure that the" +
                        " target list contains the specified intervals/contigs", ex);
            }
            logger.info("Some intervals specified; restricting genotyping to the following intervals: " + genotypingIntervals.stream()
                    .map(SimpleInterval::toString)
                    .collect(Collectors.joining(", ")));
        } else {
            logger.info("No intervals specified; considering all targets for genotyping");
            genotypingIntervals = null;
        }

        /* perform genotyping */
        final TargetCoverageSexGenotypeCalculator genotyper = new TargetCoverageSexGenotypeCalculator(rawReadCounts,
                inputTargetList, contigGermlinePloidyAnnotationList, genotypingIntervals, baselineMappingErrorProbability);
        final SexGenotypeDataCollection sampleSexGenotypeCollection = genotyper.inferSexGenotypes();

        /* save results */
        try {
            sampleSexGenotypeCollection.write(new FileWriter(outputSampleGenotypesFile));
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile("Could not write inferred sample genotypes to file", ex);
        }

        return "SUCCESS";
    }

    /**
     * Generates a maximally-covering sequence dictionary from a collection of {@link Locatable} objects. This is done by
     * extracting distinct contigs and mapping each to a {@link SAMSequenceRecord} with length {@link Integer#MAX_VALUE}.
     *
     * @param intervals a list of {@link Locatable} objects
     * @return a maximally-covering {@link SAMSequenceDictionary}
     */
    private SAMSequenceDictionary generateMaximallyCoveringSAMSequenceDictionary(@Nonnull final Collection<? extends Locatable> intervals) {
        return new SAMSequenceDictionary(intervals.stream()
                .map(Locatable::getContig)
                .distinct()
                .map(contig -> new SAMSequenceRecord(contig, Integer.MAX_VALUE))
                .collect(Collectors.toList()));
    }
}
