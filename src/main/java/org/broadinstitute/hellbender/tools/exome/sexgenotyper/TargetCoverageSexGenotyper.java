package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 * A command line tool for inferring sex genotypes from a tab-separated raw target coverage file.
 *
 * <p>
 *     In addition to the raw target coverage file, the user must provide a tab-separated contig
 *     germline ploidy annotation file. For example, a basic annotation file for homo sapiens is as follows:
 *
 *     <pre>
 *         CONTIG_NAME    CONTIG_CLASS    SEX_XX    SEX_XY
 *         1              AUTOSOMAL       2          2
 *         2              AUTOSOMAL       2          2
 *                                   ...
 *                                   ...
 *                                   ...
 *         X              ALLOSOMAL       2          0
 *         Y              ALLOSOMAL       1          1
 *     </pre>
 *
 *     CONTIG_NAME column values must be same as contig names that appear in the read counts table.
 *     CONTIG_CLASS is either AUTOSOMAL or ALLOSOMAL. "SEX_XX" and "SEX_YY" are arbitrary sex genotype tags,
 *     along with their ploidies (= number of homologs). AUTOSOMAL contigs must have the same ploidy for all
 *     sexes. Every contig that appears in the read counts table must be annotated. One may include additional
 *     sex genotypes (along with contig ploidies) for specifies having more than two sexes by adding
 *     additional columns to the contig annotation file.
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
 *     erroneously mapped to contigs with 0 actual ploidy (e.g. female homo sapiens samples may
 *     have a small number of reads aligned to the Y contig). The user must specifiy the typical mapping
 *     error probability for the tool to properly account for these errors.
 * </p>
 *
 * <p>
 *     Note: setting {@link TargetCoverageSexGenotyper#baselineMappingErrorProbability} to 0 (or to
 *     an unreasonably small number) will bias genotyping toward classes that cover more
 *     targets (SEX_XX < SEX_XY).
 * </p>
 *
 * <h2>Output File Format</h2>
 * <p>
 *     The inferred genotypes will be written to a tab-separated file such as:
 *
 *     <pre>
 *         SAMPLE_NAME                SEX_GENOTYPE      SEX_XX             SEX_XY
 *         arbitrary_XX_sample_name   SEX_XX            [log likelihood]   [log likelihood]
 *         arbitrary_XY_sample_name   SEX_XY            [log likelihood]   [log likelihood]
 *                                             ...
 *     </pre>
 * </p>
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "This tool infers sample sex genotypes from raw target read counts by calculating the likelihoods of all provided " +
                "genotypes and choosing the most likely one. The required inputs are (1) a table of raw target read counts " +
                "from one or more samples, and (2) a table of annotated contigs that includes a CONTIG_NAME column, " +
                "a CONTIG_CLASS column (AUTOSOMAL, or ALLOSOMAL), and one additional column for each sex genotype " +
                "that lists the expected germline ploidy of the contig. Sex genotypes may have arbitrary names and are " +
                "identified with their given column name. All contigs that appear in the input read count table targets " +
                "must be annotated in this file. The output is a tab-separated file that includes sample names, their " +
                "inferred sex genotypes, and the log likelihood of each sex genotype.",
        oneLineSummary = "Infers sample sex genotypes from raw read counts.",
        programGroup = CopyNumberProgramGroup.class
)

public class TargetCoverageSexGenotyper extends CommandLineProgram {

    private final Logger logger = LogManager.getLogger(TargetCoverageSexGenotyper.class);

    public static final String INPUT_READ_COUNT_COLLECTION_LONG_NAME = StandardArgumentDefinitions.INPUT_LONG_NAME;
    public static final String INPUT_READ_COUNT_COLLECTION_SHORT_NAME = StandardArgumentDefinitions.INPUT_SHORT_NAME;

    public static final String OUTPUT_SEX_GENOTYPE_LONG_NAME = StandardArgumentDefinitions.OUTPUT_LONG_NAME;
    public static final String OUTPUT_SEX_GENOTYPE_SHORT_NAME = StandardArgumentDefinitions.OUTPUT_SHORT_NAME;

    public static final String INPUT_CONTIG_ANNOTS_LONG_NAME = "contigAnnotations";
    public static final String INPUT_CONTIG_ANNOTS_SHORT_NAME = "annots";

    public static final String BASELINE_MAPPING_ERROR_PROBABILITY_LONG_NAME = "baselineMappingError";
    public static final String BASELINE_MAPPING_ERROR_PROBABILITY_SHORT_NAME = "mapErr";

    @Argument(
            doc = "Input raw read count collection tab-separated file.",
            fullName = INPUT_READ_COUNT_COLLECTION_LONG_NAME,
            shortName = INPUT_READ_COUNT_COLLECTION_SHORT_NAME,
            optional = false
    )
    protected File inputRawReadCountsFile;

    @Argument(
            doc = "Output sample sex genotype tab-separated file.",
            fullName = OUTPUT_SEX_GENOTYPE_LONG_NAME,
            shortName = OUTPUT_SEX_GENOTYPE_SHORT_NAME,
            optional = false
    )
    protected File outputSampleGenotypesFile;

    @Argument(
            doc = "Input contig annotations file.",
            fullName = INPUT_CONTIG_ANNOTS_LONG_NAME,
            shortName = INPUT_CONTIG_ANNOTS_SHORT_NAME,
            optional = false
    )
    protected File inputContigAnnotsFile;

    @Argument(
            doc = "Baseline mapping error probability.",
            fullName = BASELINE_MAPPING_ERROR_PROBABILITY_LONG_NAME,
            shortName = BASELINE_MAPPING_ERROR_PROBABILITY_SHORT_NAME,
            optional = true
    )
    protected double baselineMappingErrorProbability = 1e-4;

    @Override
    protected Object doWork() {
        /* check args */
        Utils.regularReadableUserFile(inputRawReadCountsFile);
        Utils.regularReadableUserFile(inputContigAnnotsFile);
        if (baselineMappingErrorProbability <= 0 || baselineMappingErrorProbability >= 1) {
            throw new IllegalArgumentException("Baseline mapping error probability must be a positive real number in " +
                    "range (0, 1) excluding endpoints.");
        }

        /* read input data */
        ReadCountCollection rawReadCounts;
        List<ContigGermlinePloidyAnnotation> contigGermlinePloidyAnnotationList;
        try {
            logger.info("Parsing raw read count collection file...");
            rawReadCounts = ReadCountCollectionUtils.parse(inputRawReadCountsFile);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read raw read count collection file");
        }
        try {
            logger.info("Parsing contig genotype ploidy annotations file...");
            contigGermlinePloidyAnnotationList = ContigGermlinePloidyAnnotationTableReader.readContigGermlinePloidyAnnotationsFromFile(inputContigAnnotsFile);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read contig genotype ploidy annotations file");
        }

        /* perform genotyping */
        final TargetCoverageSexGenotypeCalculator genotyper = new TargetCoverageSexGenotypeCalculator(rawReadCounts,
                contigGermlinePloidyAnnotationList, baselineMappingErrorProbability);
        final SexGenotypeDataCollection sampleSexGenotypeCollection = genotyper.inferSexGenotypes();

        /* save results */
        try {
            sampleSexGenotypeCollection.write(new FileWriter(outputSampleGenotypesFile));
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile("Could not write inferred sample genotypes to file", ex);
        }

        return "SUCCESS";
    }
}
