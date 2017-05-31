package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;

import java.io.File;
import java.util.List;

/**
 * Groups contiguous targets with the same minor allele fraction for a single sample.
 *
 *  <p>
 *     The --tumorHets file is from {@link GetBayesianHetCoverage}, preferably, but can also come from {@link org.broadinstitute.hellbender.tools.exome.GetHetCoverage}.
 *      For example,
 * </p>
 *
 * <pre>
 *    CONTIG	POSITION	REF_COUNT	ALT_COUNT
 *     1	881918	14	21
 *     1	909238	13	11
 *     1	934940	14	0
 *     1	949608	20	14
 *     1	1021415	2	10
 * </pre>
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Segment genomic data into regions of constant minor allele fraction.  Only supports one sample input.",
        oneLineSummary = "(Experimental) Segment genomic data into regions of constant minor allele fraction",
        programGroup = CopyNumberProgramGroup.class
)
public final class PerformAlleleFractionSegmentation extends CommandLineProgram {
    protected static final String INITIAL_NUM_STATES_LONG_NAME = "initialNumberOfStates";
    protected static final String INITIAL_NUM_STATES_SHORT_NAME = "initialNumStates";

    @Argument(
            doc = "Input file for tumor-sample ref/alt read counts at normal-sample heterozygous-SNP sites (output of GetHetCoverage tool).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpCountsFile;

    @Argument(
            doc = "Initial number of hidden allele-fraction states",
            fullName = INITIAL_NUM_STATES_LONG_NAME,
            shortName = INITIAL_NUM_STATES_SHORT_NAME,
            optional = false
    )
    protected int initialNumStates;

    @Argument(
            doc = "Input file for allelic-bias panel of normals.",
            fullName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_SHORT_NAME,
            optional = true
    )
    protected File allelicPoNFile;

    @Argument(
            doc = "Output file for allele-fraction segments.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File outputSegmentsFile;

    @Override
    public Object doWork() {
        final String sampleName = FilenameUtils.getBaseName(snpCountsFile.getAbsolutePath());
        final AllelicPanelOfNormals allelicPoN =
                allelicPoNFile != null ? AllelicPanelOfNormals.read(allelicPoNFile) : AllelicPanelOfNormals.EMPTY_PON;
        final AllelicCountCollection acc = new AllelicCountCollection(snpCountsFile);
        final AlleleFractionSegmenter segmenter = new AlleleFractionSegmenter(initialNumStates, acc, allelicPoN);
        final List<ModeledSegment> segments = segmenter.getModeledSegments();

        SegmentUtils.writeModeledSegmentFile(outputSegmentsFile, segments, sampleName, true);

        return "SUCCESS";
    }
}
