package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;

import java.io.File;
import java.util.List;

/**
 * Created by davidben on 5/23/16.
 */
@CommandLineProgramProperties(
        summary = "Segment genomic data into regions of constant minor allele fraction.  Only supports one sample input.",
        oneLineSummary = "Segment genomic data into regions of constant minor allele fraction",
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
        final List<ModeledSegment> segments = segmenter.findSegments();



        /* This code is for output in ACNV format so that it can be used as input to our ACNV plotting code
            during development.  When allele fraction and copy ratio HMM segmentation are combined, they will output
            a genuine ACNV modeled segment file instead of this hack.
        // for the sake of putting in a plotting script, give dummy values to ACNV stuff.  For example, give imaginary error
        // bars and a flat copy ratio
        final List<ACNVModeledSegment> outputSegments = segments.stream().map(segment -> {
            final double maf = segment.getSegmentMean();
            final PosteriorSummary copyRatioDummyPosteriorSummary = new PosteriorSummary(0, -0.1, 0.1);
            final PosteriorSummary mafPosteriorSummaryWithDummyErrorBars = new PosteriorSummary(maf, maf - 0.05, maf + 0.05);
            copyRatioDummyPosteriorSummary.setDeciles(new DecileCollection(Arrays.asList(0.25)));
            mafPosteriorSummaryWithDummyErrorBars.setDeciles(new DecileCollection(Arrays.asList(0.25)));
            return new ACNVModeledSegment(segment.getSimpleInterval(), copyRatioDummyPosteriorSummary, mafPosteriorSummaryWithDummyErrorBars);}
        ).collect(Collectors.toList());


        final ReadCountCollection dummyRcc = new ReadCountCollection(Arrays.asList(new Target("target", new SimpleInterval("1", 4, 5))),
                Arrays.asList(sampleName), new Array2DRowRealMatrix(1,1));
        final Genome genome = new Genome(dummyRcc, acc.getCounts(), sampleName);
        SegmentUtils.writeACNVModeledSegmentFile(outputSegmentsFile, outputSegments, genome);
        */

        SegmentUtils.writeModeledSegmentFile(outputSegmentsFile, segments, sampleName, true);

        return "SUCCESS";
    }
}
