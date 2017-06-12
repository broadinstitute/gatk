package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by davidben on 6/14/17.
 */
@CommandLineProgramProperties(
        summary = "Calculate pileup artifact statistics",
        oneLineSummary = "Calculate pileup artifact statistics",
        programGroup = VariantProgramGroup.class)
@BetaFeature
public class GetArtifactSummaries extends LocusWalker {

    public static final String BASES_BEFORE_AND_AFTER_LONG_NAME = "bases_before_and_after";
    public static final String BASES_BEFORE_AND_AFTER_SHORT_NAME = "bba";

    public static final String MIN_BASE_QUALITY_LONG_NAME = "min_base_quality";
    public static final String MIN_BASE_QUALITY_SHORT_NAME = "min_bq";

    private ArtifactSummary.ArtifactSummaryTableWriter writer;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table")
    private File outputTable;

    @Argument(fullName = BASES_BEFORE_AND_AFTER_LONG_NAME,
            shortName = BASES_BEFORE_AND_AFTER_SHORT_NAME,
            doc="The number of bases of context before and after the locus")
    private int basesBeforeAndAfter;

    @Argument(fullName = MIN_BASE_QUALITY_LONG_NAME,
            shortName = MIN_BASE_QUALITY_SHORT_NAME,
            doc="Minimum base quality to count substitutions")
    private byte minBaseQuality;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return makeReadFilters();
    }

    public static List<ReadFilter> makeReadFilters() {
        // The order in which we apply filters is important. Cheap filters come first so we fail fast
        List<ReadFilter> filters = new ArrayList<>();
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);
        filters.add(new MappingQualityReadFilter(60));
        filters.add(ReadFilterLibrary.PRIMARY_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);
        filters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE);
        filters.add(ReadFilterLibrary.GOOD_CIGAR);
        filters.add(new WellformedReadFilter());

        return filters;
    }

    @Override
    public void onTraversalStart() {
        try {
            writer = new ArtifactSummary.ArtifactSummaryTableWriter(outputTable);
        } catch (final IOException ex) {
            throw new UserException.BadInput("Could not open " + outputTable.getAbsolutePath() + " for writing.", ex);
        }
    }

    @Override
    public void apply(final AlignmentContext alignmentContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        if (alignmentContext.getStart() < referenceContext.getWindow().getStart() || alignmentContext.getEnd() > referenceContext.getWindow().getEnd()) {
            return;
        }
        final ArtifactSummary artifactSummary = new ArtifactSummary(alignmentContext, referenceContext, basesBeforeAndAfter, minBaseQuality);
        try {
            writer.writeRecord(artifactSummary);
        } catch (final IOException ex) {
            throw new GATKException("Problem outputting an artifact summary", ex);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        try {
            writer.close();
        } catch(final IOException ex) {
            throw new UserException.BadInput("Writing to " + outputTable.getAbsolutePath() + " failed.", ex);
        }
        return "SUCCESS";
    }
}
