package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.gvs.common.GQStateEnum;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.function.Consumer;
import java.util.stream.Stream;

public class ExtractCohortLiteEngine extends ExtractCohortEngine {
    final Logger logger;

    List<String> getFilterSetInfoTableFields() {
        return SchemaUtils.VQSLITE_YNG_FIELDS;
    }

    String getScoreFieldName() { return SchemaUtils.SCORE; }

    String getScoreKey() {
        return GATKVCFConstants.SCORE_KEY;
    }

    String getVQScoreFieldName() {
        return SchemaUtils.CALIBRATION_SENSITIVITY;
    }

    String getAlleleSpecificVQSScoreKey() {
        return GATKVCFConstants.CALIBRATION_SENSITIVITY_KEY;
    }

    String getVqScoreSNPFailureFilterName() {
        return GATKVCFConstants.CALIBRATION_SENSITIVITY_FAILURE_SNP;
    }

    String getVqScoreINDELFailureFilterName() {
        return GATKVCFConstants.CALIBRATION_SENSITIVITY_FAILURE_INDEL;
    }

    public ExtractCohortLiteEngine(final String projectID,
                                   final VCFHeader vcfHeader,
                                   final VariantAnnotatorEngine annotationEngine,
                                   final ReferenceDataSource refSource,
                                   final Map<Long, String> sampleIdToName,
                                   final String vetRangesFQDataSet,
                                   final String fqRangesExtractVetTable,
                                   final String fqRangesExtractRefTable,
                                   final GATKPath vetAvroFileName,
                                   final GATKPath refRangesAvroFileName,
                                   final List<SimpleInterval> traversalIntervals,
                                   final Long minLocation,
                                   final Long maxLocation,
                                   final String filterSetInfoTableName,
                                   final String filterSetSiteTableName,
                                   final int localSortMaxRecordsInRam,
                                   final boolean printDebugInformation,
                                   final Double vqScoreSNPThreshold,
                                   final Double vqScoreINDELThreshold,
                                   final String filterSetName,
                                   final boolean emitPLs,
                                   final boolean emitADs,
                                   final ExtractCohort.VQScoreFilteringType vqScoreFilteringType,
                                   final GQStateEnum inferredReferenceState,
                                   final boolean presortedAvroFiles,
                                   final Consumer<VariantContext> variantContextConsumer

    ) {
        super(projectID,
                vcfHeader,
                annotationEngine,
                refSource,
                sampleIdToName,
                vetRangesFQDataSet,
                fqRangesExtractVetTable,
                fqRangesExtractRefTable,
                vetAvroFileName,
                refRangesAvroFileName,
                traversalIntervals,
                minLocation,
                maxLocation,
                filterSetInfoTableName,
                filterSetSiteTableName,
                localSortMaxRecordsInRam,
                printDebugInformation,
                vqScoreSNPThreshold,
                vqScoreINDELThreshold,
                filterSetName,
                emitPLs,
                emitADs,
                vqScoreFilteringType,
                inferredReferenceState,
                presortedAvroFiles,
                variantContextConsumer);
        logger = LogManager.getLogger(ExtractCohortLiteEngine.class);
    }

    @Override
    boolean isFailingSite(final Stream<Double> vqScores, final Double vqScoreThreshold) {
        Optional<Double> minVal = vqScores
                .filter(d -> !(d.isNaN() || d.isInfinite()))
                .min(Double::compareTo);
        return minVal.isPresent() && minVal.get() > vqScoreThreshold;
    }


    @Override
    boolean isFailingGenotype(final Stream<Allele> nonRefAlleles,
                              final Map<Allele, Double> remappedVQScoreMap,
                              final Double vqScoreThreshold) {
        // get the minimum (best) calibration sensitivity for all non-Yay sites, and apply the filter
        Optional<Double> minVal =
                nonRefAlleles
                        .map(remappedVQScoreMap::get)
                        .filter(Objects::nonNull)
                        .min(Double::compareTo);

        return minVal.isPresent() && minVal.get() > vqScoreThreshold;
    }
}
