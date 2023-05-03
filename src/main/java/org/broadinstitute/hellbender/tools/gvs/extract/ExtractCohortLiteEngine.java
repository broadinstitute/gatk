package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.gvs.common.GQStateEnum;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

public class ExtractCohortLiteEngine extends ExtractCohortEngine {
    final Logger logger;

    List<String> getFilterSetInfoTableFields() {
        return SchemaUtils.VQSLITE_YNG_FIELDS;
    }
    String getVQScoreFieldName() {
        return SchemaUtils.CALIBRATION_SENSITIVITY;
    }

    String getAlleleSpecificVQSScoreKey() {
        return GATKVCFConstants.AS_VQS_SENS_KEY;
    }

    String getVqScoreSNPFailureFilterName() {
        return GATKVCFConstants.VQS_SENS_FAILURE_SNP;
    }

    String getVqScoreINDELFailureFilterName() {
        return GATKVCFConstants.VQS_SENS_FAILURE_INDEL;
    }

    public ExtractCohortLiteEngine(final String projectID,
                                   final VariantContextWriter vcfWriter,
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
                                   final ProgressMeter progressMeter,
                                   final String filterSetName,
                                   final boolean emitPLs,
                                   final boolean emitADs,
                                   final ExtractCohort.VQScoreFilteringType vqScoreFilteringType,
                                   final boolean excludeFilteredSites,
                                   final GQStateEnum inferredReferenceState,
                                   final boolean presortedAvroFiles
    ) {
        super(projectID,
                vcfWriter,
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
                progressMeter,
                filterSetName,
                emitPLs,
                emitADs,
                vqScoreFilteringType,
                excludeFilteredSites,
                inferredReferenceState,
                presortedAvroFiles);
        logger = LogManager.getLogger(ExtractCohortLiteEngine.class);
    }


    boolean isFailingSite(final List<Double> vqScoreList, final Double vqScoreThreshold) {
        Optional<Double> maxVal = vqScoreList.stream()
                .filter(d -> !(d.isNaN()||d.isInfinite()))
                .max(Double::compareTo);
        return maxVal.isPresent() && maxVal.get() > vqScoreThreshold;
    }



    boolean isFailingGenotype(final List<Allele> nonRefAlleles,
                              final LinkedHashMap<Allele, Double> remappedVQScoreMap,
                              final Double vqScoreThreshold) {
        // get the max (best) vq score (vqslod/sensitivity) for all non-Yay sites, and apply the filter
        Optional<Double> snpMax =
                nonRefAlleles.stream()
                        .map(remappedVQScoreMap::get)
                        .filter(Objects::nonNull)
                        .max(Double::compareTo);

        return snpMax.isPresent() && snpMax.get() > vqScoreThreshold;
    }
}
