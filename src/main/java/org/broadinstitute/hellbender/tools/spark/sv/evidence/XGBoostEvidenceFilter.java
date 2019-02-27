package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import biz.k11i.xgboost.Predictor;
import biz.k11i.xgboost.learner.ObjFunction;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.bed.BEDFeature;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import org.broadinstitute.hellbender.tools.spark.sv.evidence.BreakpointEvidence.*;

import java.io.*;
import java.util.*;

/**
 * A class that acts as a filter for BreakpointEvidence.
 * Features are calculated according to evidence type, overlap information, mapping quality, etc.
 * A trained classifier scores the probability the evidence overlaps a breakpoint interval, and passes evidence above
 * the specified threshold.
 */
public final class XGBoostEvidenceFilter implements Iterator<BreakpointEvidence> {
    // use fast math exp for logistic function in XGBoost?
    private static final boolean USE_FAST_MATH_EXP = true;

    private static final List<Class<?>> DEFAULT_EVIDENCE_TYPE_ORDER = Arrays.asList(
            TemplateSizeAnomaly.class, MateUnmapped.class, InterContigPair.class,
            SplitRead.class, LargeIndel.class, WeirdTemplateSize.class, SameStrandPair.class, OutiesPair.class
    );
    private static final Map<Class<?>, Integer> evidenceTypeMap = evidenceTypeOrderToImmutableMap(DEFAULT_EVIDENCE_TYPE_ORDER);
    private static final String DEFAULT_PREDICTOR_RESOURCE_PATH = "/large/sv_evidence_classifier.bin";
    private static final double DEFAULT_GOOD_GAP_OVERLAP = 0.0;
    private static final double DEFAULT_GOOD_MAPPABILITY = 1.0;
    private static final int DEFAULT_GOOD_MAPPING_QUALITY = 60;
    private static final double NON_READ_MAPPING_QUALITY = DEFAULT_GOOD_MAPPING_QUALITY; // alternatively could be Double.NaN
    private static final double NON_READ_CIGAR_LENGTHS = 0.0; // alternatively could be Double.NaN

    private final PartitionCrossingChecker partitionCrossingChecker;

    private final Predictor predictor;
    private final double thresholdProbability;
    private final ReadMetadata readMetadata;

    private final EvidenceOverlapChecker evidenceOverlapChecker;
    private final Map<BreakpointEvidence, UnscaledOverlapInfo> rawFeatureCache;

    private Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> treeItr;
    private Iterator<BreakpointEvidence> listItr;
    private final FeatureDataSource<BEDFeature> genomeGaps;
    private final FeatureDataSource<BEDFeature> umapS100Mappability;

    XGBoostEvidenceFilter(
            final Iterator<BreakpointEvidence> evidenceItr,
            final ReadMetadata readMetadata,
            final StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection params,
            final PartitionCrossingChecker partitionCrossingChecker
    ) {
        if(params.svGenomeGapsFile == null && params.runWithoutGapsAnnotation) {
            genomeGaps = null;
        } else if(params.svGenomeGapsFile != null && !params.runWithoutGapsAnnotation) {
            genomeGaps = new FeatureDataSource<>(params.svGenomeGapsFile);
        } else {
            throw new IllegalArgumentException(
                    "XGBoostEvidenceFilter requires specifying --sv-genome-gaps-file or passing --run-without-gaps-annotation (but not both)"
            );
        }
        if(params.svGenomeUmapS100File == null && params.runWithoutUmapS100Annotation) {
            umapS100Mappability = null;
        } else if(params.svGenomeUmapS100File != null && !params.runWithoutUmapS100Annotation) {
            umapS100Mappability = new FeatureDataSource<>(params.svGenomeUmapS100File);
        } else {
            throw new IllegalArgumentException(
                    "XGBoostEvidenceFilter requires specifying --sv-genome-umap-s100-file or passing --run-without-umap-s100-annotation (but not both)"
            );
        }
        predictor = loadPredictor(params.svEvidenceFilterModelFile);

        this.partitionCrossingChecker = partitionCrossingChecker;
        thresholdProbability = params.svEvidenceFilterThresholdProbability;
        this.readMetadata = readMetadata;

        evidenceOverlapChecker = new EvidenceOverlapChecker(evidenceItr, readMetadata, params.minEvidenceMapQ);
        rawFeatureCache = new HashMap<>();

        listItr = null;
        treeItr = evidenceOverlapChecker.getTreeIterator();
    }

    private static Map<Class<?>, Integer> evidenceTypeOrderToImmutableMap(final List<Class<?>> evidenceTypeOrder) {
        final HashMap<Class<?>, Integer> evidenceTypeMap = new HashMap<>();
        for(int index = 0; index < evidenceTypeOrder.size(); ++index) {
            evidenceTypeMap.put(evidenceTypeOrder.get(index), index);
        }
        return Collections.unmodifiableMap(evidenceTypeMap);
    }

    public static Predictor loadPredictor(final String modelFileLocation) {
        ObjFunction.useFastMathExp(USE_FAST_MATH_EXP);
        try(final InputStream inputStream = modelFileLocation == null ?
                resourcePathToInputStream(DEFAULT_PREDICTOR_RESOURCE_PATH) : BucketUtils.openFile(modelFileLocation)) {
            return new Predictor(inputStream);
        } catch(Exception e) {
            throw new GATKException(
                    "Unable to load predictor from classifier file "
                            + (modelFileLocation == null ? DEFAULT_PREDICTOR_RESOURCE_PATH : modelFileLocation)
                            + ": " + e.getMessage()
            );
        }
    }

    private static InputStream resourcePathToInputStream(final String resourcePath) throws IOException {
        final InputStream inputStream = XGBoostEvidenceFilter.class.getResourceAsStream(resourcePath);
        return IOUtil.hasBlockCompressedExtension(resourcePath) ?
                IOUtils.makeZippedInputStream(new BufferedInputStream(inputStream))
                : inputStream;
    }

    @Override
    public boolean hasNext() {
        if ( listItr != null && listItr.hasNext() ) {
            return true;
        }
        listItr = null;
        boolean result = false;
        while ( !result && treeItr.hasNext() ) {
            final SVIntervalTree.Entry<List<BreakpointEvidence>> entry = treeItr.next();
            final SVInterval curInterval = entry.getInterval();
            final List<BreakpointEvidence> evidenceList = entry.getValue();
            if( isValidated(entry.getValue()) || partitionCrossingChecker.onBoundary(curInterval) ) {
                // already validated (no need to mark validated again) or on partition boundary (punt for now)
                result = true;
            } else if( anyPassesFilter(evidenceList) ) {
                evidenceList.forEach(ev -> ev.setValidated(true));
                result = true;
            }

            if ( result ) {
                listItr = entry.getValue().iterator();
            }
        }
        return result;
    }

    @Override
    public BreakpointEvidence next() {
        if ( !hasNext() ) {
            throw new NoSuchElementException("No next element.");
        }
        return listItr.next();
    }

    private boolean isValidated( final List<BreakpointEvidence> evList ) {
        for ( final BreakpointEvidence ev : evList ) {
            if ( ev.isValidated() ) return true;
        }
        return false;
    }

    private boolean anyPassesFilter(final List<BreakpointEvidence> evidenceList) {
        for(final BreakpointEvidence evidence : evidenceList) {
            if(predictProbability(evidence) > thresholdProbability) {
                return true;
            }
        }
        return false;
    }

    @VisibleForTesting
    double predictProbability(final BreakpointEvidence evidence) {
        return predictor.predictSingle(getFeatures(evidence));
    }

    /**
     * Compute features vector for a piece of BreakpointEvidence
     */
    @VisibleForTesting
    EvidenceFeatures getFeatures(final BreakpointEvidence evidence) {
        // create new struct for these two, use CigarOperator to update if it's ReadEvidence
        final CigarQualityInfo cigarQualityInfo = new CigarQualityInfo(evidence);
        final double evidenceType = evidenceTypeMap.get(evidence.getClass());
        final double mappingQuality = getMappingQuality(evidence);

        // calculate these similar to BreakpointDensityFilter, but always calculate full totals, never end early.
        final CoverageScaledOverlapInfo individualOverlapInfo = getIndividualOverlapInfo(evidence);
        final CoverageScaledOverlapInfo clusterOverlapInfo = getClusterOverlapInfo(evidence);

        // calculate properties related to overlap of intervals on the reference genome
        final double referenceGapOverlap = genomeGaps == null ?
                DEFAULT_GOOD_GAP_OVERLAP
                : getGenomeIntervalsOverlap(evidence, genomeGaps, readMetadata);
        final double umapS100 = umapS100Mappability == null ?
                DEFAULT_GOOD_MAPPABILITY
                : getGenomeIntervalsOverlap(evidence, umapS100Mappability, readMetadata);

        // either templateSize is defined (for ReadEvidence) or readCount (for TemplateSizeAnomaly).
        final double templateSizeOrReadCount = getTemplateSizeOrReadCount(evidence);
        return new EvidenceFeatures(
                new double[]{
                        cigarQualityInfo.basesMatched, cigarQualityInfo.referenceLength, evidenceType, mappingQuality, templateSizeOrReadCount,
                        individualOverlapInfo.numOverlap, individualOverlapInfo.totalOverlapMappingQuality,
                        individualOverlapInfo.meanOverlapMappingQuality, individualOverlapInfo.numCoherent,
                        individualOverlapInfo.totalCoherentMappingQuality,
                        clusterOverlapInfo.numOverlap, clusterOverlapInfo.totalOverlapMappingQuality,
                        clusterOverlapInfo.meanOverlapMappingQuality, clusterOverlapInfo.numCoherent,
                        clusterOverlapInfo.totalCoherentMappingQuality,
                        referenceGapOverlap, umapS100
                }
        );
    }

    /**
     * Return mapping quality for BreakpointEvidence for the purpose of only describing this evidence (no combination
     * with overlappers).
     * For non-ReadEvidence, depending on feature-selection strategy, return NaN or "max" mapping quality (Non-ReadEvidence
     * isn't *bad* per se, so give it a good score).
     */
    private double getMappingQuality(final BreakpointEvidence evidence) {
        // Note: return "max" mapping quality for non-ReadEvidence. Reasoning: some features depend on sum or average of
        // read qualities. Non-ReadEvidence isn't *bad* per se, so give it a good score.
        return evidence instanceof ReadEvidence ? ((ReadEvidence) evidence).getMappingQuality() : NON_READ_MAPPING_QUALITY;
    }

    /**
     * Return mapping quality for BreakpointEvidence for the purposes of calculating sum over overlappers.
     * return "max" mapping quality for non-ReadEvidence. Reasoning: returning NaN will corrupt sums. Non-ReadEvidence
     * isn't *bad* per se, so give it a good score.
     */
    private int getMappingQualityForOverlap(final BreakpointEvidence evidence) {
        // Note: return "max" mapping quality for non-ReadEvidence. Reasoning: features using this function depend on
        // sum or average of read qualities. Non-ReadEvidence isn't *bad* per se, so give it a good score.
        return evidence instanceof ReadEvidence ? ((ReadEvidence) evidence).getMappingQuality() : DEFAULT_GOOD_MAPPING_QUALITY;
    }

    private double getTemplateSizeOrReadCount(final BreakpointEvidence evidence) {
        if(evidence instanceof ReadEvidence) {
            return getTemplateSize((ReadEvidence) evidence);
        } else if(evidence instanceof TemplateSizeAnomaly) {
            return getReadCounts((TemplateSizeAnomaly) evidence);
        } else {
            throw new IllegalStateException("templateSizeOrReadCount feature is only defined for ReadEvidence and TemplateSizeAnomaly, not "
                    + evidence.getClass().getName());
        }
    }

    /** For ReadEvidence, return templateSize as percentile of library's cumulative density function */
    private double getTemplateSize(final ReadEvidence readEvidence) {

        final int templateSize = readEvidence.getTemplateSize();
        final String readGroup = readEvidence.getReadGroup();
        final String library = readMetadata.getReadGroupToLibraryMap().get(readGroup);
        final LibraryStatistics libraryStatistics = readMetadata.getLibraryStatistics(library);
        final IntHistogram.CDF templateSizeCDF = libraryStatistics.getCDF();
        final int cdfBin = Integer.min(Math.abs(templateSize), templateSizeCDF.size() - 1);
        return templateSizeCDF.getFraction(cdfBin);
    }

    /** for TemplateSizeAnomaly, return readCounts scaled by average genome coverage */
    private double getReadCounts(final TemplateSizeAnomaly templateSizeAnomaly) {
        final Integer readCounts = templateSizeAnomaly.getReadCount();
        return (double)(readCounts) / readMetadata.getCoverage();
    }

    private CoverageScaledOverlapInfo getIndividualOverlapInfo(final BreakpointEvidence evidence) {
        // Since overlap info will be needed for the same evidence in different contexts, it's fastest to calculate it
        // once, cache it, and then just retrieve the info each time it's needed.
        if(!rawFeatureCache.containsKey(evidence)) {
            cacheOverlapInfo(evidence);
        }
        final UnscaledOverlapInfo evidenceFeatureCache = rawFeatureCache.get(evidence);
        // Calculate the coverage scaled overlap info
        return new CoverageScaledOverlapInfo(
                evidenceFeatureCache.numOverlap, evidenceFeatureCache.numCoherent,
                evidenceFeatureCache.totalOverlapMappingQuality, evidenceFeatureCache.totalCoherentMappingQuality,
                evidenceFeatureCache.meanOverlapMappingQuality, readMetadata.getCoverage()
        );
    }

    private CoverageScaledOverlapInfo getClusterOverlapInfo(final BreakpointEvidence evidence) {
        int clusterNumOverlap = 0;
        int clusterNumCoherent = 0;
        int clusterOverlapMappingQuality = 0;
        int clusterCoherentMappingQuality = 0;
        double clusterMeanOverlapMappingQuality = 0.0;
        for (final Iterator<BreakpointEvidence> overlapperItr = evidenceOverlapChecker.overlappers(evidence); overlapperItr.hasNext(); ) {
            final BreakpointEvidence overlapper = overlapperItr.next();
            if (overlapper.equals(evidence)) {
                continue; // don't count self-overlap in cluster features
            }
            if(!rawFeatureCache.containsKey(overlapper)) {
                cacheOverlapInfo(overlapper);
            }
            final UnscaledOverlapInfo overlapperFeatureCache = rawFeatureCache.get(overlapper);
            clusterNumOverlap = Math.max(clusterNumOverlap, overlapperFeatureCache.numOverlap);
            clusterNumCoherent = Math.max(clusterNumCoherent, overlapperFeatureCache.numCoherent);
            clusterOverlapMappingQuality = Math.max(clusterOverlapMappingQuality, overlapperFeatureCache.totalOverlapMappingQuality);
            clusterCoherentMappingQuality = Math.max(clusterCoherentMappingQuality, overlapperFeatureCache.totalCoherentMappingQuality);
            clusterMeanOverlapMappingQuality = Math.max(clusterMeanOverlapMappingQuality, overlapperFeatureCache.meanOverlapMappingQuality);
        }

        return new CoverageScaledOverlapInfo(clusterNumOverlap, clusterNumCoherent, clusterOverlapMappingQuality,
                clusterCoherentMappingQuality, clusterMeanOverlapMappingQuality, readMetadata.getCoverage());
    }

    /**
     * For given BreakpointEvidence, calculate number of overlappers, coherent pieces of BreakpointEvidence, and
     * the sums of their mapping qualities. Cache this information so that it can be looked up when computing max or
     * average of these values over the set of overlapping evidence (as in getClusterOverlapInfo).
     */
    private void cacheOverlapInfo(final BreakpointEvidence evidence) {
        int numOverlap = 0;
        int totalOverlapMappingQuality = 0;
        int numCoherent = 0;
        int totalCoherentMappingQuality = 0;
        for(final EvidenceOverlapChecker.OverlapAndCoherenceIterator overlapperItr
                = evidenceOverlapChecker.overlappersWithCoherence(evidence);
            overlapperItr.hasNext();) {
            final ImmutablePair<BreakpointEvidence, Boolean> itrResults = overlapperItr.next();
            final BreakpointEvidence overlapper = itrResults.left;
            if(overlapper.equals(evidence)) {
                continue; // don't count self-overlap
            }
            ++numOverlap;
            final int mappingQuality = getMappingQualityForOverlap(overlapper);
            totalOverlapMappingQuality += mappingQuality;

            final boolean isCoherent = itrResults.right;
            if(isCoherent) {
                ++numCoherent;
                totalCoherentMappingQuality += mappingQuality;
            }
        }
        rawFeatureCache.put(evidence,
                new UnscaledOverlapInfo(numOverlap, numCoherent, totalOverlapMappingQuality, totalCoherentMappingQuality));
    }

    private static class UnscaledOverlapInfo {
        final int numOverlap;
        final int numCoherent;
        final int totalOverlapMappingQuality;
        final int totalCoherentMappingQuality;
        final double meanOverlapMappingQuality;
        UnscaledOverlapInfo(final int numOverlap, final int numCoherent, final int totalOverlapMappingQuality,
                            final int totalCoherentMappingQuality) {
            this.numOverlap = numOverlap;
            this.numCoherent = numCoherent;
            this.totalOverlapMappingQuality = totalOverlapMappingQuality;
            this.totalCoherentMappingQuality = totalCoherentMappingQuality;
            this.meanOverlapMappingQuality = ((double)this.totalOverlapMappingQuality) / numOverlap;
        }
    }

    /**
     * Class that takes raw overlap info and automatically scales it to meanGenomeCoverage, storing the result for later retrieval
     */
    private static class CoverageScaledOverlapInfo {
        final double numOverlap;
        final double totalOverlapMappingQuality;
        final double meanOverlapMappingQuality;
        final double numCoherent;
        final double totalCoherentMappingQuality;

        CoverageScaledOverlapInfo(final int numOverlap, final int numCoherent, final int totalOverlapMappingQuality,
                                  final int totalCoherentMappingQuality, final double meanOverlapMappingQuality,
                                  final double coverage) {
            this.numOverlap = ((double)numOverlap) / coverage;
            this.totalOverlapMappingQuality = ((double) totalOverlapMappingQuality) / coverage;
            this.numCoherent = ((double)numCoherent) / coverage;
            this.totalCoherentMappingQuality = ((double) totalCoherentMappingQuality) / coverage;
            this.meanOverlapMappingQuality = meanOverlapMappingQuality;
        }
    }

    private static class CigarQualityInfo {
        final double basesMatched;
        final double referenceLength;

        CigarQualityInfo(final BreakpointEvidence evidence) {
            if(evidence instanceof ReadEvidence) {
                int numMatched = 0;
                int refLength = 0;
                final String cigarString = ((ReadEvidence) evidence).getCigarString();
                for (final CigarElement element : TextCigarCodec.decode(cigarString).getCigarElements()) {
                    final CigarOperator op = element.getOperator();
                    if (op.consumesReferenceBases()) {
                        refLength += element.getLength();
                        if (op.consumesReadBases()) {
                            numMatched += element.getLength();
                        }
                    }
                }
                basesMatched = numMatched;
                referenceLength = refLength;
            }
            else {
                basesMatched = NON_READ_CIGAR_LENGTHS;
                referenceLength = NON_READ_CIGAR_LENGTHS;
            }
        }
    }

    /**
     * Calculate fractional overlap of BreakpointEvidence with genome tract data.
     * Separately sum the number of base pairs in genomeIntervals that overlap evidence (allowing base pairs to
     * count multiple times if there is overlap in genomeIntervals) then divide by length of evidence. returned
     * value will be double >= 0, but may be larger than 1 if any genomeIntervals overlap each other.
     */
    private static double getGenomeIntervalsOverlap(final BreakpointEvidence evidence,
                                                    final FeatureDataSource<BEDFeature> genomeIntervals,
                                                    final ReadMetadata readMetadata) {
        final SVInterval location = evidence.getLocation();
        final SimpleInterval simpleInterval = new SimpleInterval(readMetadata.getContigName(location.getContig()),
                location.getStart(),
                location.getEnd() - 1); // "end - 1" because SimpleIntervals are closed on both ends
        int overlap = 0;
        for(final Iterator<BEDFeature> overlapperItr = genomeIntervals.query(simpleInterval); overlapperItr.hasNext();) {
            final BEDFeature overlapper = overlapperItr.next();
            // " + 1" because genome tract data is semi-closed, but BEDFeature is fully closed
            final int overlapLength = Math.min(simpleInterval.getEnd(), overlapper.getEnd()) + 1
                    - Math.max(simpleInterval.getStart(), overlapper.getStart());
            overlap += overlapLength;
        }
        return overlap / (double)simpleInterval.size();
    }
}
