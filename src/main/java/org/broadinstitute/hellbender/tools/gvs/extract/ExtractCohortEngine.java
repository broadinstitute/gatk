package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import io.opencensus.trace.Span;
import org.apache.avro.generic.GenericRecord;
import org.apache.commons.lang.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.GQStateEnum;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.SampleList;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.gvs.common.VariantBitSet;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.*;
import org.broadinstitute.hellbender.utils.localsort.AvroSortingCollectionCodec;
import org.broadinstitute.hellbender.utils.localsort.SortingCollection;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.function.Consumer;
import java.util.stream.Collectors;

public class ExtractCohortEngine {
    private static final Logger logger = LogManager.getLogger(ExtractCohortEngine.class);

    private final VariantContextWriter vcfWriter;

    private final boolean printDebugInformation;
    private final int localSortMaxRecordsInRam;
    private final List<SimpleInterval> traversalIntervals;
    private final Long minLocation;
    private final Long maxLocation;
    private final TableReference filterSetInfoTableRef;
    private final TableReference filterSetSiteTableRef;
    private final ReferenceDataSource refSource;
    private final Double vqsLodSNPThreshold;
    private final Double vqsLodINDELThreshold;
    private final ExtractCohort.VQSLODFilteringType VQSLODFilteringType;
    private final boolean excludeFilteredSites;

    private final ProgressMeter progressMeter;
    private final String projectID;
    private final boolean emitPLs;
    private final boolean emitADs;

    /** List of sample names seen in the variant data from BigQuery. */
    private final Set<String> sampleNames;
    private final SortedSet<Long> sampleIdsToExtract;
    private final BitSet sampleIdsToExtractBitSet;

    private final Map<Long, String> sampleIdToName;

    private final ReferenceConfidenceVariantContextMerger variantContextMerger;
    private final VariantAnnotatorEngine annotationEngine;

    private int totalNumberOfVariants = 0;
    private int totalNumberOfSites = 0;
    private int totalRangeRecords = 0;
    private int totalIrrelevantRangeRecords = 0;
    private long totalEstimatedBytesScanned = 0;

    private final GATKPath cohortAvroFileName;

    private final String vetRangesFQDataSet;

    private final String fqRangesExtractVetTable;
    private final String fqRangesExtractRefTable;

    private final GATKPath vetAvroFileName;
    private final GATKPath refRangesAvroFileName;
    private final String filterSetName;

    private final GQStateEnum inferredReferenceState;
    private final InferredReferenceRecord inferredReferenceRecord;

    private final boolean presortedAvroFiles;


    public ExtractCohortEngine(final String projectID,
                               final VariantContextWriter vcfWriter,
                               final VCFHeader vcfHeader,
                               final VariantAnnotatorEngine annotationEngine,
                               final ReferenceDataSource refSource,
                               final Map<Long, String> sampleIdToName,
                               final String cohortTableName,
                               final GATKPath cohortAvroFileName,
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
                               final Double vqsLodSNPThreshold,
                               final Double vqsLodINDELThreshold,
                               final ProgressMeter progressMeter,
                               final String filterSetName,
                               final boolean emitPLs,
                               final boolean emitADs,
                               final ExtractCohort.VQSLODFilteringType VQSLODFilteringType,
                               final boolean excludeFilteredSites,
                               final GQStateEnum inferredReferenceState,
                               final boolean presortedAvroFiles
    ) {
        this.localSortMaxRecordsInRam = localSortMaxRecordsInRam;

        this.projectID = projectID;
        this.vcfWriter = vcfWriter;
        this.refSource = refSource;
        this.sampleIdToName = sampleIdToName;
        this.sampleNames = new HashSet<>(sampleIdToName.values());
        this.sampleIdsToExtract = new TreeSet<>(this.sampleIdToName.keySet());

        long maxSampleId = sampleIdsToExtract.last();

        // Constraint due to easily available BitSet implementations, could
        // be extended in future work
        if (maxSampleId > Integer.MAX_VALUE) {
            throw new GATKException("Sample Ids > " + Integer.MAX_VALUE + " are not supported");
        }

        this.sampleIdsToExtractBitSet = new BitSet((int) maxSampleId + 1);
        for(Long id : sampleIdsToExtract) {
            sampleIdsToExtractBitSet.set(id.intValue());
        }

        this.emitPLs = emitPLs;
        this.emitADs = emitADs;

        this.vetRangesFQDataSet = vetRangesFQDataSet;
        this.fqRangesExtractVetTable = fqRangesExtractVetTable;
        this.fqRangesExtractRefTable = fqRangesExtractRefTable;

        this.cohortAvroFileName = cohortAvroFileName;
        this.vetAvroFileName = vetAvroFileName;
        this.refRangesAvroFileName = refRangesAvroFileName;

        this.traversalIntervals = traversalIntervals;
        this.minLocation = minLocation;
        this.maxLocation = maxLocation;

        this.printDebugInformation = printDebugInformation;
        this.vqsLodSNPThreshold = vqsLodSNPThreshold;
        this.vqsLodINDELThreshold = vqsLodINDELThreshold;
        this.VQSLODFilteringType = VQSLODFilteringType;
        this.excludeFilteredSites = excludeFilteredSites;

        this.filterSetInfoTableRef = VQSLODFilteringType.equals(ExtractCohort.VQSLODFilteringType.NONE) ? null : new TableReference(filterSetInfoTableName, SchemaUtils.YNG_FIELDS);
        this.filterSetSiteTableRef = VQSLODFilteringType.equals(ExtractCohort.VQSLODFilteringType.NONE) ? null : new TableReference(filterSetSiteTableName, SchemaUtils.FILTER_SET_SITE_FIELDS);

        this.progressMeter = progressMeter;

        this.filterSetName = filterSetName;

        this.annotationEngine = annotationEngine;
        this.variantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine, vcfHeader);

        this.inferredReferenceState = inferredReferenceState;

        this.presortedAvroFiles = presortedAvroFiles;

        this.inferredReferenceRecord = new InferredReferenceRecord(inferredReferenceState);

    }

    int getTotalNumberOfVariants() { return totalNumberOfVariants; }
    int getTotalNumberOfSites() { return totalNumberOfSites; }
    long getTotalEstimatedBytesScanned() { return totalEstimatedBytesScanned; }

    private void processBytesScanned(StorageAPIAvroReader reader) {
        long bytes = reader.getEstimatedTotalBytesScanned();
        totalEstimatedBytesScanned += bytes;
    }

    public void traverse() {
        //First allele here is the ref, followed by the alts associated with that ref. We need this because at this point the alleles haven't been joined and remapped to one reference allele.
        final HashMap<Long, HashMap<Allele, HashMap<Allele, Double>>> fullVqsLodMap = new HashMap<>();
        final HashMap<Long, HashMap<Allele, HashMap<Allele, String>>> fullYngMap = new HashMap<>();
        final HashMap<Long, List<String>> siteFilterMap = new HashMap<>();

        String rowRestriction = null;
        if (minLocation != null && maxLocation != null) {
            rowRestriction = "location >= " + minLocation + " AND location <= " + maxLocation;
        }
        final String rowRestrictionWithFilterSetName = rowRestriction + " AND " + SchemaUtils.FILTER_SET_NAME + " = '" + filterSetName + "'";

        boolean noVqslodFilteringRequested = VQSLODFilteringType.equals(ExtractCohort.VQSLODFilteringType.NONE);
        if (!noVqslodFilteringRequested) {
            // ensure vqslod filters are defined. this really shouldn't ever happen, said the engineer.
            if (vqsLodSNPThreshold == null || vqsLodINDELThreshold == null) {
                throw new UserException("Vqslod filtering thresholds for SNPs and INDELs must be defined.");
            }

            // get filter info (vqslod & yng values)
            try (StorageAPIAvroReader reader = new StorageAPIAvroReader(filterSetInfoTableRef, rowRestrictionWithFilterSetName, projectID)) {

                for ( final GenericRecord queryRow : reader ) {
                    final ExtractCohortFilterRecord filterRow = new ExtractCohortFilterRecord( queryRow );

                    final long location = filterRow.getLocation();
                    final Double vqslod = filterRow.getVqslod();
                    final String yng = filterRow.getYng();
                    final Allele ref = Allele.create(filterRow.getRefAllele(), true);
                    final Allele alt = Allele.create(filterRow.getAltAllele(), false);
                    fullVqsLodMap.putIfAbsent(location, new HashMap<>());
                    fullVqsLodMap.get(location).putIfAbsent(ref, new HashMap<>());
                    fullVqsLodMap.get(location).get(ref).put(alt, vqslod);
                    fullYngMap.putIfAbsent(location, new HashMap<>());
                    fullYngMap.get(location).putIfAbsent(ref, new HashMap<>());
                    fullYngMap.get(location).get(ref).put(alt, yng);
                }
                processBytesScanned(reader);
            }
        }

        // load site-level filter data into data structure
        if (filterSetSiteTableRef != null) {
            try (StorageAPIAvroReader reader = new StorageAPIAvroReader(filterSetSiteTableRef, rowRestrictionWithFilterSetName, projectID)) {
                for ( final GenericRecord queryRow : reader ) {
                    long location = Long.parseLong(queryRow.get(SchemaUtils.LOCATION_FIELD_NAME).toString());
                    List<String> filters = Arrays.asList(queryRow.get(SchemaUtils.FILTERS).toString().split(","));
                    siteFilterMap.put(location, filters);
                }
                processBytesScanned(reader);
            }
        }

        if (printDebugInformation) {
            logger.debug("using storage api with local sort");
        }
        logger.debug("Initializing Reader");

        if (!SchemaUtils.decodeContig(minLocation).equals(SchemaUtils.decodeContig(maxLocation))) {
            throw new GATKException("Can not process cross-contig boundaries for Ranges implementation");
        }

        SortedSet<Long> sampleIdsToExtract = new TreeSet<>(this.sampleIdToName.keySet());
        if (fqRangesExtractVetTable != null) {
            createVariantsFromUnsortedExtractTableBigQueryRanges(fqRangesExtractVetTable, fqRangesExtractRefTable, sampleIdsToExtract, minLocation, maxLocation, fullVqsLodMap, fullYngMap, siteFilterMap, noVqslodFilteringRequested);
        } else if (vetRangesFQDataSet != null) {
            createVariantsFromUnsortedBigQueryRanges(vetRangesFQDataSet, sampleIdsToExtract, minLocation, maxLocation, fullVqsLodMap, fullYngMap, siteFilterMap, noVqslodFilteringRequested);
        } else {
            createVariantsFromUnsortedAvroRanges(vetAvroFileName, refRangesAvroFileName, sampleIdsToExtract, minLocation, maxLocation, fullVqsLodMap, fullYngMap, siteFilterMap, noVqslodFilteringRequested, presortedAvroFiles);
        }

        logger.debug("Finished Initializing Reader");

        logger.info("Processed " + totalRangeRecords + " range records and rejected " + totalIrrelevantRangeRecords + " irrelevant ones ");
    }


    public static SortingCollection<GenericRecord> getAvroSortingCollection(org.apache.avro.Schema schema, int localSortMaxRecordsInRam) {
        final SortingCollection.Codec<GenericRecord> sortingCollectionCodec = new AvroSortingCollectionCodec(schema);
        final Comparator<GenericRecord> sortingCollectionComparator = new Comparator<GenericRecord>() {
            @Override
            public int compare( GenericRecord o1, GenericRecord o2 ) {
                final long firstPosition = (Long) o1.get(SchemaUtils.LOCATION_FIELD_NAME);
                final long secondPosition = (Long) o2.get(SchemaUtils.LOCATION_FIELD_NAME);

                final int result = Long.compare(firstPosition, secondPosition);
                if (result != 0) {
                    return result;
                } else {
                    final long firstSample = (Long) o1.get(SchemaUtils.SAMPLE_ID_FIELD_NAME);
                    final long secondSample = (Long) o2.get(SchemaUtils.SAMPLE_ID_FIELD_NAME);
                    return Long.compare(firstSample, secondSample);
                }
            }
        };
        return SortingCollection.newInstance(GenericRecord.class, sortingCollectionCodec, sortingCollectionComparator, localSortMaxRecordsInRam, true);
    }

    private SortingCollection<GenericRecord> addToVetSortingCollection(final SortingCollection<GenericRecord> sortingCollection,
                                                                       final Iterable<GenericRecord> avroReader,
                                                                       final VariantBitSet vbs) {
        int recordsProcessed = 0;
        long startTime = System.currentTimeMillis();

        // NOTE: if OverlapDetector takes too long, try using RegionChecker from tws_sv_local_assembler
        // need to manually add the upstream padding to the set of intervals
        List<SimpleInterval> overlapIntervals = new ArrayList<>();
        overlapIntervals.addAll(traversalIntervals);
        overlapIntervals.add(new SimpleInterval(SchemaUtils.decodeContig(minLocation),
                                                SchemaUtils.decodePosition(minLocation) - IngestConstants.MAX_DELETION_SIZE + 1,
                                                SchemaUtils.decodePosition(minLocation)));
        final OverlapDetector<SimpleInterval> intervalsOverlapDetector = OverlapDetector.create(overlapIntervals);

        for (final GenericRecord queryRow : avroReader) {
            long location = (Long) queryRow.get(SchemaUtils.LOCATION_FIELD_NAME);
            int position = SchemaUtils.decodePosition(location);
            SimpleInterval simpleInverval = new SimpleInterval(SchemaUtils.decodeContig(location), position, position);

            if (intervalsOverlapDetector.overlapsAny(simpleInverval)) {
                vbs.setVariant(location);
                sortingCollection.add(queryRow);
                if (++recordsProcessed % 1000000 == 0) {
                    long endTime = System.currentTimeMillis();
                    logger.info("Processed " + recordsProcessed + " VET records in " + (endTime - startTime) + " ms");
                    startTime = endTime;
                }
            }
        }

        long endTime = System.currentTimeMillis();
        logger.info("Processed " + recordsProcessed + " VET records in " + (endTime - startTime) + " ms");

        sortingCollection.printTempFileStats();
        return sortingCollection;
    }

    private SortingCollection<GenericRecord> addToRefSortingCollection(final SortingCollection<GenericRecord> sortingCollection, final Iterable<GenericRecord> avroReader, final VariantBitSet vbs) {
        int recordsProcessed = 0;
        long startTime = System.currentTimeMillis();

        for (final GenericRecord queryRow : avroReader) {
            long location = (Long) queryRow.get(SchemaUtils.LOCATION_FIELD_NAME);
            int length = ((Long) queryRow.get(SchemaUtils.LENGTH_FIELD_NAME)).intValue();

            if (vbs.containsVariant(location, location + length) ) {
                sortingCollection.add(queryRow);
            }

            if (++recordsProcessed % 1000000 == 0) {
                long endTime = System.currentTimeMillis();
                logger.info("Processed " + recordsProcessed + " Reference Ranges records in " + (endTime - startTime) + " ms");
                startTime = endTime;
            }
        }

        long endTime = System.currentTimeMillis();
        logger.info("Processed " + recordsProcessed + " Reference Ranges records in " + (endTime - startTime) + " ms");

        sortingCollection.printTempFileStats();
        return sortingCollection;
    }

    private ExtractCohortRecord mergeSampleRecord(ExtractCohortRecord r1, ExtractCohortRecord r2) {
        if (printDebugInformation) {
            logger.info("In SampleRecord Merge Logic for " + r1 + " and " + r2);
        }

        final String r1State = r1.getState();
        final String r2State = r2.getState();

        // Valid states are m, 1-6 (ref), v, *
        // For now, just handle cases where '*' should be dropped in favor anything besides missing...
        if (r1State.equals("*") && !r2State.equals("m")) {
            return r2;
        } else if (r2State.equals("*") && !r1State.equals("m")) {
            return r1;
        } else {
            return r2;
        }
    }

    protected VariantContext processSampleRecordsForLocation(final long location,
                                                             final Iterable<ExtractCohortRecord> sampleRecordsAtPosition,
                                                             final HashMap<Long, HashMap<Allele, HashMap<Allele, Double>>> fullVqsLodMap,
                                                             final HashMap<Long, HashMap<Allele, HashMap<Allele, String>>> fullYngMap,
                                                             final boolean noVqslodFilteringRequested,
                                                             final HashMap<Long, List<String>> siteFilterMap,
                                                             final ExtractCohort.VQSLODFilteringType VQSLODFilteringType) {

        final List<VariantContext> unmergedCalls = new ArrayList<>();
        final List<ReferenceGenotypeInfo> refCalls = new ArrayList<>();

        int maxSampleIdToExtract = sampleIdsToExtract.last().intValue();
        final BitSet samplesSeen = new BitSet(maxSampleIdToExtract + 1);

        boolean currentPositionHasVariant = false;
        final int currentPosition = SchemaUtils.decodePosition(location);
        final String contig = SchemaUtils.decodeContig(location);
        final Allele refAllele = getReferenceAllele(refSource, location, 1);
        int numRecordsAtPosition = 0;

        final HashMap<Allele, HashMap<Allele, Double>> vqsLodMap;
        final HashMap<Allele, HashMap<Allele, String>> yngMap;

        // TODO: optimize in the case where noFilteringRequested == true, no need to populate this

        // If there's no yng/vqslod for this site, then we'll treat these as NAYs because VQSR dropped them (they have no alt reads).
        if (fullVqsLodMap.get(SchemaUtils.encodeLocation(contig, currentPosition)) == null) {
            vqsLodMap = new HashMap<>();
            vqsLodMap.put(refAllele, new HashMap<>());
            yngMap = new HashMap<>();
            yngMap.put(refAllele, new HashMap<>());
        } else {
            vqsLodMap = fullVqsLodMap.get(SchemaUtils.encodeLocation(contig, currentPosition));
            yngMap = fullYngMap.get(SchemaUtils.encodeLocation(contig, currentPosition));
        }

        for ( final ExtractCohortRecord sampleRecord : sampleRecordsAtPosition ) {
            final String sampleName = sampleIdToName.get(sampleRecord.getSampleId());

            if (sampleName == null) {
                throw new GATKException("Unable to translate sample id " + sampleRecord.getSampleId() + " to sample name");
            }

            // PERF: BOTTLENECK (~13%)
            // Note: we've already confirmed that max sample id is an int
            samplesSeen.set(sampleRecord.getSampleId().intValue());
            ++numRecordsAtPosition;

            if ( printDebugInformation ) {
                logger.info("\t" + contig + ":" + currentPosition + ": found record for sample " + sampleName + ": " + sampleRecord);
            }

            switch (sampleRecord.getState()) {
                case "v":   // Variant
                    ++totalNumberOfVariants;
                    VariantContext vc = createVariantContextFromSampleRecord(sampleRecord, vqsLodMap, yngMap, VQSLODFilteringType);
                    unmergedCalls.add(vc);

                    currentPositionHasVariant = true;
                    break;
                case "0":   // Non Variant Block with GQ < 10
                    // Reference calls with GQ 0 should be rendered as no-call (#271)
                    // Nothing to do here -- just needed to mark the sample as seen so it doesn't get put in the high confidence ref band
                    break;
                case "1":  // Non Variant Block with 10 <=  GQ < 20
                    refCalls.add(new ReferenceGenotypeInfo(sampleName, 10));
                    break;
                case "2":  // Non Variant Block with 20 <= GQ < 30
                    refCalls.add(new ReferenceGenotypeInfo(sampleName, 20));
                    break;
                case "3":  // Non Variant Block with 30 <= GQ < 40
                    refCalls.add(new ReferenceGenotypeInfo(sampleName, 30));
                    break;
                case "4":  // Non Variant Block with 40 <= GQ < 50
                    refCalls.add(new ReferenceGenotypeInfo(sampleName, 40));
                    break;
                case "5":  // Non Variant Block with 50 <= GQ < 60
                    refCalls.add(new ReferenceGenotypeInfo(sampleName, 50));
                    break;
                case "6":  // Non Variant Block with 60 <= GQ (usually omitted from tables)
                    refCalls.add(new ReferenceGenotypeInfo(sampleName, 60));
                    break;
                case "*":   // Spanning Deletion
                    VariantContext vcStar = createVariantContextFromSampleRecord(sampleRecord, vqsLodMap, yngMap, VQSLODFilteringType);
                    unmergedCalls.add(vcStar);
                    break;
                case "m":   // Missing
                    // Nothing to do here -- just needed to mark the sample as seen so it doesn't get put in the high confidence ref band
                    break;
                default:
                    throw new GATKException("Unrecognized state: " + sampleRecord.getState());
            }

        }

        if ( printDebugInformation ) {
            logger.info(contig + ":" + currentPosition + ": processed " + numRecordsAtPosition + " total sample records");
        }

        return finalizeCurrentVariant(unmergedCalls, refCalls, samplesSeen, currentPositionHasVariant, location, contig, currentPosition, refAllele, vqsLodMap, yngMap, noVqslodFilteringRequested, siteFilterMap);
    }

    VariantContext finalizeCurrentVariant(final List<VariantContext> unmergedVariantCalls,
                                        final List<ReferenceGenotypeInfo> referenceCalls,
                                        final BitSet samplesSeen,
                                        final boolean currentPositionHasVariant,
                                        final long location,
                                        final String contig,
                                        final long start,
                                        final Allele refAllele,
                                        final HashMap<Allele, HashMap<Allele, Double>> vqsLodMap,
                                        final HashMap<Allele, HashMap<Allele, String>> yngMap,
                                        final boolean noVqslodFilteringRequested,
                                        final HashMap<Long, List<String>> siteFilterMap) {
        // If there were no variants at this site, we don't emit a record and there's nothing to do here
        if ( ! currentPositionHasVariant ) {
            return null;
        }

        VariantContext mergedVC = variantContextMerger.merge(
                unmergedVariantCalls,
                new SimpleInterval(contig, (int) start, (int) start),
                refAllele.getBases()[0],
                true,
                false,
                true);


        // Reference Sites -- first create a single VC Builder
        final VariantContextBuilder vcWithRef = new VariantContextBuilder(mergedVC);

        // create alleles (same for all genotypes)
        List<Allele> gtAlleles = Arrays.asList(mergedVC.getReference(), mergedVC.getReference());

        final GenotypesContext genotypes = GenotypesContext.copy(vcWithRef.getGenotypes());
        GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        // add known ref genotypes
        for ( final ReferenceGenotypeInfo info : referenceCalls ) {
            genotypeBuilder.reset(false);
            genotypeBuilder.name(info.getSampleName());
            genotypeBuilder.alleles(gtAlleles);
            genotypeBuilder.GQ(info.getGQ());
            genotypes.add(genotypeBuilder.make());
        }

        // Find samples for inferred state and synthesize

        // mutates in place, so rename to be clear.
        // we don't use this again as "samplesSeen" so we don't need a full copy
        BitSet samplesNotEncountered = samplesSeen;
        samplesNotEncountered.xor(sampleIdsToExtractBitSet);

        // Iterate through the samples not encountered
        for (int sampleId = samplesNotEncountered.nextSetBit(0); sampleId >= 0; sampleId = samplesNotEncountered.nextSetBit(sampleId+1)) {
            genotypeBuilder.reset(false);
            genotypeBuilder.name(sampleIdToName.get(Long.valueOf(sampleId)));
            genotypeBuilder.alleles(gtAlleles);
            genotypeBuilder.GQ(inferredReferenceState.getReferenceGQ());
            genotypes.add(genotypeBuilder.make());
        }

        vcWithRef.genotypes(genotypes);
        mergedVC = vcWithRef.make();

        ReferenceContext referenceContext = new ReferenceContext(refSource, new SimpleInterval(mergedVC));

        VariantContext genotypedVC = annotationEngine.annotateContext(mergedVC, new FeatureContext(), referenceContext, null, a -> true);

        // apply VQSLod-based filters
        VariantContext filteredVC =
                noVqslodFilteringRequested ? genotypedVC : filterSiteByAlleleSpecificVQSLOD(genotypedVC, vqsLodMap, yngMap, VQSLODFilteringType);

        // apply SiteQC-based filters, if they exist
        if ( siteFilterMap.containsKey(location) ) {
            final VariantContextBuilder sfBuilder = new VariantContextBuilder(filteredVC);

            Set<String> newFilters = new HashSet<>();

            // include existing filters, if any
            if (sfBuilder.getFilters() != null) {
                newFilters.addAll(sfBuilder.getFilters());
            }

            newFilters.addAll(siteFilterMap.get(location));
            filteredVC = sfBuilder.filters(newFilters).make();
        }

        // clean up extra annotations
        final VariantContext finalVC = removeAnnotations(filteredVC);

        return finalVC;
    }

    private void writeToVcf(VariantContext vc) {
        if ( vc != null ) {
            // Add the variant contexts that aren't filtered or add everything if we aren't excluding anything
            if (vc.isNotFiltered() || !excludeFilteredSites) {
                vcfWriter.add(vc);
            }
            progressMeter.update(vc);
        }
    }

    private VariantContext filterSiteByAlleleSpecificVQSLOD(VariantContext mergedVC, HashMap<Allele, HashMap<Allele, Double>> vqsLodMap, HashMap<Allele, HashMap<Allele, String>> yngMap, ExtractCohort.VQSLODFilteringType VQSLODFilteringType) {
        final LinkedHashMap<Allele, Double> remappedVqsLodMap = remapAllelesInMap(mergedVC, vqsLodMap, Double.NaN);
        final LinkedHashMap<Allele, String> remappedYngMap = remapAllelesInMap(mergedVC, yngMap, VCFConstants.EMPTY_INFO_FIELD);

        final LinkedHashMap<Allele, Double> relevantVqsLodMap = new LinkedHashMap<>();
        mergedVC.getAlternateAlleles().forEach(key -> Optional.ofNullable(remappedVqsLodMap.get(key)).ifPresent(value -> relevantVqsLodMap.put(key, value)));
        final LinkedHashMap<Allele, String> relevantYngMap = new LinkedHashMap<>();
        mergedVC.getAlternateAlleles().forEach(key -> Optional.ofNullable(remappedYngMap.get(key)).ifPresent(value -> relevantYngMap.put(key, value)));

        final VariantContextBuilder builder = new VariantContextBuilder(mergedVC);

        builder.attribute(GATKVCFConstants.AS_VQS_LOD_KEY, relevantVqsLodMap.values().stream().map(val -> val.equals(Double.NaN) ? VCFConstants.EMPTY_INFO_FIELD : val.toString()).collect(Collectors.toList()));
        builder.attribute(GATKVCFConstants.AS_YNG_STATUS_KEY, new ArrayList<>(relevantYngMap.values()));

        if (VQSLODFilteringType.equals(ExtractCohort.VQSLODFilteringType.SITES)) { // Note that these filters are not used with Genotype VQSLOD Filtering
            int refLength = mergedVC.getReference().length();

            // if there are any Yays, the site is PASS
            if (remappedYngMap.containsValue("Y")) {
                builder.passFilters();
            } else if (remappedYngMap.containsValue("N")) {
                builder.filter(GATKVCFConstants.NAY_FROM_YNG);
            } else {
                // if it doesn't trigger any of the filters below, we assume it passes.
                builder.passFilters();
                if (remappedYngMap.containsValue("G")) {
                    Optional<Double> snpMax = relevantVqsLodMap.entrySet().stream()
                            .filter(entry -> entry.getKey().length() == refLength)
                            .map(entry -> entry.getValue())
                            .filter(d -> !(d.isNaN()||d.isInfinite()))
                            .max(Double::compareTo);
                    if (snpMax.isPresent() && snpMax.get() < vqsLodSNPThreshold) {
                        builder.filter(GATKVCFConstants.VQSR_FAILURE_SNP);
                    }

                    Optional<Double> indelMax = relevantVqsLodMap.entrySet().stream()
                            .filter(entry -> entry.getKey().length() != refLength)
                            .map(entry -> entry.getValue())
                            .filter(d -> !(d.isNaN()||d.isInfinite()))
                            .max(Double::compareTo);

                    if (indelMax.isPresent() && indelMax.get() < vqsLodINDELThreshold) {
                        builder.filter(GATKVCFConstants.VQSR_FAILURE_INDEL);
                    }
                } else {
                    // per-conversation with Laura, if there is no information we let the site pass (ie no data does not imply failure)
                }
            }
        }

        // TODO: add in other annotations we need in output (like AF, etc?)
        final VariantContext filteredVC = builder.make();
        return filteredVC;
    }

    protected static VariantContext removeAnnotations(VariantContext vc) {

        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        List<String> rmAnnotationList = new ArrayList<>(Arrays.asList(
                GATKVCFConstants.STRAND_ODDS_RATIO_KEY,
                GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY,
                GATKVCFConstants.FISHER_STRAND_KEY));

        builder.rmAttributes(rmAnnotationList);
        return builder.make();
    }

    private <T> LinkedHashMap<Allele, T> remapAllelesInMap(VariantContext vc, HashMap<Allele, HashMap<Allele, T>> datamap, T emptyVal) {
        return remapAllelesInMap(vc.getReference(), vc.getAlternateAlleles(), vc.getContig(), vc.getStart(), datamap, emptyVal);
    }
        /*
     * Alleles from the filtering table need to be remapped to use the same ref allele that the will exist in the joined variant context.
     * This method changes the alleles in the datamap to match the representation that's in the vc.
     */
    private <T> LinkedHashMap<Allele, T> remapAllelesInMap(Allele ref, List<Allele> alternateAlleles, String contig, int start, HashMap<Allele, HashMap<Allele, T>> datamap, T emptyVal) {
        // create ordered results map
        LinkedHashMap<Allele, T> results = new LinkedHashMap<>();
        alternateAlleles.stream().forEachOrdered(allele -> results.put(allele, emptyVal));


        datamap.entrySet().stream().forEachOrdered(entry -> {
            if (entry.getKey().equals(ref)) {
                // reorder
                entry.getValue().entrySet().stream().forEach(altMapEntry -> results.put(altMapEntry.getKey(), altMapEntry.getValue()));
            } else {
                // remap
                int refLength = entry.getKey().length();
                List<Allele> allAlleles = new ArrayList<>();
                allAlleles.add(entry.getKey());
                allAlleles.addAll(entry.getValue().keySet());
                VariantContextBuilder vcb = new VariantContextBuilder("unused", contig, start, start+refLength-1, allAlleles);
                VariantContext newvc = vcb.make();

                //If the length of the reference from the filtering table is longer than the reference in the variantContext, then that allele is not present in the extracted samples and we don't need the data
                if (refLength < ref.length()) {
                    Map<Allele, Allele> alleleMapping = GATKVariantContextUtils.createAlleleMapping(ref, newvc);
                    alleleMapping.entrySet().stream().forEach(mapped -> results.put(mapped.getValue(), entry.getValue().get(mapped.getKey())));
                }
            }
        });

        return results;
    }


    // vqsLogMap and yngMap are in/out parameters for this method. i.e. they are modified by this method
    private VariantContext createVariantContextFromSampleRecord(final ExtractCohortRecord sampleRecord, HashMap<Allele, HashMap<Allele, Double>> vqsLodMap, HashMap<Allele, HashMap<Allele, String>> yngMap, ExtractCohort.VQSLODFilteringType VQSLODFilteringType) {
        final VariantContextBuilder builder = new VariantContextBuilder();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        final String contig = sampleRecord.getContig();
        final long startPosition = sampleRecord.getStart();
        final String sampleName = sampleIdToName.get(sampleRecord.getSampleId());

        if (sampleName == null) {
            throw new GATKException("Unable to translate sample id " + sampleRecord.getSampleId() + " to sample name");
        }

        builder.chr(contig);
        builder.start(startPosition);

        final List<Allele> alleles = new ArrayList<>();
        Allele ref = Allele.create(sampleRecord.getRefAllele(), true);
        alleles.add(ref);
        List<Allele> altAlleles = Arrays.stream(sampleRecord.getAltAllele().split(SchemaUtils.MULTIVALUE_FIELD_DELIMITER))
                .map(altAllele -> Allele.create(altAllele, false)).collect(Collectors.toList());

        // NOTE: gnarly needs this it seems? Should it be dropped now that we wont be using gnarly?
        altAlleles.add(Allele.NON_REF_ALLELE);

        alleles.addAll(altAlleles);
        builder.alleles(alleles);

        builder.stop(startPosition + alleles.get(0).length() - 1);

        genotypeBuilder.name(sampleName);
        vqsLodMap.putIfAbsent(ref, new HashMap<>());
        yngMap.putIfAbsent(ref, new HashMap<>());

        // need to re-prepend the leading "|" to AS_QUALapprox for use in gnarly
        if ( sampleRecord.getAsQUALApprox() != null ) {
            builder.attribute(SchemaUtils.AS_QUALapprox, "|" + sampleRecord.getAsQUALApprox());
        }

        if ( sampleRecord.getQUALApprox() != null ) {
            builder.attribute(SchemaUtils.QUALapprox, sampleRecord.getQUALApprox());
        }

        final String callGT = sampleRecord.getCallGT();
        if ("./.".equals(callGT)) {
            genotypeBuilder.alleles(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));
        } else {
            final List<Allele> genotypeAlleles =
                    Arrays.stream(callGT.split("[/|]"))
                            .map(Integer::parseInt)
                            .map(alleleIndex -> alleles.get(alleleIndex))
                            .collect(Collectors.toList());

            if (VQSLODFilteringType.equals(ExtractCohort.VQSLODFilteringType.GENOTYPE)) {
                final List<Allele> nonRefAlleles =
                        genotypeAlleles.stream()
                                .filter(a -> a.isNonReference())
                                .distinct()
                                .collect(Collectors.toList());

                final LinkedHashMap<Allele, Double> remappedVqsLodMap = remapAllelesInMap(ref, nonRefAlleles, contig, (int) startPosition, vqsLodMap, Double.NaN);
                final LinkedHashMap<Allele, String> remappedYngMap = remapAllelesInMap(ref, nonRefAlleles, contig, (int) startPosition, yngMap, VCFConstants.EMPTY_INFO_FIELD);

                // see https://github.com/broadinstitute/dsp-spec-ops/issues/291 for rationale
                // take "worst" outcome for yng/vqslod, evaluate each allele separately
                // if any allele is "N"ay, the genotype is filtered
                // if any allele is "Y"ay and the rest are "G"rey, the genotype is passed
                // if all alleles are "G"ray, the VQSLod is evaluated
                boolean anyNays = nonRefAlleles.stream().map(a -> remappedYngMap.get(a)).anyMatch(v -> "N".equals(v));
                boolean anyYays = nonRefAlleles.stream().map(a -> remappedYngMap.get(a)).anyMatch(v -> "Y".equals(v));

                // if there are any "N"s, the genotype is filtered
                if (anyNays) {
                    genotypeBuilder.filter(GATKVCFConstants.NAY_FROM_YNG);
                } else if (anyYays) {
                    // the genotype is passed, nothing to do here as non-filtered is the default
                } else {
                    // get the max (best) vqslod for all SNP non-Yay sites, and apply the filter
                    Optional<Double> snpMax =
                            nonRefAlleles.stream().filter(a -> a.length() == ref.length()).map(a -> remappedVqsLodMap.get(a)).filter(Objects::nonNull).max(Double::compareTo);

                    if (snpMax.isPresent() && snpMax.get() < vqsLodSNPThreshold) {
                        genotypeBuilder.filter(GATKVCFConstants.VQSR_FAILURE_SNP);
                    }

                    // get the max (best) vqslod for all INDEL non-Yay sites
                    Optional<Double> indelMax =
                            nonRefAlleles.stream().filter(a -> a.length() != ref.length()).map(a -> remappedVqsLodMap.get(a)).filter(Objects::nonNull).max(Double::compareTo);

                    if (indelMax.isPresent() && indelMax.get() < vqsLodINDELThreshold) {
                        genotypeBuilder.filter(GATKVCFConstants.VQSR_FAILURE_INDEL);
                    }

                }
            }
            genotypeBuilder.alleles(genotypeAlleles);
        }

        final String callGQ = sampleRecord.getCallGQ();
        if ( callGQ != null ) {
            genotypeBuilder.GQ(Integer.parseInt(callGQ));
        }

        final String callPL = sampleRecord.getCallPL();
        if ( this.emitPLs && callPL != null ) {
            genotypeBuilder.PL(Arrays.stream(callPL.split(SchemaUtils.MULTIVALUE_FIELD_DELIMITER)).mapToInt(Integer::parseInt).toArray());
        }

        final String callAD = sampleRecord.getCallAD();
        if ( this.emitADs && callAD != null ) {
            genotypeBuilder.AD(Arrays.stream(callAD.split(SchemaUtils.MULTIVALUE_FIELD_DELIMITER)).mapToInt(Integer::parseInt).toArray());
        }

        final String callRGQ = sampleRecord.getCallRGQ();
        if ( callRGQ != null ) {
            genotypeBuilder.attribute(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, callRGQ);
        }
        builder.genotypes(genotypeBuilder.make());

        return builder.make();
    }

    private SortingCollection<GenericRecord> createSortedVetCollectionFromBigQuery(final String projectID,
                                                                                          final String fqDatasetName,
                                                                                          final Set<Long> sampleIdsToExtract,
                                                                                          final Long minLocation,
                                                                                          final Long maxLocation,
                                                                                          final int localSortMaxRecordsInRam,
                                                                                          final VariantBitSet vbs
                                                                                          ) {
        SortingCollection<GenericRecord> sortedVet = null;

        Map<Integer, LinkedList<Set<Long>>> tableMap = SampleList.mapSampleIdsToTableIndexes(sampleIdsToExtract);
        for (int tableIndex : tableMap.keySet() ) {
            TableReference vetTableRef =
                    new TableReference(fqDatasetName + ".vet_" + String.format("%03d", tableIndex), SchemaUtils.EXTRACT_VET_FIELDS);

            // TODO: Comment as to why (specific sample list clause)
            for (Set<Long> chunkSampleIds : tableMap.get(tableIndex)) {
                String sampleRestriction = " AND sample_id IN (" + StringUtils.join(chunkSampleIds, ",") + ")";

                // We need to look upstream MAX_DELETION_SIZE bases in case there is a deletion that begins before
                // the requested range, but spans into our processing range.  We don't use a "length" or end position
                // because it would break the clustering indexing
                final String vetRowRestriction =
                        "location >= " + (minLocation - IngestConstants.MAX_DELETION_SIZE + 1)+ " AND location <= " + maxLocation + sampleRestriction;
                try (StorageAPIAvroReader vetReader = new StorageAPIAvroReader(vetTableRef, vetRowRestriction, projectID)) {
                    if (sortedVet == null) {
                        sortedVet = getAvroSortingCollection(vetReader.getSchema(), localSortMaxRecordsInRam);
                    }

                    addToVetSortingCollection(sortedVet, vetReader, vbs);
                    processBytesScanned(vetReader);
                }
            }
        }
        return sortedVet;
    }

    private SortingCollection<GenericRecord> createSortedReferenceRangeCollectionFromBigQuery(final String projectID,
                                                                                                     final String fqDatasetName,
                                                                                                     final Set<Long> sampleIdsToExtract,
                                                                                                     final Long minLocation,
                                                                                                     final Long maxLocation,
                                                                                                     final int localSortMaxRecordsInRam,
                                                                                                     final VariantBitSet vbs
    ) {
        Map<Integer, LinkedList<Set<Long>>> tableMap = SampleList.mapSampleIdsToTableIndexes(sampleIdsToExtract);

        SortingCollection<GenericRecord> sortedReferenceRange = null;
        for (int tableIndex : tableMap.keySet() ) {
            TableReference refTableRef =
                    new TableReference(fqDatasetName + ".ref_ranges_" + String.format("%03d", tableIndex), SchemaUtils.EXTRACT_REF_FIELDS);

            for (Set<Long> chunkSampleIds : tableMap.get(tableIndex)) {
                String sampleRestriction = " AND sample_id IN (" + StringUtils.join(chunkSampleIds, ",") + ")";

                // NOTE: MUST be written as location >= minLocation - MAX_REFERENCE_BLOCK_BASES +1 to not break cluster pruning in BigQuery
                final String refRowRestriction =
                        "location >= " + (minLocation - IngestConstants.MAX_REFERENCE_BLOCK_BASES + 1) + " AND location <= " + maxLocation + sampleRestriction;

                try (StorageAPIAvroReader refReader = new StorageAPIAvroReader(refTableRef, refRowRestriction, projectID)) {
                    if (sortedReferenceRange == null) {
                        sortedReferenceRange = getAvroSortingCollection(refReader.getSchema(), localSortMaxRecordsInRam);
                    }
                    addToRefSortingCollection(sortedReferenceRange, refReader, vbs);
                    processBytesScanned(refReader);
                }
            }
        }

        return sortedReferenceRange;
    }



    private void createVariantsFromUnsortedBigQueryRanges(
            final String fqDatasetName,
            final SortedSet<Long> sampleIdsToExtract,
            final Long minLocation,
            final Long maxLocation,
            final HashMap<Long, HashMap<Allele, HashMap<Allele, Double>>> fullVqsLodMap,
            final HashMap<Long, HashMap<Allele, HashMap<Allele, String>>> fullYngMap,
            final HashMap<Long, List<String>> siteFilterMap,
            final boolean noVqslodFilteringRequested) {

        // We could handle this by making a map of BitSets or something, but it seems unnecessary to support this
        if (!SchemaUtils.decodeContig(minLocation).equals(SchemaUtils.decodeContig(maxLocation))) {
            throw new GATKException("Can not process cross-contig boundaries");
        }

        VariantBitSet vbs = new VariantBitSet(minLocation, maxLocation);

        SortingCollection<GenericRecord> sortedVet = createSortedVetCollectionFromBigQuery(projectID,
                fqDatasetName,
                sampleIdsToExtract,
                minLocation,
                maxLocation,
                localSortMaxRecordsInRam,
                vbs);


        SortingCollection<GenericRecord> sortedReferenceRange = createSortedReferenceRangeCollectionFromBigQuery(projectID,
                fqDatasetName,
                sampleIdsToExtract,
                minLocation,
                maxLocation,
                localSortMaxRecordsInRam,
                vbs);

        createVariantsFromSortedRanges(sampleIdsToExtract, sortedVet, sortedReferenceRange, fullVqsLodMap, fullYngMap, siteFilterMap, noVqslodFilteringRequested);
    }

    //
    // BEGIN REF RANGES COHORT EXTACT
    //
    private void createVariantsFromUnsortedExtractTableBigQueryRanges(
            final String fqVetTable,
            final String fqRefTable,
            final SortedSet<Long> sampleIdsToExtract,
            final Long minLocation,
            final Long maxLocation,
            final HashMap<Long, HashMap<Allele, HashMap<Allele, Double>>> fullVqsLodMap,
            final HashMap<Long, HashMap<Allele, HashMap<Allele, String>>> fullYngMap,
            final HashMap<Long, List<String>> siteFilterMap,
            final boolean noVqslodFilteringRequested) {

        // We could handle this by making a map of BitSets or something, but it seems unnecessary to support this
        if (!SchemaUtils.decodeContig(minLocation).equals(SchemaUtils.decodeContig(maxLocation))) {
            throw new GATKException("Can not process cross-contig boundaries");
        }

        VariantBitSet vbs = new VariantBitSet(minLocation, maxLocation);

        SortingCollection<GenericRecord> sortedVet = createSortedVetCollectionFromExtractTableBigQuery(projectID,
                fqVetTable,
                minLocation,
                maxLocation,
                localSortMaxRecordsInRam,
                vbs);


        SortingCollection<GenericRecord> sortedReferenceRange = createSortedReferenceRangeCollectionFromExtractTableBigQuery(projectID,
                fqRefTable,
                minLocation,
                maxLocation,
                localSortMaxRecordsInRam,
                vbs);

        createVariantsFromSortedRanges(sampleIdsToExtract, sortedVet, sortedReferenceRange, fullVqsLodMap, fullYngMap, siteFilterMap, noVqslodFilteringRequested);
    }

    private SortingCollection<GenericRecord> createSortedVetCollectionFromExtractTableBigQuery(final String projectID,
                                                                                   final String fqVetTable,
                                                                                   final Long minLocation,
                                                                                   final Long maxLocation,
                                                                                   final int localSortMaxRecordsInRam,
                                                                                   final VariantBitSet vbs
    ) {

        TableReference tableRef =
                new TableReference(fqVetTable, SchemaUtils.EXTRACT_VET_FIELDS);

        // We need to look upstream MAX_DELETION_SIZE bases in case there is a deletion that begins before
        // the requested range, but spans into our processing range.  We don't use a "length" or end position
        // because it would break the clustering indexing
        final String vetRowRestriction =
                "location >= " + (minLocation - IngestConstants.MAX_DELETION_SIZE + 1)+ " AND location <= " + maxLocation;
        try (StorageAPIAvroReader vetReader = new StorageAPIAvroReader(tableRef, vetRowRestriction, projectID)) {
            SortingCollection<GenericRecord> sortedVet = getAvroSortingCollection(vetReader.getSchema(), localSortMaxRecordsInRam);
            addToVetSortingCollection(sortedVet, vetReader, vbs);
            processBytesScanned(vetReader);
            return sortedVet;
        }
    }

    private SortingCollection<GenericRecord> createSortedReferenceRangeCollectionFromExtractTableBigQuery(final String projectID,
                                                                                              final String fqRefTable,
                                                                                              final Long minLocation,
                                                                                              final Long maxLocation,
                                                                                              final int localSortMaxRecordsInRam,
                                                                                              final VariantBitSet vbs
    ) {

        TableReference tableRef = new TableReference(fqRefTable, SchemaUtils.EXTRACT_REF_FIELDS);

        // NOTE: MUST be written as location >= minLocation - MAX_REFERENCE_BLOCK_BASES +1 to not break cluster pruning in BigQuery
        final String refRowRestriction =
                "location >= " + (minLocation - IngestConstants.MAX_REFERENCE_BLOCK_BASES + 1) + " AND location <= " + maxLocation;

        try (StorageAPIAvroReader refReader = new StorageAPIAvroReader(tableRef, refRowRestriction, projectID)) {
            SortingCollection<GenericRecord> sortedReferenceRange = getAvroSortingCollection(refReader.getSchema(), localSortMaxRecordsInRam);
            addToRefSortingCollection(sortedReferenceRange, refReader, vbs);
            processBytesScanned(refReader);
            return sortedReferenceRange;
        }
    }

    //
    // END REF RANGES COHORT EXTRACT
    //
    private void createVariantsFromUnsortedAvroRanges(
            final GATKPath vetAvroFileName,
            final GATKPath refRangesAvroFileName,
            final SortedSet<Long> sampleIdsToExtract,
            final Long minLocation,
            final Long maxLocation,
            final HashMap<Long, HashMap<Allele, HashMap<Allele, Double>>> fullVqsLodMap,
            final HashMap<Long, HashMap<Allele, HashMap<Allele, String>>> fullYngMap,
            final HashMap<Long, List<String>> siteFilterMap,
            final boolean noVqslodFilteringRequested,
            final boolean presortedAvroFiles) {

        final AvroFileReader vetReader = new AvroFileReader(vetAvroFileName);
        final AvroFileReader refRangesReader = new AvroFileReader(refRangesAvroFileName);

        Iterable<GenericRecord> sortedVet;
        Iterable<GenericRecord> sortedReferenceRange;

        if (presortedAvroFiles) {
            sortedVet = vetReader;
            sortedReferenceRange = refRangesReader;
        } else {
            VariantBitSet vbs = new VariantBitSet(minLocation, maxLocation);

            SortingCollection<GenericRecord> localSortedVet = getAvroSortingCollection(vetReader.getSchema(), localSortMaxRecordsInRam);
            addToVetSortingCollection(localSortedVet, vetReader, vbs);

            SortingCollection<GenericRecord> localSortedReferenceRange = getAvroSortingCollection(refRangesReader.getSchema(), localSortMaxRecordsInRam);
            addToRefSortingCollection(localSortedReferenceRange, refRangesReader, vbs);

            sortedVet = localSortedVet;
            sortedReferenceRange = localSortedReferenceRange;
        }

        createVariantsFromSortedRanges(sampleIdsToExtract, sortedVet, sortedReferenceRange, fullVqsLodMap, fullYngMap, siteFilterMap, noVqslodFilteringRequested);

    }


    void createVariantsFromSortedRanges(final SortedSet<Long> sampleIdsToExtract,
                                        final Iterable<GenericRecord> sortedVet,
                                        Iterable<GenericRecord> sortedReferenceRange,
                                        final HashMap<Long, HashMap<Allele, HashMap<Allele, Double>>> fullVqsLodMap,
                                        final HashMap<Long, HashMap<Allele, HashMap<Allele, String>>> fullYngMap,
                                        final HashMap<Long, List<String>> siteFilterMap,
                                        final boolean noVqslodFilteringRequested) {
        createVariantsFromSortedRanges(sampleIdsToExtract, sortedVet, sortedReferenceRange, fullVqsLodMap, fullYngMap, siteFilterMap, noVqslodFilteringRequested, (vc) -> writeToVcf(vc));
    }

    void createVariantsFromSortedRanges(final SortedSet<Long> sampleIdsToExtract,
                                        final Iterable<GenericRecord> sortedVet,
                                        Iterable<GenericRecord> sortedReferenceRange,
                                        final HashMap<Long, HashMap<Allele, HashMap<Allele, Double>>> fullVqsLodMap,
                                        final HashMap<Long, HashMap<Allele, HashMap<Allele, String>>> fullYngMap,
                                        final HashMap<Long, List<String>> siteFilterMap,
                                        final boolean noVqslodFilteringRequested,
                                        Consumer<VariantContext> variantContextConsumer
    ) {

        long maxSampleId = sampleIdsToExtract.stream().max(Long::compare).orElseThrow(
                () -> new GATKException("Unable to calculate max sample id, sample list may be empty")
        );

        long minSampleId = sampleIdsToExtract.stream().min(Long::compare).orElseThrow(
                () -> new GATKException("Unable to calculate min sample id, sample list may be empty")
        );

        Long lastSample = null;
        Long lastPosition = null;

        final Map<Long, ExtractCohortRecord> currentPositionRecords = new HashMap<>(sampleIdToName.size() * 2);
        final Map<Long, TreeSet<ReferenceRecord>> referenceCache = new HashMap<>(sampleIdToName.size());

        // Initialize cache
        for ( Long sampleId: sampleIdsToExtract) {
            referenceCache.put(sampleId, new TreeSet<>());
        }

        // NOTE: if OverlapDetector takes too long, try using RegionChecker from tws_sv_local_assembler
        final OverlapDetector<SimpleInterval> intervalsOverlapDetector = OverlapDetector.create(traversalIntervals);

        Iterator<GenericRecord> sortedReferenceRangeIterator = sortedReferenceRange.iterator();

        for (final GenericRecord sortedRow : sortedVet) {
            final ExtractCohortRecord vetRow = new ExtractCohortRecord(sortedRow);
            long variantLocation = vetRow.getLocation();
            long variantSample = vetRow.getSampleId();

            // it's possible this variant is actually an upstream deletion before our region of interest,
            // so make sure this is not the case before we the actual variant record and discard
            if (variantLocation < minLocation) {
                handlePotentialSpanningDeletion(vetRow, referenceCache);
                continue;
            }

            // new position, fill in remainder of last position data
            // before continuing in to new position
            if (lastPosition != null && variantLocation != lastPosition) {
//                logger.info("skipped to new position " + v_position + " from " + lastPosition + " with last sample of " + lastSample);

                // if the last VET was the last sample we 're done
                if (lastSample != maxSampleId) {
                    processReferenceData(currentPositionRecords, sortedReferenceRangeIterator, referenceCache, lastPosition, lastSample + 1, maxSampleId, sampleIdsToExtract);
                }

                ++totalNumberOfSites;
                VariantContext vc = processSampleRecordsForLocation(lastPosition, currentPositionRecords.values(), fullVqsLodMap, fullYngMap, noVqslodFilteringRequested, siteFilterMap, VQSLODFilteringType);
                variantContextConsumer.accept(vc);
                currentPositionRecords.clear();

                lastSample = null;
            }

            long startingSample = (lastSample == null) ? minSampleId : lastSample + 1;
            processReferenceData(currentPositionRecords, sortedReferenceRangeIterator, referenceCache, variantLocation, startingSample, variantSample - 1, sampleIdsToExtract);

            // handle the actual variant record
            currentPositionRecords.merge(variantSample, vetRow, this::mergeSampleRecord);

            // if the variant record was a deletion, fabricate a spanning deletion row(s) for the cache
            // so that a future request for reference state at a position underlying the deletion is
            // properly handled
            handlePotentialSpanningDeletion(vetRow, referenceCache);

            lastPosition = variantLocation;
            lastSample = variantSample;
        }

        // finish writing out reference rows
        if (lastSample != maxSampleId) {
            processReferenceData(currentPositionRecords, sortedReferenceRangeIterator, referenceCache, lastPosition, lastSample + 1, maxSampleId, sampleIdsToExtract);
        }

        if (!currentPositionRecords.isEmpty()) {
            ++totalNumberOfSites;
            VariantContext vc = processSampleRecordsForLocation(lastPosition, currentPositionRecords.values(), fullVqsLodMap, fullYngMap, noVqslodFilteringRequested, siteFilterMap, VQSLODFilteringType);
            variantContextConsumer.accept(vc);
        }
    }

    private void handlePotentialSpanningDeletion(ExtractCohortRecord vetRow, Map<Long, TreeSet<ReferenceRecord>> referenceCache) {
        long position = vetRow.getLocation();
        long sample = vetRow.getSampleId();

        int refLength = vetRow.getRefAllele().length();

        List<Integer> orderedDeletionLengths =
                Arrays.stream(vetRow.getAltAllele().split(","))
                        .filter(x -> !x.equals(Allele.NO_CALL_STRING))
                        .map(x -> refLength - x.length())
                        .filter(x -> x > 0)
                        .sorted()
                        .collect(Collectors.toList());

        if (orderedDeletionLengths.size() > 0) {

            // if there is a single deletion, just encode it with the right zygosity
            if (orderedDeletionLengths.size() == 1) {
                boolean isHet
                        = Arrays.stream(vetRow.getCallGT().replace("|", "/").split("/")).filter(x -> !x.equals(Allele.NO_CALL_STRING)).distinct().count() == 2;

                mergeSpanningDeletionIntoCache(referenceCache.get(sample),
                        new SpanningDeletionRecord(position + 1, sample, orderedDeletionLengths.get(0), isHet ? "0/1" : "1/1", vetRow.getCallGQ())
                );
            } else if (orderedDeletionLengths.size() == 2) {
                // we need to add two different spanning deletions, one for the hom region and one for the het region
                int shorter = orderedDeletionLengths.get(0); // 5
                int longer = orderedDeletionLengths.get(1);  // 5
                int lengthHomRegion = shorter;
                int lengthHetRegion = longer - shorter;

                mergeSpanningDeletionIntoCache(referenceCache.get(sample),
                        new SpanningDeletionRecord(position + 1, sample, lengthHomRegion, "1/1", vetRow.getCallGQ())
                );

                // in the case of a single homozygous deletion, the length of the het region is zero so don't add it
                if (lengthHetRegion > 0) {
                    mergeSpanningDeletionIntoCache(referenceCache.get(sample),
                            new SpanningDeletionRecord(position + lengthHomRegion + 1, sample, lengthHetRegion, "0/1", vetRow.getCallGQ())
                    );
                }
            } else {
                // there are only three options here -- no deletions (do nothing) one deletion or two deletions
                throw new GATKException("Found > 2 different deletions in a single het genotype!!" + vetRow);
            }
        }
    }

    private void mergeSpanningDeletionIntoCache(TreeSet<ReferenceRecord> cache, SpanningDeletionRecord sdr) {
        // extract out any ReferenceRecords that overlap this record
        List<ReferenceRecord> recs = cache.stream().filter(x -> x.overlaps(sdr)).collect(Collectors.toList());
        // make most common case efficient
        if (recs.size() == 0) {
            cache.add(sdr);
        } else if (recs.size() == 1) {
            if (! (recs.get(0) instanceof SpanningDeletionRecord)) {
                ReferenceRecord existing = recs.get(0);

                // if this is an upstream reference record, just truncate it with the spanning deletion record
                if (existing.getStart() < sdr.getStart()) {
                    cache.remove(existing);
                    cache.add(new ReferenceRecord(existing.getLocation(), existing.getSampleId(), (int) (sdr.getLocation() - sdr.getLocation()), existing.getState()));
                } else {
                    throw new GATKException("Merging SpanningDeletionRecord " + sdr + " but found existing downstream non-SpanningDeletionRecord " + recs.get(0));
                }
            } else {

                // There are two kinds of overlap we have to deal with.  Essentially we want the later record to override the
                // earlier record in a region of overlap.
                //
                // (1) simple overlap -- intervals overlap at one end.  In this case we need to produce TWO records by just
                // truncating the earlier one so they don't overlap
                // (2) interior overlap -- one interval is completely inside the other.  In this case we need to create THREE
                // records.  The first two records are the same as the simple overlap case, but then we need to generate a
                // third record where the earlier records "resumes" after the end of the later record
                SpanningDeletionRecord existing = (SpanningDeletionRecord) recs.get(0);

                SpanningDeletionRecord earlier = (existing.getLocation() < sdr.getLocation()) ? existing : sdr;
                SpanningDeletionRecord later = (existing.getLocation() >= sdr.getLocation()) ? existing : sdr;

                // remove the existing row
                cache.remove(existing);

                // and replace it with a truncated earlier record, and the later record
                SpanningDeletionRecord s1 = new SpanningDeletionRecord(earlier.getLocation(), earlier.getSampleId(), (int) (later.getLocation() - earlier.getLocation()), earlier.getGT(), earlier.getGQ());
                cache.add(s1);
                cache.add(later);

                // and if the later record ends before the end of the earlier record emit the remainder (case 2 above)
                if (later.getEnd() < earlier.getEnd()) {
                    long newStart = later.getEndLocation() + 1;
                    int newLength = (int) (earlier.getEndLocation() - newStart + 1);
                    SpanningDeletionRecord s2 = new SpanningDeletionRecord(newStart, earlier.getSampleId(), newLength, earlier.getGT(), earlier.getGQ());
                    cache.add(s2);
                }
            }
        } else {
            throw new GATKException("Found > 2 reference records that overlap for a diploid sample! " + recs + " and " + sdr);
        }

    }
    private void processReferenceData(Map<Long, ExtractCohortRecord> currentPositionRecords, Iterator<GenericRecord> sortedReferenceRangeIterator, Map<Long, TreeSet<ReferenceRecord>> referenceCache, long location, long fromSampleId, long toSampleId, SortedSet<Long> sampleIdsToExtract) {
        // in the case where there are two adjacent samples with variants, this method is called where from is greater than to
        // this is ok, there is just no reference data to process but subSet will throw an exception so we handle it with this if block
        if (toSampleId >= fromSampleId) {
            SortedSet<Long> samples = sampleIdsToExtract.subSet(fromSampleId, toSampleId + 1); // subset is start-inclusive, end-exclusive

            for (Long s : samples) {
                ExtractCohortRecord e = processReferenceData(sortedReferenceRangeIterator, referenceCache, location, s);
                currentPositionRecords.merge(s, e, this::mergeSampleRecord);
            }
        }
    }

    private ExtractCohortRecord processReferenceData(Iterator<GenericRecord> sortedReferenceRangeIterator, Map<Long, TreeSet<ReferenceRecord>> referenceCache, long location, long sampleId) {
        ReferenceRecord r = processReferenceDataFromCache(referenceCache, location, sampleId);

        if (r == null) {
            r = processReferenceDataFromStream(sortedReferenceRangeIterator, referenceCache, location, sampleId);
        }

        // handle spanning deletion record
        if (r instanceof SpanningDeletionRecord) {
            SpanningDeletionRecord sdr = (SpanningDeletionRecord) r;
            Allele refAllele = getReferenceAllele(refSource, location, 1); // note this is location of the query/parameter, not the start of the deletion
            return new ExtractCohortRecord(location, sampleId, sdr.getState(), refAllele.getBaseString(), Allele.SPAN_DEL_STRING, sdr.getGT(), sdr.getGQ());
        } else {
            return new ExtractCohortRecord(location, sampleId, r.getState());
        }
    }

    // TODO: use this as a return record, singleton
    // TODO: so much hackery to instantiate this thing...
    public class InferredReferenceRecord extends ReferenceRecord {
        public InferredReferenceRecord(GQStateEnum inferredState) {
            super(SchemaUtils.chromAdjustment + 1, -1, -1, inferredState.getValue());
        }
    }

    private ReferenceRecord processReferenceDataFromCache(Map<Long, TreeSet<ReferenceRecord>> referenceCache, long location, long sampleId) {

        // PERF: bottleneck (iterator and remove) 6%
        Iterator<ReferenceRecord> iter = referenceCache.get(sampleId).iterator();
        while (iter.hasNext()) {
            ReferenceRecord row = iter.next();

            // if it ends before the current location, it's no longer relevant
            if (row.getEndLocation() < location) {
                iter.remove();
                continue;
            }

            // if it overlaps our location, use it!  Don't worry about removing it, if it
            // becomes irrelevant the next time we process it we will remove it
            if (row.getLocation() <= location && row.getEndLocation() >= location) {
                return row;
            }

            // completely after position, inferred state
            if (row.getLocation() > location) {
                return inferredReferenceRecord;
            }
        }

        // still here means QUEUE should be empty and we return a null state to indicate
        // no results were found.
        if (!referenceCache.get(sampleId).isEmpty()) {
            throw new GATKException("should have an empty queue: " + referenceCache.get(sampleId));
        } else {
            return null;
        }
    }

    private ReferenceRecord processReferenceDataFromStream(Iterator<GenericRecord> sortedReferenceRangeIterator, Map<Long, TreeSet<ReferenceRecord>> referenceCache, long location, long sampleId) {
        while(sortedReferenceRangeIterator.hasNext()) {
            final ReferenceRecord refRow = new ReferenceRecord(sortedReferenceRangeIterator.next());
            totalRangeRecords++;

            // skip irrelevant entries
            if (refRow.getEndLocation() < location) {
                totalIrrelevantRangeRecords++;
                continue;
            }

            // is overlapping
            if (refRow.getLocation() <= location) {

                // always add to the appropriate cache for use downstream
                //mergeIntoReferenceCache(referenceCache.get(refRow.getSampleId()), refRow);
                referenceCache.get(refRow.getSampleId()).add(refRow);

                // if this is for the requested sample, return the value
                if (refRow.getSampleId() == sampleId) {
                    return refRow;
                }
            }

            // we are now past the position, put this one entry in the cache and return the inferred state
            if (refRow.getLocation() > location) {
                //mergeIntoReferenceCache(referenceCache.get(refRow.getSampleId()), refRow);
                referenceCache.get(refRow.getSampleId()).add(refRow);
                return inferredReferenceRecord;
            }
        }

        // if we are still here... use the inferred state
        return inferredReferenceRecord;
    }

    // Turns out some data has overlapping reference blocks.  Technically this is invalid,
    // but it still exists in the wild and in WARP scientific test data.
    private void mergeIntoReferenceCache(TreeSet<ReferenceRecord> cache, ReferenceRecord rec) {
        if (cache.size() == 0) {
            cache.add(rec);
        } else {

            // iterate through cache and look for an overlap
            // NOTE: it might be possible to generalize this logic with similar logic
            // from merging in SpanningDeletionRecords but keeping it separate for now
            for (ReferenceRecord existing : cache) {
                if (existing.overlaps(rec) && !(existing instanceof SpanningDeletionRecord)) {
                    ReferenceRecord earlier = (existing.getLocation() < rec.getLocation()) ? existing : rec;
                    ReferenceRecord later = (existing.getLocation() >= rec.getLocation()) ? existing : rec;

                    // remove the existing row
                    cache.remove(existing);

                    // and replace it with a truncated earlier record, and the later record
                    ReferenceRecord s1 = new ReferenceRecord(earlier.getLocation(), earlier.getSampleId(), (int) (later.getLocation() - earlier.getLocation()), earlier.getState());
                    cache.add(s1);
                    cache.add(later);
                }
            }

        }
    }

    private static Allele getReferenceAllele(ReferenceDataSource refSource, long location, int length) {
        String contig = SchemaUtils.decodeContig(location);
        long position = SchemaUtils.decodePosition(location);
        return Allele.create(refSource.queryAndPrefetch(contig, position, position + length - 1).getBaseString(), true);
    }
}
