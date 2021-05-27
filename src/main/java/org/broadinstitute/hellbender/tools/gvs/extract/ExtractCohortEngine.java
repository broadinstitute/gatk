package org.broadinstitute.hellbender.tools.gvs.extract;

import com.google.common.collect.Sets;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.avro.generic.GenericRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.gnarlyGenotyper.GnarlyGenotyperEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.*;
import org.broadinstitute.hellbender.utils.localsort.AvroSortingCollectionCodec;
import org.broadinstitute.hellbender.utils.localsort.SortingCollection;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.util.*;
import java.util.stream.Collectors;

public class ExtractCohortEngine {
    private static final Logger logger = LogManager.getLogger(ExtractCohortEngine.class);

    private final VariantContextWriter vcfWriter;

    private final boolean printDebugInformation;
    private final int localSortMaxRecordsInRam;
    private final TableReference cohortTableRef;
    private final List<SimpleInterval> traversalIntervals;
    private final Long minLocation;
    private final Long maxLocation;
    private final TableReference filterSetInfoTableRef;
    private final TableReference filterSetSiteTableRef;
    private final ReferenceDataSource refSource;
    private Double vqsLodSNPThreshold;
    private Double vqsLodINDELThreshold;
    private boolean performGenotypeVQSLODFiltering;
    private boolean excludeFilteredSites;

    private final ProgressMeter progressMeter;
    private final String projectID;
    private final CommonCode.ModeEnum mode;

    /** List of sample names seen in the variant data from BigQuery. */
    private Set<String> sampleNames;
    private final ReferenceConfidenceVariantContextMerger variantContextMerger;
    private final boolean disableGnarlyGenotyper;
    private final GnarlyGenotyperEngine gnarlyGenotyper;

    private int totalNumberOfVariants = 0;
    private int totalNumberOfSites = 0;

    private final String cohortAvroFileName;
    private final String filterSetName;

    /**
     * The conf threshold above which variants are not included in the position tables.
     * This value is used to construct the genotype information of those missing samples
     * when they are merged together into a {@link VariantContext} object
     */
    public static int MISSING_CONF_THRESHOLD = 60;


    public ExtractCohortEngine(final String projectID,
                               final VariantContextWriter vcfWriter,
                               final VCFHeader vcfHeader,
                               final VariantAnnotatorEngine annotationEngine,
                               final ReferenceDataSource refSource,
                               final Set<String> sampleNames,
                               final CommonCode.ModeEnum mode,
                               final String cohortTableName,
                               final String cohortAvroFileName,
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
                               final boolean disableGnarlyGenotyper,
                               final boolean performGenotypeVQSLODFiltering,
                               final boolean excludeFilteredSites
    ) {
        this.localSortMaxRecordsInRam = localSortMaxRecordsInRam;

        this.projectID = projectID;
        this.vcfWriter = vcfWriter;
        this.refSource = refSource;
        this.sampleNames = sampleNames;
        this.mode = mode;

        this.cohortTableRef = cohortTableName == null || "".equals(cohortTableName) ? null :
                new TableReference(cohortTableName, emitPLs ? SchemaUtils.COHORT_FIELDS : SchemaUtils.COHORT_FIELDS_NO_PL);

        this.cohortAvroFileName = cohortAvroFileName;
        this.traversalIntervals = traversalIntervals;
        this.minLocation = minLocation;
        this.maxLocation = maxLocation;
        this.filterSetInfoTableRef = filterSetInfoTableName == null || "".equals(filterSetInfoTableName) ? null : new TableReference(filterSetInfoTableName, SchemaUtils.YNG_FIELDS);
        this.filterSetSiteTableRef = filterSetSiteTableName == null || "".equals(filterSetSiteTableName) ? null : new TableReference(filterSetSiteTableName, SchemaUtils.FILTER_SET_SITE_FIELDS);

        this.printDebugInformation = printDebugInformation;
        this.vqsLodSNPThreshold = vqsLodSNPThreshold;
        this.vqsLodINDELThreshold = vqsLodINDELThreshold;
        this.performGenotypeVQSLODFiltering = performGenotypeVQSLODFiltering;
        this.excludeFilteredSites = excludeFilteredSites;

        this.progressMeter = progressMeter;

        this.filterSetName = filterSetName;

        this.variantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine, vcfHeader);
        this.disableGnarlyGenotyper = disableGnarlyGenotyper;
        this.gnarlyGenotyper = disableGnarlyGenotyper ? null : new GnarlyGenotyperEngine(false, 30, false, emitPLs, true);

    }

    private final static double INDEL_QUAL_THRESHOLD = GenotypeCalculationArgumentCollection.DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING - 10 * Math.log10(HomoSapiensConstants.INDEL_HETEROZYGOSITY);
    private final static double SNP_QUAL_THRESHOLD = GenotypeCalculationArgumentCollection.DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING - 10 * Math.log10(HomoSapiensConstants.SNP_HETEROZYGOSITY);

    int getTotalNumberOfVariants() { return totalNumberOfVariants; }
    int getTotalNumberOfSites() { return totalNumberOfSites; }

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

        boolean noVqslodFilteringRequested = (filterSetInfoTableRef == null);
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
            }
        }

        if (printDebugInformation) {
            logger.debug("using storage api with local sort");
        }
        logger.debug("Initializing Reader");
        if (cohortTableRef != null){
            final StorageAPIAvroReader storageAPIAvroReader = new StorageAPIAvroReader(cohortTableRef, rowRestriction, projectID);
            createVariantsFromUnsortedResult(storageAPIAvroReader, fullVqsLodMap, fullYngMap, siteFilterMap, noVqslodFilteringRequested);
        }
        else {
            final AvroFileReader avroFileReader = new AvroFileReader(cohortAvroFileName);
            createVariantsFromUnsortedResult(avroFileReader, fullVqsLodMap, fullYngMap, siteFilterMap, noVqslodFilteringRequested);
        }
        logger.debug("Finished Initializing Reader");
    }


    public SortingCollection<GenericRecord> getAvroSortingCollection(org.apache.avro.Schema schema, int localSortMaxRecordsInRam) {
        final SortingCollection.Codec<GenericRecord> sortingCollectionCodec = new AvroSortingCollectionCodec(schema);
        final Comparator<GenericRecord> sortingCollectionComparator = new Comparator<GenericRecord>() {
            @Override
            public int compare( GenericRecord o1, GenericRecord o2 ) {
                final long firstPosition = Long.parseLong(o1.get(SchemaUtils.LOCATION_FIELD_NAME).toString());
                final long secondPosition = Long.parseLong(o2.get(SchemaUtils.LOCATION_FIELD_NAME).toString());

                return Long.compare(firstPosition, secondPosition);
            }
        };
        return SortingCollection.newInstance(GenericRecord.class, sortingCollectionCodec, sortingCollectionComparator, localSortMaxRecordsInRam, true);
    }


    private void createVariantsFromUnsortedResult(final GATKAvroReader avroReader,
                                                  final HashMap<Long, HashMap<Allele, HashMap<Allele, Double>>> fullVqsLodMap,
                                                  final HashMap<Long, HashMap<Allele, HashMap<Allele, String>>> fullYngMap,
                                                  final HashMap<Long, List<String>> siteFilterMap,
                                                  final boolean noVqslodFilteringRequested) {

        final org.apache.avro.Schema schema = avroReader.getSchema();

        SortingCollection<GenericRecord> sortingCollection =  getAvroSortingCollection(schema, localSortMaxRecordsInRam);

        int recordsProcessed = 0;
        long startTime = System.currentTimeMillis();

        for ( final GenericRecord queryRow : avroReader ) {

            sortingCollection.add(queryRow);
            if (recordsProcessed++ % 1000000 == 0) {
                long endTime = System.currentTimeMillis();
                logger.info("Processed " + recordsProcessed + " from BigQuery Read API in " + (endTime-startTime) + " ms");
                startTime = endTime;
            }
        }

        sortingCollection.printTempFileStats();

        final Map<String, ExtractCohortRecord> currentPositionRecords = new HashMap<>(sampleNames.size() * 2);

        long currentLocation = -1;

        // NOTE: if OverlapDetector takes too long, try using RegionChecker from tws_sv_local_assembler
        final OverlapDetector<SimpleInterval> intervalsOverlapDetector = OverlapDetector.create(traversalIntervals);

        for ( final GenericRecord sortedRow : sortingCollection ) {
            final ExtractCohortRecord cohortRow = new ExtractCohortRecord( sortedRow );

            if ( intervalsOverlapDetector.overlapsAny(cohortRow) ) {
                final long location = cohortRow.getLocation();
                final String sampleName = cohortRow.getSampleName();

                if (location != currentLocation && currentLocation != -1) {
                    ++totalNumberOfSites;
                    processSampleRecordsForLocation(currentLocation, currentPositionRecords.values(), fullVqsLodMap, fullYngMap, noVqslodFilteringRequested, siteFilterMap, performGenotypeVQSLODFiltering);

                    currentPositionRecords.clear();
                }

                currentPositionRecords.merge(sampleName, cohortRow, this::mergeSampleRecord);
                currentLocation = location;
            }
        }

        if ( ! currentPositionRecords.isEmpty() ) {
            ++totalNumberOfSites;
            processSampleRecordsForLocation(currentLocation, currentPositionRecords.values(), fullVqsLodMap, fullYngMap, noVqslodFilteringRequested, siteFilterMap, performGenotypeVQSLODFiltering);
        }
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

    private void processSampleRecordsForLocation(final long location,
                                                 final Iterable<ExtractCohortRecord> sampleRecordsAtPosition,
                                                 final HashMap<Long, HashMap<Allele, HashMap<Allele, Double>>> fullVqsLodMap,
                                                 final HashMap<Long, HashMap<Allele, HashMap<Allele, String>>> fullYngMap,
                                                 final boolean noVqslodFilteringRequested,
                                                 final HashMap<Long, List<String>> siteFilterMap,
                                                 final boolean performGenotypeVQSLODFiltering) {

        final List<VariantContext> unmergedCalls = new ArrayList<>();
        final Set<String> currentPositionSamplesSeen = new HashSet<>();
        boolean currentPositionHasVariant = false;
        final int currentPosition = SchemaUtils.decodePosition(location);
        final String contig = SchemaUtils.decodeContig(location);
        final Allele refAllele = Allele.create(refSource.queryAndPrefetch(contig, currentPosition, currentPosition).getBaseString(), true);
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
            final String sampleName = sampleRecord.getSampleName();
            currentPositionSamplesSeen.add(sampleName);
            ++numRecordsAtPosition;

            if ( printDebugInformation ) {
                logger.info("\t" + contig + ":" + currentPosition + ": found record for sample " + sampleName + ": " + sampleRecord);
            }

            switch (sampleRecord.getState()) {
                case "v":   // Variant
                    ++totalNumberOfVariants;
                    VariantContext vc = createVariantContextFromSampleRecord(sampleRecord, vqsLodMap, yngMap, performGenotypeVQSLODFiltering);
                    unmergedCalls.add(vc);

                    currentPositionHasVariant = true;
                    break;
                case "0":   // Non Variant Block with GQ < 10
                    // Reference calls with GQ 0 should be rendered as no-call (#271)
                    // Nothing to do here -- just needed to mark the sample as seen so it doesn't get put in the high confidence ref band
                    break;
                case "1":  // Non Variant Block with 10 <=  GQ < 20
                    unmergedCalls.add(createRefSiteVariantContextWithGQ(sampleName, contig, currentPosition, refAllele, 10));
                    break;
                case "2":  // Non Variant Block with 20 <= GQ < 30
                    unmergedCalls.add(createRefSiteVariantContextWithGQ(sampleName, contig, currentPosition, refAllele, 20));
                    break;
                case "3":  // Non Variant Block with 30 <= GQ < 40
                    unmergedCalls.add(createRefSiteVariantContextWithGQ(sampleName, contig, currentPosition, refAllele, 30));
                    break;
                case "4":  // Non Variant Block with 40 <= GQ < 50
                    unmergedCalls.add(createRefSiteVariantContextWithGQ(sampleName, contig, currentPosition, refAllele, 40));
                    break;
                case "5":  // Non Variant Block with 50 <= GQ < 60
                    unmergedCalls.add(createRefSiteVariantContextWithGQ(sampleName, contig, currentPosition, refAllele, 50));
                    break;
                case "6":  // Non Variant Block with 60 <= GQ (usually omitted from tables)
                    unmergedCalls.add(createRefSiteVariantContextWithGQ(sampleName, contig, currentPosition, refAllele, 60));
                    break;
                case "*":   // Spanning Deletion - do nothing. just mark the sample as seen
                    break;
                case "m":   // Missing
                    // Nothing to do here -- just needed to mark the sample as seen so it doesn't get put in the high confidence ref band
                    break;
                case "u":   // unknown GQ used for array data
                    unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele));
                    break;
                default:
                    throw new GATKException("Unrecognized state: " + sampleRecord.getState());
            }

        }

        if ( printDebugInformation ) {
            logger.info(contig + ":" + currentPosition + ": processed " + numRecordsAtPosition + " total sample records");
        }

        finalizeCurrentVariant(unmergedCalls, currentPositionSamplesSeen, currentPositionHasVariant, location, contig, currentPosition, refAllele, vqsLodMap, yngMap, noVqslodFilteringRequested, siteFilterMap);
    }

    private void finalizeCurrentVariant(final List<VariantContext> unmergedCalls,
                                        final Set<String> currentVariantSamplesSeen,
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
            return;
        }

        // Find samples for dropped state and synthesize. If not arrays use GQ 60
        final Set<String> samplesNotEncountered = Sets.difference(sampleNames, currentVariantSamplesSeen);
        for ( final String missingSample : samplesNotEncountered ) {
            unmergedCalls.add(createRefSiteVariantContextWithGQ(missingSample, contig, start, refAllele, MISSING_CONF_THRESHOLD));
        }

        // we only need to retain the NonRefSymbolicAllele if we are using gnarly
        boolean removeNonRefSymbolicAllele = this.disableGnarlyGenotyper;

        final VariantContext mergedVC = variantContextMerger.merge(
                unmergedCalls,
                new SimpleInterval(contig, (int) start, (int) start),
                refAllele.getBases()[0],
                removeNonRefSymbolicAllele,
                false,
                true);

        final VariantContext genotypedVC = this.disableGnarlyGenotyper ? mergedVC : gnarlyGenotyper.finalizeGenotype(mergedVC);

        // Gnarly will indicate dropping a site by returning a null from the finalizeGenotype method
        if (!this.disableGnarlyGenotyper && genotypedVC == null) {
            return;
        }

        // apply VQSLod-based filters
        VariantContext filteredVC =
                noVqslodFilteringRequested ? genotypedVC : filterSiteByAlleleSpecificVQSLOD(genotypedVC, vqsLodMap, yngMap, performGenotypeVQSLODFiltering);

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

        if ( finalVC != null ) {
            // Add the variant contexts that aren't filtered or add everything if we aren't excluding anything
            if (finalVC.isNotFiltered() || !excludeFilteredSites) {
                vcfWriter.add(finalVC);
            }
            progressMeter.update(finalVC);
        }
    }

    private VariantContext filterSiteByAlleleSpecificVQSLOD(VariantContext mergedVC, HashMap<Allele, HashMap<Allele, Double>> vqsLodMap, HashMap<Allele, HashMap<Allele, String>> yngMap, boolean onlyAnnotate) {
        final LinkedHashMap<Allele, Double> remappedVqsLodMap = remapAllelesInMap(mergedVC, vqsLodMap, Double.NaN);
        final LinkedHashMap<Allele, String> remappedYngMap = remapAllelesInMap(mergedVC, yngMap, VCFConstants.EMPTY_INFO_FIELD);

        final LinkedHashMap<Allele, Double> relevantVqsLodMap = new LinkedHashMap<>();
        mergedVC.getAlternateAlleles().forEach(key -> Optional.ofNullable(remappedVqsLodMap.get(key)).ifPresent(value -> relevantVqsLodMap.put(key, value)));
        final LinkedHashMap<Allele, String> relevantYngMap = new LinkedHashMap<>();
        mergedVC.getAlternateAlleles().forEach(key -> Optional.ofNullable(remappedYngMap.get(key)).ifPresent(value -> relevantYngMap.put(key, value)));

        final VariantContextBuilder builder = new VariantContextBuilder(mergedVC);

        builder.attribute(GATKVCFConstants.AS_VQS_LOD_KEY, relevantVqsLodMap.values().stream().map(val -> val.equals(Double.NaN) ? VCFConstants.EMPTY_INFO_FIELD : val.toString()).collect(Collectors.toList()));
        builder.attribute(GATKVCFConstants.AS_YNG_STATUS_KEY, new ArrayList<>(relevantYngMap.values()));

        if (!onlyAnnotate) {
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

    protected VariantContext removeAnnotations(VariantContext vc) {

        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        List<String> rmAnnotationList = new ArrayList<>(Arrays.asList(GATKVCFConstants.STRAND_ODDS_RATIO_KEY,
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
    private VariantContext createVariantContextFromSampleRecord(final ExtractCohortRecord sampleRecord, HashMap<Allele, HashMap<Allele, Double>> vqsLodMap, HashMap<Allele, HashMap<Allele, String>> yngMap, boolean performGenotypeVQSLODFiltering) {
        final VariantContextBuilder builder = new VariantContextBuilder();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        final String contig = sampleRecord.getContig();
        final long startPosition = sampleRecord.getStart();
        final String sample = sampleRecord.getSampleName();

        builder.chr(contig);
        builder.start(startPosition);

        final List<Allele> alleles = new ArrayList<>();
        Allele ref = Allele.create(sampleRecord.getRefAllele(), true);
        alleles.add(ref);
        List<Allele> altAlleles = Arrays.stream(sampleRecord.getAltAllele().split(SchemaUtils.MULTIVALUE_FIELD_DELIMITER))
                .map(altAllele -> Allele.create(altAllele, false)).collect(Collectors.toList());

        // NOTE: gnarly needs this it seems?
        altAlleles.add(Allele.NON_REF_ALLELE);

        alleles.addAll(altAlleles);
        builder.alleles(alleles);

        builder.stop(startPosition + alleles.get(0).length() - 1);

        genotypeBuilder.name(sample);
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

            if (performGenotypeVQSLODFiltering) {
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
        if ( callPL != null ) {
            genotypeBuilder.PL(Arrays.stream(callPL.split(SchemaUtils.MULTIVALUE_FIELD_DELIMITER)).mapToInt(Integer::parseInt).toArray());
        }

        final String callRGQ = sampleRecord.getCallRGQ();
        if ( callRGQ != null ) {
            genotypeBuilder.attribute(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, callRGQ);
        }
        builder.genotypes(genotypeBuilder.make());

        return builder.make();
    }

    private VariantContext createRefSiteVariantContextWithGQ(final String sample, final String contig, final long start, final Allele refAllele, final Integer gq) {
        final VariantContextBuilder builder = new VariantContextBuilder("", contig, start, start, new ArrayList<>(Arrays.asList(refAllele, Allele.NON_REF_ALLELE)));
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample, new ArrayList<>(Arrays.asList(refAllele, refAllele)));

        builder.attribute(VCFConstants.END_KEY, Long.toString(start));
        if (gq != null) {
            genotypeBuilder.GQ(gq);
        }
        builder.genotypes(genotypeBuilder.make());
        return builder.make();
    }
    private VariantContext createRefSiteVariantContext(final String sample, final String contig, final long start, final Allele refAllele) {
        return createRefSiteVariantContextWithGQ(sample, contig, start, refAllele, null);
    }
}
