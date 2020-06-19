package org.broadinstitute.hellbender.tools.variantdb;

import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.avro.generic.GenericRecord;
import org.apache.commons.lang.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.*;
import org.broadinstitute.hellbender.utils.localsort.AvroSortingCollectionCodec;
import org.broadinstitute.hellbender.utils.localsort.SortingCollection;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;
import static org.broadinstitute.hellbender.tools.variantdb.ExtractCohortBQ.*;

public class ExtractCohortEngine {
    private static final Logger logger = LogManager.getLogger(ExtractCohortEngine.class);

    private final VariantContextWriter vcfWriter;

    private final boolean printDebugInformation;
    private final int localSortMaxRecordsInRam;
    private final TableReference cohortTableRef;
    private final TableReference filteringTableRef;
    private final ReferenceDataSource refSource;
    private double vqsLodSNPThreshold = 0;
    private double vqsLodINDELThreshold = 0;

    private final ProgressMeter progressMeter;
    private final String projectID;
    private final CommonCode.ModeEnum mode;

    /** List of sample names seen in the variant data from BigQuery. */
    private Set<String> sampleNames;
    private final ReferenceConfidenceVariantContextMerger variantContextMerger;
    private final ExtractCohort.QueryMode queryMode;

    private int totalNumberOfVariants = 0;
    private int totalNumberOfSites = 0;

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
                               final String filteringTableName,
                               final int localSortMaxRecordsInRam,
                               final boolean printDebugInformation,
                               final double vqsLodSNPThreshold,
                               final double vqsLodINDELThreshold,
                               final ProgressMeter progressMeter,
                               final ExtractCohort.QueryMode queryMode) {
        this.localSortMaxRecordsInRam = localSortMaxRecordsInRam;

        this.projectID = projectID;
        this.vcfWriter = vcfWriter;
        this.refSource = refSource;
        this.sampleNames = sampleNames;
        this.mode = mode;
        this.cohortTableRef = mode.equals(CommonCode.ModeEnum.ARRAYS) ? new TableReference(cohortTableName, SchemaUtils.ARRAY_COHORT_FIELDS) : new TableReference(cohortTableName, SchemaUtils.COHORT_FIELDS);
        this.filteringTableRef = mode.equals(CommonCode.ModeEnum.ARRAYS) ? null : new TableReference(filteringTableName, SchemaUtils.YNG_FIELDS);
        this.printDebugInformation = printDebugInformation;
        this.vqsLodSNPThreshold = vqsLodSNPThreshold;
        this.vqsLodINDELThreshold = vqsLodINDELThreshold;
        this.progressMeter = progressMeter;
        this.queryMode = queryMode;

        this.variantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine, vcfHeader);

    }

    int getTotalNumberOfVariants() { return totalNumberOfVariants; }
    int getTotalNumberOfSites() { return totalNumberOfSites; }

    public void traverse() {
        switch (queryMode) {
            case LOCAL_SORT:
                if (printDebugInformation) {
                    logger.debug("using storage api with local sort");
                }
                final StorageAPIAvroReader storageAPIAvroReader = new StorageAPIAvroReader(cohortTableRef);
                createVariantsFromUngroupedTableResult(storageAPIAvroReader);
                break;
            case QUERY:
                if (printDebugInformation) {
                    logger.debug("using query api with order by");
                }
                // create the query string
                String q = "SELECT " + StringUtils.join(SchemaUtils.ARRAY_COHORT_FIELDS,",") + " FROM " + cohortTableRef.getFQTableName() + " ORDER BY " + SchemaUtils.LOCATION_FIELD_NAME;
                TableResult tr = BigQueryUtils.executeQuery(BigQueryUtils.getBigQueryEndPoint(), cohortTableRef.tableProject, cohortTableRef.tableDataset, q);
                createVariantsFromSortedTableResults(new QueryAPIRowReader(tr));
                break;
        }
    }

    private void createVariantsFromSortedTableResults(final QueryAPIRowReader reader) {

//        final Set<String> columnNames = new HashSet<>();
//        if ( schema.getField(POSITION_FIELD_NAME) == null || schema.getField(VALUES_ARRAY_FIELD_NAME) == null ) {
//            throw new UserException("Records must contain position and values columns");
//        }
//        schema.getField(VALUES_ARRAY_FIELD_NAME).schema().getElementType().getFields().forEach(field -> columnNames.add(field.name()));
//        validateSchema(columnNames);

        List<GenericRecord> sampleRecords = new ArrayList<>();
        long currentLocation = -1;

        while (reader.hasNext()) {
            ++totalNumberOfSites;
            FieldValueList row = reader.next();

            final long location = row.get(SchemaUtils.LOCATION_FIELD_NAME).getLongValue();
            if (currentLocation < 0) {
                currentLocation = location;
            }

            if (currentLocation != location) {
                if ( printDebugInformation ) {
                    logger.info(currentLocation + ": processing records");
                }
                // TODO this should start a thread or something - i.e. scatter
                processSampleRecordsForLocation(currentLocation, sampleRecords, new HashSet<>(SchemaUtils.ARRAY_COHORT_FIELDS));
                currentLocation = location;
                sampleRecords = new ArrayList<>();
            }

            sampleRecords.add(new QueryRecord(row));

        }
        processSampleRecordsForLocation(currentLocation, sampleRecords, new HashSet<>(SchemaUtils.ARRAY_COHORT_FIELDS));

    }

   



    private void createVariantsFromUngroupedTableResult(final GATKAvroReader avroReader) {

        final org.apache.avro.Schema schema = avroReader.getSchema();

        final Set<String> columnNames = new HashSet<>();
//        if ( schema.getField(SchemaUtils.POSITION_FIELD_NAME) == null ) {
//            throw new UserException("Records must contain a position column");
//        }
        schema.getFields().forEach(field -> columnNames.add(field.name()));
//        validateSchema(columnNames);

        SortingCollection<GenericRecord> sortingCollection =  getAvroSortingCollection(schema, localSortMaxRecordsInRam);

        for ( final GenericRecord queryRow : avroReader ) {
            sortingCollection.add(queryRow);
        }

        sortingCollection.printTempFileStats();

        final List<GenericRecord> currentPositionRecords = new ArrayList<>(sampleNames.size() * 2);
        long currentLocation = -1;

        for ( final GenericRecord sortedRow : sortingCollection ) {
            final long location = Long.parseLong(sortedRow.get(SchemaUtils.LOCATION_FIELD_NAME).toString());

            if ( location != currentLocation && currentLocation != -1 ) {
                ++totalNumberOfSites;
                processSampleRecordsForLocation(currentLocation, currentPositionRecords, columnNames);

                currentPositionRecords.clear();
            }

            currentPositionRecords.add(sortedRow);
            currentLocation = location;
        }

        if ( ! currentPositionRecords.isEmpty() ) {
            ++totalNumberOfSites;
            processSampleRecordsForLocation(currentLocation, currentPositionRecords, columnNames);
        }
    }

    private void processSampleRecordsForLocation(final long location, final Iterable<GenericRecord> sampleRecordsAtPosition, final Set<String> columnNames) {
        final List<VariantContext> unmergedCalls = new ArrayList<>();
        final Set<String> currentPositionSamplesSeen = new HashSet<>();
        boolean currentPositionHasVariant = false;
        final long currentPosition = SchemaUtils.decodePosition(location);
        final String contig = SchemaUtils.decodeContig(location);
        final Allele refAllele = Allele.create(refSource.queryAndPrefetch(contig, currentPosition, currentPosition).getBaseString(), true);
        int numRecordsAtPosition = 0;

        final HashMap<Allele, HashMap<Allele, Double>> vqsLodMap = new HashMap<>();
        final HashMap<Allele, HashMap<Allele, String>> yngMap = new HashMap<>();

        for ( final GenericRecord sampleRecord : sampleRecordsAtPosition ) {
            final String sampleName = sampleRecord.get(SchemaUtils.SAMPLE_NAME_FIELD_NAME).toString();
            currentPositionSamplesSeen.add(sampleName);
            ++numRecordsAtPosition;

            if ( printDebugInformation ) {
                logger.info("\t" + contig + ":" + currentPosition + ": found record for sample " + sampleName + ": " + sampleRecord);
            }

            switch (sampleRecord.get(SchemaUtils.STATE_FIELD_NAME).toString()) {
                case "v":   // Variant
                    ++totalNumberOfVariants;
                    unmergedCalls.add(createVariantContextFromSampleRecord(sampleRecord, columnNames, contig, currentPosition, sampleName, vqsLodMap, yngMap));
                    currentPositionHasVariant = true;
                    break;
                case "0":   // Non Variant Block with GQ < 10
                    unmergedCalls.add(createRefSiteVariantContextWithGQ(sampleName, contig, currentPosition, refAllele, 0));
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
                    throw new GATKException("Unrecognized state: " + sampleRecord.get(SchemaUtils.STATE_FIELD_NAME).toString());
            }

        }

        if ( printDebugInformation ) {
            logger.info(contig + ":" + currentPosition + ": processed " + numRecordsAtPosition + " total sample records");
        }

        finalizeCurrentVariant(unmergedCalls, currentPositionSamplesSeen, currentPositionHasVariant, contig, currentPosition, refAllele, vqsLodMap, yngMap);
    }

    private void finalizeCurrentVariant(final List<VariantContext> unmergedCalls, final Set<String> currentVariantSamplesSeen, final boolean currentPositionHasVariant, final String contig, final long start, final Allele refAllele, HashMap<Allele, HashMap<Allele, Double>> vqsLodMap, HashMap<Allele, HashMap<Allele, String>> yngMap) {
        // If there were no variants at this site, we don't emit a record and there's nothing to do here
        if ( ! currentPositionHasVariant ) {
            return;
        }

        // Find samples for dropped state and synthesize. If not arrays use GQ 60
        final Set<String> samplesNotEncountered = Sets.difference(sampleNames, currentVariantSamplesSeen);
        for ( final String missingSample : samplesNotEncountered ) {
            if (mode.equals(CommonCode.ModeEnum.ARRAYS)) {
                unmergedCalls.add(createRefSiteVariantContext(missingSample, contig, start, refAllele));

            } else {
                unmergedCalls.add(createRefSiteVariantContextWithGQ(missingSample, contig, start, refAllele, MISSING_CONF_THRESHOLD));
            }
        }

        final VariantContext mergedVC = variantContextMerger.merge(
                unmergedCalls,
                new SimpleInterval(contig, (int) start, (int) start),
                refAllele.getBases()[0],
                true,
                false,
                true);


        final VariantContext finalVC = mode.equals(CommonCode.ModeEnum.ARRAYS) ? mergedVC : filterVariants(mergedVC);
//        final VariantContext annotatedVC = enableVariantAnnotator ?
//                variantAnnotator.annotateContext(finalizedVC, new FeatureContext(), null, null, a -> true): finalVC;

//        if ( annotatedVC != null ) {
//            vcfWriter.add(annotatedVC);
//            progressMeter.update(annotatedVC);
//        }


        if ( finalVC != null ) {
            vcfWriter.add(finalVC);
            progressMeter.update(finalVC);
        } else {
            // TODO should i print a warning here?
            vcfWriter.add(mergedVC);
            progressMeter.update(mergedVC);
        }
    }

    private VariantContext filterVariants(VariantContext mergedVC) {
        return mergedVC;
        //TODO not needed for arrays
//            LinkedHashMap<Allele, Double> remappedVqsLodMap = remapAllelesInMap(mergedVC, vqsLodMap, Double.NaN);
//            LinkedHashMap<Allele, String> remappedYngMap = remapAllelesInMap(mergedVC, yngMap, VCFConstants.EMPTY_INFO_FIELD);
//
//            final VariantContextBuilder builder = new VariantContextBuilder(mergedVC);
//
//            builder.attribute(GATKVCFConstants.AS_VQS_LOD_KEY, remappedVqsLodMap.values().stream().map(val -> val.equals(Double.NaN) ? VCFConstants.EMPTY_INFO_FIELD : val.toString()).collect(Collectors.toList()));
//            builder.attribute(GATKVCFConstants.AS_YNG_STATUS_KEY, remappedYngMap.values());
//
//            int refLength = mergedVC.getReference().length();
//
//            // if there are any Yays, the site is PASS
//            if (remappedYngMap.values().contains("Y")) {
//                builder.filter("PASS");
//            } else if (remappedYngMap.values().contains("N")) {
//                // TODO: do we want to remove this variant?
//                builder.filter("NAY");
//            } else {
//                if (remappedYngMap.values().contains("G")) {
//                    // TODO change the initial query to include the filtername from the tranches tables
//                    Optional<Double> snpMax = remappedVqsLodMap.entrySet().stream().filter(entry -> entry.getKey().length() == refLength).map(entry -> entry.getValue().equals(Double.NaN) ? 0.0 : entry.getValue()).max(Double::compareTo);
//                    if (snpMax.isPresent() && snpMax.get() < vqsLodSNPThreshold) {
//                        // TODO: add in sensitivities
//                        builder.filter("VQSRTrancheSNP");
//                    }
//                    Optional<Double> indelMax = remappedVqsLodMap.entrySet().stream().filter(entry -> entry.getKey().length() != refLength).map(entry -> entry.getValue().equals(Double.NaN) ? 0.0 : entry.getValue()).max(Double::compareTo);
//                    if (indelMax.isPresent() && indelMax.get() < vqsLodINDELThreshold) {
//                        // TODO: add in sensitivities
//                        builder.filter("VQSRTrancheINDEL");
//                    }
//                }
//                // TODO: what if there is nothing in the YNG table?
//                // this shouldn't happen
//            }
//
//            final VariantContext filteredVC = builder.make();
        }

    // vqsLogMap and yngMap are in/out parameters for this method. i.e. they are modified by this method
    private VariantContext createVariantContextFromSampleRecord(final GenericRecord sampleRecord, final Set<String> columnNames, final String contig, final long startPosition, final String sample, HashMap<Allele, HashMap<Allele, Double>> vqsLodMap, HashMap<Allele, HashMap<Allele, String>> yngMap) {
        final VariantContextBuilder builder = new VariantContextBuilder();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        builder.chr(contig);
        builder.start(startPosition);

        final List<Allele> alleles = new ArrayList<>();
        Allele ref = Allele.create(sampleRecord.get(SchemaUtils.REF_ALLELE_FIELD_NAME).toString(), true);
        alleles.add(ref);
        List<Allele> altAlleles = Arrays.stream(sampleRecord.get(SchemaUtils.ALT_ALLELE_FIELD_NAME).toString().split(SchemaUtils.MULTIVALUE_FIELD_DELIMITER))
                .map(altAllele -> Allele.create(altAllele, false)).collect(Collectors.toList());
        alleles.addAll(altAlleles);
        builder.alleles(alleles);

        builder.stop(startPosition + alleles.get(0).length() - 1);

        genotypeBuilder.name(sample);

        vqsLodMap.putIfAbsent(ref, new HashMap<>());
        yngMap.putIfAbsent(ref, new HashMap<>());


        for ( final String columnName : columnNames ) {
            if ( SchemaUtils.REQUIRED_FIELDS.contains(columnName) ||
                    columnName.equals(SchemaUtils.LOCATION_FIELD_NAME) ) {
                continue;
            }

            final Object columnValue = sampleRecord.get(columnName);
            if ( columnValue == null ) {
                continue;
            }
            final String columnValueString = columnValue.toString();

            if ( columnName.startsWith(SchemaUtils.GENOTYPE_FIELD_PREFIX) ) {
                final String genotypeAttributeName = columnName.substring(SchemaUtils.GENOTYPE_FIELD_PREFIX.length());

                if ( genotypeAttributeName.equals(VCFConstants.GENOTYPE_KEY) ) {
                    final List<Allele> genotypeAlleles =
                            Arrays.stream(columnValueString.split("[/|]"))
                                    .map(Integer::parseInt)
                                    .map(alleleIndex -> alleles.get(alleleIndex))
                                    .collect(Collectors.toList());
                    genotypeBuilder.alleles(genotypeAlleles);
                } else if ( genotypeAttributeName.equals(VCFConstants.GENOTYPE_QUALITY_KEY) ) {
                    genotypeBuilder.GQ(Integer.parseInt(columnValueString));
                } else if ( genotypeAttributeName.equals(VCFConstants.GENOTYPE_PL_KEY) ) {
                    genotypeBuilder.PL(Arrays.stream(columnValueString.split(SchemaUtils.MULTIVALUE_FIELD_DELIMITER)).mapToInt(Integer::parseInt).toArray());
                } else if ( genotypeAttributeName.equals(VCFConstants.DEPTH_KEY) ) {
                    genotypeBuilder.DP(Integer.parseInt(columnValueString));
                } else if ( genotypeAttributeName.equals(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY) ) {
                    genotypeBuilder.attribute(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, Integer.parseInt(columnValueString));
                } else if ( genotypeAttributeName.equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS) ) {
                    genotypeBuilder.AD(Arrays.stream(columnValueString.split(SchemaUtils.MULTIVALUE_FIELD_DELIMITER)).mapToInt(Integer::parseInt).toArray());
                } else if ( genotypeAttributeName.equals(GATKVCFConstants.AS_VQS_LOD_KEY) ) {
                    HashMap<Allele, Double> innerVqslodMap = vqsLodMap.get(ref);
                    double[] vqslodValues = Arrays.stream(columnValueString.split(SchemaUtils.MULTIVALUE_FIELD_DELIMITER)).mapToDouble(Double::parseDouble).toArray();
                    // should we use put or putIfAbsent - this may insert the same value multiple times. same with YNG below
                    new IndexRange(0, vqslodValues.length).forEach(i -> innerVqslodMap.put(altAlleles.get(i), vqslodValues[i]));
                } else if ( genotypeAttributeName.equals("YNG_STATUS") ) {
                    HashMap<Allele, String> innerYNGMap = yngMap.get(ref);
                    String[] yngValues = columnValueString.split(SchemaUtils.MULTIVALUE_FIELD_DELIMITER);
                    new IndexRange(0, yngValues.length).forEach(i -> innerYNGMap.put(altAlleles.get(i), yngValues[i]));
                } else {
                    genotypeBuilder.attribute(genotypeAttributeName, columnValueString);
                }
            } else {
                builder.attribute(columnName, columnValueString);
            }
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
