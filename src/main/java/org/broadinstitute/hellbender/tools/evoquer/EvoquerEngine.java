package org.broadinstitute.hellbender.tools.evoquer;

import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericRecord;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.gnarlyGenotyper.GnarlyGenotyperEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;
import java.util.stream.Collectors;

/**
 * EvoquerEngine ("EvokerEngine"):
 * Extract Variants Out of big QUERy Engine.
 *
 * Performs the work for {@link Evoquer} in a manner that is compartmentalized and reusable.
 *
 * Created by jonn on 4/17/19.
 */
class EvoquerEngine {
    private static final Logger logger = LogManager.getLogger(EvoquerEngine.class);


    public static final String SAMPLE_TABLE_NAME = "sample_list";

    public static final String POSITION_FIELD_NAME = "position";
    public static final String VALUES_ARRAY_FIELD_NAME = "values";

    public static final String SAMPLE_FIELD_NAME = "sample";
    public static final String STATE_FIELD_NAME = "state";
    public static final String REF_ALLELE_FIELD_NAME = "ref";
    public static final String ALT_ALLELE_FIELD_NAME = "alt";

    public static final String GENOTYPE_FIELD_PREFIX = "call_";
    public static final String MULTIVALUE_FIELD_DELIMITER = ",";

    public static final ImmutableSet<String> REQUIRED_FIELDS = ImmutableSet.of(
            SAMPLE_FIELD_NAME,
            STATE_FIELD_NAME,
            REF_ALLELE_FIELD_NAME,
            ALT_ALLELE_FIELD_NAME
    );

    /**
     * The conf threshold above which variants are not included in the position tables.
     * This value is used to construct the genotype information of those missing samples
     * when they are merged together into a {@link VariantContext} object
     */
    public static final int MISSING_CONF_THRESHOLD = 60;

    private final boolean precomputedResultsMode;

    private final VariantContextWriter vcfWriter;

    private final VCFHeader vcfHeader;

    private final ReferenceDataSource refSource;

    private final ReferenceConfidenceVariantContextMerger variantContextMerger;

    private final GnarlyGenotyperEngine gnarlyGenotyper;

    private final String projectID;

    /** Set of sample names seen in the variant data from BigQuery. */
    private final Set<String> sampleNames = new HashSet<>();

    private final Set<Long> seenPositions = new HashSet<>();

    /**
     * Map between contig name and the BigQuery table containing position data from that contig.
     */
    private final Map<String, String> contigToPositionTableMap;

    /**
     * Map between contig name and the BigQuery table containing variant data from that contig.
     */
    private final Map<String, String> contigToVariantTableMap;

    /**
     * Map between contig name and the BigQuery table containing the list of samples
     */
    private final Map<String, String> contigToSampleTableMap;
    
    private final int queryRecordLimit;

    private final boolean runQueryOnly;

    private final boolean disableGnarlyGenotyper;

    private final boolean runQueryInBatchMode;

    private final boolean printDebugInformation;

    private final ProgressMeter progressMeter;

    private int totalNumberOfVariants = 0;
    private int totalNumberOfSites = 0;

    EvoquerEngine( final String projectID,
                   final Map<String, Evoquer.EvoquerDataset> datasetMap,
                   final int queryRecordLimit,
                   final VariantContextWriter vcfWriter,
                   final Set<VCFHeaderLine> toolDefaultVCFHeaderLines,
                   final VariantAnnotatorEngine annotationEngine,
                   final ReferenceDataSource refSource,
                   final boolean runQueryOnly,
                   final boolean disableGnarlyGenotyper,
                   final boolean keepAllSitesInGnarlyGenotyper,
                   final boolean runQueryInBatchMode,
                   final boolean printDebugInformation,
                   final ProgressMeter progressMeter ) {

        // We were given a dataset map, so we're going to do live queries against BigQuery
        this.precomputedResultsMode = false;

        this.vcfWriter = vcfWriter;
        this.refSource = refSource;
        this.projectID = projectID;
        this.queryRecordLimit = queryRecordLimit;
        this.runQueryOnly = runQueryOnly;
        this.disableGnarlyGenotyper = disableGnarlyGenotyper;
        this.runQueryInBatchMode = runQueryInBatchMode;
        this.printDebugInformation = printDebugInformation;
        this.progressMeter = progressMeter;

        final Map<String, String> tmpContigToPositionTableMap = new HashMap<>();
        final Map<String, String> tmpContigToVariantTableMap = new HashMap<>();
        final Map<String, String> tmpContigToSampleTableMap = new HashMap<>();
        for ( final Map.Entry<String, Evoquer.EvoquerDataset> datasetEntry : datasetMap.entrySet() ) {
            tmpContigToPositionTableMap.put(datasetEntry.getKey(), datasetEntry.getValue().getDatasetName() + "." + datasetEntry.getValue().getPetTableName());
            tmpContigToVariantTableMap.put(datasetEntry.getKey(), datasetEntry.getValue().getDatasetName() + "." + datasetEntry.getValue().getVetTableName());
            tmpContigToSampleTableMap.put(datasetEntry.getKey(), datasetEntry.getValue().getDatasetName() + "." + SAMPLE_TABLE_NAME);
        }
        contigToPositionTableMap = Collections.unmodifiableMap(tmpContigToPositionTableMap);
        contigToVariantTableMap = Collections.unmodifiableMap(tmpContigToVariantTableMap);
        contigToSampleTableMap = Collections.unmodifiableMap(tmpContigToSampleTableMap);

        // Get the samples used in the dataset:
        populateSampleNames();
        this.vcfHeader = generateVcfHeader(toolDefaultVCFHeaderLines, refSource.getSequenceDictionary());

        this.variantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine, vcfHeader);
        this.gnarlyGenotyper = new GnarlyGenotyperEngine(keepAllSitesInGnarlyGenotyper, GenotypeCalculationArgumentCollection.DEFAULT_MAX_ALTERNATE_ALLELES, false, false);
    }

    EvoquerEngine( final List<String> sampleNames,
                   final VariantContextWriter vcfWriter,
                   final Set<VCFHeaderLine> toolDefaultVCFHeaderLines,
                   final VariantAnnotatorEngine annotationEngine,
                   final ReferenceDataSource refSource,
                   final boolean disableGnarlyGenotyper,
                   final boolean keepAllSitesInGnarlyGenotyper,
                   final boolean printDebugInformation,
                   final ProgressMeter progressMeter ) {

        // We weren't given a dataset map, so we're not going to do any live queries against BigQuery,
        // and will instead expect to be given URIs to precomputed results.
        this.precomputedResultsMode = true;

        this.vcfWriter = vcfWriter;
        this.refSource = refSource;
        this.disableGnarlyGenotyper = disableGnarlyGenotyper;
        this.printDebugInformation = printDebugInformation;
        this.progressMeter = progressMeter;

        this.projectID = null;
        this.contigToPositionTableMap = null;
        this.contigToVariantTableMap = null;
        this.contigToSampleTableMap = null;
        this.queryRecordLimit = 0;
        this.runQueryOnly = false;
        this.runQueryInBatchMode = false;

        this.sampleNames.addAll(sampleNames);
        this.vcfHeader = generateVcfHeader(toolDefaultVCFHeaderLines, refSource.getSequenceDictionary());
        this.variantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine, vcfHeader);
        this.gnarlyGenotyper = new GnarlyGenotyperEngine(keepAllSitesInGnarlyGenotyper, GenotypeCalculationArgumentCollection.DEFAULT_MAX_ALTERNATE_ALLELES, false, false);
    }

    /**
     * Connects to the BigQuery table for the given interval and pulls out the information on the samples that
     * contain variants.
     * @param interval {@link SimpleInterval} over which to query the BigQuery table.
     */
    void evokeInterval(final SimpleInterval interval) {
        if ( precomputedResultsMode ) {
            throw new GATKException("Cannot do live queries in precomputedResultsMode");
        }

        if ( contigToPositionTableMap.containsKey(interval.getContig()) ) {
            // Get the query string:
            final String variantQueryString = getVariantQueryString(interval);

            // Execute the query:
            final BigQueryUtils.StorageAPIAvroReader storageAPIAvroReader = BigQueryUtils.executeQueryWithStorageAPI(variantQueryString, runQueryInBatchMode);

            createVariantsFromTableResult(storageAPIAvroReader, interval.getContig());
        }
        else {
            logger.warn("Contig missing from contigPositionExpandedTableMap, ignoring interval: " + interval.toString());
        }
    }

    void evokeAvroResult(final String avroResultFile, final String contig) {
        if ( ! precomputedResultsMode ) {
            throw new GATKException("Must be in precomputed results mode to accept precomputed avro inputs");
        }

        final GATKAvroReader avroReader = new GCSAvroReader(avroResultFile);
        createVariantsFromTableResult(avroReader, contig);
    }

    /**
     * Generates a {@link VCFHeader} object based on the VariantContext objects queried from the BigQuery backend.
     * If no objects have been queried, this will return a default {@link VCFHeader}.
     * @param defaultHeaderLines The default header lines to be added to the top of the VCF header.
     * @param sequenceDictionary The SequenceDictionary of the reference on which the variants are based.
     * @return A {@link VCFHeader} object representing the header for all variants that have been queried from the BigQuery backend.
     */
    private VCFHeader generateVcfHeader(final Set<VCFHeaderLine> defaultHeaderLines,
                                       final SAMSequenceDictionary sequenceDictionary) {
        final Set<VCFHeaderLine> headerLines = new HashSet<>();

        headerLines.addAll( getEvoquerVcfHeaderLines() );
        headerLines.addAll( defaultHeaderLines );

        final VCFHeader header = new VCFHeader(headerLines, sampleNames);
        header.setSequenceDictionary(sequenceDictionary);

        return header;
    }

    VCFHeader getHeader() {
        return vcfHeader;
    }

    int getTotalNumberOfVariants() { return totalNumberOfVariants; }
    int getTotalNumberOfSites() { return totalNumberOfSites; }

    //==================================================================================================================
    // Private Instance Methods:

    private Set<VCFHeaderLine> getEvoquerVcfHeaderLines() {
        final Set<VCFHeaderLine> headerLines = new HashSet<>();

        // TODO: Get a list of all possible values here so that we can make sure they're in the VCF Header!

        // Add standard VCF fields first:
        VCFStandardHeaderLines.addStandardInfoLines( headerLines, true,
                VCFConstants.STRAND_BIAS_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.RMS_MAPPING_QUALITY_KEY,
                VCFConstants.ALLELE_COUNT_KEY,
                VCFConstants.ALLELE_FREQUENCY_KEY,
                VCFConstants.ALLELE_NUMBER_KEY,
                VCFConstants.END_KEY
        );

        VCFStandardHeaderLines.addStandardFormatLines(headerLines, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS
        );

        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
        
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_MAP_QUAL_RANK_SUM_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_QUAL_APPROX_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_READ_POS_RANK_SUM_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_SB_TABLE_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.SB_TABLE_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_VARIANT_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VARIANT_DEPTH_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_STRAND_ODDS_RATIO_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.STRAND_ODDS_RATIO_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_FISHER_STRAND_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.FISHER_STRAND_KEY));


        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.QUAL_BY_DEPTH_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EXCESS_HET_KEY));

        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.SB_TABLE_KEY));
        
        return headerLines;
    }

    private String getPositionTableForContig(final String contig ) {
        return contigToPositionTableMap.get(contig);
    }

    private String getVariantTableForContig(final String contig ) {
        return contigToVariantTableMap.get(contig);
    }

    private String getTableQualifier() {
        return projectID;
    }

    private String getFQTableName( final String tableName ) {
        return getTableQualifier() + "." + tableName;
    }

    /**
     * Get the fully-qualified table name corresponding to the table in BigQuery that contains the position
     * data specified in the given {@code interval}.
     *
     * Uses {@link #projectID} for the project of the BigQuery table.
     * Assumes the tables have dataset information in them.
     *
     * @param interval The {@link SimpleInterval} for which to get the corresponding table in BigQuery.
     * @return The name of the table corresponding to the given {@code interval}, or {@code null} if no such table exists.
     */
    private String getFQPositionTable(final SimpleInterval interval) {
        return getFQTableName(getPositionTableForContig( interval.getContig() ));
    }

    /**
     * Get the fully-qualified table name corresponding to the table in BigQuery that contains the variant
     * data specified in the given {@code interval}.
     *
     * Uses {@link #projectID} for the project of the BigQuery table.
     * Assumes the tables have dataset information in them.
     *
     * @param interval The {@link SimpleInterval} for which to get the corresponding table in BigQuery.
     * @return The name of the table corresponding to the given {@code interval}, or {@code null} if no such table exists.
     */
    private String getFQVariantTable(final SimpleInterval interval) {
        return getFQTableName(getVariantTableForContig( interval.getContig() ));
    }

    private void createVariantsFromTableResult( final GATKAvroReader avroReader, final String contig) {

        final org.apache.avro.Schema schema = avroReader.getSchema();

        final Set<String> columnNames = new HashSet<>();
        if ( schema.getField(POSITION_FIELD_NAME) == null || schema.getField(VALUES_ARRAY_FIELD_NAME) == null ) {
            throw new UserException("Records must contain position and values columns");
        }
        schema.getField(VALUES_ARRAY_FIELD_NAME).schema().getElementType().getFields().forEach(field -> columnNames.add(field.name()));
        validateSchema(columnNames);

        for ( final GenericRecord row : avroReader ) {
            if ( runQueryOnly ) {
                continue;
            }

            ++totalNumberOfSites;

            final long currentPosition = Long.parseLong(row.get(POSITION_FIELD_NAME).toString());
            final List<VariantContext> unmergedCalls = new ArrayList<>();
            final Set<String> currentPositionSamplesSeen = new HashSet<>();
            boolean currentPositionHasVariant = false;
            final Allele refAllele = Allele.create(refSource.queryAndPrefetch(contig, currentPosition, currentPosition).getBaseString(), true);

            if (!seenPositions.contains(currentPosition)) {
                seenPositions.add(currentPosition);
            } else {
                logger.warn("Duplicate variant at site " + currentPosition + " will be ignored.");
                break;
            }

            if ( printDebugInformation ) {
                logger.info(contig + ":" + currentPosition + ": found record: " + row);
            }

            final GenericData.Array<?> sampleArray = (GenericData.Array) row.get(VALUES_ARRAY_FIELD_NAME);

            for ( final Object sampleArrayEntry : sampleArray ) {
                final GenericData.Record sampleRecord = (GenericData.Record)sampleArrayEntry;
                final String sampleName = sampleRecord.get(SAMPLE_FIELD_NAME).toString();
                currentPositionSamplesSeen.add(sampleName);

                if ( printDebugInformation ) {
                    logger.info("\t" + contig + ":" + currentPosition + ": found struct for sample " + sampleName + ": " + sampleRecord);
                }

                switch (sampleRecord.get(STATE_FIELD_NAME).toString()) {
                    case "v":   // Variant
                        ++totalNumberOfVariants;
                        unmergedCalls.add(createVariantContextFromSampleRecord(sampleRecord, columnNames, contig, currentPosition, sampleName));
                        currentPositionHasVariant = true;
                        break;
                    case "0":   // Non Variant Block with GQ < 10
                        unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 0));
                        break;
                    case "1":  // Non Variant Block with 10 <=  GQ < 20
                        unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 10));
                        break;
                    case "2":  // Non Variant Block with 20 <= GQ < 30
                        unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 20));
                        break;
                    case "3":  // Non Variant Block with 30 <= GQ < 40
                        unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 30));
                        break;
                    case "4":  // Non Variant Block with 40 <= GQ < 50
                        unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 40));
                        break;
                    case "5":  // Non Variant Block with 50 <= GQ < 60
                        unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 50));
                        break;
                    case "6":  // Non Variant Block with 60 <= GQ (usually omitted from tables)
                        unmergedCalls.add(createRefSiteVariantContext(sampleName, contig, currentPosition, refAllele, 60));
                        break;
                    case "*":   // Spanning Deletion
                        unmergedCalls.add(createVariantContextForSpanningDelete(sampleName, contig, currentPosition, refAllele));
                        break;
                    case "m":   // Missing
                        // Nothing to do here -- just needed to mark the sample as seen so it doesn't get put in the high confidence ref band
                        break;
                    default:
                        throw new GATKException("Unrecognized state: " + sampleRecord.get(STATE_FIELD_NAME).toString());
                }
            }

            finalizeCurrentVariant(unmergedCalls, currentPositionSamplesSeen, currentPositionHasVariant, contig, currentPosition, refAllele);
        }
    }

    private void validateSchema(final Set<String> columnNames) {
        for ( final String requiredField : REQUIRED_FIELDS ) {
            if ( ! columnNames.contains(requiredField) ) {
                throw new UserException("Missing required column: " + requiredField +
                    ". Actual columns encountered were: " + columnNames);
            }
        }
    }

    private VariantContext createVariantContextFromSampleRecord(final GenericData.Record sampleRecord, final Set<String> columnNames, final String contig, final long startPosition, final String sample) {
        final VariantContextBuilder builder = new VariantContextBuilder();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        builder.chr(contig);
        builder.start(startPosition);

        final List<Allele> alleles = new ArrayList<>();
        alleles.add(Allele.create(sampleRecord.get(REF_ALLELE_FIELD_NAME).toString(), true));
        Arrays.stream(sampleRecord.get(ALT_ALLELE_FIELD_NAME).toString().split(MULTIVALUE_FIELD_DELIMITER))
                .forEach(altAllele -> alleles.add(Allele.create(altAllele, false)));
        builder.alleles(alleles);

        builder.stop(startPosition + alleles.get(0).length() - 1);

        genotypeBuilder.name(sample);

        for ( final String columnName : columnNames ) {
            if ( REQUIRED_FIELDS.contains(columnName) ) {
                continue;
            }

            final Object columnValue = sampleRecord.get(columnName);
            if ( columnValue == null ) {
                continue;
            }
            final String columnValueString = columnValue.toString();

            if ( columnName.startsWith(GENOTYPE_FIELD_PREFIX) ) {
                final String genotypeAttributeName = columnName.substring(GENOTYPE_FIELD_PREFIX.length());

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
                    genotypeBuilder.PL(Arrays.stream(columnValueString.split(MULTIVALUE_FIELD_DELIMITER)).mapToInt(Integer::parseInt).toArray());
                } else if ( genotypeAttributeName.equals(VCFConstants.DEPTH_KEY) ) {
                    genotypeBuilder.DP(Integer.parseInt(columnValueString));
                } else if ( genotypeAttributeName.equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS) ) {
                    genotypeBuilder.AD(Arrays.stream(columnValueString.split(MULTIVALUE_FIELD_DELIMITER)).mapToInt(Integer::parseInt).toArray());
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

    private VariantContext createVariantContextForSpanningDelete(final String sample, final String contig, final long start, final Allele refAllele) {
        final VariantContextBuilder builder = new VariantContextBuilder();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        builder.chr(contig);
        builder.start(start);

        final List<Allele> alleles = new ArrayList<>();
        alleles.add(refAllele);
        alleles.add(Allele.SPAN_DEL);
        alleles.add(Allele.NON_REF_ALLELE);
        builder.alleles(alleles);

        builder.stop(start);

        genotypeBuilder.name(sample);

        builder.attribute(VCFConstants.END_KEY, Long.toString(start));

        // is this correct? why in the build for ref is the ref added twice?
        final List<Allele> genotypeAlleles = new ArrayList<>();
        genotypeAlleles.add(refAllele);
        genotypeAlleles.add(Allele.SPAN_DEL);
        genotypeBuilder.alleles(genotypeAlleles);

        // is there a gq for a spanning deletion?
//        genotypeBuilder.GQ(gq);

        builder.genotypes(genotypeBuilder.make());

        return builder.make();

    }

    private VariantContext createRefSiteVariantContext(final String sample, final String contig, final long start, final Allele refAllele, final int gq) {
        final VariantContextBuilder builder = new VariantContextBuilder();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        builder.chr(contig);
        builder.start(start);

        final List<Allele> alleles = new ArrayList<>();
        alleles.add(refAllele);
        alleles.add(Allele.NON_REF_ALLELE);
        builder.alleles(alleles);

        builder.stop(start);

        genotypeBuilder.name(sample);

        builder.attribute(VCFConstants.END_KEY, Long.toString(start));

        final List<Allele> genotypeAlleles = new ArrayList<>();
        genotypeAlleles.add(refAllele);
        genotypeAlleles.add(refAllele);
        genotypeBuilder.alleles(genotypeAlleles);

        genotypeBuilder.GQ(gq);

        builder.genotypes(genotypeBuilder.make());
        
        return builder.make();
    }

    private void finalizeCurrentVariant(final List<VariantContext> unmergedCalls, final Set<String> currentVariantSamplesSeen, final boolean currentPositionHasVariant, final String contig, final long start, final Allele refAllele) {
        // If there were no variants at this site, we don't emit a record and there's nothing to do here
        if ( ! currentPositionHasVariant ) {
            return;
        }

        // Find missing samples and synthesize GQ 60s
        final Set<String> samplesNotEncountered = Sets.difference(sampleNames, currentVariantSamplesSeen);
        for ( final String missingSample : samplesNotEncountered ) {
            unmergedCalls.add(createRefSiteVariantContext(missingSample, contig, start, refAllele, MISSING_CONF_THRESHOLD));
        }

        // Note that we remove NON_REF in our variantContextMerger if the GnarlyGenotyper is disabled, but
        // keep NON_REF around if the GnarlyGenotyper is enabled. This is because the GnarlyGenotyper expects
        // the NON_REF allele to still be present, and will give incorrect results if it's not.
        final VariantContext mergedVC = variantContextMerger.merge(unmergedCalls, new SimpleInterval(contig, (int) start, (int) start), refAllele.getBases()[0], disableGnarlyGenotyper, false);

        final VariantContext finalizedVC = disableGnarlyGenotyper ? mergedVC : gnarlyGenotyper.finalizeGenotype(mergedVC);

        if ( finalizedVC != null ) { // GnarlyGenotyper returns null for variants it refuses to output
            vcfWriter.add(finalizedVC);
            progressMeter.update(finalizedVC);
        }
        else {
            logger.warn(String.format("GnarlyGenotyper returned null for site %s:%s", contig, start));
            progressMeter.update(mergedVC);
        }
    }

    // HACK to deal with malformed allele-specific annotations in our test dataset
    private boolean alleleSpecificAnnotationsAreMalformed( final VariantContext vc ) {
        // AS_RAW_MQ  AS_RAW_MQRankSum  AS_QUALapprox  AS_RAW_ReadPosRankSum  AS_SB_TABLE  AS_VarDP
        final int AS_RAW_MQ_values = StringUtils.countMatches(vc.getAttributeAsString("AS_RAW_MQ", ""), "|") + 1;
        final int AS_RAW_MQRankSum_values = StringUtils.countMatches(vc.getAttributeAsString("AS_RAW_MQRankSum", ""), "|") + 1;
        final int AS_QUALapprox_values = StringUtils.countMatches(vc.getAttributeAsString("AS_QUALapprox", ""), "|") + 1;
        final int AS_RAW_ReadPosRankSum_values = StringUtils.countMatches(vc.getAttributeAsString("AS_RAW_ReadPosRankSum", ""), "|") + 1;
        final int AS_SB_TABLE_values = StringUtils.countMatches(vc.getAttributeAsString("AS_SB_TABLE", ""), "|") + 1;
        final int AS_VarDP_values = StringUtils.countMatches(vc.getAttributeAsString("AS_VarDP", ""), "|") + 1;

        if ( AS_RAW_MQ_values != vc.getNAlleles() ||
             AS_RAW_MQRankSum_values != vc.getNAlleles() ||
             AS_QUALapprox_values != vc.getNAlleles() - 1 ||
             AS_RAW_ReadPosRankSum_values != vc.getNAlleles() ||
             AS_SB_TABLE_values != vc.getNAlleles() ||
             AS_VarDP_values != vc.getNAlleles() ) {

            return true;
        }

        return false;
    }

    private void populateSampleNames() {
        // TODO: For now, use the sample list from an arbitrary dataset's sample_list table.
        // TODO: Eventually, may need to crosscheck sample lists across datasets
        final String sampleTableName = contigToSampleTableMap.entrySet().iterator().next().getValue();

        // Get the query string:
        final String sampleListQueryString = getSampleListQueryString(sampleTableName);
        
        // Execute the query:
        final TableResult result = BigQueryUtils.executeQuery(sampleListQueryString);

        // Show our pretty results:
        if ( printDebugInformation ) {
            logger.info("Sample names returned:");
            final String prettyQueryResults = BigQueryUtils.getResultDataPrettyString(result);
            logger.info("\n" + prettyQueryResults);
        }

        // Add our samples to our map:
        for ( final FieldValueList row : result.iterateAll() ) {
            sampleNames.add( row.get(0).getStringValue() );
        }
    }

    private String getVariantQueryString( final SimpleInterval interval ) {
        String limitString = "";
        if ( queryRecordLimit > 0 ) {
            limitString = "LIMIT " + queryRecordLimit;
        }
        
        return String.format(
                "SELECT position, ARRAY_AGG(STRUCT( sample, state, ref, alt, AS_RAW_MQ, AS_RAW_MQRankSum, AS_QUALapprox, AS_RAW_ReadPosRankSum, AS_SB_TABLE, AS_VarDP, call_GT, call_AD, call_DP, call_GQ, call_PGT, call_PID, call_PL  )) AS values\n" +
                "FROM `%s` AS pet\n" +
                "LEFT OUTER JOIN `%s` AS vet\n" +
                "USING (position, sample)\n" +
                "WHERE (position >= %d AND position <= %d)\n" +
                "GROUP BY position\n" +
                limitString,
                getFQPositionTable(interval),
                getFQVariantTable(interval),
                interval.getStart(),
                interval.getEnd());
    }

    private String getOptimizedVariantQueryString(final SimpleInterval interval) {
        String limitString = "";
        if (queryRecordLimit > 0) {
            limitString = "LIMIT " + queryRecordLimit;
        }

        return String.format(
                "WITH new_pet AS (SELECT * FROM `%s` WHERE position in (SELECT DISTINCT position FROM `%s` WHERE position >= %d AND position <= %d AND state = 'v'))\n" +
                        "SELECT new_pet.position, ARRAY_AGG(STRUCT( new_pet.sample, state, ref, alt, AS_RAW_MQ, AS_RAW_MQRankSum, AS_QUALapprox, AS_RAW_ReadPosRankSum, AS_SB_TABLE, AS_VarDP, call_GT, call_AD, call_DP, call_GQ, call_PGT, call_PID, call_PL  )) AS values\n" +
                        "FROM new_pet\n" +
                        "LEFT OUTER JOIN `%s` AS vet\n" +
                        "USING (position, sample)\n" +
                        "GROUP BY position\n" +
                        limitString,
                getFQPositionTable(interval),
                getFQPositionTable(interval),
                interval.getStart(),
                interval.getEnd(),
                getFQVariantTable(interval));
    }

    private String getSampleListQueryString(final String sampleTableName) {
        return "SELECT sample FROM `" + getFQTableName(sampleTableName)+ "`";
    }
}
