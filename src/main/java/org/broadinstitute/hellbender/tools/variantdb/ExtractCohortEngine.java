package org.broadinstitute.hellbender.tools.variantdb;

import com.google.common.collect.Sets;
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
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.GATKAvroReader;
import org.broadinstitute.hellbender.utils.bigquery.StorageAPIAvroReader;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.localsort.AvroSortingCollectionCodec;
import org.broadinstitute.hellbender.utils.localsort.EvoquerSortingCollection;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;

public class ExtractCohortEngine {
    private static final Logger logger = LogManager.getLogger(ExtractCohortEngine.class);

    private final VariantContextWriter vcfWriter;

    private final boolean printDebugInformation;
    private final int localSortMaxRecordsInRam;
//    private final TableReference sampleTableRef;
    private final TableReference cohortTableRef;
    private final TableReference filteringTableRef;
    private final ReferenceDataSource refSource;
    private double vqsLodSNPThreshold = 0;
    private double vqsLodINDELThreshold = 0;

    private final ProgressMeter progressMeter;
    private final String projectID;

    /** List of sample names seen in the variant data from BigQuery. */
    private Set<String> sampleNames;
    private final ReferenceConfidenceVariantContextMerger variantContextMerger;


    private int totalNumberOfVariants = 0;
    private int totalNumberOfSites = 0;

    /**
     * The conf threshold above which variants are not included in the position tables.
     * This value is used to construct the genotype information of those missing samples
     * when they are merged together into a {@link VariantContext} object
     */
    public static final int MISSING_CONF_THRESHOLD = 60;


    public ExtractCohortEngine(final String projectID,
                               final VariantContextWriter vcfWriter,
                               final VCFHeader vcfHeader,
                               final VariantAnnotatorEngine annotationEngine,
                               final ReferenceDataSource refSource,
                               final Set<String> sampleNames,
                               final String cohortTableName,
                               final String filteringTableName,
                               final int localSortMaxRecordsInRam,
                               final boolean printDebugInformation,
                               final double vqsLodSNPThreshold,
                               final double vqsLodINDELThreshold,
                               final ProgressMeter progressMeter) {
        this.localSortMaxRecordsInRam = localSortMaxRecordsInRam;

        this.projectID = projectID;
        this.vcfWriter = vcfWriter;
        this.refSource = refSource;
        this.sampleNames = sampleNames;
        this.cohortTableRef = new TableReference(cohortTableName, SchemaConstants.COHORT_FIELDS);
        this.filteringTableRef = filteringTableName == null ? null : new TableReference(filteringTableName, SchemaConstants.YNG_FIELDS);
        this.printDebugInformation = printDebugInformation;
        this.vqsLodSNPThreshold = vqsLodSNPThreshold;
        this.vqsLodINDELThreshold = vqsLodINDELThreshold;
        this.progressMeter = progressMeter;

        this.variantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine, vcfHeader);

    }

    public void traverse() {
        final StorageAPIAvroReader storageAPIAvroReader = new StorageAPIAvroReader(cohortTableRef);
        createVariantsFromUngroupedTableResult(storageAPIAvroReader);
    }

    private EvoquerSortingCollection<GenericRecord> getEvoquerSortingCollection(org.apache.avro.Schema schema) {
        final EvoquerSortingCollection.Codec<GenericRecord> sortingCollectionCodec = new AvroSortingCollectionCodec(schema);
        final Comparator<GenericRecord> sortingCollectionComparator = new Comparator<GenericRecord>() {
            @Override
            public int compare( GenericRecord o1, GenericRecord o2 ) {
                final long firstPosition = Long.parseLong(o1.get(SchemaConstants.POSITION_FIELD_NAME).toString());
                final long secondPosition = Long.parseLong(o2.get(SchemaConstants.POSITION_FIELD_NAME).toString());

                return Long.compare(firstPosition, secondPosition);
            }
        };
        return EvoquerSortingCollection.newInstance(GenericRecord.class, sortingCollectionCodec, sortingCollectionComparator, localSortMaxRecordsInRam);
    }



    private void createVariantsFromUngroupedTableResult(final GATKAvroReader avroReader) {

        final org.apache.avro.Schema schema = avroReader.getSchema();

        final Set<String> columnNames = new HashSet<>();
        if ( schema.getField(SchemaConstants.POSITION_FIELD_NAME) == null ) {
            throw new UserException("Records must contain a position column");
        }
        schema.getFields().forEach(field -> columnNames.add(field.name()));
//        validateSchema(columnNames);

        EvoquerSortingCollection<GenericRecord> sortingCollection =  getEvoquerSortingCollection(schema);

        for ( final GenericRecord queryRow : avroReader ) {
            sortingCollection.add(queryRow);
        }

        sortingCollection.printTempFileStats();

        final List<GenericRecord> currentPositionRecords = new ArrayList<>(sampleNames.size() * 2);
        long currentPosition = -1;
        String currentContig = "";

        for ( final GenericRecord sortedRow : sortingCollection ) {
            final long rowPosition = Long.parseLong(sortedRow.get(SchemaConstants.POSITION_FIELD_NAME).toString());
            currentContig = sortedRow.get(SchemaConstants.CHROM_FIELD_NAME).toString();

            if ( rowPosition != currentPosition && currentPosition != -1 ) {
                ++totalNumberOfSites;
                processSampleRecordsForPosition(currentPosition, currentContig, currentPositionRecords, columnNames);

                currentPositionRecords.clear();
            }

            currentPositionRecords.add(sortedRow);
            currentPosition = rowPosition;
        }

        if ( ! currentPositionRecords.isEmpty() ) {
            ++totalNumberOfSites;
            processSampleRecordsForPosition(currentPosition, currentContig, currentPositionRecords, columnNames);
        }
    }

    private void processSampleRecordsForPosition(final long currentPosition, final String contig, final Iterable<GenericRecord> sampleRecordsAtPosition, final Set<String> columnNames) {
        final List<VariantContext> unmergedCalls = new ArrayList<>();
        final Set<String> currentPositionSamplesSeen = new HashSet<>();
        boolean currentPositionHasVariant = false;
        final Allele refAllele = Allele.create(refSource.queryAndPrefetch(contig, currentPosition, currentPosition).getBaseString(), true);
        int numRecordsAtPosition = 0;

        final HashMap<Allele, HashMap<Allele, Double>> vqsLodMap = new HashMap<>();
        final HashMap<Allele, HashMap<Allele, String>> yngMap = new HashMap<>();

        for ( final GenericRecord sampleRecord : sampleRecordsAtPosition ) {
            final String sampleName = sampleRecord.get(SchemaConstants.SAMPLE_FIELD_NAME).toString();
            currentPositionSamplesSeen.add(sampleName);
            ++numRecordsAtPosition;

            if ( printDebugInformation ) {
                logger.info("\t" + contig + ":" + currentPosition + ": found record for sample " + sampleName + ": " + sampleRecord);
            }

            switch (sampleRecord.get(SchemaConstants.STATE_FIELD_NAME).toString()) {
                case "v":   // Variant
                    ++totalNumberOfVariants;
                    unmergedCalls.add(createVariantContextFromSampleRecord(sampleRecord, columnNames, contig, currentPosition, sampleName, vqsLodMap, yngMap));
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
                case "*":   // Spanning Deletion: TODO handle correctly
                    break;
                case "m":   // Missing
                    // Nothing to do here -- just needed to mark the sample as seen so it doesn't get put in the high confidence ref band
                    break;
                default:
                    throw new GATKException("Unrecognized state: " + sampleRecord.get(SchemaConstants.STATE_FIELD_NAME).toString());
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

        // Find missing samples and synthesize GQ 60s
        final Set<String> samplesNotEncountered = Sets.difference(sampleNames, currentVariantSamplesSeen);
        for ( final String missingSample : samplesNotEncountered ) {
            unmergedCalls.add(createRefSiteVariantContext(missingSample, contig, start, refAllele, MISSING_CONF_THRESHOLD));
        }

        // Note that we remove NON_REF in our variantContextMerger if the GnarlyGenotyper is disabled, but
        // keep NON_REF around if the GnarlyGenotyper is enabled. This is because the GnarlyGenotyper expects
        // the NON_REF allele to still be present, and will give incorrect results if it's not.
        final VariantContext mergedVC = variantContextMerger.merge(unmergedCalls, new SimpleInterval(contig, (int) start, (int) start), refAllele.getBases()[0], true, false, true);


        final VariantContext finalVC = filteringTableRef != null ? mergedVC : filterVariants(mergedVC);
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
        Allele ref = Allele.create(sampleRecord.get(SchemaConstants.REF_ALLELE_FIELD_NAME).toString(), true);
        alleles.add(ref);
        List<Allele> altAlleles = Arrays.stream(sampleRecord.get(SchemaConstants.ALT_ALLELE_FIELD_NAME).toString().split(SchemaConstants.MULTIVALUE_FIELD_DELIMITER))
                .map(altAllele -> Allele.create(altAllele, false)).collect(Collectors.toList());
        alleles.addAll(altAlleles);
        builder.alleles(alleles);

        builder.stop(startPosition + alleles.get(0).length() - 1);

        genotypeBuilder.name(sample);

        vqsLodMap.putIfAbsent(ref, new HashMap<>());
        yngMap.putIfAbsent(ref, new HashMap<>());


        for ( final String columnName : columnNames ) {
            if ( SchemaConstants.REQUIRED_FIELDS.contains(columnName) ||
                    columnName.equals(SchemaConstants.POSITION_FIELD_NAME) ||
                    columnName.equals(SchemaConstants.CHROM_FIELD_NAME) ) {
                continue;
            }

            final Object columnValue = sampleRecord.get(columnName);
            if ( columnValue == null ) {
                continue;
            }
            final String columnValueString = columnValue.toString();

            if ( columnName.startsWith(SchemaConstants.GENOTYPE_FIELD_PREFIX) ) {
                final String genotypeAttributeName = columnName.substring(SchemaConstants.GENOTYPE_FIELD_PREFIX.length());

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
                    genotypeBuilder.PL(Arrays.stream(columnValueString.split(SchemaConstants.MULTIVALUE_FIELD_DELIMITER)).mapToInt(Integer::parseInt).toArray());
                } else if ( genotypeAttributeName.equals(VCFConstants.DEPTH_KEY) ) {
                    genotypeBuilder.DP(Integer.parseInt(columnValueString));
                } else if ( genotypeAttributeName.equals(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY) ) {
                    genotypeBuilder.attribute(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, Integer.parseInt(columnValueString));
                } else if ( genotypeAttributeName.equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS) ) {
                    genotypeBuilder.AD(Arrays.stream(columnValueString.split(SchemaConstants.MULTIVALUE_FIELD_DELIMITER)).mapToInt(Integer::parseInt).toArray());
                } else if ( genotypeAttributeName.equals(GATKVCFConstants.AS_VQS_LOD_KEY) ) {
                    HashMap<Allele, Double> innerVqslodMap = vqsLodMap.get(ref);
                    double[] vqslodValues = Arrays.stream(columnValueString.split(SchemaConstants.MULTIVALUE_FIELD_DELIMITER)).mapToDouble(Double::parseDouble).toArray();
                    // should we use put or putIfAbsent - this may insert the same value multiple times. same with YNG below
                    new IndexRange(0, vqslodValues.length).forEach(i -> innerVqslodMap.put(altAlleles.get(i), vqslodValues[i]));
                } else if ( genotypeAttributeName.equals("YNG_STATUS") ) {
                    HashMap<Allele, String> innerYNGMap = yngMap.get(ref);
                    String[] yngValues = columnValueString.split(SchemaConstants.MULTIVALUE_FIELD_DELIMITER);
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
}
