package org.broadinstitute.hellbender.tools.variantdb.arrays;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.avro.generic.GenericRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.variantdb.arrays.BasicArrayData.ArrayGenotype;
import org.broadinstitute.hellbender.tools.variantdb.CommonCode;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.tools.variantdb.arrays.tables.ProbeInfo;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.*;
import org.broadinstitute.hellbender.utils.localsort.SortingCollection;

import java.text.DecimalFormat;
import java.util.*;
import static org.broadinstitute.hellbender.tools.variantdb.arrays.ExtractCohortBQ.*;


public class ArrayExtractCohortEngine {
    private final DecimalFormat df = new DecimalFormat();
    private final String DOT = ".";
    
    private static final Logger logger = LogManager.getLogger(ArrayExtractCohortEngine.class);

    private final VariantContextWriter vcfWriter;

    private final boolean useCompressedData;
    private final boolean printDebugInformation;
    private final int localSortMaxRecordsInRam;
    private final TableReference cohortTableRef;
    private final ReferenceDataSource refSource;

    private final ProgressMeter progressMeter;
    private final String projectID;

    /** List of sample names seen in the variant data from BigQuery. */
    private final Map<Integer, String> sampleIdMap;
    private final Set<String> sampleNames;

    private final Map<Long, ProbeInfo> probeIdMap;
    private final ReferenceConfidenceVariantContextMerger variantContextMerger;

    private int totalNumberOfVariants = 0;
    private int totalNumberOfSites = 0;

    public ArrayExtractCohortEngine(final String projectID,
                                    final VariantContextWriter vcfWriter,
                                    final VCFHeader vcfHeader,
                                    final VariantAnnotatorEngine annotationEngine,
                                    final ReferenceDataSource refSource,
                                    final Map<Integer, String> sampleIdMap,
                                    final Map<Long, ProbeInfo> probeIdMap,
                                    final String cohortTableName,
                                    final int localSortMaxRecordsInRam,
                                    final boolean useCompressedData,
                                    final boolean printDebugInformation,
                                    final ProgressMeter progressMeter) {

        this.df.setMaximumFractionDigits(3);
        this.df.setGroupingSize(0);
                                
        this.localSortMaxRecordsInRam = localSortMaxRecordsInRam;

        this.projectID = projectID;
        this.vcfWriter = vcfWriter;
        this.refSource = refSource;
        this.sampleIdMap = sampleIdMap;
        this.sampleNames = new HashSet<>(sampleIdMap.values());

        this.probeIdMap = probeIdMap;

        this.cohortTableRef = new TableReference(cohortTableName, useCompressedData? SchemaUtils.RAW_ARRAY_COHORT_FIELDS_COMPRESSED:SchemaUtils.RAW_ARRAY_COHORT_FIELDS_UNCOMPRESSED);

        this.useCompressedData = useCompressedData;
        this.printDebugInformation = printDebugInformation;
        this.progressMeter = progressMeter;

        // TODO: what is the right variant context merger for arrays?
        this.variantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine, vcfHeader);

    }

    int getTotalNumberOfVariants() { return totalNumberOfVariants; }
    int getTotalNumberOfSites() { return totalNumberOfSites; }

    public void traverse() {
        if (printDebugInformation) {
            logger.debug("using storage api with local sort");
        }
        final StorageAPIAvroReader storageAPIAvroReader = new StorageAPIAvroReader(cohortTableRef);
        createVariantsFromUngroupedTableResult(storageAPIAvroReader);
    }


    private void createVariantsFromUngroupedTableResult(final GATKAvroReader avroReader) {

        // stream out the data and sort locally
        final org.apache.avro.Schema schema = avroReader.getSchema();
        final Set<String> columnNames = new HashSet<>();
        schema.getFields().forEach(field -> columnNames.add(field.name()));

        Comparator<GenericRecord> comparator = this.useCompressedData ? COMPRESSED_PROBE_ID_COMPARATOR : UNCOMPRESSED_PROBE_ID_COMPARATOR;

        SortingCollection<GenericRecord> sortingCollection =  getAvroProbeIdSortingCollection(schema, localSortMaxRecordsInRam, comparator);
        for ( final GenericRecord queryRow : avroReader ) {
            sortingCollection.add(queryRow);
        }

        sortingCollection.printTempFileStats();

        // iterate through records and process them
        final List<GenericRecord> currentPositionRecords = new ArrayList<>(sampleIdMap.size() * 2);
        long currentProbeId = -1;

        for ( final GenericRecord sortedRow : sortingCollection ) {
            long probeId;
            if (useCompressedData) {
                final long bits = (Long) sortedRow.get(SchemaUtils.BASIC_ARRAY_DATA_FIELD_NAME);
                BasicArrayData data = new BasicArrayData(bits);
                probeId = data.probeId;
            } else {
                probeId = (Long) sortedRow.get("probe_id");
            }

            if ( probeId != currentProbeId && currentProbeId != -1 ) {
                ++totalNumberOfSites;
                processSampleRecordsForLocation(currentProbeId, currentPositionRecords, columnNames);
                currentPositionRecords.clear();
            }

            currentPositionRecords.add(sortedRow);
            currentProbeId = probeId;
        }

        if ( ! currentPositionRecords.isEmpty() ) {
            ++totalNumberOfSites;
            processSampleRecordsForLocation(currentProbeId, currentPositionRecords, columnNames);
        }
    }

    private void processSampleRecordsForLocation(final long probeId, final Iterable<GenericRecord> sampleRecordsAtPosition, final Set<String> columnNames) {
        final List<VariantContext> unmergedCalls = new ArrayList<>();
        final Set<String> currentPositionSamplesSeen = new HashSet<>();
        boolean currentPositionHasVariant = false;

        final ProbeInfo probeInfo = probeIdMap.get(probeId);
        if (probeInfo == null) {
            throw new RuntimeException("Unable to find probeInfo for " + probeId);
        }

        final String contig = probeInfo.contig;
        final long position = probeInfo.position;
        final Allele refAllele = Allele.create(refSource.queryAndPrefetch(contig, position, position).getBaseString(), true);

        int numRecordsAtPosition = 0;

        for ( final GenericRecord sampleRecord : sampleRecordsAtPosition ) {
            final long sampleId;
            if (useCompressedData) {
                final long bits = (Long) sampleRecord.get(SchemaUtils.BASIC_ARRAY_DATA_FIELD_NAME);
                BasicArrayData data = new BasicArrayData(bits);
                sampleId = data.sampleId;
            } else {
                sampleId = (Long) sampleRecord.get(SchemaUtils.SAMPLE_ID_FIELD_NAME);

                // TODO: hack to test roundtrip

            }

            // TODO: handle missing values
            String sampleName = sampleIdMap.get((int) sampleId);            
            currentPositionSamplesSeen.add(sampleName);

            ++numRecordsAtPosition;

            if ( printDebugInformation ) {
                logger.info("\t" + contig + ":" + position + ": found record for sample " + sampleName + ": " + sampleRecord);
            }

            ++totalNumberOfVariants;
            unmergedCalls.add(createVariantContextFromSampleRecord(probeInfo, sampleRecord, columnNames, contig, position, sampleName));

        }

        if ( printDebugInformation ) {
            logger.info(contig + ":" + position + ": processed " + numRecordsAtPosition + " total sample records");
        }

        finalizeCurrentVariant(unmergedCalls, currentPositionSamplesSeen, contig, position, refAllele);
    }

    private void finalizeCurrentVariant(final List<VariantContext> unmergedCalls, final Set<String> currentVariantSamplesSeen, final String contig, final long start, final Allele refAllele) {

        // TODO: this is where we infer missing data points... once we know what we want to drop
        // final Set<String> samplesNotEncountered = Sets.difference(sampleNames, currentVariantSamplesSeen);
        // for ( final String missingSample : samplesNotEncountered ) {
        //         unmergedCalls.add(createRefSiteVariantContext(missingSample, contig, start, refAllele));
        // }

        final VariantContext mergedVC = variantContextMerger.merge(
                unmergedCalls,
                new SimpleInterval(contig, (int) start, (int) start),
                refAllele.getBases()[0],
                true,
                false,
                true);


        final VariantContext finalVC = mergedVC;

        // TODO: this was commented out... probably need to re-enable
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

    private String formatFloatForVcf(final Float value) {
        if (value == null || Double.isNaN(value)) {
            return DOT;
        }
        return df.format(value);
    }

    private Float getNullableFloatFromDouble(Object d) {
        return d == null ? null : (float)  ((Double) d).doubleValue();
    }

    private VariantContext createVariantContextFromSampleRecord(final ProbeInfo probeInfo, final GenericRecord sampleRecord, final Set<String> columnNames, final String contig, final long startPosition, final String sample) {
        final VariantContextBuilder builder = new VariantContextBuilder();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        builder.chr(contig);
        builder.start(startPosition);
        builder.id(probeInfo.name);
        
        final List<Allele> alleles = new ArrayList<>();
        Allele ref = Allele.create(probeInfo.ref, true);        
        alleles.add(ref);

        Allele alleleA = Allele.create(probeInfo.alleleA, false);
        Allele alleleB = Allele.create(probeInfo.alleleB, false);

        boolean alleleAisRef = probeInfo.ref.equals(probeInfo.alleleA);
        boolean alleleBisRef = probeInfo.ref.equals(probeInfo.alleleB);

        if (alleleAisRef) {
            alleleA = ref;
        } else {
            alleles.add(alleleA);
        }

        if (alleleBisRef) {
            alleleB = ref;
        } else {
            alleles.add(alleleB);
        }

        builder.alleles(alleles);
        builder.stop(startPosition + alleles.get(0).length() - 1);

        Float normx;
        Float normy;
        Float baf;
        Float lrr;
        List<Allele> genotypeAlleles = new ArrayList<Allele>();

        if (this.useCompressedData) {
            final BasicArrayData basicData = new BasicArrayData((Long) sampleRecord.get(SchemaUtils.BASIC_ARRAY_DATA_FIELD_NAME));
            Object rd = sampleRecord.get(SchemaUtils.RAW_ARRAY_DATA_FIELD_NAME);

            final RawArrayData rawData = new RawArrayData((Long) rd);
            normx = rawData.normx;
            normy = rawData.normy;
            lrr = rawData.lrr;
            baf = rawData.baf;

            if (basicData.genotype == ArrayGenotype.AA) {
                genotypeAlleles.add(alleleA);
                genotypeAlleles.add(alleleA);
            } else if (basicData.genotype == ArrayGenotype.AB) {
                genotypeAlleles.add(alleleA);
                genotypeAlleles.add(alleleB);
            } else if (basicData.genotype == ArrayGenotype.BB) {
                genotypeAlleles.add(alleleB);
                genotypeAlleles.add(alleleB);
            } else {
                genotypeAlleles.add(Allele.NO_CALL);
                genotypeAlleles.add(Allele.NO_CALL);
            }
        } else {
            Object gt = sampleRecord.get("GT_encoded");
            ArrayGenotype agt;
            if ("AA".equals(gt.toString())) {
                genotypeAlleles.add(alleleA);
                genotypeAlleles.add(alleleA);
                agt =  ArrayGenotype.AA;
            } else if ("AB".equals(gt.toString())) {
                genotypeAlleles.add(alleleA);
                genotypeAlleles.add(alleleB);
                agt =  ArrayGenotype.AB;
            } else if ("BB".equals(gt.toString())) {
                genotypeAlleles.add(alleleB);
                genotypeAlleles.add(alleleB);
                agt =  ArrayGenotype.BB;
            } else if (".".equals(gt.toString())) {
                genotypeAlleles.add(Allele.NO_CALL);
                genotypeAlleles.add(Allele.NO_CALL);
                agt =  ArrayGenotype.NO_CALL;
            } else {
                System.out.println("Processing getnotype " + gt.toString());
                throw new RuntimeException();
            }

            // TODO: constantize
            try {
                normx = getNullableFloatFromDouble(sampleRecord.get("NORMX"));
                normy = getNullableFloatFromDouble(sampleRecord.get("NORMY"));            
                baf = getNullableFloatFromDouble(sampleRecord.get("BAF"));
                lrr = getNullableFloatFromDouble(sampleRecord.get("LRR"));

                // Hack to pack and unpack data
                BasicArrayData b = new BasicArrayData(0, (int) probeInfo.probeId, agt);
                RawArrayData d = new RawArrayData(normx, normy, lrr, baf);

                long bits = d.encode();
                RawArrayData d2 = new RawArrayData(bits);
                normx = d2.normx;
                normy = d2.normy;
                baf = d2.baf;
                lrr = d2.lrr;

            } catch (NullPointerException npe) {
                System.out.println("NPE on " + sampleRecord);
                System.out.println("NPE on BAF " + sampleRecord.get("BAF"));
                System.out.println("NPE on LRR " +sampleRecord.get("LRR"));
                throw npe;
            }


        }
        genotypeBuilder.alleles(genotypeAlleles);

        genotypeBuilder.attribute(CommonCode.NORMX, formatFloatForVcf(normx));
        genotypeBuilder.attribute(CommonCode.NORMY, formatFloatForVcf(normy));
        genotypeBuilder.attribute(CommonCode.BAF, formatFloatForVcf(baf));
        genotypeBuilder.attribute(CommonCode.LRR, formatFloatForVcf(lrr));      

        genotypeBuilder.name(sample);

        builder.genotypes(genotypeBuilder.make());

        try {
            VariantContext vc = builder.make();
            return vc;
        } catch (Exception e) {
            System.out.println("Error: "+ e.getMessage() + " processing " + sampleRecord + " and ref: " +ref + " PI: " + probeInfo.alleleA + "/" +probeInfo.alleleB + " with ga " + genotypeAlleles + " and alleles " + alleles);
            throw e;
        }
    }
}
