package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.avro.generic.GenericRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.FisherStrand;
import org.broadinstitute.hellbender.tools.walkers.annotator.QualByDepth;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandOddsRatio;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StrandBiasTest;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.bigquery.GATKAvroReader;
import org.broadinstitute.hellbender.utils.bigquery.StorageAPIAvroReader;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.localsort.AvroSortingCollection;
import org.broadinstitute.hellbender.utils.localsort.SortingCollection;

import java.util.*;

import static org.broadinstitute.hellbender.tools.variantdb.SchemaUtils.LOCATION_FIELD_NAME;

public class ExtractFeaturesEngine {

    private static final Logger logger = LogManager.getLogger(ExtractFeaturesEngine.class);

    private final VariantContextWriter vcfWriter;

    private final boolean printDebugInformation;
    private final int localSortMaxRecordsInRam;
    private final boolean trainingSitesOnly;

    private final ReferenceDataSource refSource;

    private final ReferenceConfidenceVariantContextMerger variantContextMerger;
    private final String projectID;

    private final TableReference altAlleleTable;
    private final TableReference sampleListTable;
    private final ProgressMeter progressMeter;
    private final Long minLocation;
    private final Long maxLocation;
    private final boolean useBatchQueries;

//    /** Set of sample names seen in the variant data from BigQuery. */
//    private final Set<String> sampleNames = new HashSet<>();

    public ExtractFeaturesEngine(final String projectID,
                               final VariantContextWriter vcfWriter,
                               final VCFHeader vcfHeader,
                               final VariantAnnotatorEngine annotationEngine,
                               final ReferenceDataSource refSource,
                               final boolean trainingSitesOnly,
                               final String fqAltAlleleTable,
                               final TableReference sampleListTable,
                               final Long minLocation,
                               final Long maxLocation,
                               final int localSortMaxRecordsInRam,
                               final boolean printDebugInformation,
                               final boolean useBatchQueries,
                               final ProgressMeter progressMeter) {
        this.localSortMaxRecordsInRam = localSortMaxRecordsInRam;

        this.projectID = projectID;
        this.vcfWriter = vcfWriter;
        this.refSource = refSource;
        this.trainingSitesOnly = trainingSitesOnly;
        this.altAlleleTable = new TableReference(fqAltAlleleTable, SchemaUtils.ALT_ALLELE_FIELDS);
        this.sampleListTable = sampleListTable;
        this.printDebugInformation = printDebugInformation;
        this.useBatchQueries = useBatchQueries;
        this.progressMeter = progressMeter;
        this.minLocation = minLocation;
        this.maxLocation = maxLocation;

        this.variantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine, vcfHeader);

    }
    public void traverse() {
        final String featureQueryString = ExtractFeaturesBQ.getVQSRFeatureExtractQueryString(altAlleleTable, sampleListTable, minLocation, maxLocation, trainingSitesOnly);
        logger.info(featureQueryString);
        final StorageAPIAvroReader storageAPIAvroReader = BigQueryUtils.executeQueryWithStorageAPI(featureQueryString, SchemaUtils.FEATURE_EXTRACT_FIELDS, projectID, useBatchQueries);

        createVQSRInputFromTableResult(storageAPIAvroReader);
    }

    private void createVQSRInputFromTableResult(final GATKAvroReader avroReader) {
        final org.apache.avro.Schema schema = avroReader.getSchema();
        final Set<String> columnNames = new HashSet<>();
        if ( schema.getField(LOCATION_FIELD_NAME) == null ) {
            throw new UserException("Records must contain a location column");
        }
        schema.getFields().forEach(field -> columnNames.add(field.name()));

        // TODO: this hardcodes a list of required fields... which this case doesn't require!
        // should genericize/parameterize so we can also validate the schema here
        // validateSchema(columnNames);

        SortingCollection<GenericRecord> sortingCollection = AvroSortingCollection.getAvroSortingCollection(schema, localSortMaxRecordsInRam, SchemaUtils.LOCATION_COMPARATOR);
        for ( final GenericRecord queryRow : avroReader ) {
            sortingCollection.add(queryRow);
        }
        sortingCollection.printTempFileStats();
        for ( final GenericRecord row : sortingCollection ) {
            processVQSRRecordForPosition(row);
        }
    }

    private void processVQSRRecordForPosition(GenericRecord rec) {
        final long location = Long.parseLong(rec.get(LOCATION_FIELD_NAME).toString());

        String contig = SchemaUtils.decodeContig(location);
        int position = SchemaUtils.decodePosition(location);
        // I don't understand why the other modes iterate through a list of all columns, and then
        //  switch based on the column name rather than just getting the columns desired directly (like I'm going to do here).
        // It might be because those field names are just prefixed versions of standard VCF fields?

        // TODO: de-python names (no  _)
        String ref = rec.get("ref").toString();
        String allele = rec.get("allele").toString();

        if (allele == null || allele.equals("")) {
            logger.warn("SEVERE WARNING: skipping " + contig + ":" + position + "(location="+location+") because it has a null alternate allele!");
            return;
        }

        // Numbers are returned as Long (sci notation)
        Double qual = Double.valueOf(rec.get(SchemaUtils.RAW_QUAL).toString());

        Object o_raw_ref_ad = rec.get("ref_ad");
        Double raw_ref_ad = (o_raw_ref_ad==null)?0:Double.valueOf(o_raw_ref_ad.toString());

        Object o_AS_MQRankSum = rec.get(SchemaUtils.AS_MQRankSum);
        Float AS_MQRankSum = (o_AS_MQRankSum==null)?null:Float.parseFloat(o_AS_MQRankSum.toString());

        Object o_AS_ReadPosRankSum = rec.get(SchemaUtils.AS_ReadPosRankSum);
        Float AS_ReadPosRankSum = (o_AS_ReadPosRankSum==null)?null:Float.parseFloat(o_AS_ReadPosRankSum.toString());

        Double raw_mq = Double.valueOf(rec.get(SchemaUtils.RAW_MQ).toString());
        Double raw_ad = Double.valueOf(rec.get(SchemaUtils.RAW_AD).toString());
        Double raw_ad_gt_1 = Double.valueOf(rec.get("RAW_AD_GT_1").toString());

        // TODO: KCIBUL QUESTION -- if we skip this... we won't have YNG Info @ extraction time?
        // if (raw_ad == 0) {
        //     logger.info("skipping " + contig + ":" + position + "(location="+location+") because it has no alternate reads!");
        //     return;
        // }

        int sb_ref_plus = Double.valueOf(rec.get("SB_REF_PLUS").toString()).intValue();
        int sb_ref_minus = Double.valueOf(rec.get("SB_REF_MINUS").toString()).intValue();
        int sb_alt_plus = Double.valueOf(rec.get("SB_ALT_PLUS").toString()).intValue();
        int sb_alt_minus = Double.valueOf(rec.get("SB_ALT_MINUS").toString()).intValue();


//        logger.info("processing " + contig + ":" + position);



        // NOTE: if VQSR required a merged VCF (e.g. multiple alleles on a  given row) we have to do some merging here...

        final VariantContextBuilder builder = new VariantContextBuilder();

        builder.chr(contig);
        builder.start(position);

        final List<Allele> alleles = new ArrayList<>();
        alleles.add(Allele.create(ref, true));
        alleles.add(Allele.create(allele, false));
        builder.alleles(alleles);
        builder.stop(position + alleles.get(0).length() - 1);


        double qd_depth = (raw_ref_ad + raw_ad_gt_1);
        double as_qd = QualByDepth.fixTooHighQD( (qual /  qd_depth) );

        final int[][] refAltTable = new int[][]{new int[]{sb_ref_plus, sb_ref_minus}, new int[]{ sb_alt_plus, sb_alt_minus}};

        double fs = QualityUtils.phredScaleErrorRate(Math.max(FisherStrand.pValueForContingencyTable(refAltTable), AS_StrandBiasTest.MIN_PVALUE));
        double sor = StrandOddsRatio.calculateSOR(refAltTable);

        double mq = Math.sqrt( raw_mq / raw_ad);

        builder.attribute("AS_QD", String.format("%.2f", as_qd) );
        builder.attribute("AS_FS", String.format("%.3f", fs));
        builder.attribute("AS_MQ", String.format("%.2f", mq) );
        builder.attribute("AS_MQRankSum", AS_MQRankSum==null?".":String.format("%.3f", AS_MQRankSum) );
        builder.attribute("AS_ReadPosRankSum", AS_ReadPosRankSum==null?".":String.format("%.3f", AS_ReadPosRankSum));
        builder.attribute("AS_SOR", String.format("%.3f", sor));

        // check out 478765 -- we need to "merge" different variant  contexts  at the same position.
        // I think if we just set the right "base" annotations it's possible the standard annotation processing
        // merging will take care of it?
        // Also -- is it possible to write a JS UDF like areEqual(ref1, alt2, ref2, alt2)?
        VariantContext vc = builder.make();
        vcfWriter.add(vc);
        progressMeter.update(vc);


    }
}
