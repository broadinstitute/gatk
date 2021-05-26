package org.broadinstitute.hellbender.tools.gvs.filtering;

import htsjdk.samtools.util.OverlapDetector;
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
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.ExcessHet;
import org.broadinstitute.hellbender.tools.walkers.annotator.FisherStrand;
import org.broadinstitute.hellbender.tools.walkers.annotator.QualByDepth;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandOddsRatio;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StrandBiasTest;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.utils.GenotypeCounts;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.bigquery.GATKAvroReader;
import org.broadinstitute.hellbender.utils.bigquery.StorageAPIAvroReader;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.localsort.AvroSortingCollection;
import org.broadinstitute.hellbender.utils.localsort.SortingCollection;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.util.*;

import static org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils.LOCATION_FIELD_NAME;

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
    private final List<SimpleInterval> traversalIntervals;
    private final Long minLocation;
    private final Long maxLocation;
    private final boolean useBatchQueries;
    private final int numSamples;
    private final int hqGenotypeGQThreshold;
    private final int hqGenotypeDepthThreshold;
    private final double hqGenotypeABThreshold;
    private final int excessAllelesThreshold;
    private final List<String> queryLabels;

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
                               final List<SimpleInterval> traversalIntervals,
                               final Long minLocation,
                               final Long maxLocation,
                               final int localSortMaxRecordsInRam,
                               final boolean printDebugInformation,
                               final boolean useBatchQueries,
                               final ProgressMeter progressMeter,
                               final int numSamples,
                               final int hqGenotypeGQThreshold,
                               final int hqGenotypeDepthThreshold,
                               final double hqGenotypeABThreshold,
                               final int excessAllelesThreshold,
                               final List<String> queryLabels
    ) {

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
        this.traversalIntervals = traversalIntervals;
        this.minLocation = minLocation;
        this.maxLocation = maxLocation;
        this.numSamples = numSamples;
        this.hqGenotypeGQThreshold = hqGenotypeGQThreshold;
        this.hqGenotypeDepthThreshold = hqGenotypeDepthThreshold;
        this.hqGenotypeABThreshold = hqGenotypeABThreshold;
        this.excessAllelesThreshold = excessAllelesThreshold;
        this.queryLabels = queryLabels;

        this.variantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine, vcfHeader);
    }

    // taken from GnarlyGenotypingEngine
    private final static double INDEL_QUAL_THRESHOLD = GenotypeCalculationArgumentCollection.DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING - 10 * Math.log10(HomoSapiensConstants.INDEL_HETEROZYGOSITY);
    private final static double SNP_QUAL_THRESHOLD = GenotypeCalculationArgumentCollection.DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING - 10 * Math.log10(HomoSapiensConstants.SNP_HETEROZYGOSITY);

    public void traverse() {

        final String featureQueryString = ExtractFeaturesBQ.getVQSRFeatureExtractQueryString(altAlleleTable,
                                                                                             sampleListTable,
                                                                                             minLocation,
                                                                                             maxLocation,
                                                                                             trainingSitesOnly,
                                                                                             hqGenotypeGQThreshold,
                                                                                             hqGenotypeDepthThreshold,
                                                                                             hqGenotypeABThreshold);

        final String userDefinedFunctions = ExtractFeaturesBQ.getVQSRFeatureExtractUserDefinedFunctionsString();
        Map<String, String> cleanQueryLabels = createQueryLabels(queryLabels);

        final StorageAPIAvroReader storageAPIAvroReader = BigQueryUtils.executeQueryWithStorageAPI(
                featureQueryString,
                SchemaUtils.FEATURE_EXTRACT_FIELDS,
                projectID,
                userDefinedFunctions,
                useBatchQueries,
                cleanQueryLabels);

        createVQSRInputFromTableResult(storageAPIAvroReader);
    }

    static Map<String, String>  createQueryLabels(List<String> labelStringList) {
        // static labels are added to the query to make tracking this workflow easier downstream
        Map<String, String> labelForQuery = new HashMap<String, String>();
        labelForQuery.put("gvs_tool_name", "extract-features");
        labelForQuery.put("gvs_query_name", "extract-features");
        // add additional key value pair labels

        // Each resource can have multiple labels, up to a maximum of 64. -- labelStringList cannot be >62
        // The key portion of a label must be unique. However, you can use the same key with multiple resources.

        if ( labelStringList.size() > 62 ) {
            throw new UserException("Only 62 unique label keys are allowed per resource.");
            // BQ allows for 64, and we add two above
        }

        for (String labelMapString: labelStringList) {
            if ( !labelMapString.contains("=") ) {
                throw new UserException("All labels must contain a key and value in key=value format");
            }
            String[] label = labelMapString.split("=");
            String labelKey = label[0];
            String labelValue = label[1];
            // validate the label key and label value according to GCP Requirements for labels
            // Each label must be a key-value pair.
            // Keys have a minimum length of 1 character and a maximum length of 63 characters, and cannot be empty.
            // Values can be empty, and have a maximum length of 63 characters.
            // Keys and values can contain only lowercase letters, numeric characters, underscores, and dashes.
            // All characters must use UTF-8 encoding, and international characters are allowed.
            // Keys must start with a lowercase letter or international character.
            if ( labelKey.length() > 63 || labelKey.length() < 1 || labelKey.isEmpty()) {
                throw new UserException("Label key length must be between 1 and 63 characters");
            }
            if ( labelValue.length() > 63 ) {
                throw new UserException("Label value length must be less than 63 characters");
            }
            String labelRegex = "[a-z0-9_-]+$";
            if ( !labelKey.matches(labelRegex) || !labelValue.matches(labelRegex)  ) {
                throw new UserException("Label keys and values can contain only" +
                    " lowercase letters, numeric characters, underscores, and dashes");
            }
            labelForQuery.put(labelKey, labelValue);
        }

        return labelForQuery;
    }

    private void createVQSRInputFromTableResult(final GATKAvroReader avroReader) {
        final org.apache.avro.Schema schema = avroReader.getSchema();
        if ( schema.getField(LOCATION_FIELD_NAME) == null ) {
            throw new UserException("Records must contain a location column");
        }

        // TODO: this hardcodes a list of required fields... which this case doesn't require!
        // should genericize/parameterize so we can also validate the schema here
        // validateSchema(columnNames);

        SortingCollection<GenericRecord> sortingCollection = AvroSortingCollection.getAvroSortingCollection(schema, localSortMaxRecordsInRam, SchemaUtils.LOCATION_COMPARATOR);
        for ( final GenericRecord queryRow : avroReader ) {
            sortingCollection.add(queryRow);
        }
        sortingCollection.printTempFileStats();


        // NOTE: if OverlapDetector takes too long, try using RegionChecker from tws_sv_local_assembler
        final OverlapDetector<SimpleInterval> intervalsOverlapDetector = OverlapDetector.create(traversalIntervals);

        for ( final GenericRecord genericRow : sortingCollection ) {
            final ExtractFeaturesRecord row = new ExtractFeaturesRecord(genericRow);
            if ( intervalsOverlapDetector.overlapsAny(row) ) {
                processVQSRRecordForPosition(row);
            }
        }
    }

    private void processVQSRRecordForPosition(ExtractFeaturesRecord rec) {
        final long location = rec.getLocation();

        String contig = rec.getContig();
        int position = rec.getStart();
        // I don't understand why the other modes iterate through a list of all columns, and then
        //  switch based on the column name rather than just getting the columns desired directly (like I'm going to do here).
        // It might be because those field names are just prefixed versions of standard VCF fields?

        String ref = rec.getRef();
        String allele = rec.getAllele();

        if (allele == null || allele.equals("")) {
            logger.warn("SEVERE WARNING: skipping " + contig + ":" + position + "(location="+location+") because it has a null alternate allele!");
            return;
        }

        // Numbers are returned as Long (sci notation)
        Double qual = rec.getRawQual();

        Double raw_ref_ad = rec.getRefAD();  // if null, will return 0

        Float AS_MQRankSum = rec.getAsMQRankSum();

        Float AS_ReadPosRankSum = rec.getAsReadPosRankSum();

        Double raw_mq = rec.getRawMQ();
        Double raw_ad = rec.getRawAD();
        Double raw_ad_gt_1 = rec.getRawADGT1();


        int sb_ref_plus = rec.getSbRefPlus();
        int sb_ref_minus = rec.getSbRefMinus();
        int sb_alt_plus = rec.getSbAltPlus();
        int sb_alt_minus = rec.getSbAltMinus();

//        logger.info("processing " + contig + ":" + position);

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

        double qualapprox = rec.getQualApprox();

        builder.attribute(GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY, String.format("%.2f", as_qd) );
        builder.attribute(GATKVCFConstants.AS_FISHER_STRAND_KEY, String.format("%.3f", fs));
        builder.attribute(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY, String.format("%.2f", mq) );
        builder.attribute(GATKVCFConstants.AS_MAP_QUAL_RANK_SUM_KEY, AS_MQRankSum==null?".":String.format("%.3f", AS_MQRankSum) );
        builder.attribute(GATKVCFConstants.AS_READ_POS_RANK_SUM_KEY, AS_ReadPosRankSum==null?".":String.format("%.3f", AS_ReadPosRankSum));
        builder.attribute(GATKVCFConstants.AS_STRAND_ODDS_RATIO_KEY, String.format("%.3f", sor));
        builder.attribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY, String.format("%.3f", qualapprox));

//        From the warp JointGenotyping pipeline
//        # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
//        # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
        double excess_het_threshold = 54.69;

        int hets = rec.getNumHetSamples();
        int homvars = rec.getNumHomvarSamples();

        int samplesMinusVariants = numSamples - (hets + homvars);
        GenotypeCounts gcApprox = new GenotypeCounts(samplesMinusVariants, hets, homvars);
        double excessHetApprox = ExcessHet.calculateEH(gcApprox, numSamples).getRight();
        if (excessHetApprox > excess_het_threshold) {
            builder.filter(GATKVCFConstants.EXCESS_HET_KEY);
        }
        builder.attribute(GATKVCFConstants.EXCESS_HET_KEY, String.format("%.3f", excessHetApprox));

        if (rec.getDistinctAlleles() > excessAllelesThreshold) {
            builder.filter(GATKVCFConstants.EXCESS_ALLELES);
        }

        if (rec.getHqGenotypeSamples() < 1) {
            builder.filter(GATKVCFConstants.NO_HQ_GENOTYPES);
        }

        boolean isSNP = rec.getNumSnpAlleles() > 0;
        if (qualapprox < (isSNP ? SNP_QUAL_THRESHOLD : INDEL_QUAL_THRESHOLD)) {
            builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);
        }

        VariantContext vc = builder.make();
        vcfWriter.add(vc);
        progressMeter.update(vc);
    }
}
