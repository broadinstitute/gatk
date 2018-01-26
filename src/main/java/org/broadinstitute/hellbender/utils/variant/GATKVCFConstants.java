package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.Allele;

import java.util.Arrays;
import java.util.List;

/**
 * This class contains any constants (primarily FORMAT/INFO keys) in VCF files used by the GATK.
 * Note that VCF-standard constants are in VCFConstants, in htsjdk.  Keys in header lines should
 * have matching entries in GATKVCFHeaderLines
 */
public final class GATKVCFConstants {

    public static final String CONTIG_ID_KEY =                      "ID";
    public static final String CONTIG_LENGTH_KEY =                  "length";
    public static final String ASSEMBLY_NAME_KEY =                  "assembly";

    //INFO keys
    public static final String ALLELE_SPECIFIC_PREFIX =             "AS_";
    public static final String AS_FILTER_STATUS_KEY =               "AS_FilterStatus";
    public static final String RAW_RMS_MAPPING_QUALITY_KEY =        "RAW_MQ";
    public static final String AS_RMS_MAPPING_QUALITY_KEY =         "AS_MQ";
    public static final String AS_RAW_RMS_MAPPING_QUALITY_KEY =     "AS_RAW_MQ";
    public static final String AS_CULPRIT_KEY =                     "AS_culprit";
    public static final String AS_VQS_LOD_KEY =                     "AS_VQSLOD";
    public static final String ALLELE_BALANCE_HET_KEY =             "ABHet";
    public static final String ALLELE_BALANCE_HOM_KEY =             "ABHom";
    public static final String ORIGINAL_AC_KEY =                    "AC_Orig"; //SelectVariants
    public static final String ORIGINAL_AF_KEY =                    "AF_Orig"; //SelectVariants
    public static final String ORIGINAL_AN_KEY =                    "AN_Orig"; //SelectVariants
    public static final String BASE_COUNTS_KEY =                    "BaseCounts";
    public static final String BASE_QUAL_RANK_SUM_KEY =             "BaseQRankSum";
    public static final String AS_BASE_QUAL_RANK_SUM_KEY =          "AS_BaseQRankSum";
    public static final String AS_RAW_BASE_QUAL_RANK_SUM_KEY =      "AS_RAW_BaseQRankSum";
    public static final String GENOTYPE_AND_VALIDATE_STATUS_KEY =   "callStatus";
    public static final String CLIPPING_RANK_SUM_KEY =              "ClippingRankSum";
    public static final String CULPRIT_KEY =                        "culprit";
    public static final String SPANNING_DELETIONS_KEY =             "Dels";
    public static final String ORIGINAL_DP_KEY =                    "DP_Orig"; //SelectVariants
    public static final String DOWNSAMPLED_KEY =                    "DS";
    public static final String EVENT_COUNT_IN_HAPLOTYPE_KEY =       "ECNT"; //M2
    public static final String FISHER_STRAND_KEY =                  "FS";
    public static final String AS_FISHER_STRAND_KEY =               "AS_FS";
    public static final String AS_SB_TABLE_KEY =                    "AS_SB_TABLE";
    public static final String GC_CONTENT_KEY =                     "GC";
    public static final String GQ_MEAN_KEY =                        "GQ_MEAN";
    public static final String GQ_STDEV_KEY =                       "GQ_STDDEV";
    public static final String HAPLOTYPE_SCORE_KEY =                "HaplotypeScore";
    public static final String HI_CONF_DENOVO_KEY =                 "hiConfDeNovo";
    public static final String HOMOPOLYMER_RUN_KEY =                "HRun";
    public static final String HARDY_WEINBERG_KEY =                 "HW";
    public static final String AVG_INTERVAL_DP_KEY =                "IDP"; //DiagnoseTargets
    public static final String INTERVAL_GC_CONTENT_KEY =            "IGC";
    public static final String INBREEDING_COEFFICIENT_KEY =         "InbreedingCoeff";
    public static final String AS_INBREEDING_COEFFICIENT_KEY =      "AS_InbreedingCoeff";
    public static final String EXCESS_HET_KEY =                     "ExcessHet";
    public static final String AS_HETEROZYGOSITY_KEY =              "AS_InbreedingCoeff";
    public static final String LIKELIHOOD_RANK_SUM_KEY =            "LikelihoodRankSum";
    public static final String LO_CONF_DENOVO_KEY =                 "loConfDeNovo";
    public static final String LOW_MQ_KEY =                         "LowMQ";
    public static final String MLE_ALLELE_COUNT_KEY =               "MLEAC";
    public static final String MLE_ALLELE_FREQUENCY_KEY =           "MLEAF";
    public static final String MLE_PER_SAMPLE_ALLELE_COUNT_KEY =    "MLPSAC";
    public static final String MLE_PER_SAMPLE_ALLELE_FRACTION_KEY = "MLPSAF";
    public static final String MAP_QUAL_RANK_SUM_KEY =              "MQRankSum";
    public static final String RAW_MAP_QUAL_RANK_SUM_KEY =          "RAW_MQRankSum";
    public static final String AS_MAP_QUAL_RANK_SUM_KEY =           "AS_MQRankSum";
    public static final String AS_RAW_MAP_QUAL_RANK_SUM_KEY =       "AS_RAW_MQRankSum";
    public static final String MENDEL_VIOLATION_LR_KEY =            "MVLR";
    public static final String NOCALL_CHROM_KEY =                   "NCC";
    public static final String NUMBER_OF_DISCOVERED_ALLELES_KEY =   "NDA";
    public static final String NEGATIVE_LABEL_KEY =                 "NEGATIVE_TRAIN_SITE";
    public static final String NON_DIPLOID_RATIO_KEY =              "OND";
    public static final String ORIGINAL_CONTIG_KEY =                "OriginalChr"; //LiftoverVariants
    public static final String ORIGINAL_START_KEY =                 "OriginalStart"; //LiftoverVariants
    public static final String N_BASE_COUNT_KEY =                   "PercentNBase";
    public static final String RBP_INCONSISTENT_KEY =               "PhasingInconsistent"; //ReadBackedPhasing
    public static final String GENOTYPE_PRIOR_KEY =                 "PG";
    public static final String POSITIVE_LABEL_KEY =                 "POSITIVE_TRAIN_SITE";
    public static final String QUAL_BY_DEPTH_KEY =                  "QD";
    public static final String AS_QUAL_BY_DEPTH_KEY =               "AS_QD";
    public static final String AS_QUAL_KEY =                        "AS_QUAL";
    public static final String READ_POS_RANK_SUM_KEY =              "ReadPosRankSum";
    public static final String AS_READ_POS_RANK_SUM_KEY =           "AS_ReadPosRankSum";
    public static final String AS_RAW_READ_POS_RANK_SUM_KEY =       "AS_RAW_ReadPosRankSum";
    public static final String REFSAMPLE_DEPTH_KEY =                "REFDEPTH";
    public static final String REPEATS_PER_ALLELE_KEY =             "RPA";
    public static final String REPEAT_UNIT_KEY =                    "RU";
    public static final String SAMPLE_LIST_KEY =                    "Samples";
    public static final String STRAND_ODDS_RATIO_KEY =              "SOR";
    public static final String AS_STRAND_ODDS_RATIO_KEY =           "AS_SOR";
    public static final String STR_PRESENT_KEY =                    "STR";
    public static final String TRANSMISSION_DISEQUILIBRIUM_KEY =    "TDT";
    public static final String VARIANT_TYPE_KEY =                   "VariantType";
    public static final String VQS_LOD_KEY =                        "VQSLOD";
    public static final String CNN_1D_KEY =                         "CNN_1D";
    public static final String F1R2_KEY =                           "F1R2";
    public static final String F2R1_KEY =                           "F2R1";

    // Mutect2-specific INFO keys
    public static final String TUMOR_LOD_KEY =                      "TLOD";
    public static final String NORMAL_LOD_KEY =                     "NLOD";
    public static final String IN_PON_VCF_ATTRIBUTE =               "IN_PON";
    public static final String NORMAL_ARTIFACT_LOD_ATTRIBUTE =      "N_ART_LOD";
    public static final String POPULATION_AF_VCF_ATTRIBUTE =        "POP_AF";
    public static final String GERMLINE_POSTERIORS_VCF_ATTRIBUTE =  "P_GERMLINE";
    public static final List<String> STANDARD_MUTECT_INFO_FIELDS = Arrays.asList(NORMAL_LOD_KEY, TUMOR_LOD_KEY, NORMAL_ARTIFACT_LOD_ATTRIBUTE,
            EVENT_COUNT_IN_HAPLOTYPE_KEY, IN_PON_VCF_ATTRIBUTE, POPULATION_AF_VCF_ATTRIBUTE, GERMLINE_POSTERIORS_VCF_ATTRIBUTE);

    //FORMAT keys
    public static final String ALLELE_BALANCE_KEY =                 "AB";
    public static final String ALLELE_FRACTION_KEY =                "AF"; //M2
    public static final String PL_FOR_ALL_SNP_ALLELES_KEY =         "APL";
    public static final String RBP_HAPLOTYPE_KEY =                  "HP"; //ReadBackedPhasing
    public static final String AVG_INTERVAL_DP_BY_SAMPLE_KEY =      "IDP"; //DiagnoseTargets
    public static final String JOINT_LIKELIHOOD_TAG_NAME =          "JL"; //FamilyLikelihoodsUtils
    public static final String JOINT_POSTERIOR_TAG_NAME =           "JP"; //FamilyLikelihoodsUtils
    public static final String LOW_COVERAGE_LOCI =                  "LL"; //DiagnoseTargets
    public final static String MIN_DP_FORMAT_KEY =                  "MIN_DP";
    public static final String MAPPING_QUALITY_ZERO_BY_SAMPLE_KEY = "MQ0";
    public static final String HAPLOTYPE_CALLER_PHASING_GT_KEY =    "PGT";
    public static final String HAPLOTYPE_CALLER_PHASING_ID_KEY =    "PID";
    public static final String PHRED_SCALED_POSTERIORS_KEY =        "PP"; //FamilyLikelihoodsUtils / PosteriorLikelihoodsUtils
    public static final String REFERENCE_GENOTYPE_QUALITY =         "RGQ";
    public static final String STRAND_COUNT_BY_SAMPLE_KEY =         "SAC";
    public static final String STRAND_BIAS_BY_SAMPLE_KEY =          "SB";
    public final static String TRANSMISSION_PROBABILITY_KEY =       "TP"; //PhaseByTransmission
    public static final String ZERO_COVERAGE_LOCI =                 "ZL"; //DiagnoseTargets

    //FILTERS
    /* Note that many filters used throughout GATK (most notably in VariantRecalibration) are dynamic,
       their names (or descriptions) depend on some threshold.  Those filters are not included here
     */
    public static final String CLUSTERED_EVENTS_FILTER_NAME =                 "clustered_events"; //M2
    public static final String GERMLINE_RISK_FILTER_NAME =                    "germline_risk"; //M2
    public static final String LOW_QUAL_FILTER_NAME =                         "LowQual";
    public static final String PON_FILTER_NAME =                              "panel_of_normals"; //M2
    public static final String STR_CONTRACTION_FILTER_NAME =                  "str_contraction"; //M2
    public static final String TUMOR_LOD_FILTER_NAME =                        "t_lod"; //M2
    public static final String MULTIALLELIC_FILTER_NAME =                     "multiallelic"; //M2
    public static final String STRAND_ARTIFACT_FILTER_NAME =                  "strand_artifact"; // M2
    public static final String DUPLICATED_EVIDENCE_FILTER_NAME =              "duplicate_evidence";
    public final static String ARTIFACT_IN_NORMAL_FILTER_NAME =                "artifact_in_normal";
    public final static String MEDIAN_BASE_QUALITY_FILTER_NAME =                "base_quality";
    public final static String MEDIAN_MAPPING_QUALITY_FILTER_NAME =             "mapping_quality";
    public final static String MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME = "fragment_length";
    public final static String READ_POSITION_FILTER_NAME =                      "read_position";
    public final static String CONTAMINATION_FILTER_NAME =                      "contamination";

    public static final List<String> MUTECT_FILTER_NAMES = Arrays.asList(STR_CONTRACTION_FILTER_NAME,
            PON_FILTER_NAME, CLUSTERED_EVENTS_FILTER_NAME, TUMOR_LOD_FILTER_NAME, GERMLINE_RISK_FILTER_NAME,
            MULTIALLELIC_FILTER_NAME, STRAND_ARTIFACT_FILTER_NAME, ARTIFACT_IN_NORMAL_FILTER_NAME,
            MEDIAN_BASE_QUALITY_FILTER_NAME, MEDIAN_MAPPING_QUALITY_FILTER_NAME,
            MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME,
            READ_POSITION_FILTER_NAME, CONTAMINATION_FILTER_NAME, DUPLICATED_EVIDENCE_FILTER_NAME);

    // Symbolic alleles
    public final static String SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG = "ALT";
    public final static String NON_REF_SYMBOLIC_ALLELE_NAME = "NON_REF";
    public final static String SPANNING_DELETION_SYMBOLIC_ALLELE_NAME_DEPRECATED = "*:DEL";
    public final static Allele SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED = Allele.create("<" + SPANNING_DELETION_SYMBOLIC_ALLELE_NAME_DEPRECATED + ">", false); // represents any possible spanning deletion allele at this si


    public static boolean isSpanningDeletion(final Allele allele){
        return allele.equals(Allele.SPAN_DEL) || allele.equals(SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED);
    }
}
