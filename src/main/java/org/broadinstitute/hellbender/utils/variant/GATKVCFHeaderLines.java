package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashMap;
import java.util.Map;

import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.*;

/**
 * This class contains the {@link VCFHeaderLine} definitions for the annotation keys in {@link GATKVCFConstants}.
 * VCF-standard header lines are in {@link VCFStandardHeaderLines}, in htsjdk
 */
public class GATKVCFHeaderLines {

    public static VCFInfoHeaderLine getInfoLine(final String id) { return infoLines.get(id); }
    public static VCFFormatHeaderLine getFormatLine(final String id) { return formatLines.get(id); }
    public static VCFFilterHeaderLine getFilterLine(final String id) { return filterLines.get(id); }

    private static final Map<String, VCFInfoHeaderLine> infoLines = new HashMap<>(60);
    private static final Map<String, VCFFormatHeaderLine> formatLines = new HashMap<>(25);
    private static final Map<String, VCFFilterHeaderLine> filterLines = new HashMap<>(2);

    private static void addFormatLine(final VCFFormatHeaderLine line) {
        Utils.nonNull(line);
        formatLines.put(line.getID(), line);
    }

    private static void addInfoLine(final VCFInfoHeaderLine line) {
        Utils.nonNull(line);
        infoLines.put(line.getID(), line);
    }

    private static void addFilterLine(final VCFFilterHeaderLine line) {
        Utils.nonNull(line);
        filterLines.put(line.getID(), line);
    }

    static {
        addFilterLine(new VCFFilterHeaderLine(LOW_QUAL_FILTER_NAME, "Low quality"));

        // M2-related filters
        addFilterLine(new VCFFilterHeaderLine(ALT_ALLELE_IN_NORMAL_FILTER_NAME, "Evidence seen in the normal sample"));
        addFilterLine(new VCFFilterHeaderLine(CLUSTERED_EVENTS_FILTER_NAME, "Clustered events observed in the tumor"));
        addFilterLine(new VCFFilterHeaderLine(GERMLINE_RISK_FILTER_NAME, "Evidence indicates this site is germline, not somatic"));
        addFilterLine(new VCFFilterHeaderLine(HOMOLOGOUS_MAPPING_EVENT_FILTER_NAME, "More than three events were observed in the tumor"));
        addFilterLine(new VCFFilterHeaderLine(MULTI_EVENT_ALT_ALLELE_IN_NORMAL_FILTER_NAME, "Multiple events observed in tumor and normal"));
        addFilterLine(new VCFFilterHeaderLine(PON_FILTER_NAME, "Seen in at least 2 samples in the panel of normals"));
        addFilterLine(new VCFFilterHeaderLine(TUMOR_LOD_FILTER_NAME, "Tumor does not meet likelihood threshold"));
        addFilterLine(new VCFFilterHeaderLine(STR_CONTRACTION_FILTER_NAME, "Site filtered due to contraction of short tandem repeat region"));

        addFormatLine(new VCFFormatHeaderLine(ALLELE_BALANCE_KEY, 1, VCFHeaderLineType.Float, "Allele balance for each het genotype"));
        addFormatLine(new VCFFormatHeaderLine(MAPPING_QUALITY_ZERO_BY_SAMPLE_KEY, 1, VCFHeaderLineType.Integer, "Number of Mapping Quality Zero Reads per sample"));
        addFormatLine(new VCFFormatHeaderLine(MLE_PER_SAMPLE_ALLELE_COUNT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Maximum likelihood expectation (MLE) for the alternate allele count, in the same order as listed, for each individual sample"));
        addFormatLine(new VCFFormatHeaderLine(MLE_PER_SAMPLE_ALLELE_FRACTION_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Maximum likelihood expectation (MLE) for the alternate allele fraction, in the same order as listed, for each individual sample"));
        addFormatLine(new VCFFormatHeaderLine(STRAND_COUNT_BY_SAMPLE_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Number of reads on the forward and reverse strand supporting each allele (including reference)"));
        addFormatLine(new VCFFormatHeaderLine(STRAND_BIAS_BY_SAMPLE_KEY, 4, VCFHeaderLineType.Integer, "Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias."));
        addFormatLine(new VCFFormatHeaderLine(MLE_PER_SAMPLE_ALLELE_COUNT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Maximum likelihood expectation (MLE) for the alternate allele count, in the same order as listed, for each individual sample"));
        addFormatLine(new VCFFormatHeaderLine(MLE_PER_SAMPLE_ALLELE_FRACTION_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Maximum likelihood expectation (MLE) for the alternate allele fraction, in the same order as listed, for each individual sample"));
        addFormatLine(new VCFFormatHeaderLine(PL_FOR_ALL_SNP_ALLELES_KEY, 10, VCFHeaderLineType.Integer, "Phred-scaled genotype likelihoods for all 4 possible bases regardless of whether there is statistical evidence for them. Ordering is always PL for AA AC CC GA GC GG TA TC TG TT."));
        addFormatLine(new VCFFormatHeaderLine(HAPLOTYPE_CALLER_PHASING_ID_KEY, 1, VCFHeaderLineType.String, "Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group"));
        addFormatLine(new VCFFormatHeaderLine(HAPLOTYPE_CALLER_PHASING_GT_KEY, 1, VCFHeaderLineType.String, "Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another"));
        addFormatLine(new VCFFormatHeaderLine(MIN_DP_FORMAT_KEY, 1, VCFHeaderLineType.Integer, "Minimum DP observed within the GVCF block"));
        addFormatLine(new VCFFormatHeaderLine(REFERENCE_GENOTYPE_QUALITY, 1, VCFHeaderLineType.Integer, "Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)"));
        addFormatLine(new VCFFormatHeaderLine(TRANSMISSION_PROBABILITY_KEY, 1, VCFHeaderLineType.Integer, "Phred score of the genotype combination and phase given that the genotypes are correct"));
        addFormatLine(new VCFFormatHeaderLine(RBP_HAPLOTYPE_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Read-backed phasing haplotype identifiers"));
        addFormatLine(new VCFFormatHeaderLine(AVG_INTERVAL_DP_BY_SAMPLE_KEY, 1, VCFHeaderLineType.Float, "Average sample depth across the interval. Sum of the sample specific depth in all loci divided by interval size."));
        addFormatLine(new VCFFormatHeaderLine(LOW_COVERAGE_LOCI, 1, VCFHeaderLineType.Integer, "Number of loci for this sample, in this interval with low coverage (below the minimum coverage) but not zero."));
        addFormatLine(new VCFFormatHeaderLine(ZERO_COVERAGE_LOCI, 1, VCFHeaderLineType.Integer, "Number of loci for this sample, in this interval with zero coverage."));
        addFormatLine(new VCFFormatHeaderLine(PHRED_SCALED_POSTERIORS_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Phred-scaled Posterior Genotype Probabilities"));
        addFormatLine(new VCFFormatHeaderLine(JOINT_LIKELIHOOD_TAG_NAME, 1, VCFHeaderLineType.Integer, "Phred-scaled joint likelihood of the genotype combination (before applying family priors)"));
        addFormatLine(new VCFFormatHeaderLine(JOINT_POSTERIOR_TAG_NAME, 1, VCFHeaderLineType.Integer, "Phred-scaled joint posterior probability of the genotype combination (after applying family priors)"));

        // M2-related info lines
        addFormatLine(new VCFFormatHeaderLine(ALLELE_FRACTION_KEY, 1, VCFHeaderLineType.Float, "Allele fraction of the event in the tumor"));
        addFormatLine(new VCFFormatHeaderLine(OXOG_ALT_F1R2_KEY, 1, VCFHeaderLineType.Integer, "Count of reads in F1R2 pair orientation supporting the alternate allele"));
        addFormatLine(new VCFFormatHeaderLine(OXOG_ALT_F2R1_KEY, 1, VCFHeaderLineType.Integer, "Count of reads in F2R1 pair orientation supporting the alternate allele"));
        addFormatLine(new VCFFormatHeaderLine(OXOG_REF_F1R2_KEY, 1, VCFHeaderLineType.Integer, "Count of reads in F1R2 pair orientation supporting the reference allele"));
        addFormatLine(new VCFFormatHeaderLine(OXOG_REF_F2R1_KEY, 1, VCFHeaderLineType.Integer, "Count of reads in F2R1 pair orientation supporting the reference allele"));
        addFormatLine(new VCFFormatHeaderLine(OXOG_FRACTION_KEY, 1, VCFHeaderLineType.Float, "Fraction of alt reads indicating OxoG error"));
        addFormatLine(new VCFFormatHeaderLine(QUALITY_SCORE_SUM_KEY, 1, VCFHeaderLineType.Integer, "Sum of base quality scores for each allele"));


        addInfoLine(new VCFInfoHeaderLine(MLE_ALLELE_COUNT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed"));
        addInfoLine(new VCFInfoHeaderLine(MLE_ALLELE_FREQUENCY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed"));
        addInfoLine(new VCFInfoHeaderLine(DOWNSAMPLED_KEY, 0, VCFHeaderLineType.Flag, "Were any of the samples downsampled?"));
        addInfoLine(new VCFInfoHeaderLine(ALLELE_BALANCE_HET_KEY, 1, VCFHeaderLineType.Float, "Allele Balance for heterozygous calls (ref/(ref+alt))"));
        addInfoLine(new VCFInfoHeaderLine(ALLELE_BALANCE_HOM_KEY, 1, VCFHeaderLineType.Float, "Allele Balance for homozygous calls (A/(A+O)) where A is the allele (ref or alt) and O is anything other"));
        addInfoLine(new VCFInfoHeaderLine(NON_DIPLOID_RATIO_KEY, 1, VCFHeaderLineType.Float, "Overall non-diploid ratio (alleles/(alleles+non-alleles))"));
        addInfoLine(new VCFInfoHeaderLine(BASE_COUNTS_KEY, 4, VCFHeaderLineType.Integer, "Counts of each base"));
        addInfoLine(new VCFInfoHeaderLine(LOW_MQ_KEY, 3, VCFHeaderLineType.Float, "3-tuple: <fraction of reads with MQ=0>,<fraction of reads with MQ<=10>,<total number of reads>"));
        addInfoLine(new VCFInfoHeaderLine(N_BASE_COUNT_KEY, 1, VCFHeaderLineType.Float, "Percentage of N bases in the pileup"));
        addInfoLine(new VCFInfoHeaderLine(BASE_QUAL_RANK_SUM_KEY, 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities"));
        addInfoLine(new VCFInfoHeaderLine(CLIPPING_RANK_SUM_KEY, 1, VCFHeaderLineType.Float, "Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases"));
        addInfoLine(new VCFInfoHeaderLine(FISHER_STRAND_KEY, 1, VCFHeaderLineType.Float, "Phred-scaled p-value using Fisher's exact test to detect strand bias"));
        addInfoLine(new VCFInfoHeaderLine(GC_CONTENT_KEY, 1, VCFHeaderLineType.Float, "GC content around the variant (see docs for window size details)"));
        addInfoLine(new VCFInfoHeaderLine(NOCALL_CHROM_KEY, 1, VCFHeaderLineType.Integer, "Number of no-called samples"));
        addInfoLine(new VCFInfoHeaderLine(GQ_MEAN_KEY, 1, VCFHeaderLineType.Float, "Mean of all GQ values"));
        addInfoLine(new VCFInfoHeaderLine(GQ_STDEV_KEY, 1, VCFHeaderLineType.Float, "Standard deviation of all GQ values"));
        addInfoLine(new VCFInfoHeaderLine(HAPLOTYPE_SCORE_KEY, 1, VCFHeaderLineType.Float, "Consistency of the site with at most two segregating haplotypes"));
        addInfoLine(new VCFInfoHeaderLine(HARDY_WEINBERG_KEY, 1, VCFHeaderLineType.Float, "Phred-scaled p-value for Hardy-Weinberg violation"));
        addInfoLine(new VCFInfoHeaderLine(HOMOPOLYMER_RUN_KEY, 1, VCFHeaderLineType.Integer, "Largest Contiguous Homopolymer Run of Variant Allele In Either Direction"));
        addInfoLine(new VCFInfoHeaderLine(INBREEDING_COEFFICIENT_KEY, 1, VCFHeaderLineType.Float, "Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation"));
        addInfoLine(new VCFInfoHeaderLine(EXCESS_HET_KEY, 1, VCFHeaderLineType.Float, "Phred-scaled p-value for exact test of excess heterozygosity"));
        addInfoLine(new VCFInfoHeaderLine(LIKELIHOOD_RANK_SUM_KEY, 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt Vs. Ref haplotype likelihoods"));
        addInfoLine(new VCFInfoHeaderLine(MAP_QUAL_RANK_SUM_KEY, 1, VCFHeaderLineType.Float, "Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities"));
        addInfoLine(new VCFInfoHeaderLine(MENDEL_VIOLATION_LR_KEY, 1, VCFHeaderLineType.Float, "Mendelian violation likelihood ratio: L[MV] - L[No MV]"));
        addInfoLine(new VCFInfoHeaderLine(HI_CONF_DENOVO_KEY, 1, VCFHeaderLineType.String, "High confidence possible de novo mutation (GQ >= 20 for all trio members)=[comma-delimited list of child samples]"));
        addInfoLine(new VCFInfoHeaderLine(LO_CONF_DENOVO_KEY, 1, VCFHeaderLineType.String, "Low confidence possible de novo mutation (GQ >= 10 for child, GQ > 0 for parents)=[comma-delimited list of child samples]"));
        addInfoLine(new VCFInfoHeaderLine(QUAL_BY_DEPTH_KEY, 1, VCFHeaderLineType.Float, "Variant Confidence/Quality by Depth"));
        addInfoLine(new VCFInfoHeaderLine(READ_POS_RANK_SUM_KEY, 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias"));
        addInfoLine(new VCFInfoHeaderLine(SAMPLE_LIST_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "List of polymorphic samples"));
        addInfoLine(new VCFInfoHeaderLine(SPANNING_DELETIONS_KEY, 1, VCFHeaderLineType.Float, "Fraction of Reads Containing Spanning Deletions"));
        addInfoLine(new VCFInfoHeaderLine(STRAND_ODDS_RATIO_KEY, 1, VCFHeaderLineType.Float, "Symmetric Odds Ratio of 2x2 contingency table to detect strand bias"));
        addInfoLine(new VCFInfoHeaderLine(STR_PRESENT_KEY, 0, VCFHeaderLineType.Flag, "Variant is a short tandem repeat"));
        addInfoLine(new VCFInfoHeaderLine(REPEAT_UNIT_KEY, 1, VCFHeaderLineType.String, "Tandem repeat unit (bases)"));
        addInfoLine(new VCFInfoHeaderLine(REPEATS_PER_ALLELE_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Number of times tandem repeat unit is repeated, for each allele (including reference)"));
        addInfoLine(new VCFInfoHeaderLine(TRANSMISSION_DISEQUILIBRIUM_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Test statistic from Wittkowski transmission disequilibrium test."));
        addInfoLine(new VCFInfoHeaderLine(VARIANT_TYPE_KEY, 1, VCFHeaderLineType.String, "Variant type description"));
        addInfoLine(new VCFInfoHeaderLine(NUMBER_OF_DISCOVERED_ALLELES_KEY, 1, VCFHeaderLineType.Integer, "Number of alternate alleles discovered (but not necessarily genotyped) at this site"));
        addInfoLine(new VCFInfoHeaderLine(REFSAMPLE_DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Total reference sample depth"));
        addInfoLine(new VCFInfoHeaderLine(ORIGINAL_AC_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Original AC"));
        addInfoLine(new VCFInfoHeaderLine(ORIGINAL_AF_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Original AF"));
        addInfoLine(new VCFInfoHeaderLine(ORIGINAL_AN_KEY, 1, VCFHeaderLineType.Integer, "Original AN"));
        addInfoLine(new VCFInfoHeaderLine(ORIGINAL_DP_KEY, 1, VCFHeaderLineType.Integer, "Original DP"));
        addInfoLine(new VCFInfoHeaderLine(ORIGINAL_CONTIG_KEY, 1, VCFHeaderLineType.String, "Original contig name for the record"));
        addInfoLine(new VCFInfoHeaderLine(ORIGINAL_START_KEY, 1, VCFHeaderLineType.Integer, "Original start position for the record"));
        addInfoLine(new VCFInfoHeaderLine(VQS_LOD_KEY, 1, VCFHeaderLineType.Float, "Log odds ratio of being a true variant versus being false under the trained gaussian mixture model"));
        addInfoLine(new VCFInfoHeaderLine(CULPRIT_KEY, 1, VCFHeaderLineType.String, "The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out"));
        addInfoLine(new VCFInfoHeaderLine(POSITIVE_LABEL_KEY, 1, VCFHeaderLineType.Flag, "This variant was used to build the positive training set of good variants"));
        addInfoLine(new VCFInfoHeaderLine(NEGATIVE_LABEL_KEY, 1, VCFHeaderLineType.Flag, "This variant was used to build the negative training set of bad variants"));
        addInfoLine(new VCFInfoHeaderLine(RBP_INCONSISTENT_KEY, 0, VCFHeaderLineType.Flag, "Are the reads significantly haplotype-inconsistent?"));
        addInfoLine(new VCFInfoHeaderLine(GENOTYPE_AND_VALIDATE_STATUS_KEY, 1, VCFHeaderLineType.String, "Value from the validation VCF"));
        addInfoLine(new VCFInfoHeaderLine(AVG_INTERVAL_DP_KEY, 1, VCFHeaderLineType.Float, "Average depth across the interval. Sum of the depth in a loci divided by interval size."));
        addInfoLine(new VCFInfoHeaderLine(INTERVAL_GC_CONTENT_KEY, 1, VCFHeaderLineType.Float, "GC Content of the interval"));
        addInfoLine(new VCFInfoHeaderLine(GENOTYPE_PRIOR_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Genotype Likelihood Prior"));

        // M2-related info lines
        addInfoLine(new VCFInfoHeaderLine(EVENT_COUNT_IN_HAPLOTYPE_KEY, 1, VCFHeaderLineType.String, "Number of events in this haplotype"));
        addInfoLine(new VCFInfoHeaderLine(EVENT_DISTANCE_MAX_KEY, 1, VCFHeaderLineType.Integer, "Maximum distance between events in this active region"));
        addInfoLine(new VCFInfoHeaderLine(EVENT_DISTANCE_MIN_KEY, 1, VCFHeaderLineType.Integer, "Minimum distance between events in this active region"));
        addInfoLine(new VCFInfoHeaderLine(HAPLOTYPE_COUNT_KEY, 1, VCFHeaderLineType.String, "Number of haplotypes that support this variant"));
        addInfoLine(new VCFInfoHeaderLine(NORMAL_LOD_KEY, 1, VCFHeaderLineType.String, "Normal LOD score"));
        addInfoLine(new VCFInfoHeaderLine(PANEL_OF_NORMALS_COUNT_KEY, 1, VCFHeaderLineType.String, "Count from Panel of Normals"));
        addInfoLine(new VCFInfoHeaderLine(TUMOR_LOD_KEY, 1, VCFHeaderLineType.String, "Tumor LOD score"));
    }
}
