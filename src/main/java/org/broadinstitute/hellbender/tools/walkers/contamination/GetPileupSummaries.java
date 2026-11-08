package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * <p>Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results can be used with {@link CalculateContamination}.</p>
 *
 * <p>
 * The tool requires a <i>common</i> germline variant sites VCF, e.g. derived from the gnomAD resource, with population allele frequencies (AF) in the INFO field.
 * This resource must contain only biallelic SNPs and can be an eight-column sites-only VCF.
 * The tool ignores the filter status of the variant calls in this germline resource.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 *     In particular, the mutect_resources.wdl script prepares a suitable resource from a larger dataset. An example excerpt is shown.
 * </p>
 *
 * <pre>
 * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
 * chr6	29942512	.	G	C	2974860	VQSRTrancheSNP99.80to99.90	AF=0.063
 * chr6	29942517	.	C	A	2975860	VQSRTrancheSNP99.80to99.90	AF=0.062
 * chr6	29942525	.	G	C	2975600	VQSRTrancheSNP99.60to99.80	AF=0.063
 * chr6	29942547	rs114945359	G	C	15667700	PASS	AF=0.077
 * </pre>
 *
 * <h3>Usage examples</h3>
 *
 * <pre>
 * gatk GetPileupSummaries \
 *   -I tumor.bam \
 *   -V common_biallelic.vcf.gz \
 *   -L common_biallelic.vcf.gz \
 *   -O pileups.table
 * </pre>
 *
 * <pre>
 * gatk GetPileupSummaries \
 *   -I normal.bam \
 *   -V common_biallelic.vcf.gz \
 *   -L common_biallelic.vcf.gz \
 *   -O pileups.table
 * </pre>
 *
 * Although the sites (-L) and variants (-V) resources will often be identical, this need not be the case.  For example,
 * <pre>
 * gatk GetPileupSummaries \
 *   -I normal.bam \
 *   -V gnomad.vcf.gz \
 *   -L common_snps.interval_list \
 *   -O pileups.table
 * </pre>
 * attempts to get pileups at a list of common snps and emits output for those sites that are present in gnomAD, using the
 * allele frequencies from gnomAD.  Note that the sites may be a subset of the variants, the variants may be a subset of the sites,
 * or they may overlap partially.  In all cases pileup summaries are emitted for the overlap and nowhere else.  The most common use
 * case in which sites and variants differ is when the variants resources is a large file and the sites is an interval list subset from that file.
 *
 * <p>
 * GetPileupSummaries tabulates results into six columns as shown below.
 * The alt_count and allele_frequency correspond to the ALT allele in the germline resource.
 * </p>
 *
 * <pre>
 * contig	position	ref_count	alt_count	other_alt_count	allele_frequency
 * chr6	29942512	9	0	0	0.063
 * chr6	29942517	13	1	0	0.062
 * chr6	29942525	13	7	0	0.063
 * chr6	29942547	36	0	0	0.077
 * </pre>
 *
 * <p>
 * Note the default maximum population AF ({@code --maximum-population-allele-frequency} or {@code -max-af})
 * is set to 0.2, which limits the sites the tool considers to those in the variants resource file that have
 * AF of 0.2 or less. Likewise, the default minimum population AF ({@code --minimum-population-allele-frequency}
 * or {@code -min-af}) is set to 0.01, which limits the sites the tool considers to those in the variants resource
 * file that have AF of 0.01 or more.  {@link CalculateContamination} accepts the same arguments, so one can collect
 * a wider range of sites here and restrict the contamination estimate to the usual range there.
 * </p>
 *
 * <p>
 * Finally for those using mappers other than bwa mem or dragen-os, {@code --minimum-mapping-quality} threshold is
 * set to 50, which limits the usable reads that tool considers for generating pileups. Certain mappers are known to
 * assign scores less that this threshold even for the unique mappings. If you observe all empty results in your
 * summary file please adjust the {@code --minimum-mapping-quality} parameter according to your input files.
 * </p>
 *
 */
@CommandLineProgramProperties(
        summary = "Tabulates pileup metrics for inferring contamination",
        oneLineSummary = "Tabulates pileup metrics for inferring contamination",
        programGroup = CoverageAnalysisProgramGroup.class)
@DocumentedFeature
public class GetPileupSummaries extends LocusWalker {

    public static final String MAX_SITE_AF_LONG_NAME = "maximum-population-allele-frequency";
    public static final String MIN_SITE_AF_LONG_NAME = "minimum-population-allele-frequency";
    public static final String MAX_SITE_AF_SHORT_NAME = "max-af";
    public static final String MIN_SITE_AF_SHORT_NAME = "min-af";

    private static final double DEFAULT_MIN_POPULATION_AF = 0.01;
    private static final double DEFAULT_MAX_POPULATION_AF = 0.2;
    private static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 50;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table", optional=false)
    private File outputTable;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants and allele frequencies")
    public FeatureInput<VariantContext> variants;

    @Argument(fullName = MIN_SITE_AF_LONG_NAME,
            shortName = MIN_SITE_AF_SHORT_NAME,
            doc = "Minimum population allele frequency of sites to consider.  A low value increases accuracy at the expense of speed.", optional = true)
    private double minPopulationAlleleFrequency = DEFAULT_MIN_POPULATION_AF;

    @Argument(fullName = MAX_SITE_AF_LONG_NAME,
            shortName = MAX_SITE_AF_SHORT_NAME,
            doc = "Maximum population allele frequency of sites to consider.", optional = true)
    private double maxPopulationAlleleFrequency = DEFAULT_MAX_POPULATION_AF;

    private PileupSummaryCollector collector;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return false;
    }

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public boolean requiresFeatures() {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = new ArrayList<>();
        filters.add(new MappingQualityReadFilter(DEFAULT_MINIMUM_MAPPING_QUALITY));
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.PRIMARY_LINE);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);
        filters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE);
        filters.add(ReadFilterLibrary.GOOD_CIGAR);
        filters.add(new WellformedReadFilter());
        return filters;
    }

    @Override
    public void onTraversalStart() {
        PileupSummaryCollector.validateHeader((VCFHeader) getHeaderForFeatures(variants));

        // the engine has already filtered the reads, and the table counts the reads of every sample under the first sample's name
        collector = new PileupSummaryCollector(variants, minPopulationAlleleFrequency, maxPopulationAlleleFrequency);
        final String sampleName = ReadUtils.getSamplesFromHeader(getHeaderForReads()).stream().findFirst().get();
        collector.addTable(outputTable, sampleName, read -> true);
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        collector.process(alignmentContext, featureContext);
    }

    @Override
    public Object onTraversalSuccess() {
        collector.validateAlleleFrequencies();
        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        if (collector != null) {
            collector.close();
        }
    }
}
