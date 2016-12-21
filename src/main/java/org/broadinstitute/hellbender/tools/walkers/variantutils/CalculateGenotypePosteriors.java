package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.samples.*;
import org.broadinstitute.hellbender.utils.variant.*;

import java.io.File;
import java.util.*;

/**
 * Calculate genotype posterior probabilities given panel data
 *
 * <p>
 * Given a VCF with genotype likelihoods from the HaplotypeCaller, UnifiedGenotyper, or another source which provides
 * <b>unbiased</b> genotype likelihoods, calculate the posterior genotype state and probability given allele frequency
 * information from both the samples themselves and input VCFs describing allele frequencies in related populations.</p>
 *
 * <p>The AF field will not be used in this calculation as it does not provide a way to estimate the confidence interval
 * or uncertainty around the allele frequency, while AN provides this necessary information. This uncertainty is
 * modeled by a Dirichlet distribution: that is, the frequency is known up to a Dirichlet distribution with
 * parameters AC1+q,AC2+q,...,(AN-AC1-AC2-...)+q, where "q" is the global frequency prior (typically q << 1). The
 * genotype priors applied then follow a Dirichlet-Multinomial distribution, where 2 alleles per sample are drawn
 * independently. This assumption of independent draws is the assumption Hardy-Weinberg Equilibrium. Thus, HWE is
 * imposed on the likelihoods as a result of CalculateGenotypePosteriors.</p>
 *
 * <h3>Input</h3>
 * <p>
 *     <ul>
 *         <li>A VCF with genotype likelihoods, and optionally genotypes, AC/AN fields, or MLEAC/AN fields</li>
 *         <li>(Optional) A PED pedigree file containing the description of the individuals relationships.</li>
 *     </ul>
 * </p>
 *
 * <p>
 * A collection of VCFs to use for informing allele frequency priors. Each VCF must have one of
 * </p>
 * <ul>
 *     <li>AC field and AN field</li>
 *     <li>MLEAC field and AN field</li>
 *     <li>genotypes</li>
 * </ul>
 * </p>
 *
 * <h3>Output</h3>
 * <p>A new VCF with:</p>
 * <ul>
 *     <li>Genotype posteriors added to the genotype fields ("PP")</li>
 *     <li>Genotypes and GQ assigned according to these posteriors</li>
 *     <li>Per-site genotype priors added to the INFO field ("PG")</li>
 *     <li>(Optional) Per-site, per-trio joint likelihoods (JL) and joint posteriors (JL) given as Phred-scaled probability
 *  of all genotypes in the trio being correct based on the PLs for JL and the PPs for JP. These annotations are added to
 *  the genotype fields.</li>
 * </ul>
 *
 * <h3>Notes</h3>
 * <p>
 * Using the default behavior, priors will only be applied for each variants (provided each variant has at least 10
 * called samples.) SNP sites in the input callset that have a SNP at the matching site in the supporting VCF will have
 * priors applied based on the AC from the supporting samples and the input callset (unless the --ignoreInputSamples
 * flag is used). If the site is not called in the supporting VCF, priors will be applied using the discovered AC from
 * the input samples (unless the --discoveredACpriorsOff flag is used). Flat priors are applied for any non-SNP sites in
 * the input callset.
 * </p>
 *
 * <h3>Usage examples</h3>
 * <h4>Inform the genotype assignment of NA12878 using the 1000G Euro panel</h4>
 * <pre>
 * ./gatk-launch \
 *   CalculateGenotypePosteriors \
 *   -V NA12878.wgs.HC.vcf \
 *   -supporting 1000G_EUR.genotypes.combined.vcf \
 *   -O NA12878.wgs.HC.posteriors.vcf
 * </pre>
 *
 * <h4>Refine the genotypes of a large panel based on the discovered allele frequency</h4>
 * <pre>
 * ./gatk-launch \
 *   CalculateGenotypePosteriors \
 *   -V input.vcf \
 *   -O output.withPosteriors.vcf
 * </pre>
 *
 * <h4>Apply frequency and HWE-based priors to the genotypes of a family without including the family allele counts
 * in the allele frequency estimates the genotypes of a large panel based on the discovered allele frequency</h4>
 * <pre>
 * ./gatk-launch \
 *   CalculateGenotypePosteriors \
 *   -V input.vcf \
 *   -O output.withPosteriors.vcf \
 *   --ignoreInputSamples
 * </pre>
 *
 * <h4>Calculate the posterior genotypes of a callset, and impose that a variant *not seen* in the external panel
 * is tantamount to being AC=0, AN=100 within that panel</h4>
 * <pre>
 * ./gatk-launch \
 *   CalculateGenotypePosteriors \
 *   -supporting external.panel.vcf \
 *   -V input.vcf \
 *   -O output.withPosteriors.vcf \
 *   --numRefSamplesIfNoCall 100
 * </pre>
 *
 * <h4>Apply only family priors to a callset</h4>
 * <pre>
 * ./gatk-launch \
 *   CalculateGenotypePosteriors \
 *   -V input.vcf \
 *   --skipPopulationPriors \
 *   -ped family.ped \
 *   -O output.withPosteriors.vcf
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Given a VCF with genotype likelihoods from the HaplotypeCaller, UnifiedGenotyper, or another source which provides" +
                  "unbiased genotype likelihoods, calculate the posterior genotype state and likelihood given allele frequency" +
                  "information from both the samples themselves and input VCFs describing allele frequencies in related populations.",
        oneLineSummary = "Calculate genotype posterior probabilities given panel data",
        programGroup = VariantProgramGroup.class
)
public final class CalculateGenotypePosteriors extends VariantWalker {

    private static final Logger logger = LogManager.getLogger(CalculateGenotypePosteriors.class);

    /**
     * Supporting external panels. Allele counts from these panels (taken from AC,AN or MLEAC,AN or raw genotypes) will
     * be used to inform the frequency distribution underlying the genotype priors. These files must be VCF 4.2 spec or later.
     */
    @Argument(fullName="supporting", shortName = "supporting", doc="Other callsets to use in generating genotype posteriors", optional=true)
    public List<FeatureInput<VariantContext>> supportVariants = new ArrayList<>();

    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String out = null;

    /**
     * The global prior of a variant site -- i.e. the expected allele frequency distribution knowing only that N alleles
     * exist, and having observed none of them. This is the "typical" 1/x trend, modeled here as not varying
     * across alleles. The calculation for this parameter is (Effective population size) * (steady state mutation rate)
     *
     */
     @Argument(fullName="globalPrior",shortName="G",doc="The global Dirichlet prior parameters for the allele frequency",optional=true)
     public double globalPrior = HomoSapiensConstants.SNP_HETEROZYGOSITY;

    /**
     * The mutation prior -- i.e. the probability that a new mutation occurs. Sensitivity analysis on known de novo 
     * mutations suggests a default value of 10^-6.
     *
     */
    @Argument(fullName="deNovoPrior",shortName="DNP",doc="The de novo mutation prior",optional=true)
    public double deNovoPrior = 1e-6;

    /**
     * When a variant is not seen in a panel, whether to infer (and with what effective strength) that only reference
     * alleles were ascertained at that site. E.g. "If not seen in 1000Genomes, treat it as AC=0, AN=2000". This is
     * applied across all external panels, so if numRefIsMissing = 10, and the variant is absent in two panels, this
     * confers evidence of AC=0,AN=20
     */
    @Argument(fullName="numRefSamplesIfNoCall",shortName="nrs",doc="The number of homozygous reference to infer were " +
            "seen at a position where an \"other callset\" contains no site or genotype information",optional=true)
    public int numRefIfMissing = 0;

    /**
     * Rather than looking for the MLEAC field first, and then falling back to AC; first look for the AC field and then
     * fall back to MLEAC or raw genotypes
     */
    @Argument(fullName="defaultToAC",shortName="useAC",doc="Use the AC field as opposed to MLEAC. Does nothing if VCF lacks MLEAC field",optional=true)
    public boolean defaultToAC = false;

    /**
     * Do not use the [MLE] allele count from the input samples (the ones for which you're calculating posteriors)
     * in the site frequency distribution; only use the AC and AN calculated from external sources.
     */
    @Argument(fullName="ignoreInputSamples",shortName="ext",doc="Use external information only; do not inform genotype priors by " +
          "the discovered allele frequency in the callset whose posteriors are being calculated. Useful for callsets containing " +
          "related individuals.",optional=true)
    public boolean ignoreInputSamples = false;

    /**
     * Calculate priors for missing external variants from sample data -- default behavior is to apply flat priors
     */
    @Argument(fullName="discoveredACpriorsOff",shortName="useACoff",doc="Do not use discovered allele count in the input callset " +
            "for variants that do not appear in the external callset. ", optional=true)
    public boolean useACoff = false;

    /**
     * Skip application of population-based priors
     */
    @Argument(fullName="skipPopulationPriors",shortName="skipPop",doc="Skip application of population-based priors", optional=true)
    public boolean skipPopulationPriors = false;

    /**
     * Skip application of family-based priors.  Note: if pedigree file is absent, family-based priors will be skipped.
     */
    @Argument(fullName="skipFamilyPriors",shortName="skipFam",doc="Skip application of family-based priors", optional=true)
    public boolean skipFamilyPriors = false;

    /**
     * <p>Reads PED file-formatted tabular text files describing meta-data about the samples being
     * processed in the GATK.</p>
     *
     * <ul>
     *  <li>see <a href="http://www.broadinstitute.org/mpg/tagger/faq.html">http://www.broadinstitute.org/mpg/tagger/faq.html</a></li>
     *  <li>see <a href="http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped">http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped</a></li>
     * </ul>
     *
     * <p>The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:</p>
     *
     * <ul>
     *  <li>Family ID</li>
     *  <li>Individual ID</li>
     *  <li>Paternal ID</li>
     *  <li>Maternal ID</li>
     *  <li>Sex (1=male; 2=female; other=unknown)</li>
     *  <li>Phenotype</li>
     * </ul>
     *
     *  <p>The IDs are alphanumeric: the combination of family and individual ID should uniquely identify a person.
     *  A PED file must have 1 and only 1 phenotype in the sixth column. The phenotype can be either a
     *  quantitative trait or an affection status column: GATK will automatically detect which type
     *  (i.e. based on whether a value other than 0, 1, 2 or the missing genotype code is observed).</p>
     *
     *  <p>If an individual's sex is unknown, then any character other than 1 or 2 can be used.</p>
     *
     *  <p>You can add a comment to a PED or MAP file by starting the line with a # character. The rest of that
     *  line will be ignored. Do not start any family IDs with this character therefore.</p>
     *
     *  <p>Affection status should be coded:</p>
     *
     * <ul>
     *  <li>-9 missing</li>
     *   <li>0 missing</li>
     *   <li>1 unaffected</li>
     *   <li>2 affected</li>
     * </ul>
     *
     * <p>If any value outside of -9,0,1,2 is detected than the samples are assumed
     * to phenotype values are interpreted as string phenotype values.  In this case -9 uniquely
     * represents the missing value.</p>
     *
     * <p>Genotypes (column 7 onwards) cannot be specified to the GATK.</p>
     *
     * <p>For example, here are two individuals (one row = one person):</p>
     *
     * <pre>
     *   FAM001  1  0 0  1  2
     *   FAM001  2  0 0  1  2
     * </pre>
     *
     * <p>Each -ped argument can be tagged with NO_FAMILY_ID, NO_PARENTS, NO_SEX, NO_PHENOTYPE to
     * tell the GATK PED parser that the corresponding fields are missing from the ped file.</p>
     *
     * <p>Note that most GATK walkers do not use pedigree information.  Walkers that require pedigree
     * data should clearly indicate so in their arguments and will throw errors if required pedigree
     * information is missing.</p>
     */
    @Argument(fullName="pedigree", shortName="ped", doc="Pedigree file for samples", optional=true)
    private File pedigreeFile = null;

    private FamilyLikelihoods famUtils;
    private SampleDB sampleDB = null;

    private VariantContextWriter vcfWriter;

    @Override
    public void onTraversalStart() {
        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(out).setOutputFileType(VariantContextWriterBuilder.OutputType.VCF);
        if (hasReference()){
            vcfWriter = builder.setReferenceDictionary(getBestAvailableSequenceDictionary()).setOption(Options.INDEX_ON_THE_FLY).build();
        } else {
            vcfWriter = builder.unsetOption(Options.INDEX_ON_THE_FLY).build();
            logger.info("Can't make an index for output file " + out + " because a reference dictionary is required for creating Tribble indices on the fly");
        }

        sampleDB = initializeSampleDB();

        // Get list of samples to include in the output
        final Map<String, VCFHeader> vcfHeaders = Collections.singletonMap(getDrivingVariantsFeatureInput().getName(), getHeaderForVariants());
        final Set<String> vcfSamples = VcfUtils.getSortedSampleSet(vcfHeaders, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        //Get the trios from the families passed as ped
        if (!skipFamilyPriors){
            final Set<Trio> trios = sampleDB.getTrios();
            if(trios.isEmpty()) {
                logger.info("No PED file passed or no *non-skipped* trios found in PED file. Skipping family priors.");
                skipFamilyPriors = true;
            }
        }

        final VCFHeader header = vcfHeaders.values().iterator().next();
        if ( ! header.hasGenotypingData() ) {
            throw new UserException("VCF has no genotypes");
        }

        if ( header.hasInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY) ) {
            final VCFInfoHeaderLine mleLine = header.getInfoHeaderLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY);
            if ( mleLine.getCountType() != VCFHeaderLineCount.A ) {
                throw new UserException("VCF does not have a properly formatted MLEAC field: the count type should be \"A\"");
            }

            if ( mleLine.getType() != VCFHeaderLineType.Integer ) {
                throw new UserException("VCF does not have a properly formatted MLEAC field: the field type should be \"Integer\"");
            }
        }

        // Initialize VCF header
        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfHeaders.values(), true);
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.GENOTYPE_PRIOR_KEY));
        if (!skipFamilyPriors) {
            headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.JOINT_LIKELIHOOD_TAG_NAME));
            headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.JOINT_POSTERIOR_TAG_NAME));
        }
        headerLines.add(new VCFHeaderLine("source", "CalculateGenotypePosteriors"));

        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));

        final Map<String,Set<Sample>> families = sampleDB.getFamilies(vcfSamples);
        famUtils = new FamilyLikelihoods(sampleDB, deNovoPrior, vcfSamples, families);
    }

    /**
     * Entry-point function to initialize the samples database from input data
     */
    private SampleDB initializeSampleDB() {
        final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        if (pedigreeFile != null) {
            sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));
        }
        return sampleDBBuilder.getFinalSampleDB();
    }

    @Override
    public void apply(final VariantContext variant,
                      final ReadsContext readsContext,
                      final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        final Collection<VariantContext> vcs = featureContext.getValues(getDrivingVariantsFeatureInput());

        final Collection<VariantContext> otherVCs = featureContext.getValues(supportVariants);

        final int missing = supportVariants.size() - otherVCs.size();

        for ( final VariantContext vc : vcs ) {
            final VariantContext vc_familyPriors;
            final VariantContext vc_bothPriors;

            //do family priors first (if applicable)
            final VariantContextBuilder builder = new VariantContextBuilder(vc);
            //only compute family priors for biallelelic sites
            if (!skipFamilyPriors && vc.isBiallelic()){
                final GenotypesContext gc = famUtils.calculatePosteriorGLs(vc);
                builder.genotypes(gc);
            }
            VariantContextUtils.calculateChromosomeCounts(builder, false);
            vc_familyPriors = builder.make();

            if (!skipPopulationPriors) {
                vc_bothPriors = PosteriorProbabilitiesUtils.calculatePosteriorProbs(vc_familyPriors, otherVCs, missing * numRefIfMissing, globalPrior, !ignoreInputSamples, defaultToAC, useACoff);
            } else {
                final VariantContextBuilder builder2 = new VariantContextBuilder(vc_familyPriors);
                VariantContextUtils.calculateChromosomeCounts(builder, false);
                vc_bothPriors = builder2.make();
            }
            vcfWriter.add(vc_bothPriors);
        }
    }

    @Override
    public void closeTool(){
        vcfWriter.close();
    }
}

