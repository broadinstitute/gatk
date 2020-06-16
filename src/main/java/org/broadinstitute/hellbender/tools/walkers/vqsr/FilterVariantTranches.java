package org.broadinstitute.hellbender.tools.walkers.vqsr;

import java.util.*;
import java.util.stream.Collectors;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

/**
 * Apply tranche filtering to VCF based on scores from an annotation in the INFO field.
 * The annotation can come from the {@link CNNScoreVariants} tool (CNNLOD), VQSR (VQSLOD),
 * or any other variant scoring tool which adds numeric annotations in a VCF's INFO field.
 *
 * Tranches are specified in percent sensitivity to the variants in the resource files.
 * For example, if you specify INDEL tranches 98.0 and 99.0 using the CNN_2D score
 * the filtered VCF will contain 2 filter tranches for INDELS: CNN_2D_INDEL_Tranche_98.00_99.00
 * and CNN_2D_INDEL_Tranche_99.00_100.00. Variants that scored better than the 98th percentile of variants in the
 * resources pass through the filter and will have `PASS` in the filter field. We expect variants in the tranche
 * CNN_2D_INDEL_Tranche_99.00_100.00 to be more sensitive, but less precise than CNN_2D_INDEL_Tranche_98.00_99.00,
 * because variants in CNN_2D_INDEL_Tranche_99.00_100.00 have lower scores than variants in the tranche
 * CNN_2D_INDEL_Tranche_98.00_99.00.
 *
 * The default tranche filtering threshold for SNPs is 99.95 and for INDELs it is 99.4.
 * These thresholds maximize the F1 score (the harmonic mean of sensitivity and precision)
 * for whole genome human data but may need to be tweaked for different datasets.
 *
 *
 * <h3>Inputs</h3>
 * <ul>
 *      <li>The input variants to tranche filter.</li>
 *      <li>resource A VCF containing known SNP and or INDEL sites. Can be supplied as many times as necessary </li>
 *      <li>info-key The key from the INFO field of the VCF which contains the values that will be used to filter.</li>
 *      <li>tranche List of percent sensitivities to the known sites at which we will filter.  Must be between 0 and 100.</li>
 * </ul>
 *
 * If you want to remove existing filters from your VCF add the argument `--invalidate-previous-filters`.
 *
 *
 * <h3>Outputs</h3>
 * <ul>
 * <li>A tranche filtered VCF.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 *
 * <h4>Apply tranche filters based on the scores in the info field with key CNN_1D.</h4>
 * <pre>
 * gatk FilterVariantTranches \
 *   -V input.vcf.gz \
 *   --resource hapmap.vcf \
 *   --resource mills.vcf \
 *   --info-key CNN_1D \
 *   --snp-tranche 99.95 \
 *   --indel-tranche 99.4 \
 *   -O filtered.vcf
 * </pre>
 *
 * <h4>Apply tranche filters based on the scores in the info field with key CNN_2D
 * and remove any existing filters from the VCF.</h4>
 * <pre>
 * gatk FilterVariantTranches \
 *   -V input.vcf.gz \
 *   --resource hapmap.vcf \
 *   --resource mills.vcf \
 *   --info-key CNN_2D \
 *   --snp-tranche 99.95 \
 *   --indel-tranche 99.4 \
 *   --invalidate-previous-filters \
 *   -O filtered.vcf
 * </pre>
 *
 * <h4>Apply several tranche filters based on the scores in the info field with key CNN_2D.</h4>
 * This will result in a VCF with filters: CNN_2D_SNP_Tranche_99.90_99.95, CNN_2D_SNP_Tranche_99.95_100.00
 * CNN_2D_INDEL_Tranche_99.00_99.50, CNN_2D_INDEL_Tranche_99.50_100.00.
 * The interpretation is that `PASS` variants are the best, variants in the tranche `CNN_2D_INDEL_Tranche_99.00_99.50`
 * are ok and variants in the tranche `CNN_2D_INDEL_Tranche_99.50_100.00` are the worst.
 * <pre>
 * gatk FilterVariantTranches \
 *   -V input.vcf.gz \
 *   --resource hapmap.vcf \
 *   --resource mills.vcf \
 *   --info-key CNN_2D \
 *   --snp-tranche 99.9 --snp-tranche 99.95 \
 *   --indel-tranche 99.0 --indel-tranche 99.4 \
 *   -O filtered.vcf
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Apply tranche filtering based on a truth VCF of known common sites of variation and a score from VCF INFO field",
        oneLineSummary = "Apply tranche filtering",
        programGroup = VariantFilteringProgramGroup.class
)

public class FilterVariantTranches extends TwoPassVariantWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF file")
    private GATKPath outputVcf = null;

    @Argument(fullName="snp-tranche",
            shortName="snp-tranche",
            doc="The level(s) of sensitivity to SNPs in the resource VCFs at which to filter SNPs. " +
                    "Higher numbers mean more desired sensitivity and thus less stringent filtering." +
                    "Specified in percents, i.e. 99.9 for 99.9 percent and 1.0 for 1 percent.",
            optional=true)
    private List<Double> snpTranches = new ArrayList<>(Arrays.asList(99.95));

    @Argument(fullName="indel-tranche",
            shortName="indel-tranche",
            doc="The level(s) of sensitivity to indels in the resource VCFs at which to filter indels. " +
                    "Higher numbers mean more desired sensitivity and thus less stringent filtering." +
                    "Specified in percents, i.e. 99.9 for 99.9 percent and 1.0 for 1 percent.",
            optional=true)
    private List<Double> indelTranches = new ArrayList<>(Arrays.asList(99.4));

    @Argument(fullName=StandardArgumentDefinitions.RESOURCE_LONG_NAME,
            doc="A list of validated VCFs with known sites of common variation",
            optional=false)
    private List<FeatureInput<VariantContext>> resources = new ArrayList<>();

    @Argument(fullName = "info-key", shortName = "info-key", doc = "The key must be in the INFO field of the input VCF.")
    private String infoKey = GATKVCFConstants.CNN_2D_KEY;

    @Argument(fullName = StandardArgumentDefinitions.INVALIDATE_PREVIOUS_FILTERS_LONG_NAME,
            doc = "Remove all filters that already exist in the VCF.",
            optional=true)
    private boolean removeOldFilters = false;

    private VariantContextWriter vcfWriter;
    private List<Double> resourceSNPScores = new ArrayList<>();
    private List<Double> snpCutoffs = new ArrayList<>();
    private List<Double> resourceIndelScores = new ArrayList<>();
    private List<Double> indelCutoffs = new ArrayList<>();

    private int scoredSnps = 0;
    private int filteredSnps = 0;
    private int scoredIndels = 0;
    private int filteredIndels = 0;

    private static String SNPString = "SNP";
    private static String INDELString = "INDEL";

    @Override
    public void onTraversalStart() {
        snpTranches = validateTranches(snpTranches);
        indelTranches = validateTranches(indelTranches);
        vcfWriter = createVCFWriter(outputVcf);
        writeVCFHeader(vcfWriter);
    }

    @Override
    public void firstPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (!variant.hasAttribute(infoKey)){
            return;
        } else if (variant.isSNP()){
            scoredSnps++;
        } else if (variant.isIndel()){
            scoredIndels++;
        }

        for (FeatureInput<VariantContext> featureSource : resources) {
            for (VariantContext v : featureContext.getValues(featureSource)) {
                for (final Allele a : variant.getAlternateAlleles()) {
                    try {
                        if ((variant.getStart() == v.getStart()) && GATKVariantContextUtils.isAlleleInList(variant.getReference(), a, v.getReference(), v.getAlternateAlleles())) {
                            if (variant.isSNP()) {
                                resourceSNPScores.add(Double.parseDouble((String) variant.getAttribute(infoKey)));
                                return;
                            } else {
                                resourceIndelScores.add(Double.parseDouble((String) variant.getAttribute(infoKey)));
                                return;
                            }
                        }
                    } catch (IllegalStateException e) {
                        throw new UserException.BadInput(String.format("The provided variant file(s) have inconsistent references " +
                                "for the same position(s) at %s:%d, %s in input vs. %s in resource", v.getContig(), v.getStart(), variant.getReference(), v.getReference()));
                    }
                }
            }
        }
    }

    @Override
    public void afterFirstPass() {
        logger.info(String.format("Found %d SNPs and %d indels with INFO score key:%s.", scoredSnps, scoredIndels, infoKey));
        logger.info(String.format("Found %d SNPs and %d indels in the resources.", resourceSNPScores.size(), resourceIndelScores.size()));

        if (scoredSnps == 0 && scoredIndels == 0) {
            throw new UserException.BadInput("VCF contains no variants or no variants with INFO score key \"" + infoKey + "\"");
        }

        if (resourceSNPScores.size() == 0 && resourceIndelScores.size() == 0) {
            throw new UserException.BadInput("Neither SNP nor indel resource contains variants overlapping input.  Filtering cannot be performed.");
        }

        if ((scoredSnps > 0 && resourceSNPScores.size() == 0)) {
            throw new UserException.BadInput("SNPs are present in input VCF, but cannot be filtered because no overlapping SNPs were found in the resources.");
        }

        if ((scoredIndels > 0 && resourceIndelScores.size() == 0)) {
            throw new UserException.BadInput("Indels are present in input VCF, but cannot be filtered because no overlapping indels were found in the resources.");
        }

        Collections.sort(resourceSNPScores, Collections.reverseOrder());
        Collections.sort(resourceIndelScores, Collections.reverseOrder());

        if (resourceSNPScores.size() != 0) {
            for (double t : snpTranches) {
                int snpIndex = (int) ((t / 100.0) * (double) (resourceSNPScores.size() - 1));
                snpCutoffs.add(resourceSNPScores.get(snpIndex));
            }
        }

        if (resourceIndelScores.size() != 0) {
            for (double t : indelTranches) {
                int indelIndex = (int) ((t / 100.0) * (double) (resourceIndelScores.size() - 1));
                indelCutoffs.add(resourceIndelScores.get(indelIndex));
            }
        }

    }

    @Override
    protected void secondPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final VariantContextBuilder builder = new VariantContextBuilder(variant);

        if (removeOldFilters) {
            builder.unfiltered();
        }

        if (variant.hasAttribute(infoKey)) {
            final double score = Double.parseDouble((String) variant.getAttribute(infoKey));
            if (variant.isSNP() && snpCutoffs.size() != 0 && isTrancheFiltered(score, snpCutoffs)) {
                builder.filter(filterStringFromScore(SNPString, score, snpTranches, snpCutoffs));
                filteredSnps++;
            } else if (variant.isIndel() && indelCutoffs.size() != 0 && isTrancheFiltered(score, indelCutoffs)) {
                builder.filter(filterStringFromScore(INDELString, score, indelTranches, indelCutoffs));
                filteredIndels++;
            }
        }

        if (builder.getFilters() == null || builder.getFilters().size() == 0){
            builder.passFilters();
        }
        
        vcfWriter.add(builder.make());
    }

    @Override
    public void closeTool() {
        logger.info(String.format("Filtered %d SNPs out of %d and filtered %d indels out of %d with INFO score: %s.",
                filteredSnps, scoredSnps, filteredIndels, scoredIndels, infoKey));

        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

    private void writeVCFHeader(VariantContextWriter vcfWriter) {
        // setup the header fields
        final VCFHeader inputHeader = getHeaderForVariants();
        Set<VCFHeaderLine> hInfo = new LinkedHashSet<VCFHeaderLine>();
        hInfo.addAll(inputHeader.getMetaDataInSortedOrder());

        boolean hasInfoKey = hInfo.stream().anyMatch(
                x -> x instanceof VCFInfoHeaderLine && ((VCFInfoHeaderLine) x).getID().equals(infoKey));
        if (!hasInfoKey){
            throw new UserException(String.format("Input VCF does not contain a header line for specified info key:%s", infoKey));
        }

        if (removeOldFilters){
            hInfo.removeIf(x -> x instanceof VCFFilterHeaderLine);
        }

        addTrancheHeaderFields(SNPString, snpTranches, hInfo);
        addTrancheHeaderFields(INDELString, indelTranches, hInfo);

        final TreeSet<String> samples = new TreeSet<>();
        samples.addAll(inputHeader.getGenotypeSamples());
        hInfo.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    private void addTrancheHeaderFields(String type, List<Double> tranches, Set<VCFHeaderLine> hInfo){
        if( tranches.size() >= 2 ) {
            for(int i = 0; i < tranches.size() - 1; i++) {
                String filterKey = filterKeyFromTranches(type, infoKey, tranches.get(i), tranches.get(i+1));
                String filterDescription = filterDescriptionFromTranches(type, infoKey, tranches.get(i), tranches.get(i+1));
                hInfo.add(new VCFFilterHeaderLine(filterKey, filterDescription));
            }
        }
        String filterKey = filterKeyFromTranches(type, infoKey, tranches.get(tranches.size()-1), 100.0);
        String filterDescription = filterDescriptionFromTranches(type, infoKey, tranches.get(tranches.size()-1), 100.0);
        hInfo.add(new VCFFilterHeaderLine(filterKey, filterDescription));
    }

    private String filterKeyFromTranches(String type, String infoKey, double t1, double t2){
        return String.format("%s_%s_Tranche_%.2f_%.2f", infoKey, type, t1, t2);
    }
    
    private String filterDescriptionFromTranches(String type, String infoKey, double t1, double t2){
        return String.format("%s truth resource sensitivity between %.2f and %.2f for info key %s", type, t1, t2, infoKey);
    }

    private boolean isTrancheFiltered(double score, List<Double> cutoffs) {
        return score <= cutoffs.get(0);
    }

    private String filterStringFromScore(String type, double score, List<Double> tranches, List<Double> cutoffs){
        for (int i = 0; i < cutoffs.size(); i++){
            if (score > cutoffs.get(i) && i == 0){
                throw new GATKException("Trying to add a filter to a passing variant.");
            } else if (score > cutoffs.get(i)){
                return filterKeyFromTranches(type, infoKey, tranches.get(i-1), tranches.get(i));
            }
        }
        return filterKeyFromTranches(type, infoKey, tranches.get(tranches.size()-1), 100.0);
    }

    private List<Double> validateTranches(List<Double> tranches){
        if (tranches.size() < 1 || tranches.stream().anyMatch(d -> d < 0 || d >= 100.0)){
            throw new CommandLineException("At least 1 tranche value must be given and all tranches must be greater than 0 and less than 100.");
        }
        List<Double> newTranches = tranches.stream().distinct().collect(Collectors.toList());
        newTranches.sort(Double::compareTo);
        return newTranches;
    }

}
