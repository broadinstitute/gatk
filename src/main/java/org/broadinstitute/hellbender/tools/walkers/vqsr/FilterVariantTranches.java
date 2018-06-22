package org.broadinstitute.hellbender.tools.walkers.vqsr;

import java.util.*;
import java.io.File;
import java.util.stream.Collectors;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import picard.cmdline.programgroups.VariantFilteringProgramGroup;

/**
 * Apply tranche filtering to VCF based on scores from an annotation in the INFO field.
 *
 * <h3>Inputs</h3>
 * <ul>
 *      <li>The input variants to tranche filter.</li>
 *      <li>resource A VCF containing known SNP and or INDEL sites. Can be supplied as many times as necessary </li>
 *      <li>info-key The key from the INFO field of the VCF which contains the values that will be used to filter.</li>
 *      <li>tranche List of percent sensitivities to the known sites at which we will filter.  Must be between 0 and 100.</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 * <li>A tranche filtered VCF.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 *
 * <h4>Apply tranche filters based on CNN_1D scores</h4>
 * <pre>
 * gatk FilterVariantTranches \
 *   -V input.vcf.gz \
 *   --resource hapmap.vcf \
 *   --resource mills.vcf \
 *   --info-key CNN_1D \
 *   --tranche 99.9 --tranche 99.0 --tranche 95 \
 *   -O filtered.vcf
 * </pre>
 *
 */
@DocumentedFeature
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Apply tranche filtering based on a truth VCF of known common sites of variation and a score from VCF INFO field",
        oneLineSummary = "Apply tranche filtering",
        programGroup = VariantFilteringProgramGroup.class
)

public class FilterVariantTranches extends TwoPassVariantWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF file")
    private String outputVcf = null;

    @Argument(fullName="tranche",
            shortName="t",
            doc="The levels of truth sensitivity at which to slice the data. (in percents, i.e. 99.9 for 99.9 percent and 1.0 for 1 percent)",
            optional=true)
    private List<Double> tranches = new ArrayList<>(Arrays.asList(99.9, 99.99));

    @Argument(fullName="resource",
            shortName = "resource",
            doc="A list of validated VCFs with known sites of common variation",
            optional=false)
    private List<FeatureInput<VariantContext>> resources = new ArrayList<>();

    @Argument(fullName = "info-key", shortName = "info-key", doc = "The key must be in the INFO field of the input VCF.")
    private String infoKey = GATKVCFConstants.CNN_2D_KEY;

    private VariantContextWriter vcfWriter;
    private List<Double> snpScores = new ArrayList<>();
    private List<Double> snpCutoffs = new ArrayList<>();
    private List<Double> indelScores = new ArrayList<>();
    private List<Double> indelCutoffs = new ArrayList<>();

    private int scoredSnps = 0;
    private int filteredSnps = 0;
    private int scoredIndels = 0;
    private int filteredIndels = 0;

    @Override
    public void onTraversalStart() {
        if (tranches.size() < 1 || tranches.stream().anyMatch(d -> d < 0 || d >= 100.0)){
            throw new GATKException("At least 1 tranche value must be given and all tranches must be greater than 0 and less than 100.");
        }
        tranches = tranches.stream().distinct().collect(Collectors.toList());
        tranches.sort(Double::compareTo);
        vcfWriter = createVCFWriter(new File(outputVcf));
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
                if (variant.isSNP()){
                    if(variant.getAlternateAlleles().stream().anyMatch(v.getAlternateAlleles()::contains)) {
                        snpScores.add(Double.parseDouble((String) variant.getAttribute(infoKey)));
                        return;
                    }
                } else if (variant.isIndel()){
                    if(variant.getAlternateAlleles().stream().anyMatch(v.getAlternateAlleles()::contains)){
                        indelScores.add(Double.parseDouble((String)variant.getAttribute(infoKey)));
                        return;
                    }
                }
            }
        }
    }

    @Override
    public void afterFirstPass() {
        logger.info(String.format("Found %d SNPs %d indels with INFO score key:%s.", scoredSnps, scoredIndels, infoKey));
        logger.info(String.format("Found %d SNPs %d indels in the resources.", snpScores.size(), indelScores.size()));

        if (scoredSnps == 0 || scoredIndels == 0 || snpScores.size() == 0 || indelScores.size() == 0){
            throw new GATKException("VCF must contain SNPs and indels with scores and resources must contain matching SNPs and indels.");
        }

        Collections.sort(snpScores, Collections.reverseOrder());
        Collections.sort(indelScores, Collections.reverseOrder());

        for(double t : tranches) {
            int snpIndex = (int)((t/100.0)*(double)(snpScores.size()-1));
            snpCutoffs.add(snpScores.get(snpIndex));
            int indelIndex = (int)((t/100.0)*(double)(indelScores.size()-1));
            indelCutoffs.add(indelScores.get(indelIndex));
        }

    }

    @Override
    protected void secondPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final VariantContextBuilder builder = new VariantContextBuilder(variant);

        if (variant.hasAttribute(infoKey)) {
            final double score = Double.parseDouble((String) variant.getAttribute(infoKey));
            if (variant.isSNP() && isTrancheFiltered(score, snpCutoffs)) {
                builder.filter(filterStringFromScore(score, snpCutoffs));
                filteredSnps++;
            } else if (variant.isIndel() && isTrancheFiltered(score, indelCutoffs)) {
                builder.filter(filterStringFromScore(score, indelCutoffs));
                filteredIndels++;
            }
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
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();
        final Set<VCFHeaderLine> hInfo = new HashSet<>(inputHeaders);

        if( tranches.size() >= 2 ) {
            for(int i = 0; i < tranches.size() - 1; i++) {
                String filterKey = filterKeyFromTranches(infoKey, tranches.get(i), tranches.get(i+1));
                String filterDescription = filterDescriptionFromTranches(infoKey, tranches.get(i), tranches.get(i+1));
                hInfo.add(new VCFFilterHeaderLine(filterKey, filterDescription));
            }
        }
        String filterKey = filterKeyFromTranches(infoKey, tranches.get(tranches.size()-1), 100.0);
        String filterDescription = filterDescriptionFromTranches(infoKey, tranches.get(tranches.size()-1), 100.0);
        hInfo.add(new VCFFilterHeaderLine(filterKey, filterDescription));
        final TreeSet<String> samples = new TreeSet<>();
        samples.addAll(inputHeader.getGenotypeSamples());
        hInfo.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    private String filterKeyFromTranches(String infoKey, double t1, double t2){
        return String.format("%s_Tranche_%.2f_%.2f", infoKey, t1, t2);
    }
    
    private String filterDescriptionFromTranches(String infoKey, double t1, double t2){
        return String.format("Truth sensitivity between %.2f and %.2f for info key %s", t1, t2, infoKey);
    }

    private boolean isTrancheFiltered(double score, List<Double> cutoffs) {
        return score <= cutoffs.get(0);
    }

    private String filterStringFromScore(double score, List<Double> cutoffs){
        for (int i = 0; i < cutoffs.size(); i++){
            if (score > cutoffs.get(i) && i == 0){
                throw new GATKException("Trying to add a filter to a passing variant.");
            } else if (score > cutoffs.get(i)){
                return filterKeyFromTranches(infoKey, tranches.get(i-1), tranches.get(i));
            }
        }
        return filterKeyFromTranches(infoKey, tranches.get(tranches.size()-1), 100.0);
    }

}
