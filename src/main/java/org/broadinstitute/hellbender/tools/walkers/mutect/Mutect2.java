package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;

/**
 * Call somatic SNPs and indels via local re-assembly of haplotypes
 *
 * <p>Mutect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original Mutect (<a href='http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html'>Cibulskis et al., 2013</a>) with the assembly-based machinery of HaplotypeCaller.</p>
 *
 * <p>The basic operation of Mutect2 proceeds similarly to that of the <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php">HaplotypeCaller</a>   </p>
 *
 * <h3>Differences from HaplotypeCaller</h3>
 * <p>While the HaplotypeCaller relies on a ploidy assumption (diploid by default) to inform its genotype likelihood and
 * variant quality calculations, Mutect2 allows for a varying allelic fraction for each variant, as is often seen in tumors with purity less
 * than 100%, multiple subclones, and/or copy number variation (either local or aneuploidy). Mutect2 also differs from the HaplotypeCaller in that it does apply some hard filters
 * to variants before producing output.</p>
 *
 * <h3>Usage examples</h3>
 * <p>These are example commands that show how to run Mutect2 for typical use cases. Square brackets ("[ ]")
 * indicate optional arguments. Note that parameter values shown here may not be the latest recommended; see the
 * Best Practices documentation for detailed recommendations. </p>
 *
 * <br />
 * <h4>Tumor/Normal variant calling</h4>
 * <pre>
 *   java
 *     -jar GenomeAnalysisTK.jar \
 *     -T Mutect2 \
 *     -R reference.fasta \
 *     -I tumor.bam \
 *     -I normal.bam \
 *     -tumor tumorSampleName \ // as in the BAM header
 *     -normal normalSampleName \ // as in the BAM header
 *     [--dbsnp dbSNP.vcf] \
 *     [--cosmic COSMIC.vcf] \
 *     [-L targets.interval_list] \
 *     -o output.vcf
 * </pre>
 *
 * <h4>Normal-only calling for panel of normals creation</h4>
 * <pre>
 *   java
 *     -jar GenomeAnalysisTK.jar
 *     -T Mutect2
 *     -R reference.fasta
 *     -I:tumor normal1.bam \
 *     [--dbsnp dbSNP.vcf] \
 *     [--cosmic COSMIC.vcf] \
 *     --artifact_detection_mode \
 *     [-L targets.interval_list] \
 *     -o output.normal1.vcf
 * </pre>
 * <br />
 * For full PON creation, call each of your normals separately in artifact detection mode. Then use CombineVariants to
 * output only sites where a variant was seen in at least two samples:
 * <pre>
 * java -jar GenomeAnalysisTK.jar
 *     -T CombineVariants
 *     -R reference.fasta
 *     -V output.normal1.vcf -V output.normal2.vcf [-V output.normal2.vcf ...] \
 *     -minN 2 \
 *     --setKey "null" \
 *     --filteredAreUncalled \
 *     --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
 *     [-L targets.interval_list] \
 *     -o Mutect2_PON.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>Mutect2 currently only supports the calling of a single tumor-normal pair at a time</li>
 * </ul>
 *
 */
@CommandLineProgramProperties(
        summary = "Call somatic SNPs and indels via local re-assembly of haplotypes",
        oneLineSummary = "Call somatic SNPs and indels via local re-assembly of haplotypes",
        programGroup = VariantProgramGroup.class
)
public final class Mutect2 extends AssemblyRegionWalker {

    @ArgumentCollection
    protected M2ArgumentCollection MTAC = new M2ArgumentCollection();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to which variants should be written")
    public String outputVCF = null;

    private VariantContextWriter vcfWriter;

    private Mutect2Engine m2Engine;

    @Override
    protected int defaultReadShardSize() { return 5000; }

    @Override
    protected int defaultReadShardPadding() { return 100; }

    @Override
    protected int defaultMinAssemblyRegionSize() { return 50; }

    @Override
    protected int defaultMaxAssemblyRegionSize() { return 300; }

    @Override
    protected int defaultAssemblyRegionPadding() { return 100; }

    @Override
    protected int defaultMaxReadsPerAlignmentStart() { return 50; }

    @Override
    protected double defaultActiveProbThreshold() { return 0.002; }

    @Override
    protected int defaultMaxProbPropagationDistance() { return 50; }

    @Override
    public void onTraversalStart() {
        m2Engine = new Mutect2Engine(MTAC, getHeaderForReads(), referenceArguments.getReferenceFileName());
        vcfWriter = GATKVariantContextUtils.createVCFWriter(new File(outputVCF), getHeaderForReads().getSequenceDictionary(), false);
        m2Engine.writeHeader(vcfWriter);
    }

    @Override
    public AssemblyRegionEvaluator assemblyRegionEvaluator() { return m2Engine; }

    @Override
    public void apply(final AssemblyRegion region, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        m2Engine.callRegion(region, featureContext).stream()
                // Only include calls that start within the current read shard (as opposed to the padded regions around it).
                // This is critical to avoid duplicating events that span shard boundaries!
                .filter(call -> getCurrentReadShardBounds().contains(call))
                .forEach(vcfWriter::add);
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }

        if ( m2Engine != null ) {
            m2Engine.shutdown();
        }
    }
}
