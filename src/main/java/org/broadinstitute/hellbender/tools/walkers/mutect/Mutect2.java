package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReadOrientationArtifact;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReferenceBases;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.downsampling.MutectDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;

import java.io.File;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * <p>Call somatic short variants via local assembly of haplotypes.
 * Short variants include single nucleotide (SNV) and insertion and deletion (indel) variants.
 * The caller combines the DREAM challenge-winning somatic genotyping engine of the original MuTect
 * (<a href='http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html'>Cibulskis et al., 2013</a>) with the
 * assembly-based machinery of <a href="https://www.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php">HaplotypeCaller</a>.
 * </p>
 * <p>
 *     This tool is featured in the <i>Somatic Short Mutation calling Best Practice Workflow</i>.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 *     Although we present the tool for somatic calling, it may apply to other contexts, such as mitochondrial variant calling.
 * </p>
 *
 * <h3>Usage examples</h3>
 * <p>Example commands show how to run Mutect2 for typical scenarios. The two modes are (i) <i>somatic mode</i> where a tumor sample is matched with a normal sample in analysis and
 * (ii) <i>tumor-only mode</i> where a single sample's alignment data undergoes analysis. </p>
 *
 * <h4>(i) Tumor with matched normal</h4>
 * <p>
 *     Given a matched normal, Mutect2 is designed to call somatic variants only. The tool includes logic to skip
 *     emitting variants that are clearly present in the germline based on provided evidence, e.g. in the matched normal.
 *     This is done at an early stage to avoid spending computational resources on germline events. If the variant's germline status is
 *     borderline, then Mutect2 will emit the variant to the callset for subsequent filtering and review.
 * </p>
 *
 * <pre>
 * gatk Mutect2 \
 *   -R reference.fa \
 *   -I tumor.bam \
 *   -tumor tumor_sample_name \
 *   -I normal.bam \
 *   -normal normal_sample_name \
 *   --germline-resource af-only-gnomad.vcf.gz \
 *   --af-of-alleles-not-in-resource 0.00003125 \
 *   --panel-of-normals pon.vcf.gz \
 *   -O somatic.vcf.gz
 * </pre>
 *
 * <p>The --af-of-alleles-not-in-resource argument value should match expectations for alleles not found in the provided germline resource.
 * Note the tool does not require a germline resource nor a panel of normals (PoN) to run.
 * The tool prefilters sites for the matched normal and the PoN. For the germline resource, the tool prefilters on the allele.
 * Below is an excerpt of a known variants resource with population allele frequencies</p>
 * <pre>
 *     #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
 *      1       10067   .       T       TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC      30.35   PASS    AC=3;AF=7.384E-5
 *      1       10108   .       CAACCCT C       46514.32        PASS    AC=6;AF=1.525E-4
 *      1       10109   .       AACCCTAACCCT    AAACCCT,*       89837.27        PASS    AC=48,5;AF=0.001223,1.273E-4
 *      1       10114   .       TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTA  *,CAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTA,T      36728.97        PASS    AC=55,9,1;AF=0.001373,2.246E-4,2.496E-5
 *      1       10119   .       CT      C,*     251.23  PASS    AC=5,1;AF=1.249E-4,2.498E-5
 *      1       10120   .       TA      CA,*    14928.74        PASS    AC=10,6;AF=2.5E-4,1.5E-4
 *      1       10128   .       ACCCTAACCCTAACCCTAAC    A,*     285.71  PASS    AC=3,1;AF=7.58E-5,2.527E-5
 *      1       10131   .       CT      C,*     378.93  PASS    AC=7,5;AF=1.765E-4,1.261E-4
 *      1       10132   .       TAACCC  *,T     18025.11        PASS    AC=12,2;AF=3.03E-4,5.049E-5
 * </pre>
 *
 * <h4>(ii) Tumor-only mode</h4>
 * <p>This mode runs on a single sample, e.g. single tumor or single normal sample.
 * To create a PoN, call on each normal sample in this mode, then use {@link CreateSomaticPanelOfNormals} to generate the PoN. </p>
 *
 * <pre>
 *  gatk Mutect2 \
 *   -R reference.fa \
 *   -I sample.bam \
 *   -tumor sample_name \
 *   -O single_sample.vcf.gz
 * </pre>
 *
 *
 * <h3>Further points of interest</h3>
 * <p>
 *     Additional parameters that factor towards filtering, including normal-artifact-lod (default threshold 0.0) and
 *     tumor-lod (default threshold 5.3), are available in {@link FilterMutectCalls}. While the tool calculates
 *     normal-lod assuming a diploid genotype, it calculates
 *     normal-artifact-lod with the same approach it uses for tumor-lod, i.e. with a variable ploidy assumption.
 * </p>
 *
 * <dl>
 *     <dd>- If the normal artifact log odds becomes large, then FilterMutectCalls applies the artifact-in-normal filter.
 *     For matched normal samples with tumor contamination, consider increasing the normal-artifact-lod threshold.</dd>
 *     <dd>- The tumor log odds, which is calculated independently of any matched normal, determines whether to filter a tumor
 *     variant. Variants with tumor LODs exceeding the threshold pass filtering.</dd>
 * </dl>
 *
 * <p>
 *     If a variant is absent from a given germline resource, then the value for --af-of-alleles-not-in-resource applies.
 *     For example, gnomAD's 16,000 samples (~32,000 homologs per locus) becomes a probability of one in 32,000 or less.
 *     Thus, an allele's absence from the germline resource becomes evidence that it is not a germline variant.
 * </p>
 *
 * <h3>Caveats</h3>
 * <p>
 *     Although GATK4 Mutect2 accomodates varying coverage depths, further optimization of parameters
 *     may improve calling for extreme high depths, e.g. 1000X.
 * </p>
 */
 @CommandLineProgramProperties(
         summary = "Call somatic SNVs and indels via local assembly of haplotypes",
         oneLineSummary = "Call somatic SNVs and indels via local assembly of haplotypes",
         programGroup = ShortVariantDiscoveryProgramGroup.class
 )
@DocumentedFeature
public final class Mutect2 extends AssemblyRegionWalker {

    @ArgumentCollection
    protected M2ArgumentCollection MTAC = new M2ArgumentCollection();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to which variants should be written")
    public File outputVCF;

    private VariantContextWriter vcfWriter;

    private Mutect2Engine m2Engine;

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
    protected boolean includeReadsWithDeletionsInIsActivePileups() { return true; }

    @Override
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Mutect2Engine.makeStandardMutect2ReadFilters();
    }

    @Override
    public ReadTransformer makePostReadFilterTransformer() {
        return super.makePostReadFilterTransformer().andThen(Mutect2Engine.makeStandardMutect2PostFilterReadTransformer(referenceArguments.getReferencePath(), !MTAC.dontClipITRArtifacts));
    }

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Mutect2Engine.getStandardMutect2AnnotationGroups();
    }

    @Override
    protected ReadsDownsampler createDownsampler() {
        return new MutectDownsampler(maxReadsPerAlignmentStart, MTAC.maxSuspiciousReadsPerAlignmentStart, MTAC.downsamplingStride);
    }

    @Override
    public AssemblyRegionEvaluator assemblyRegionEvaluator() { return m2Engine; }

    @Override
    public void onTraversalStart() {
        VariantAnnotatorEngine annotatorEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), null, Collections.emptyList(), false);
        m2Engine = new Mutect2Engine(MTAC, createOutputBamIndex, createOutputBamMD5, getHeaderForReads(), referenceArguments.getReferenceFileName(), annotatorEngine);
        vcfWriter = createVCFWriter(outputVCF);
        m2Engine.writeHeader(vcfWriter, getDefaultToolVCFHeaderLines());
    }

    @Override
    public Collection<Annotation> makeVariantAnnotations(){
        final Collection<Annotation> annotations = super.makeVariantAnnotations();

        if (MTAC.artifactPriorTable != null){
            // Enable the annotations associated with the read orientation model
            annotations.add(new ReadOrientationArtifact(MTAC.artifactPriorTable));
            annotations.add(new ReferenceBases());
        }
        return annotations;
    }

    @Override
    public Object onTraversalSuccess() {
        return "SUCCESS";
    }

    @Override
    public void apply(final AssemblyRegion region, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        m2Engine.callRegion(region, referenceContext, featureContext).forEach(vcfWriter::add);
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
        if (m2Engine != null) {
            m2Engine.shutdown();
        }
    }
}
