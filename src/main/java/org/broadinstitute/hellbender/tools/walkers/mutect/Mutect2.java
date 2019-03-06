package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.downsampling.MutectDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.variant.writers.SomaticGVCFWriter;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * <p>Call somatic short mutations via local assembly of haplotypes.
 * Short mutations include single nucleotide (SNA) and insertion and deletion (indel) alterations.
 * The caller uses a Bayesian somatic genotyping model that differs from the original MuTect by
 * <a href='http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html'>Cibulskis et al., 2013</a>
 * and uses the assembly-based machinery of
 * <a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php">HaplotypeCaller</a>.
 * Of note, Mutect2 v4.1.0.0 onwards enables joint analysis of multiple samples.</p>
 *
 * <p>This tool is featured in the <i>Somatic Short Mutation calling Best Practice Workflow</i>.
 * See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 * step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 * for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 * <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 * For pipelines with example data, see the <a href="https://github.com/gatk-workflows/gatk4-somatic-snvs-indels">gatk-workflows repository</a>.
 * Although we present the tool for somatic calling, it may apply to other contexts, such as mitochondrial variant calling and detection of somatic mosaicism.</p>
 *
 * <p>Starting with v4.1.0.0 Mutect2 accomodates extreme high depths, e.g. 20,000X. See the following articles for details on this and additional applications.
 * <ul>
 *   <li><a href="https://software.broadinstitute.org/gatk/blog?id=23400">Blog#23400</a> details general improvements to
 *   Mutect2 v4.1.0.0.</li>
 *   <li><a href="https://software.broadinstitute.org/gatk/blog?id=23598">Blog#23598</a> details Mutect2 mitochondrial
 *   mode.</li>
 *   <li><a href="https://software.broadinstitute.org/gatk/blog?id=XXX">Blog#XXX</a> (link to come) details use of Mutect2 in extremely
 *   low allele fraction variant detection.</li>
 * </ul></p>
 *
 * <h3>Usage examples</h3>
 * <p>Example commands show how to run Mutect2 for typical scenarios. The three modes are (i) tumor-normal mode where a
 * tumor sample is matched with a normal sample in analysis, (ii) tumor-only mode where a single sample's alignment
 * data undergoes analysis, and (iii) mitochondrial mode where sensitive calling at high depths is desirable.
 * <ul>
 *     <li>As of v4.1, there is no longer a need to specify the tumor sample name with -tumor. You need only specify
 *     the normal sample name with <nobr>-normal,</nobr> if you include a normal.</li>
 *     <li>Starting with v4.0.4.0, GATK recommends the default setting of --af-of-alleles-not-in-resource, which the
 *     tool dynamically adjusts for different modes. tumor-only calling sets the default to 5e-8, tumor-normal calling
 *     sets it to 1e-6 and mitochondrial mode sets it to 4e-3. For previous versions, the default was 0.001, the
 *     average heterozygosity of humans. For other organisms, change --af-of-alleles-not-in-resource to
 *     1/(ploidy*samples in resource).</li>
 * </ul></p>
 *
 * <h4>(i) Tumor with matched normal</h4>
 * <p>Given a matched normal, Mutect2 is designed to call somatic variants only. The tool includes logic to skip
 * emitting variants that are clearly present in the germline based on provided evidence, e.g. in the matched normal.
 * This is done at an early stage to avoid spending computational resources on germline events. If the variant's
 * germline status is borderline, then Mutect2 will emit the variant to the callset for subsequent filtering by
 * <a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_FilterMutectCalls.php">FilterMutectCalls</a>
 * and review.</p>
 *
 * <pre>
 *     gatk Mutect2 \
 *     -R reference.fa \
 *     -I tumor.bam \
 *     -I normal.bam \
 *     -normal normal_sample_name \
 *     --germline-resource af-only-gnomad.vcf.gz \
 *     --panel-of-normals pon.vcf.gz \
 *     -O somatic.vcf.gz
 * </pre>
 *
 * <p> As of v4.1 Mutect2 supports joint calling of multiple tumor and normal samples from the same individual. The
 * only difference is that -I and <nobr>-normal</nobr> must be specified for the extra samples.</p>
 *
 * <pre>
 *     gatk Mutect2 \
 *     -R reference.fa \
 *     -I tumor1.bam \
 *     -I tumor2.bam \
 *     -I normal1.bam \
 *     -I normal2.bam \
 *     -normal normal1_sample_name \
 *     -normal normal2_sample_name \
 *     --germline-resource af-only-gnomad.vcf.gz \
 *     --panel-of-normals pon.vcf.gz \
 *     -O somatic.vcf.gz
 * </pre>
 *
 * <h4>(ii) Tumor-only mode</h4>
 * <p>This mode runs on a single type of sample, e.g. the tumor or the normal.
 * To create a PoN, call on each normal sample in this mode, then use
 * <a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_CreateSomaticPanelOfNormals.php">CreateSomaticPanelOfNormals</a>
 * to generate the PoN.</p>
 *
 * <pre>
 *  gatk Mutect2 \
 *   -R reference.fa \
 *   -I sample.bam \
 *   -O single_sample.vcf.gz
 * </pre>
 *
 * <p>To call mutations on a tumor sample, call in this mode using a PoN and germline resource. After FilterMutectCalls
 * filtering, consider additional filtering by functional significance with
 * <a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_funcotator_Funcotator.php">Funcotator</a>.</p>
 *
 * <pre>
 *  gatk Mutect2 \
 *  -R reference.fa \
 *  -I sample.bam \
 *  --germline-resource af-only-gnomad.vcf.gz \
 *  --panel-of-normals pon.vcf.gz \
 *  -O single_sample.vcf.gz
 * </pre>
 *
 * <h4>(iii) Mitochondrial mode</h4>
 * <p>Mutect2 automatically sets parameters appropriately for calling on mitochondria with the <nobr>--mitochondria</nobr> flag.
 * Specifically, the mode sets <nobr>–-initial-tumor-lod</nobr> to 0, <nobr>–-tumor-lod-to-emit</nobr> to 0, <nobr>--af-of-alleles-not-in-resource</nobr> to
 * 4e-3, and the advanced parameter <nobr>--pruning-lod-threshold</nobr> to -4.</p>
 *
 * <pre>
 *  gatk Mutect2 \
 *  -R reference.fa \
 *  -L chrM \
 *  --mitochondria \
 *  --median-autosomal-coverage 30 \
 *  -I mitochondria.bam \
 *  -O mitochondria.vcf.gz
 * </pre>
 *
 * <p>Setting the advanced option --median-autosomal-coverage argument (default 0) activates a recommended filter against
 * likely erroneously mapped  <a href="https://en.wikipedia.org/wiki/NUMT">NuMTs (nuclear mitochondrial DNA segments)</a>.
 * For the value, provide the median coverage expected in autosomal regions with coverage.  The mode accepts only a
 * single sample, which can be provided in multiple files.</p>
 *
 * <h3>Notes</h3>
 * <ol>
 *     <li>Mutect2 does not require a germline resource nor a panel of normals (PoN) to run. The tool prefilters sites
 * for the matched normal and the PoN.
 *     <li>For the germline resource, the tool prefilters on the allele. If a variant is absent from a given germline
 *     resource, then the value for --af-of-alleles-not-in-resource applies such that an allele's absence from the
 *     germline resource becomes evidence that it is not a germline variant. Below is an excerpt of a known variants
 *     resource with population allele frequencies.
 *
 *     <pre>
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
 *     </pre>
 *     </li>
 *
 * <li>Additional parameters that factor towards filtering, including normal-artifact-lod (default threshold 0.0) and
 * tumor-lod (default threshold 5.3), are available in
 * <a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_FilterMutectCalls.php">FilterMutectCalls</a>.
 * While the tool calculates normal-lod assuming a diploid genotype, it calculates
 * normal-artifact-lod with the same approach it uses for tumor-lod, i.e. with a variable ploidy assumption.
 *
 * <ul>
 *     <li>If the normal artifact log odds becomes large, then FilterMutectCalls applies the artifact-in-normal filter.
 *     For matched normal samples with tumor contamination, consider increasing the normal-artifact-lod threshold.</li>
 *     <li>The tumor log odds, which is calculated independently of any matched normal, determines whether to filter a tumor
 *     variant. Variants with tumor LODs exceeding the threshold pass filtering.</li>
 * </ul></p>
 *</li>
 * </ol>
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
        if (m2Engine.emitReferenceConfidence()) {
            logger.warn("Note that the Mutect2 reference confidence mode is in BETA -- the likelihoods model and output format are subject to change in subsequent versions.");
            if ( MTAC.emitReferenceConfidence == ReferenceConfidenceMode.GVCF ) {
                try {
                    vcfWriter = new SomaticGVCFWriter(vcfWriter, new ArrayList<Number>(MTAC.GVCFGQBands));
                } catch ( IllegalArgumentException e ) {
                    throw new CommandLineException.BadArgumentValue("GQBands", "are malformed: " + e.getMessage());
                }
            }
        }
        m2Engine.writeHeader(vcfWriter, getDefaultToolVCFHeaderLines());
    }

    @Override
    public Collection<Annotation> makeVariantAnnotations(){
        final Collection<Annotation> annotations = super.makeVariantAnnotations();

        if (!MTAC.artifactPriorTables.isEmpty()){
            // Enable the annotations associated with the read orientation model
            annotations.add(new ReadOrientationArtifact(MTAC.artifactPriorTables));
            annotations.add(new ReferenceBases());
        }
        if (MTAC.autosomalCoverage > 0) {
            annotations.add(new PolymorphicNuMT(MTAC.autosomalCoverage));
        }
        if (MTAC.mitochondria) {
            annotations.add(new OriginalAlignment());
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
