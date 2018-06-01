package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.FileNotFoundException;
import java.util.*;
import java.util.ArrayList;
import java.util.List;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import java.nio.file.Path;


/**
 * Call germline SNPs and indels via local re-assembly of haplotypes
 *
 * <p>The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper.</p>
 *
 * <p>In the GVCF workflow used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to generate an intermediate GVCF (not to be used in final analysis), which can then be used in GenotypeGVCFs for joint genotyping of multiple samples in a very efficient way. The GVCF workflow enables rapid incremental processing of samples as they roll off the sequencer, as well as scaling to very large cohort sizes (e.g. the 92K exomes of ExAC).</p>
 *
 * <p>In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. Note however that the algorithms used to calculate variant likelihoods is not well suited to extreme allele frequencies (relative to ploidy) so its use is not recommended for somatic (cancer) variant discovery. For that purpose, use Mutect2 instead.</p>
 *
 * <p>Finally, HaplotypeCaller is also able to correctly handle the splice junctions that make RNAseq a challenge for most variant callers,
 * on the condition that the input read data has previously been processed according to our recommendations as documented <a href='https://software.broadinstitute.org/gatk/documentation/article?id=4067'>here</a>.</p>
 *
 * <h3>How HaplotypeCaller works</h3>
 *
 * <br />
 * <h4><a href='https://software.broadinstitute.org/gatk/documentation/article?id=4147'>1. Define active regions </a></h4>
 *
 * <p>The program determines which regions of the genome it needs to operate on (active regions), based on the presence of
 * evidence for variation.
 *
 * <br />
 * <h4><a href='https://software.broadinstitute.org/gatk/documentation/article?id=4146'>2. Determine haplotypes by assembly of the active region </a></h4>
 *
 * <p>For each active region, the program builds a De Bruijn-like graph to reassemble the active region and identifies
 * what are the possible haplotypes present in the data. The program then realigns each haplotype against the reference
 * haplotype using the Smith-Waterman algorithm in order to identify potentially variant sites. </p>
 *
 * <br />
 * <h4><a href='https://software.broadinstitute.org/gatk/documentation/article?id=4441'>3. Determine likelihoods of the haplotypes given the read data </a></h4>
 *
 * <p>For each active region, the program performs a pairwise alignment of each read against each haplotype using the
 * PairHMM algorithm. This produces a matrix of likelihoods of haplotypes given the read data. These likelihoods are
 * then marginalized to obtain the likelihoods of alleles for each potentially variant site given the read data.   </p>
 *
 * <br />
 * <h4><a href='https://software.broadinstitute.org/gatk/documentation/article?id=4442'>4. Assign sample genotypes </a></h4>
 *
 * <p>For each potentially variant site, the program applies Bayes' rule, using the likelihoods of alleles given the
 * read data to calculate the likelihoods of each genotype per sample given the read data observed for that
 * sample. The most likely genotype is then assigned to the sample.    </p>
 *
 * <h3>Input</h3>
 * <p>
 * Input bam file(s) from which to make variant calls
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Either a VCF or GVCF file with raw, unfiltered SNP and indel calls. Regular VCFs must be filtered either by variant
 * recalibration (Best Practice) or hard-filtering before use in downstream analyses. If using the GVCF workflow, the
 * output is a GVCF file that must first be run through GenotypeGVCFs and then filtering before further analysis.
 * </p>
 *
 * <h3>Usage examples</h3>
 *
 * <p>These are example commands that show how to run HaplotypeCaller for typical use cases. Have a look at the <a href='https://software.broadinstitute.org/gatk/documentation/article?id=3893'>method documentation</a> for the basic GVCF workflow. </p>
 *
 * <br />
 * <h4>Single-sample GVCF calling (outputs intermediate GVCF)</h4>
 * <pre>
 * gatk --java-options "-Xmx4g" HaplotypeCaller  \
 *   -R Homo_sapiens_assembly38.fasta \
 *   -I input.bam \
 *   -O output.g.vcf.gz \
 *   -ERC GVCF
 * </pre>
 *
 * <h4>Single-sample GVCF calling with <a href='https://software.broadinstitute.org/gatk/documentation/article?id=9622'>allele-specific annotations</a></h4>
 * <pre>
 * gatk --java-options "-Xmx4g" HaplotypeCaller  \
 *   -R Homo_sapiens_assembly38.fasta \
 *   -I input.bam \
 *   -O output.g.vcf.gz \
 *   -ERC GVCF \
 *   -G Standard \
 *   -G AS_Standard
 * </pre>
 *
 * <h4>Variant calling with <a href='https://software.broadinstitute.org/gatk/documentation/article?id=5484'>bamout</a> to show realigned reads</h4>
 * <pre>
 * gatk --java-options "-Xmx4g" HaplotypeCaller  \
 *   -R Homo_sapiens_assembly38.fasta \
 *   -I input.bam \
 *   -O output.vcf.gz \
 *   -bamout bamout.bam
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>We have not yet fully tested the interaction between the GVCF-based calling or the multisample calling and the
 * RNAseq-specific functionalities. Use those in combination at your own risk.</li>
 * </ul>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool is able to handle many non-diploid use cases; the desired ploidy can be specified using the -ploidy
 * argument. Note however that very high ploidies (such as are encountered in large pooled experiments) may cause
 * performance challenges including excessive slowness. We are working on resolving these limitations.</p>
 *
 * <h3>Additional Notes</h3>
 * <ul>
 *     <li>When working with PCR-free data, be sure to set `-pcr_indel_model NONE` (see argument below).</li>
 *     <li>When running in `-ERC GVCF` or `-ERC BP_RESOLUTION` modes, the confidence threshold
 *     is automatically set to 0. This cannot be overridden by the command line. The threshold can be set manually
 *     to the desired level in the next step of the workflow (GenotypeGVCFs)</li>
 *     <li>We recommend using a list of intervals to speed up analysis. See <a href='https://software.broadinstitute.org/gatk/documentation/article?id=4133'>this document</a> for details.</li>
 * </ul>
 *
 */
@CommandLineProgramProperties(
        summary = "Call germline SNPs and indels via local re-assembly of haplotypes",
        oneLineSummary = "Call germline SNPs and indels via local re-assembly of haplotypes",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public final class HaplotypeCaller extends AssemblyRegionWalker {

    //NOTE: many of these settings are referenced by HaplotypeCallerSpark
    public static final int DEFAULT_MIN_ASSEMBLY_REGION_SIZE = 50;
    public static final int DEFAULT_MAX_ASSEMBLY_REGION_SIZE = 300;
    public static final int DEFAULT_ASSEMBLY_REGION_PADDING = 100;
    public static final int DEFAULT_MAX_READS_PER_ALIGNMENT = 50;
    public static final double DEFAULT_ACTIVE_PROB_THRESHOLD = 0.002;
    public static final int DEFAULT_MAX_PROB_PROPAGATION_DISTANCE = 50;
    @ArgumentCollection
    private HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();

    /**
     * A raw, unfiltered, highly sensitive callset in VCF format.
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to which variants should be written")
    public String outputVCF = null;

    private VariantContextWriter vcfWriter;

    private HaplotypeCallerEngine hcEngine;

    @Override
    protected int defaultMinAssemblyRegionSize() { return DEFAULT_MIN_ASSEMBLY_REGION_SIZE; }

    @Override
    protected int defaultMaxAssemblyRegionSize() { return DEFAULT_MAX_ASSEMBLY_REGION_SIZE; }

    @Override
    protected int defaultAssemblyRegionPadding() { return DEFAULT_ASSEMBLY_REGION_PADDING; }

    @Override
    protected int defaultMaxReadsPerAlignmentStart() { return DEFAULT_MAX_READS_PER_ALIGNMENT; }

    @Override
    protected double defaultActiveProbThreshold() { return DEFAULT_ACTIVE_PROB_THRESHOLD; }

    @Override
    protected int defaultMaxProbPropagationDistance() { return DEFAULT_MAX_PROB_PROPAGATION_DISTANCE; }

    @Override
    protected boolean includeReadsWithDeletionsInIsActivePileups() { return true; }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return HaplotypeCallerEngine.makeStandardHCReadFilters();
    }

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() { return HaplotypeCallerEngine.getStandardHaplotypeCallerAnnotationGroups();}

    @Override
    public boolean useVariantAnnotations() { return true;}

    /**
     * If we are in reference confidence mode we want to filter the annotations as there are certain annotations in the standard
     * HaplotypeCaller set which are no longer relevant, thus we filter them out before constructing the
     * VariantAnnotationEngine because the user args will have been parsed by that point.
     *
     * @see GATKTool#makeVariantAnnotations()
     * @return a collection of annotation arguments with alterations depending on hcArgs.emitReferenceConfidence
     */
    @Override
    public Collection<Annotation> makeVariantAnnotations() {
        final boolean confidenceMode = hcArgs.emitReferenceConfidence != ReferenceConfidenceMode.NONE;
        final Collection<Annotation> annotations = super.makeVariantAnnotations();
        return confidenceMode? HaplotypeCallerEngine.filterReferenceConfidenceAnnotations(annotations): annotations;
    }

    @Override
    public AssemblyRegionEvaluator assemblyRegionEvaluator() {
        return hcEngine;
    }

    @Override
    public void onTraversalStart() {
        if (hcArgs.emitReferenceConfidence == ReferenceConfidenceMode.GVCF && hcArgs.maxMnpDistance > 0) {
            throw new CommandLineException.BadArgumentValue("Non-zero maxMnpDistance is incompatible with GVCF mode.");
        }
        final ReferenceSequenceFile referenceReader = getReferenceReader(referenceArguments);
        final VariantAnnotatorEngine variantAnnotatorEngine = new VariantAnnotatorEngine(makeVariantAnnotations(),
                hcArgs.dbsnp.dbsnp, hcArgs.comps,  hcArgs.emitReferenceConfidence != ReferenceConfidenceMode.NONE);
        hcEngine = new HaplotypeCallerEngine(hcArgs, createOutputBamIndex, createOutputBamMD5, getHeaderForReads(), referenceReader, variantAnnotatorEngine);

        // The HC engine will make the right kind (VCF or GVCF) of writer for us
        final SAMSequenceDictionary sequenceDictionary = getHeaderForReads().getSequenceDictionary();
        vcfWriter = hcEngine.makeVCFWriter(outputVCF, sequenceDictionary, createOutputVariantIndex, createOutputVariantMD5, outputSitesOnlyVCFs);
        hcEngine.writeHeader(vcfWriter, sequenceDictionary, getDefaultToolVCFHeaderLines());
    }

    private static CachingIndexedFastaSequenceFile getReferenceReader(ReferenceInputArgumentCollection referenceArguments) {
        final CachingIndexedFastaSequenceFile referenceReader;
        final Path reference = IOUtils.getPath(referenceArguments.getReferenceFileName());
        try {
            referenceReader = new CachingIndexedFastaSequenceFile(reference);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotReadInputFile(reference, e);
        }
        return referenceReader;
    }

    @Override
    public void apply(final AssemblyRegion region, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        hcEngine.callRegion(region, featureContext).forEach(vcfWriter::add);
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }

        if ( hcEngine != null ) {
            hcEngine.shutdown();
        }
    }
}
