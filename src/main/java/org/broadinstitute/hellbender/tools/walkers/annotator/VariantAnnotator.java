package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


/**
 * Annotate variant calls with context information
 *
 * <p>
 * This tool is designed to annotate variant calls based on their context (as opposed to functional annotation).
 * Various annotation modules are available; see the "Annotation Modules" page linked in the Tool Documentation sidebar for a complete list.
 *
 * <h3>Input</h3>
 * <p>
 * A variant set to annotate and optionally one or more BAM files.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * An annotated VCF.
 * </p>
 *
 * <h3>Usage examples</h3>
 * <br />
 *
 * <h4>Annotate a VCF with dbSNP IDs and depth of coverage for each sample</h4>
 * <pre>
 *   VariantAnnotator \
 *   -R reference.fasta \
 *   -I input.bam \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -A Coverage \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 * <h4>Annotate a VCF with allele frequency by an external resource. Annotation will only occur if there is allele concordance between the resource and the input VCF </h4>
 * <pre>
 *   VariantAnnotator \
 *   -R reference.fasta \
 *   -I input.bam \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -L anotherInput.vcf \
 *   --resource foo:resource.vcf \
 *   -E foo.AF \
 *   --resource-allele-concordance
 * </pre>
 *
 * <h4>Annotate with AF and FILTER fields from an external resource </h4>
 * <pre>
 *   VariantAnnotator \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o output.vcf \
 *   --resource foo:resource.vcf \
 *   --expression foo.AF \
 *   --expression foo.FILTER
 * </pre>
 *
 * <h3>Caveat</h3>
 * <p>This tool outputs no annotations by default, all annotations/groups must be specified explicitly. </p>
 *
 * <h3>Special note on RankSumTestAnnotations</h3>
 * <p>RankSumAnnotations produced by this tool are not the same as those produced by the HaplotypeCaller. Without the
 * likelihoods, the tool resorts to a pileup heuristic to categorize reads which means that RankSumAnnotations will only
 * be present for SNP variants.</p>
 */
@CommandLineProgramProperties(summary="Tool for adding annotations to VCF files",
        oneLineSummary = "Tool for adding annotations to VCF files",
        programGroup = VariantManipulationProgramGroup.class)
@DocumentedFeature
public class VariantAnnotator extends VariantWalker {
    public static final String EXPRESSION_LONG_NAME = "expression";
    public static final String EXPRESSION_SHORT_NAME = "E";

    private  SAMSequenceDictionary sequenceDictionary;
    private VariantContextWriter vcfWriter;
    private static final int REFERENCE_PADDING = 100;

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     * as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field). Records that are
     * filtered in the comp track will be ignored. Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    @Advanced
    @Argument(fullName = StandardArgumentDefinitions.COMPARISON_LONG_NAME, shortName = StandardArgumentDefinitions.COMPARISON_SHORT_NAME, doc = "Comparison VCF file(s)", optional = true)
    public List<FeatureInput<VariantContext>> comps = new ArrayList<>();

    /**
     * An external resource VCF file or files from which to annotate.
     *
     * Use this option to add annotations from a resource file to the output.
     * For example, if you want to annotate your callset with the AC field value from a VCF file named
     * 'resource_file.vcf', you tag it with '-resource:my_resource resource_file.vcf' and you additionally specify
     * '-E my_resource.AC' (-E is short for --expression, also documented on this page). In the resulting output
     * VCF, any records for which there is a record at the same position in the resource file will be annotated with
     * 'my_resource.AC=N'. Note that if there are multiple records in the resource file that overlap the given
     * position, one is chosen randomly. Check for allele concordance if using --resource-allele-concordance, otherwise
     * the match is based on position only.
     */
    @Argument(fullName= StandardArgumentDefinitions.RESOURCE_LONG_NAME, doc="External resource VCF file", optional=true)
    public List<FeatureInput<VariantContext>> resources;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The file to which variants should be written", optional=false)
    protected File outputFile;

    /**
     * This option enables you to add annotations from one VCF to another.
     *
     * For example, if you want to annotate your callset with the AC field value from a VCF file named
     * 'resource_file.vcf', you tag it with '-resource:my_resource resource_file.vcf' (see the -resource argument, also
     * documented on this page) and you specify '-E my_resource.AC'. In the resulting output VCF, any records for
     * which there is a record at the same position in the resource file will be annotated with 'my_resource.AC=N'.
     * INFO field data, ID, ALT, and FILTER fields may be used as expression values.
     * Note that if there are multiple records in the resource file that overlap the given position, one is chosen
     * randomly.
     */
    @Argument(fullName= EXPRESSION_LONG_NAME, shortName= EXPRESSION_SHORT_NAME, doc="One or more specific expressions to apply to variant calls", optional=true)
    protected Set<String> expressionsToUse = new HashSet<>();

    /**
     * If this argument is specified, add annotations (specified by --expression) from an external resource
     * (specified by --resource) to the input VCF (specified by --variant) only if the alleles are
     * concordant between input and the resource VCFs. Otherwise, always add the annotations.
     */
    @Argument(fullName="resource-allele-concordance", shortName="rac", doc="Check for allele concordances when using an external resource VCF file", optional=true)
    protected Boolean expressionAlleleConcordance = false;

    /**
     * Bases with a quality below this threshold will not be used for annotation.
     */
    @Argument(fullName = AssemblyBasedCallerArgumentCollection.MIN_BASE_QUALITY_SCORE_LONG_NAME,  doc = "Minimum base quality required to confidently assign a read to an allele", optional = true)
    public byte minBaseQualityScore = 10;

    public List<ReadFilter> getDefaultReadFilters() {
        return Lists.newArrayList( new WellformedReadFilter(),
                ReadFilterLibrary.NOT_DUPLICATE,
                ReadFilterLibrary.PRIMARY_LINE,
                ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK,
                ReadFilterLibrary.MAPPED);
    }

    @Override
    public boolean useVariantAnnotations() { return true;}

    private VariantAnnotatorEngine annotatorEngine;
    private SampleList variantSamples;

    /**
     * Prepare the output file and the list of available features.
     */
    public void onTraversalStart() {
        // get the list of all sample names from the variant VCF, if applicable
        final List<String> samples = getHeaderForVariants().getGenotypeSamples();
        sequenceDictionary = getBestAvailableSequenceDictionary();
        variantSamples = new IndexedSampleList(samples);

        // Check that the reads have compatible samples to the variants:
        if (hasReads()) {
            final Set<String> readsSamples = getHeaderForReads().getReadGroups().stream().map(rg -> rg.getSample()).collect(Collectors.toSet());
            readsSamples.forEach(readSample -> {
                if (!samples.contains(readSample))
                    throw new UserException(String.format("Reads sample '%s' from readgroups tags does not match any sample in the variant genotypes", readSample));
            });
        }
        annotatorEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), dbsnp.dbsnp, comps, false, false);
        annotatorEngine.addExpressions(expressionsToUse, resources, expressionAlleleConcordance);

        // setup the header fields
        // note that if any of the definitions conflict with our new ones, then we want to overwrite the old ones
        Set<VCFHeaderLine> hInfo = new HashSet<>();
        hInfo.addAll(annotatorEngine.getVCFAnnotationDescriptions(false));
        hInfo.addAll(getHeaderForVariants().getMetaDataInInputOrder());
        hInfo = VcfUtils.updateHeaderContigLines(
                hInfo,
                hasReference() ? referenceArguments.getReferencePath() : null,
                hasReference() ? getReferenceDictionary() : sequenceDictionary,
                false);

        vcfWriter = createVCFWriter(outputFile);
        vcfWriter.writeHeader(new VCFHeader(hInfo, samples));
    }

    /**
     * For each site of interest, annotate based on the requested annotation types
     *
     * @param vc
     * @param readsContext Reads overlapping the current variant. Will be an empty, but non-null, context object
     *                     if there is no backing source of reads data (in which case all queries on it will return
     *                     an empty array/iterator)
     * @param refContext
     * @param fc
     */
    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {

        // if the reference is present and base is not ambiguous, we can annotate
        if (refContext.getBases().length ==0 || BaseUtils.simpleBaseToBaseIndex(refContext.getBase()) != -1 ) {
            final ReferenceContext expandedRefContext = !hasReference() ? refContext :
                    new ReferenceContext(refContext, new SimpleInterval(vc).expandWithinContig(REFERENCE_PADDING, sequenceDictionary));
            VariantContext annotatedVC = annotatorEngine.annotateContext(vc, fc, expandedRefContext, makeLikelihoods(vc, readsContext), a -> true);
            vcfWriter.add(annotatedVC);
        } else {
            vcfWriter.add(vc);
        }
    }

    private AlleleLikelihoods<GATKRead, Allele> makeLikelihoods(final VariantContext vc, final ReadsContext readsContext) {
        final List<GATKRead> reads = Utils.stream(readsContext).collect(Collectors.toList());
        final AlleleLikelihoods<GATKRead, Allele> result = new AlleleLikelihoods<>(variantSamples,
                new IndexedAlleleList<>(vc.getAlleles()), AssemblyBasedCallerUtils.splitReadsBySample(variantSamples, getHeaderForReads(), reads));

        final ReadPileup readPileup = new ReadPileup(vc, readsContext);
        final Map<String, ReadPileup> pileupsBySample = readPileup.splitBySample(getHeaderForReads(), "__UNKNOWN__");

        final Allele refAllele = vc.getReference();
        final List<Allele> altAlleles = vc.getAlternateAlleles();
        final int numAlleles = vc.getNAlleles();

        // manually fill each sample's likelihoods one read at a time, assigning probability 1 (0 in log space) to the matching allele, if present
        for (final Map.Entry<String, ReadPileup> samplePileups : pileupsBySample.entrySet()) {
            final LikelihoodMatrix<GATKRead, Allele> sampleMatrix = result.sampleMatrix(result.indexOfSample(samplePileups.getKey()));

            final MutableInt readIndex = new MutableInt(0);
            for (final PileupElement pe : samplePileups.getValue()) {
                // initialize all alleles as equally unlikely
                IntStream.range(0, numAlleles).forEach(a -> sampleMatrix.set(a, readIndex.intValue(), Double.NEGATIVE_INFINITY));

                // if present set the matching allele as likely
                final Allele pileupAllele = GATKVariantContextUtils.chooseAlleleForRead(pe, refAllele, altAlleles, minBaseQualityScore);
                if (pileupAllele != null) {
                    sampleMatrix.set(sampleMatrix.indexOfAllele(pileupAllele), readIndex.intValue(), 0);
                }


                readIndex.increment();
            }
        }

        return result;
    }

    /**
     * Make sure that the writer is closed upon completing the tool.
     */
    @Override
    public void closeTool() {
        if (vcfWriter !=null) {
            vcfWriter.close();
        }
    }
}

