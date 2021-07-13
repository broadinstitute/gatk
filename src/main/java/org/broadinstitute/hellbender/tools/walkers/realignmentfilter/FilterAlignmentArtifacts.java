package org.broadinstitute.hellbender.tools.walkers.realignmentfilter;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import com.google.common.collect.Multisets;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedOnStart;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.CountingVariantFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>Filter false positive alignment artifacts from a VCF callset.</p>
 *
 * <p>
 *     FilterAlignmentArtifacts identifies alignment artifacts, that is, apparent variants due to reads being mapped to the wrong genomic locus.
 * </p>
 * <p>
 *     Alignment artifacts can occur whenever there is sufficient sequence similarity between two or more regions in the genome
 *     to confuse the alignment algorithm.  This can occur when the aligner for whatever reason overestimate how uniquely a read
 *     maps, thereby assigning it too high of a mapping quality.  It can also occur through no fault of the aligner due to gaps in
 *     the reference, which can also hide the true position to which a read should map.  By using a good alignment algorithm
 *     (the GATK wrapper of BWA-MEM), giving it sensitive settings (which may have been impractically slow for the original
 *     bam alignment) and mapping to the best available reference we can avoid these pitfalls.  The last point is especially important:
 *     one can (and should) use a BWA-MEM index image corresponding to the best reference, regardless of the reference to which
 *     the bam was aligned.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 * </p>
 * <p>
 *     The input bam to this tool should be the same tumor bam that Mutect2 was run on.  The reference passed with the -R argument
 *     must be the reference to which the input bam was aligned.  This does not need to correspond to the reference of the BWA-MEM
 *     index image.  The latter should be derived from the best available reference, for example hg38 in humans as of February 2018.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk FilterAlignmentArtifacts \
 *   -R hg19.fasta
 *   -V somatic.vcf.gz \
 *   -I somatic_bamout.bam \
 *   --bwa-mem-index-image hg38.index_image \
 *   -O filtered.vcf.gz
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Filter alignment artifacts from a vcf callset.",
        oneLineSummary = "Filter alignment artifacts from a vcf callset.",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class FilterAlignmentArtifacts extends MultiVariantWalkerGroupedOnStart {
    public static final int DEFAULT_DISTANCE_TO_GROUP_VARIANTS = 1000;
    public static final int DEFAULT_REF_PADDING = 100;
    public static final int DEFAULT_MAX_GROUPED_SPAN = 10_000;
    private static final int MIN_UNITIG_LENGTH = 30;
    private static final int ASSEMBLY_PADDING = 50;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output filtered VCF file", optional=false)
    private final GATKPath outputVcf = null;

    public static final int DEFAULT_INDEL_START_TOLERANCE = 5;
    public static final String INDEL_START_TOLERANCE_LONG_NAME = "indel-start-tolerance";
    @Argument(fullName = INDEL_START_TOLERANCE_LONG_NAME, doc="Max distance between indel start of aligned read in the bam and the variant in the vcf", optional=true)
    private int indelStartTolerance = DEFAULT_INDEL_START_TOLERANCE;

    public static final int DEFAULT_KMER_SIZE = 21;
    public static final String KMER_SIZE_LONG_NAME = "kmer-size";
    @Argument(fullName = KMER_SIZE_LONG_NAME, doc="Kmer size for reassembly", optional=true)
    private int kmerSize = DEFAULT_KMER_SIZE;

    public static final String DONT_SKIP_ALREADY_FILTERED_VARIANTS_LONG_NAME = "dont-skip-filtered-variants";
    @Argument(fullName = DONT_SKIP_ALREADY_FILTERED_VARIANTS_LONG_NAME,
            doc="Try to realign all variants, even ones that have already been filtered.", optional=true)
    private boolean dontSkipFilteredVariants = false;

    @Argument(fullName= AssemblyBasedCallerArgumentCollection.BAM_OUTPUT_LONG_NAME, shortName= AssemblyBasedCallerArgumentCollection.BAM_OUTPUT_SHORT_NAME, doc="File to which assembled haplotypes should be written", optional = true)
    public String bamOutputPath = null;

    @Advanced
    @Argument(fullName = AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_LONG_NAME, doc = "Which Smith-Waterman implementation to use, generally FASTEST_AVAILABLE is the right choice", optional = true)
    public SmithWatermanAligner.Implementation smithWatermanImplementation = SmithWatermanAligner.Implementation.JAVA;


    @ArgumentCollection
    protected RealignmentArgumentCollection realignmentArgumentCollection = new RealignmentArgumentCollection();

    private SmithWatermanAligner smithWatermanAligner;
    private VariantContextWriter vcfWriter;
    private RealignmentEngine realignmentEngine;
    private SAMFileHeader bamHeader;
    private SampleList samplesList;
    private CachingIndexedFastaSequenceFile referenceReader;
    private ReadThreadingAssembler assemblyEngine;
    private final M2ArgumentCollection MTAC = new M2ArgumentCollection();
    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine;
    private Optional<HaplotypeBAMWriter> haplotypeBAMWriter;

    @Override
    public List<ReadFilter> getDefaultReadFilters() { return Mutect2Engine.makeStandardMutect2ReadFilters(); }

    @Override
    protected CountingVariantFilter makeVariantFilter() {
        return new CountingVariantFilter(dontSkipFilteredVariants ? VariantFilterLibrary.ALLOW_ALL_VARIANTS : VariantFilterLibrary.PASSES_FILTERS);
    }

    @Override
    public boolean requiresReads() { return true; }

    @Override
    protected int defaultDistanceToGroupVariants() { return DEFAULT_DISTANCE_TO_GROUP_VARIANTS; }

    @Override
    protected int defaultReferenceWindowPadding() { return DEFAULT_REF_PADDING; }

    @Override
    protected int defaultMaxGroupedSpan() {
        return DEFAULT_MAX_GROUPED_SPAN;
    }

    @Override
    public void onTraversalStart() {
        smithWatermanAligner = SmithWatermanAligner.getAligner(smithWatermanImplementation);
        realignmentEngine = new RealignmentEngine(realignmentArgumentCollection);
        vcfWriter = createVCFWriter(outputVcf);

        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeader.getMetaDataInSortedOrder());
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.ALIGNMENT_ARTIFACT_FILTER_NAME));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.UNITIG_SIZES_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ALIGNMENT_SCORE_DIFFERENCE_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.JOINT_ALIGNMENT_COUNT_KEY));
        headerLines.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter.writeHeader(vcfHeader);
        bamHeader = getHeaderForReads();
        samplesList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSamplesFromHeader(bamHeader)));
        referenceReader = ReferenceUtils.createReferenceReader(Utils.nonNull(referenceArguments.getReferenceSpecifier()));
        assemblyEngine = MTAC.createReadThreadingAssembler();
        likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(MTAC.likelihoodArgs, true);
        haplotypeBAMWriter = bamOutputPath == null ? Optional.empty() :
                Optional.of(new HaplotypeBAMWriter(HaplotypeBAMWriter.WriterType.ALL_POSSIBLE_HAPLOTYPES, IOUtils.getPath(bamOutputPath), true, false, getHeaderForSAMWriter()));
    }

    @Override
    public Object onTraversalSuccess() {
        return "SUCCESS";
    }

    @Override
    public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext, final List<ReadsContext> readsContexts) {


        // for now we do one variant at a time but eventually we will want to combine all reads supporting all variants
        // into a single graph.  This is non-trivial because there may be more than one phasing between variants.
        for (final VariantContext vc : variantContexts) {
            final AssemblyRegion assemblyRegion = makeAssemblyRegionFromVariantReads(readsContexts, vc);


            // TODO: give this tool M2 Assembler args to allow override default M2ArgumentCollection?
            final AssemblyResultSet assemblyResult = AssemblyBasedCallerUtils.assembleReads(assemblyRegion, Collections.emptyList(), MTAC, bamHeader, samplesList, logger, referenceReader, assemblyEngine, smithWatermanAligner, false);
            final AssemblyRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();

            final Map<String,List<GATKRead>> reads = AssemblyBasedCallerUtils.splitReadsBySample(samplesList, bamHeader, regionForGenotyping.getReads());

            final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult,samplesList,reads);
            readLikelihoods.switchToNaturalLog();
            final Map<GATKRead,GATKRead> readRealignments = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc(), smithWatermanAligner);
            readLikelihoods.changeEvidence(readRealignments);
            writeBamOutput(assemblyResult, readLikelihoods, new HashSet<>(readLikelihoods.alleles()), regionForGenotyping.getSpan());

            final LocusIteratorByState libs = new LocusIteratorByState(regionForGenotyping.getReads().iterator(), DownsamplingMethod.NONE, false, samplesList.asListOfSamples(), bamHeader, true);

            final List<byte[]> unitigs = getUnitigs(libs);

            final VariantContextBuilder vcb = new VariantContextBuilder(vc)
                    .attribute(GATKVCFConstants.UNITIG_SIZES_KEY, unitigs.stream().mapToInt(u -> u.length).toArray());

            final List<List<BwaMemAlignment>> unitigAlignments = unitigs.stream()
                    .map(realignmentEngine::realign).collect(Collectors.toList());

            final List<List<BwaMemAlignment>> jointAlignments = RealignmentEngine.findJointAlignments(unitigAlignments, realignmentArgumentCollection.maxReasonableFragmentLength);
            vcb.attribute(GATKVCFConstants.JOINT_ALIGNMENT_COUNT_KEY, jointAlignments.size());
            jointAlignments.sort(Comparator.comparingInt(FilterAlignmentArtifacts::jointAlignmentScore).reversed());

            // best mapping to another contig
            if (!jointAlignments.isEmpty() && jointAlignments.get(0).get(0).getRefId() != getReferenceDictionary().getSequenceIndex(vc.getContig())) {
                vcb.filter(GATKVCFConstants.ALIGNMENT_ARTIFACT_FILTER_NAME);
            } else if (jointAlignments.size() > 1) {

                final int totalBases = unitigs.stream().mapToInt(unitig -> unitig.length).sum();
                final int scoreDiff = jointAlignmentScore(jointAlignments.get(0)) - jointAlignmentScore(jointAlignments.get(1));
                final int mismatchDiff = totalMismatches(jointAlignments.get(1)) - totalMismatches(jointAlignments.get(0));

                vcb.attribute(GATKVCFConstants.ALIGNMENT_SCORE_DIFFERENCE_KEY, scoreDiff);

                final boolean multimapping = (double) scoreDiff / totalBases < realignmentArgumentCollection.minAlignerScoreDifferencePerBase
                        && (double) mismatchDiff / totalBases < realignmentArgumentCollection.minMismatchDifferencePerBase;

                if (multimapping) {
                    vcb.filter(GATKVCFConstants.ALIGNMENT_ARTIFACT_FILTER_NAME);
                }
            }

            vcfWriter.add(vcb.make());
        }
    }

    private AssemblyRegion makeAssemblyRegionFromVariantReads(final List<ReadsContext> readsContexts, final VariantContext vc) {
        final Set<String> variantReadNames = readsContexts.stream().flatMap(Utils::stream)
                .filter(read -> RealignmentEngine.supportsVariant(read, vc, indelStartTolerance))
                .map(GATKRead::getName)
                .collect(Collectors.toSet());

        final List<GATKRead> variantReads = readsContexts.stream().flatMap(Utils::stream)
                .filter(read -> variantReadNames.contains(read.getName()))
                .sorted(Comparator.comparingInt(GATKRead::getStart))
                .collect(Collectors.toList());

        final int firstReadStart = variantReads.stream().mapToInt(GATKRead::getStart).min().orElse(vc.getStart());
        final int lastReadEnd = variantReads.stream().mapToInt(GATKRead::getEnd).max().orElse(vc.getEnd());
        final SimpleInterval assemblyWindow = new SimpleInterval(vc.getContig(), Math.max(firstReadStart - ASSEMBLY_PADDING,1), lastReadEnd + ASSEMBLY_PADDING);

        final AssemblyRegion assemblyRegion = new AssemblyRegion(assemblyWindow, 0, bamHeader);
        assemblyRegion.addAll(variantReads);

        return assemblyRegion;
    }

    private void writeBamOutput(final AssemblyResultSet assemblyResult, final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods, final Set<Haplotype> haplotypes, Locatable callableRegion) {
        haplotypeBAMWriter.ifPresent(writer -> writer.writeReadsAlignedToHaplotypes(
                assemblyResult.getHaplotypeList(),
                assemblyResult.getPaddedReferenceLoc(),
                assemblyResult.getHaplotypeList(),
                haplotypes,
                readLikelihoods,
                callableRegion));
    }

    // TODO: what about deletions in pileup?
    private List<byte[]> getUnitigs(final LocusIteratorByState libs) {
        final List<StringBuilder> unitigBuilders = new ArrayList<>();
        int lastCoveredLocus = Integer.MIN_VALUE;
        while (libs.hasNext()) {
            final ReadPileup pileup = libs.next().getBasePileup();
            if (pileup.isEmpty()) {
                continue;
            }

            // begin new unitig if this pileup isn't contiguous with the last
            final int currentLocus = pileup.getLocation().getStart();
            if (currentLocus != lastCoveredLocus + 1) {
                unitigBuilders.add(new StringBuilder());
            }
            lastCoveredLocus = currentLocus;
            final StringBuilder currentUnitigBuilder = unitigBuilders.get(unitigBuilders.size() - 1);

            // add no bases (deletion) or consensus bases.
            final int[] baseCounts = pileup.getBaseCounts();
            final int deletionCount = (int) Utils.stream(pileup).filter(PileupElement::isDeletion).count();
            if (deletionCount < pileup.size() / 2) {
                final byte consensusBase = BaseUtils.baseIndexToSimpleBase(MathUtils.maxElementIndex(baseCounts));
                currentUnitigBuilder.append((char) consensusBase);

                // in addition to consensus base, add inserted bases if needed
                final Multiset<String> insertedBases = Utils.stream(pileup)
                        .map(PileupElement::getBasesOfImmediatelyFollowingInsertion)
                        .filter(s -> s != null)
                        .collect(Collectors.toCollection(HashMultiset::create));

                if (insertedBases.size() > pileup.size() / 2) {
                    final String consensusInsertion = Multisets.copyHighestCountFirst(insertedBases).entrySet().iterator().next().getElement();
                    currentUnitigBuilder.append(consensusInsertion);
                }
            }
        }

        return unitigBuilders.stream()
                .map(builder -> builder.toString().getBytes())
                .filter(unitig -> unitig.length > MIN_UNITIG_LENGTH)
                .collect(Collectors.toList());
    }

    private static int jointAlignmentScore(final List<BwaMemAlignment> alignments) {
        return alignments.stream().mapToInt(BwaMemAlignment::getAlignerScore).sum();
    }

    private static int totalMismatches(final List<BwaMemAlignment> alignments) {
        return alignments.stream().mapToInt(BwaMemAlignment::getNMismatches).sum();
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

}
