package org.broadinstitute.hellbender.tools.walkers;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.VariantAnnotationArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedOnStart;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.*;

/**
 * Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file
 *
 * <p>
 * CombineGVCFs is meant to be used for hierarchical merging of gVCFs that will eventually be input into GenotypeGVCFs.
 * One would use this tool when needing to genotype too large a number of individual gVCFs; One would first use
 * CombineGVCFs to combine them into a single GVCF and pass this into GenotypeGVCFs.</p>
 *
 * <h3>Input</h3>
 * <p>
 * Two or more Haplotype Caller gVCFs to combine.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined multisample gVCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * ./gatk-launch \
 *   CombineGVCFs \
 *   -R reference.fasta \
 *   --variant sample1.g.vcf \
 *   --variant sample2.g.vcf \
 *   -O cohort.g.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <p>Only gVCF files produced by HaplotypeCaller (or CombineGVCFs) can be used as input for this tool. Some other
 * programs produce files that they call gVCFs but those lack some important information (accurate genotype likelihoods
 * for every position) that GenotypeGVCFs requires for its operation.</p>
 * <p>If the gVCF files contain allele specific annotations, add -G Standard -G AS_Standard to the command line.</p>
 *
 */
@CommandLineProgramProperties(summary = "Merges one or more HaplotypeCaller GVCF files into a single GVCF with appropriate annotations", oneLineSummary = "Merges one or more HaplotypeCaller GVCF files into a single GVCF with appropriate annotations", programGroup = VariantProgramGroup.class)
@DocumentedFeature
public final class CombineGVCFs extends MultiVariantWalkerGroupedOnStart {

    private VariantAnnotatorEngine annotationEngine;
    private VariantContextWriter vcfWriter;
    private SAMSequenceDictionary sequenceDictionary;
    private ReferenceConfidenceVariantContextMerger referenceConfidenceVariantContextMerger;

    /**
     * Which groups of annotations to add to the output VCF file.
     */
    @ArgumentCollection
    public VariantAnnotationArgumentCollection variantAnnotationArgumentCollection = new VariantAnnotationArgumentCollection(
            Arrays.asList(StandardAnnotation.class.getSimpleName()),
            Collections.emptyList(),
            Collections.emptyList());

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The combined GVCF output file", optional=false)
    private File outputFile;

    @Argument(fullName="convertToBasePairResolution", shortName="bpResolution", doc = "If specified, convert banded gVCFs to all-sites gVCFs", optional=true)
    protected boolean useBpResolution = false;

    /**
     * To reduce file sizes our gVCFs group similar reference positions into bands.  However, there are cases when users will want to know that no bands
     * span across a given genomic position (e.g. when scatter-gathering jobs across a compute farm).  The option below enables users to break bands at
     * pre-defined positions.  For example, a value of 10,000 would mean that we would ensure that no bands span across chr1:10000, chr1:20000, etc.
     *
     * Note that the --convertToBasePairResolution argument is just a special case of this argument with a value of 1.
     */
    @Argument(fullName="breakBandsAtMultiplesOf", shortName="breakBandsAtMultiplesOf", doc = "If > 0, reference bands will be broken up at genomic positions that are multiples of this number", optional=true)
    protected int multipleAtWhichToBreakBands = 0;

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    // the annotation engine
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // State that gets accumulated between calls of apply()
    private final LinkedList<VariantContext> VariantContextsOverlappingCurrentMerge = new LinkedList<>();
    private final Set<String> samples = new HashSet<>();
    private SimpleInterval prevPos = null;
    private byte refAfterPrevPos;
    private ReferenceContext storedReferenceContext;


    public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext) {
        // If we need to stop at an intermediate site since the last apply, do so (caused by gvcfBlocks, contexts ending, etc...)
        if (!VariantContextsOverlappingCurrentMerge.isEmpty()) {
            Locatable last = prevPos==null ? VariantContextsOverlappingCurrentMerge.get(0) : prevPos;
            // If on a different contig, close out all the queued states on the current contig
            int end = last.getContig().equals(referenceContext.getWindow().getContig())
                    ? referenceContext.getInterval().getStart() - 1
                    : VariantContextsOverlappingCurrentMerge.stream().mapToInt(VariantContext::getEnd).max().getAsInt();

            createIntermediateVariants( new SimpleInterval(last.getContig(), last.getStart(), end));
        }

        mergeWithnewVCs(variantContexts, referenceContext);

        // Update the stored reference if it has a later stop position than the current stored reference
        if ( (storedReferenceContext == null) ||
                (!referenceContext.getWindow().contigsMatch(storedReferenceContext.getWindow()) ) ||
                (storedReferenceContext.getWindow().getEnd() < referenceContext.getWindow().getEnd())) {
            storedReferenceContext = referenceContext;
        }
    }

    /**
     * Method that calculates if there are any sites in the provided interval where we should expect the tool to create
     * a new variant context object by calling endPreviousStates() and closes them by providing appropriate reference
     * information and an empty list of new variant contexts.
     *
     */
    @VisibleForTesting
    void createIntermediateVariants(SimpleInterval intervalToClose) {
        Set<Integer> sitesToStop = new HashSet<>();
        resizeReferenceIfNeeded(intervalToClose);

        // Break up the GVCF according to the provided reference blocking scheme
        if ( multipleAtWhichToBreakBands > 0) {
            for (int i = ((intervalToClose.getStart())/multipleAtWhichToBreakBands)*multipleAtWhichToBreakBands; i <= intervalToClose.getEnd(); i+=multipleAtWhichToBreakBands) {
                sitesToStop.add(i-1); // Subtract 1 here because we want to split before this base
            }
        }

        // If any variant contexts ended (or were spanning deletions) the last context compute where we should stop them
        for (VariantContext vc : VariantContextsOverlappingCurrentMerge) {

            // Asking if the number of alleles > 2 is a shorthand for a variant being present, as we expect <non-ref>
            // symbolic alleles to be present in all VariantContext. This might also be the case if we saw a spanning
            // deletion that reads into the current site, as we would expect ReferenceConfidenceVariantContextMerger to
            // insert symbolic alleles for those spanning variants.
            if (vc.getNAlleles() > 2) {
                for (int i = vc.getStart(); i <= vc.getEnd(); i++ ) {
                    sitesToStop.add(i);
                }
            } else if (vc.getEnd() <= intervalToClose.getEnd()) {
                sitesToStop.add(vc.getEnd());
            }
        }

        List<Integer> stoppedLocs = new ArrayList<>(sitesToStop);
        stoppedLocs.sort(Comparator.naturalOrder());

        // For each stopped loc, create a fake QueuedContextState and pass it to endPreviousStats
        for (int stoppedLoc : stoppedLocs) {
            SimpleInterval loc = new SimpleInterval(intervalToClose.getContig(), stoppedLoc, stoppedLoc);
            if (( stoppedLoc <= intervalToClose.getEnd() && stoppedLoc>= intervalToClose.getStart()) && isWithinInterval(loc)) {
                byte[] refBases = Arrays.copyOfRange(storedReferenceContext.getBases(), stoppedLoc - storedReferenceContext.getWindow().getStart(), stoppedLoc - storedReferenceContext.getWindow().getStart() + 2);
                endPreviousStates(loc, refBases, Collections.emptyList(), true);
            }
        }

    }

    /**
     * Resize {@link #storedReferenceContext} to cover at least as much as intervalToClose
     * @param intervalToClose
     */
    private void resizeReferenceIfNeeded(SimpleInterval intervalToClose) {
        final int leftEdge = storedReferenceContext.getInterval().getStart() - intervalToClose.getStart();
        final int rightEdge = intervalToClose.getEnd() - storedReferenceContext.getInterval().getEnd();

        storedReferenceContext.setWindow(Math.max(1, leftEdge), Math.max(1, rightEdge));
    }


    @Override
    public void onTraversalStart() {
        // create the annotation engine
        annotationEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(variantAnnotationArgumentCollection, dbsnp.dbsnp, Collections.emptyList());

        vcfWriter = getVCFWriter();

        referenceConfidenceVariantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine);

        //now that we have all the VCF headers, initialize the annotations (this is particularly important to turn off RankSumTest dithering in integration tests)'
        sequenceDictionary = getBestAvailableSequenceDictionary();

        // optimization to prevent mods when we always just want to break bands
        if ( multipleAtWhichToBreakBands == 1 || useBpResolution) {
            useBpResolution = true;
            multipleAtWhichToBreakBands = 1;
        }
    }

    private VariantContextWriter getVCFWriter() {
        final SortedSet<String> samples = getSamplesForVariants();

        final VCFHeader inputVCFHeader = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder(), samples);

        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions());

        // add headers for annotations added by this tool
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        VariantContextWriter writer = createVCFWriter(outputFile);

        final Set<String> sampleNameSet = new IndexedSampleList(samples).asSetOfSamples();
        final VCFHeader vcfHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));
        writer.writeHeader(vcfHeader);

        return writer;
    }

    /**
     * Method which calls endPreviousStates at the appropriate places on the given a new startingStates object
     * and an OverallState object corresponding to the currently accumulated reads.
     *
     * @param variantContexts list of variant contexts with the same start position to be reduced
     * @param referenceContext ReferenceContext object overlapping the provided VariantContexts
     */
    public void mergeWithnewVCs(final List<VariantContext> variantContexts, final ReferenceContext referenceContext) {
        if ( !variantContexts.isEmpty() ) {
            if ( ! okayToSkipThisSite(variantContexts, referenceContext) ) {
                SimpleInterval loc = referenceContext.getInterval();
                endPreviousStates( new SimpleInterval(loc.getContig(),loc.getStart()-1,loc.getStart()-1),
                        Arrays.copyOfRange(referenceContext.getBases(), 1, referenceContext.getWindow().getLengthOnReference()),
                        variantContexts,
                        false);
            }
            VariantContextsOverlappingCurrentMerge.addAll(variantContexts);
            for(final VariantContext vc : VariantContextsOverlappingCurrentMerge){
                samples.addAll(vc.getSampleNames());
            }
        }
    }

    /**
     * Is it okay to skip the given position?
     *
     * @param variantContexts  the query variant contexts representing the current position
     * @param referenceContext  Reference context object overlapping the variant contexts
     * @return true if it is okay to skip this position, false otherwise
     */
    private boolean okayToSkipThisSite(List<VariantContext> variantContexts, ReferenceContext referenceContext) {
        Set<String> intersection = new HashSet<>(getSamples(variantContexts));
        intersection.retainAll(samples);

        //if there's a starting VC with a sample that's already in a current VC, don't skip this position
        return prevPos != null && referenceContext.getInterval().getStart() == prevPos.getStart() + 1 && intersection.isEmpty(); //TODO off by 1 bug
    }

    private Set<String> getSamples(List<VariantContext> variantContexts) {
        Set<String> output = new HashSet<>();
        for (final VariantContext vc : variantContexts) {
            output.addAll(vc.getSampleNames());
        }
        return output;
    }

    /**
     * Disrupt the VariantContexts so that they all stop at the given pos, write them out, and put the remainder back in the list.
     * @param pos   the position for the starting variantContexts
     * @param variantContexts the current variant contexts with the same starting position
     * @param forceOutputAtCurrentPosition  indicates whether we output a variant at the current position, independent of VCF start/end, i.e. in BP resolution mode
     */
    private void endPreviousStates(final SimpleInterval pos, final byte[] refBases, final List<VariantContext> variantContexts, boolean forceOutputAtCurrentPosition) {
        Set<String> newSamples = getSamples(variantContexts);

        final byte refBase = refBases[0];
        //if we're in BP resolution mode or a VC ends at the current position then the reference for the next output VC (refNextBase)
        // will be advanced one base
        final byte refNextBase = (forceOutputAtCurrentPosition) ? (refBases.length > 1 ? refBases[1] : (byte)'N' ): refBase;

        final List<VariantContext> stoppedVCs = new ArrayList<>(VariantContextsOverlappingCurrentMerge.size());

        for (int i = VariantContextsOverlappingCurrentMerge.size() - 1; i >= 0; i-- ) {
            final VariantContext vc = VariantContextsOverlappingCurrentMerge.get(i);
            //the VC for the previous state will be stopped if its position is previous to the current position or it we've moved to a new contig
            if ( vc.getStart() <= pos.getStart() || !vc.contigsMatch(pos)) {

                stoppedVCs.add(vc);

                // if it was ending anyways, then remove it from the future state
                // or if ending vc is the same sample as a starting VC, then remove it from the future state
                if((vc.getEnd() == pos.getStart()) || (variantContexts.size() > 0 && !forceOutputAtCurrentPosition && newSamples.containsAll(vc.getSampleNames()))) {
                    samples.removeAll(vc.getSampleNames());
                    VariantContextsOverlappingCurrentMerge.remove(i);
                }
            }
        }

        //output the stopped variantContexts if there is no previous output (state.prevPos == null) or our current position is past
        // the last write position (state.prevPos)
        //NOTE: BP resolution with have current position == state.prevPos because it gets output via a different control flow
        if ( !stoppedVCs.isEmpty() &&  (prevPos == null || IntervalUtils.isAfter(pos,prevPos,sequenceDictionary) )) {
            final SimpleInterval closingSpot = new SimpleInterval(stoppedVCs.get(0).getContig(), pos.getStart(), pos.getStart());

            // we need the specialized merge if the site contains anything other than ref blocks
            final VariantContext mergedVC;
            if ( containsTrueAltAllele(stoppedVCs) ) {
                mergedVC = referenceConfidenceVariantContextMerger.merge(stoppedVCs, closingSpot, refBase, false, false);
            } else {
                mergedVC = referenceBlockMerge(stoppedVCs, pos.getStart());
            }

            vcfWriter.add(mergedVC);
            prevPos = closingSpot;
            refAfterPrevPos = refNextBase;
        }
    }

    /**
     * Combine a list of reference block VariantContexts.
     * We can't use GATKVariantContextUtils.simpleMerge() because it is just too slow for this sort of thing.
     *
     * @param VCs   the variant contexts to merge
     * @param end   the end of this block (inclusive)
     * @return a new merged VariantContext
     */
    private VariantContext referenceBlockMerge(final List<VariantContext> VCs, final int end) {

        final VariantContext first = VCs.get(0);

        // ref allele and start
        final Allele refAllele;
        final int start;
        if ( prevPos == null || !prevPos.getContig().equals(first.getContig()) || first.getStart() >= prevPos.getStart() + 1) {
            start = first.getStart();
            refAllele = first.getReference();
        } else {
            start = prevPos.getStart() + 1;
            refAllele = Allele.create(refAfterPrevPos, true);
        }

        // attributes
        final Map<String, Object> attrs = new HashMap<>(1);
        if ( !useBpResolution && end != start ) {
            attrs.put(VCFConstants.END_KEY, Integer.toString(end));
        }

        // genotypes
        final GenotypesContext genotypes = GenotypesContext.create();
        for ( final VariantContext vc : VCs ) {
            for ( final Genotype g : vc.getGenotypes() ) {
                genotypes.add(new GenotypeBuilder(g).alleles(GATKVariantContextUtils.noCallAlleles(g.getPloidy())).make());
            }
        }

        return new VariantContextBuilder("", first.getContig(), start, end, Arrays.asList(refAllele, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)).attributes(attrs).genotypes(genotypes).make();
    }

    /**
     * Does the given list of VariantContexts contain any with an alternate allele other than <NON_REF>?
     *
     * @param VCs  list of VariantContexts
     * @return true if there are one or more variantContexts that contain a true alternate allele, false otherwise
     */
    private boolean containsTrueAltAllele(final List<VariantContext> VCs) {

        for ( final VariantContext vc : VCs ) {
            if ( vc.getNAlleles() > 2 )
                return true;
        }
        return false;
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();

        SimpleInterval interval = prevPos != null ? new SimpleInterval(prevPos.getContig(), prevPos.getStart(), VariantContextsOverlappingCurrentMerge.stream().map(VariantContext::getEnd).max(Comparator.naturalOrder()).get()) :
                storedReferenceContext.getInterval();

        createIntermediateVariants(interval);

        // there shouldn't be any state left unless the user cut in the middle of a gVCF block
        if ( !VariantContextsOverlappingCurrentMerge.isEmpty() ) {
            logger.warn("You have asked for an interval that cuts in the middle of one or more gVCF blocks. Please note that this will cause you to lose records that don't end within your interval.");
        }

        return null;
    }

    @Override
    public void closeTool(){
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
