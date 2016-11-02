package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Bytes;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.InversionType;

@CommandLineProgramProperties(
        summary        = "Tool to genotype inversions on a single sample, based on assembled contigs from previous steps in the SV pipeline, " +
                         "and outputs genotype likelihoods.",
        oneLineSummary = "Genotype single sample inversion variant on Spark.",
        programGroup   = StructuralVariationSparkProgramGroup.class)
public final class SingleDiploidSampleBiallelicInversionGenotyperSpark extends SingleDiploidSampleBiallelicSVGenotyperSpark {
    private static final long serialVersionUID = 1L;

    static final int readLength = 151;
    private static final double expectedBaseErrorRate = 1.0/readLength; // blunt guess, will change
    // LikelihoodEngineArgumentCollection.phredScaledGlobalReadMismappingRate default value, used in HC
    private static final double maximumLikelihoodDifferenceCap = 45;

    // -----------------------------------------------------------------------------------------------
    // Overrides
    // -----------------------------------------------------------------------------------------------

    @Override
    public void runTool(final JavaSparkContext ctx){

        super.runTool(ctx);
    }

    /**
     * Custom class for holding information around two inversion breakpoints.
     */
    @VisibleForTesting
    static final class InversionJunction extends SVJunction {
        private static final long serialVersionUID = 1L;

        @VisibleForTesting
        final InversionType invType;


        InversionJunction(final VariantContext vc,
                          final ReferenceMultiSource reference,
                          final Map<Long, List<LocalAssemblyContig>> assemblyID2AlignmentRecords){

            super(vc, assemblyID2AlignmentRecords);

            final int isFiveToThree = vc.getAttributeAsBoolean(GATKSVVCFHeaderLines.INV_5_TO_3, false) ? 0b1 : 0b0;
            final int bit = isFiveToThree + (vc.getAttributeAsBoolean(GATKSVVCFHeaderLines.INV_3_TO_5, false) ? 0b10 : 0b0);
            if (bit == 0b0) {
                invType = InversionType.INV_NONE;
            } else if (bit == 0b01) {
                invType = InversionType.INV_5_TO_3;
            } else if (bit == 0b10) {
                invType = InversionType.INV_3_TO_5;
            } else {
                throw new IllegalArgumentException("Seemingly broken VCF, where the inversion breakpoint is of both type 5-to-3 and 3-to-5. Site: "
                        + vc.getContig() + ":" + vc.getID());
            }
            setWhichEnd();
            constructAlleles(reference);
        }

        @VisibleForTesting
        @Override
        protected void setWhichEnd(){
            if (invType == InversionType.INV_NONE) {
                whichEnd = 0;
            } else if (invType == InversionType.INV_5_TO_3) {
                whichEnd = 1;
            } else if (invType == InversionType.INV_3_TO_5) {
                whichEnd = 2;
            } else {
                throw new IllegalArgumentException("Seemingly broken VCF, un-recognizable inversion type, Site: "
                        + vc.getContig() + vc.getStart());
            }
        }

        /**
         * TODO: confirm that the breakpoint locations set by caller is what's assumed: BP is where the ref and alt begin to differ
         * TODO: get the short window case correct.
         */
        @VisibleForTesting
        @Override
        protected Tuple2<SimpleInterval, SimpleInterval> constructReferenceWindows(){

            final List<SimpleInterval> leftAlignedBreakpointLocations = getLeftAlignedBreakpointLocations();
            final SimpleInterval fiveEndBPLoc = leftAlignedBreakpointLocations.get(0);
            final SimpleInterval threeEndBPLoc = leftAlignedBreakpointLocations.get(1);

            final String contig = fiveEndBPLoc.getContig();

            final int flank = readLength - 1; // -1, think about it
            final int extraFlanking = (maybeNullHomology==null) ? 0 : maybeNullHomology.length;
            final int ll = fiveEndBPLoc.getStart()  - flank;
            final int lr = fiveEndBPLoc.getEnd()    + (flank-1) + (invType.equals(InversionType.INV_NONE) ? 0 : extraFlanking);
            final int rl = threeEndBPLoc.getStart() - flank - (invType.equals(InversionType.INV_NONE) ? 0 : extraFlanking);
            final int rr = threeEndBPLoc.getEnd()   + (flank-1);

            return new Tuple2<>(new SimpleInterval(contig, ll, lr), new SimpleInterval(contig, rl, rr));
        }

        //TODO: confirm that the breakpoint locations set by caller is what's assumed: BP is where the ref and alt begin to differ
        // TODO: because the assembler may decide to stop short of extending to the full windows around identified breakpoint,
        //       making use of the assembled contig to construct the alt allele is not mature yet until we can make sense out of the assembly graph
        //       so now we simply use the reference bases to construct the alt allele,
        //       what could be done intermediately is to amend the allele bases with the contigs, although this way
        //       there will be more one one alt allele for one end, because that's the reason why the assembler decides keep the bubble there
        @Override
        protected final ArrayList<SVDummyAllele> constructAlternateAlleles(final ReferenceMultiSource reference){

            try {
                final int flank = readLength - 1;
                final int hom = maybeNullHomology == null ? 0 : maybeNullHomology.length;
                final int ins = maybeNullInsertedSeq == null ? 0 : maybeNullInsertedSeq.length;

                final byte[] ll = reference.getReferenceBases(null, new SimpleInterval(fiveEndBPLoc.getContig(), fiveEndBPLoc.getStart() - flank, fiveEndBPLoc.getEnd() - 1)).getBases(); // flank bases left to 5-BP
                final byte[] rr = reference.getReferenceBases(null, new SimpleInterval(threeEndBPLoc.getContig(), threeEndBPLoc.getEnd(), threeEndBPLoc.getEnd()+flank-1)).getBases();
                if (hom!=0 && ins!=0) { // both insertion and homology

                    final byte[] lr = reference.getReferenceBases(null, new SimpleInterval(fiveEndBPLoc.getContig(), fiveEndBPLoc.getEnd()+hom, fiveEndBPLoc.getEnd()+hom+flank-1)).getBases();
                    final byte[] rl = reference.getReferenceBases(null, new SimpleInterval(threeEndBPLoc.getContig(), threeEndBPLoc.getStart()-hom-flank, threeEndBPLoc.getEnd()-hom-1)).getBases();
                    SequenceUtil.reverseComplement(rl);

                    byte[] firstHalfOfLeft;
                    byte[] firstHalfOfRight;

                    // case 3a: homology + insertion in '+' strand representation of a 5-to-3 inversion
                    firstHalfOfLeft  = Bytes.concat(ll, maybeNullHomology, maybeNullInsertedSeq);
                    final SVDummyAllele leftAltAlleleHomIns = new SVDummyAllele(Bytes.concat(firstHalfOfLeft, rl), false);
                    firstHalfOfRight = Bytes.concat(maybeNullHomology, maybeNullInsertedSeq, lr);
                    SequenceUtil.reverseComplement(firstHalfOfRight);
                    final SVDummyAllele rightAltAlleleHomIns = new SVDummyAllele(Bytes.concat(firstHalfOfRight, rr), false);

                    // case 3b: insertion + homology in '+' strand representation of a 5-to-3 inversion
                    firstHalfOfLeft  = Bytes.concat(ll, maybeNullInsertedSeq, maybeNullHomology);
                    final SVDummyAllele leftAltAlleleInsHom = new SVDummyAllele(Bytes.concat(firstHalfOfLeft, rl), false);
                    firstHalfOfRight = Bytes.concat(maybeNullInsertedSeq, maybeNullHomology, lr);
                    SequenceUtil.reverseComplement(firstHalfOfRight);
                    final SVDummyAllele rightAltAlleleInsHom = new SVDummyAllele(Bytes.concat(firstHalfOfRight, rr), false);

                    return new ArrayList<>( Arrays.asList(leftAltAlleleHomIns, rightAltAlleleHomIns, leftAltAlleleInsHom, rightAltAlleleInsHom) );
                } else {
                    byte[] tobeFlipped = reference.getReferenceBases(null, new SimpleInterval(threeEndBPLoc.getContig(), threeEndBPLoc.getStart()-flank-hom, threeEndBPLoc.getStart()-hom-1)).getBases();
                    SequenceUtil.reverseComplement(tobeFlipped);

                    final SVDummyAllele leftAltAllele;
                    if (hom==0 && ins==0) leftAltAllele = new SVDummyAllele(Bytes.concat(ll,                       tobeFlipped), false);    // nothing, very clean
                    else if(hom!=0)       leftAltAllele = new SVDummyAllele(Bytes.concat(ll, maybeNullHomology,    tobeFlipped), false);    // only homology
                    else                  leftAltAllele = new SVDummyAllele(Bytes.concat(ll, maybeNullInsertedSeq, tobeFlipped), false);    // only insertion

                    tobeFlipped = reference.getReferenceBases(null, new SimpleInterval(fiveEndBPLoc.getContig(), fiveEndBPLoc.getEnd()+hom, fiveEndBPLoc.getEnd()+hom+flank-1)).getBases();
                    if (hom!=0)      tobeFlipped = Bytes.concat(maybeNullHomology, tobeFlipped);
                    else if (ins!=0) tobeFlipped = Bytes.concat(maybeNullInsertedSeq, tobeFlipped);
                    SequenceUtil.reverseComplement(tobeFlipped);
                    final SVDummyAllele rightAltAllele = new SVDummyAllele(Bytes.concat(tobeFlipped, rr), false);

                    return new ArrayList<>(Arrays.asList(leftAltAllele, rightAltAllele));
                }
            } catch (final IOException ioex){
                throw new GATKException("Cannot resolve reference for constructing ref allele for inversion junction");
            }
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            if (!super.equals(o)) return false;

            final InversionJunction that = (InversionJunction) o;

            return invType == that.invType;

        }

        @Override
        public int hashCode() {
            int result = super.hashCode();
            result = 31 * result + invType.ordinal();
            return result;
        }
    }

    /**
     * {@inheritDoc}
     */
    @VisibleForTesting
    @Override
    protected InversionJunction convertToSVJunction(final VariantContext vc,
                                                    final Broadcast<Map<Long, List<LocalAssemblyContig>>> assembly2Alignments,
                                                    final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast){
        return new InversionJunction(vc, referenceMultiSourceBroadcast.getValue(), assembly2Alignments.getValue());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected boolean readSuitableForGenotypingJunction(final SVJunction junction, final GATKRead read){
        if(read.getAttributeAsInteger("RC")!=null && read.getAttributeAsInteger("RC")==1) return true; // these are FASTQ reconstructed reads, always use
        return readResideInJunctionWindow((InversionJunction) junction, read); // TODO: read could also share k-mer with ref or alt
    }

    /**
     * {@inheritDoc}
     *
     * Normalizes the read likelihoods, see {@link ReadLikelihoods#normalizeLikelihoods(boolean, double)}.
     * filter out poorly modeled reads, see {@link ReadLikelihoods#filterPoorlyModeledReads(double)},
     * and
     * marginalizes the read likelihoods by reducing from multiple ref alleles and multiple alt alleles
     * used in the RLC (read likelihood calculation) step to a single symbolic ref allele and
     * a single symbolic alt allele (see {@link ReadLikelihoods#marginalize(Map)} for logic)
     * so that the following assumption is made:
     * (P=2, A=2) =&gt; [0/0, 0/1, 1/1] three possible genotypes for a diploid biallelic sample.
     * TODO: is this marginalization logic correct?
     *
     * TODO: Re-alignment of reads back to their best allele is not implemented yet (as it happens in HC).
     */
    @VisibleForTesting
    @Override
    protected ReadLikelihoods<SVDummyAllele> updateReads(final SVJunction inversionJunction, final ReadLikelihoods<SVDummyAllele> matrix){

        matrix.normalizeLikelihoods(false, QualityUtils.qualToErrorProbLog10(maximumLikelihoodDifferenceCap));
        matrix.filterPoorlyModeledReads(expectedBaseErrorRate);

        final Map<SVDummyAllele, List<SVDummyAllele>> symbolicOnes2RealAlleles = new HashMap<>();
        final List<SVDummyAllele> alleles = inversionJunction.getAlleles();
        final List<SVDummyAllele> vcAlleles = inversionJunction.getOriginalVC().getAlleles().stream().map(SVDummyAllele::new).collect(Collectors.toList());
        symbolicOnes2RealAlleles.put(vcAlleles.get(0), alleles.stream().filter(SVDummyAllele::isReference).collect(Collectors.toList()));
        symbolicOnes2RealAlleles.put(vcAlleles.get(1), alleles.stream().filter(SVDummyAllele::isNonReference).collect(Collectors.toList()));

        final ReadLikelihoods<SVDummyAllele> result = matrix.marginalize(symbolicOnes2RealAlleles);
        if(SingleDiploidSampleBiallelicSVGenotyperSpark.in_debug_state){
            final StringBuilder builder = new StringBuilder();
            builder.append("UPDATED_RLL_MATRIX\n").append(result.sampleMatrix(0).toString());
            inversionJunction.debugString += builder.toString();
        }
        return result;
    }

    // -----------------------------------------------------------------------------------------------
    // Utilities
    // -----------------------------------------------------------------------------------------------

    /**
     * Test if a particular read starts after or ends before the two windows spanned by the inversion junction.
     */
    @VisibleForTesting
    boolean readResideInJunctionWindow(final InversionJunction invJunction, final GATKRead read){

        if(read.isUnmapped()) return false; // TODO: really stupid!

        final Tuple2<SimpleInterval, SimpleInterval> windows = invJunction.getReferenceWindows();
        return windows._1().contains(read) || windows._2().contains(read);
    }
}
