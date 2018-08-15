package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.tools.walkers.annotator.ChromosomeCounts;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;

import java.io.File;
import java.nio.file.Path;
import java.util.*;

/**
 * Left-align indels in a variant callset
 *
 * <p>
 * This tool takes a VCF file, left-aligns the indels and trims common bases from indels,
 * leaving them with a minimum representation. The same indel can often be placed at multiple positions and still
 * represent the same haplotype. While the standard convention with VCF is to place an indel at the left-most position
 * this isn't always done, so this tool can be used to left-align them. This tool optionally splits multiallelic
 * sites into biallelics and left-aligns individual alleles. Optionally, the tool will not trim common bases from indels.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * A variant call set to left-align and trim.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A left-aligned VCF.
 * </p>
 *
 * <h3>Usage examples</h3>
 *
 * <h4>Left align and trim alleles</h4>
 * <pre>
 * gatk LeftAlignAndTrimVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -O output.vcf
 * </pre>
 *
 * <h4>Left align and don't trim alleles</h4>
 * <pre>
 * gatk LeftAlignAndTrimVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -O output.vcf \
 *   --dontTrimAlleles
 * </pre>
 *
 * <h4>Left align and trim alleles, process alleles <= 208 bases</h4>
 * <pre>
 * gatk LeftAlignAndTrimVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -O output.vcf \
 *   --maxIndelSize 208
 * </pre>
 *  * <h4>Left align and trim alleles, look back in reference at 208 bases before allele</h4>
 *  * <pre>
 *  * gatk LeftAlignAndTrimVariants \
 *  *   -R reference.fasta \
 *  *   -V input.vcf \
 *  *   -O output.vcf \
 *  *   --referenceWindow 208
 *  * </pre>
 *
 * <h4>Split multiallics into biallelics, left align and trim alleles</h4>
 * <pre>
 * gatk LeftAlignAndTrimVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -O output.vcf \
 *   --splitMultiallelics
 * </pre>
 *
 * <h4>Split multiallelics into biallics, left align but don't trim alleles, and store the original AC, AF, and AN values</h4>
 * <pre>
 * gatk LeftAlignAndTrimVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -O output.vcf \
 *   --splitMultiallelics \
 *   --dontTrimAlleles
 *   --keepOriginalAC
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "This tool takes a VCF file, left-aligns the indels and trims common bases from indels," +
            "leaving them with a minimum representation. The same indel can often be placed at multiple positions and still" +
            "represent the same haplotype. While the standard convention with VCF is to place an indel at the left-most position" +
            "this isn't always done, so this tool can be used to left-align them. This tool optionally splits multiallelic" +
            "sites into biallelics and left-aligns individual alleles. Optionally, the tool will not trim common bases from indels.",
            oneLineSummary = "Left align and trim vairants",
            programGroup = VariantManipulationProgramGroup.class
)
@DocumentedFeature
public class LeftAlignAndTrimVariants extends VariantWalker {

    /**
     * Output file to which to write left aligned variants
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written")
    public File outFile = null;
    /**
     * If this argument is set, bases common to all alleles will not be removed and will not leave their minimal representation.
     */
    @Argument(fullName="dontTrimAlleles", shortName="notrim", doc="Do not Trim alleles to remove bases common to all of them", optional=true)
    protected boolean dontTrimAlleles = false;

    /**
     * If this argument is set, split multiallelic records and left-align individual alleles.
     * If this argument is not set, multiallelic records are not attempted to left-align and will be copied as is.
     */
    @Argument(fullName="splitMultiallelics", shortName="split", doc="Split multiallelic records and left-align individual alleles", optional=true)
    protected boolean splitMultiallelics = false;

    /**
     * When subsetting a callset, this tool recalculates the AC, AF, and AN values corresponding to the contents of the
     * subset. If this flag is enabled, the original values of those annotations will be stored in new annotations called
     * AC_Orig, AF_Orig, and AN_Orig.
     */
    @Argument(fullName="keepOriginalAC", shortName="keepOriginalAC", doc="Store the original AC, AF, and AN values after subsetting", optional=true)
    private boolean keepOriginalChrCounts = false;

    /**
     * Maximum indel size to realign.  Indels larger than this will be left unadjusted.
     */
    @Argument(fullName="maxIndelSize",shortName="maxIndelSize",doc="Set maximum indel size to realign",optional=true)
    private static int maxIndelSize=200;

    /**
     * Distance in reference to look back before allele
     */
    @Argument(fullName="maxLeadingBases",shortName = "maxLeadingBases",doc="Set max reference window size to look back before allele",optional=true)
    private static int maxLeadingBases=1000;

    @Hidden
    @Argument(fullName="suppress-reference-path", optional=true,
            doc="Suppress reference path in output for test result differencing")
    private boolean suppressReferencePath = false;

    private VariantContextWriter vcfWriter = null;

    private int numRealignedVariants;

    /**
     * Set up the VCF writer, samples
     *
     */
    @Override
    public void onTraversalStart() {
        numRealignedVariants=0;
        final Map<String, VCFHeader> vcfHeaders = Collections.singletonMap(getDrivingVariantsFeatureInput().getName(), getHeaderForVariants());
        final SortedSet<String> vcfSamples = VcfUtils.getSortedSampleSet(vcfHeaders, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        // Initialize VCF header lines
        final Set<VCFHeaderLine> headerLines = createVCFHeaderLineList(vcfHeaders);

        Set<VCFHeaderLine> actualLines = null;
        SAMSequenceDictionary sequenceDictionary = null;
        if (hasReference()) {
            Path refPath = referenceArguments.getReferencePath();
            sequenceDictionary= this.getReferenceDictionary();
            actualLines = VcfUtils.updateHeaderContigLines(headerLines, refPath, sequenceDictionary, suppressReferencePath);
        }
        else {
            sequenceDictionary = getHeaderForVariants().getSequenceDictionary();
            if (null != sequenceDictionary) {
                actualLines = VcfUtils.updateHeaderContigLines(headerLines, null, sequenceDictionary, suppressReferencePath);
            }
            else {
                actualLines = headerLines;
            }
        }

        vcfWriter = createVCFWriter(outFile);
        vcfWriter.writeHeader(new VCFHeader(actualLines, vcfSamples));
    }

    /**
     * Prepare the VCF header lines
     */
    private Set<VCFHeaderLine> createVCFHeaderLineList(Map<String, VCFHeader> vcfHeaders) {

        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfHeaders.values(), true);
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        if (splitMultiallelics) {
            GATKVariantContextUtils.addChromosomeCountsToHeader(headerLines);
        }

        if (keepOriginalChrCounts) {
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AC_KEY));
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AF_KEY));
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AN_KEY));
        }
        return headerLines;
    }


    /**
     * 
     * @param vc
     * @param readsContext Reads overlapping the current variant. Will be an empty, but non-null, context object
     *                     if there is no backing source of reads data (in which case all queries on it will return
     *                     an empty array/iterator)
     * @param ref
     * @param featureContext Features spanning the current variant. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     */
    @Override
    public void apply(VariantContext vc, ReadsContext readsContext, ReferenceContext ref, FeatureContext featureContext) {
        final List<VariantContext> vcList= splitMultiallelics ? GATKVariantContextUtils.splitVariantContextToBiallelics(vc, false,
                GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, keepOriginalChrCounts) : Collections.singletonList(vc);

        for (final VariantContext v: vcList) {
            numRealignedVariants +=trimAlign(v,ref);
        }
    }

    /**
     * Print out message of how many variants were realigned
     * @return
     */
    @Override
    public Object onTraversalSuccess() {
        String msg=numRealignedVariants + (numRealignedVariants==1 ? " variant " : " variants ") + "left aligned";
        return msg;
    }

    /**
     * Close out the new variants file.
     */
    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }

    /**
     * Trim, align and write out the vc.
     *
     * @param vc                Input VC with variants to left align
     * @param ref               Reference context
     * @return                  Number of records left-aligned (0 or 1)
     */
    protected int trimAlign(final VariantContext vc, final ReferenceContext ref){

        final int refLength =  vc.getReference().length();

        // ignore if the reference length is greater than the maxIndelSize
        if ( refLength > maxIndelSize ) {
            logger.info(String.format("%s (%d) at position %s:%d; skipping that record. Set --maxIndelSize >= %d",
                    "Reference allele is too long", refLength, vc.getContig(), vc.getStart(), refLength));
            return 0;
        }

        // optionally don't trim VC
        final VariantContext v = dontTrimAlleles ? vc : GATKVariantContextUtils.trimAlleles(vc, true, true);

        // align the VC
        final VariantContext alignedV = leftAlign(v, ref);

        // write out new VC
        vcfWriter.add(alignedV);

        // number of records left aligned
        return (v.equals(alignedV)? 0 : 1);
    }

    /**
     * Main routine workhorse. By definition, it will only take biallelic vc's. Splitting into multiple alleles has to be
     * handled by calling routine.
     * @param vc                  Input VC with variants to left align
     * @param ref                 Reference context
     * @return                    new VC.
     */

    protected VariantContext leftAlign (final VariantContext vc, final ReferenceContext ref) {
        if (!vc.isIndel() || vc.isComplexIndel() ) {
            return vc;
        }

        // get the indel length
        final int indelLength;
        if ( vc.isSimpleDeletion() )
            indelLength = vc.getReference().length() - 1;
        else
            indelLength = vc.getAlternateAllele(0).length() - 1;

        int leadingBases=50;

        while (leadingBases < maxLeadingBases) {
            ref.setWindow(leadingBases,maxIndelSize);
            final byte[] refSeq = ref.getBases();

            // create an indel haplotype.
            //
            final int originalIndex = vc.getStart() - ref.getWindow().getStart() + 1;
            if (originalIndex < 0 || originalIndex >= ref.getBases().length)
                return vc;

            final byte[] originalIndel = makeHaplotype(vc, refSeq, originalIndex, indelLength);
            // create a CIGAR string to represent the event
            ArrayList<CigarElement> elements = new ArrayList<CigarElement>();
            elements.add(new CigarElement(originalIndex, CigarOperator.M));
            elements.add(new CigarElement(indelLength, vc.isSimpleDeletion() ? CigarOperator.D : CigarOperator.I));
            elements.add(new CigarElement(refSeq.length - originalIndex, CigarOperator.M));
            Cigar originalCigar = new Cigar(elements);

            // left align the CIGAR
            Cigar newCigar = AlignmentUtils.leftAlignIndel(originalCigar, refSeq, originalIndel, 0, 0, true);
            if (newCigar.equals(originalCigar)) {
                return vc;
            }
            if (newCigar.getCigarElement(0).getOperator()==CigarOperator.M) {
                int difference = originalIndex - newCigar.getCigarElement(0).getLength();
                VariantContext newVC = new VariantContextBuilder(vc).start(vc.getStart()-difference).stop(vc.getEnd()-difference).make();

                final int indelIndex = originalIndex-difference;
                final byte[] newBases = new byte[indelLength + 1];
                newBases[0] = refSeq[indelIndex-1];
                System.arraycopy((vc.isSimpleDeletion() ? refSeq : originalIndel), indelIndex, newBases, 1, indelLength);
                final Allele newAllele = Allele.create(newBases, vc.isSimpleDeletion());
                newVC = updateAllele(newVC, newAllele);
                // overwrite default return value with new left-aligned VC
                return newVC;
            }
            // if we have left aligned to the beginning of the reference window, double leadingBases length and try again
            leadingBases *= 2;
        }
        //if we have made it here we have failed to align the indel before reaching the maximum number of leading Bases
        logger.info("Indel left aligned to beginning of reference window.  Please set --maxLeadingBases to greater than " + maxLeadingBases +
            "and try again if you would like to align this indel");
        return vc;
    }

    protected VariantContext updateAllele(final VariantContext vc, final Allele newAllele) {
        // create a mapping from original allele to new allele
        HashMap<Allele, Allele> alleleMap = new HashMap<Allele, Allele>(vc.getAlleles().size());
        if ( newAllele.isReference() ) {
            alleleMap.put(vc.getReference(), newAllele);
            alleleMap.put(vc.getAlternateAllele(0), Allele.create(newAllele.getBases()[0], false));
        } else {
            alleleMap.put(vc.getReference(), Allele.create(newAllele.getBases()[0], true));
            alleleMap.put(vc.getAlternateAllele(0), newAllele);
        }

        // create new Genotype objects
        GenotypesContext newGenotypes = GenotypesContext.create(vc.getNSamples());
        for ( final Genotype genotype : vc.getGenotypes() ) {
            List<Allele> newAlleles = new ArrayList<Allele>();
            for ( Allele allele : genotype.getAlleles() ) {
                Allele newA = alleleMap.get(allele);
                if ( newA == null )
                    newA = Allele.NO_CALL;
                newAlleles.add(newA);
            }
            newGenotypes.add(new GenotypeBuilder(genotype).alleles(newAlleles).make());
        }

        return new VariantContextBuilder(vc).alleles(alleleMap.values()).genotypes(newGenotypes).make();
    }

    /**
     * Make a haplotype from a given alt allele, using bases in input reference, index of an input reference
     * @param vc                                Input VC - will use only alt allele from it
     * @param ref                               Ref bases
     * @param indexOfRef                        Index in ref where to create indel
     * @param indelLength                       Indel length
     * @return                                  Haplotype, containing the reference and the indel
     */
    private static byte[] makeHaplotype(VariantContext vc, byte[] ref, int indexOfRef, int indelLength) {
        byte[] hap = new byte[ref.length + (indelLength * (vc.isSimpleDeletion() ? -1 : 1))];

        // add the bases before the indel
        System.arraycopy(ref, 0, hap, 0, indexOfRef);
        int currentPos = indexOfRef;

        // take care of the indel
        if ( vc.isSimpleDeletion() ) {
            indexOfRef += indelLength;
        } else {
            System.arraycopy(vc.getAlternateAllele(0).getBases(), 1, hap, currentPos, indelLength);
            currentPos += indelLength;
        }

        // add the bases after the indel
        System.arraycopy(ref, indexOfRef, hap, currentPos, ref.length - indexOfRef);

        return hap;
    }

}
