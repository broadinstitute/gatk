/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2015 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.hellbender.engine.writers;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

/**
 * Genome-wide VCF writer
 *
 * User: depristo
 * Date: 6/24/13
 * Time: 2:51 PM
 */
public class GVCFWriter implements VariantContextWriter {

    //
    // Final fields initialized in constructor
    //
    /** Where we'll ultimately write our VCF records */
    final private VariantContextWriter underlyingWriter;

    final private List<HomRefBlock> GQPartitions;

    /** fields updated on the fly during GVCFWriter operation */
    int nextAvailableStart = -1;
    String contigOfNextAvailableStart = null;
    private String sampleName = null;
    private HomRefBlock currentBlock = null;
    private final int defaultPloidy;

    /**
     * Is the proposed GQ partitions well-formed?
     *
     * @param GQPartitions proposed GQ partitions
     * @return a non-null string if something is wrong (string explains issue)
     */
    protected static List<HomRefBlock> parsePartitions(final List<Integer> GQPartitions, final int defaultPloidy) {
        if ( GQPartitions == null ) throw new IllegalArgumentException("GQpartitions cannot be null");
        if ( GQPartitions.isEmpty() ) throw new IllegalArgumentException("GQpartitions cannot be empty");

        final List<HomRefBlock> result = new LinkedList<>();
        int lastThreshold = 0;
        for ( final Integer value : GQPartitions ) {
            if ( value == null ) throw new IllegalArgumentException("GQPartitions contains a null integer");
            if ( value < lastThreshold ) throw new IllegalArgumentException("GQPartitions is out of order.  Last is " + lastThreshold + " but next is " + value);
            if ( value == lastThreshold ) throw new IllegalArgumentException("GQPartitions is equal elements: Last is " + lastThreshold + " but next is " + value);
            result.add(new HomRefBlock(lastThreshold, value,defaultPloidy));
            lastThreshold = value;
        }
        result.add(new HomRefBlock(lastThreshold, Integer.MAX_VALUE,defaultPloidy));

        return result;
    }

    /**
     * Create a new GVCF writer
     *
     * Should be a non-empty list of boundaries.  For example, suppose this variable is
     *
     * [A, B, C]
     *
     * We would partition our hom-ref sites into the following bands:
     *
     * X < A
     * A <= X < B
     * B <= X < C
     * X >= C
     *
     * @param underlyingWriter the ultimate destination of the GVCF records
     * @param GQPartitions a well-formed list of GQ partitions
     * @param defaultPloidy the assumed ploidy for input variant context without one.
     */
    public GVCFWriter(final VariantContextWriter underlyingWriter, final List<Integer> GQPartitions, final int defaultPloidy) {
        if ( underlyingWriter == null ) throw new IllegalArgumentException("underlyingWriter cannot be null");
        this.underlyingWriter = underlyingWriter;
        this.GQPartitions = parsePartitions(GQPartitions,defaultPloidy);
        this.defaultPloidy = defaultPloidy;
    }

    /**
     * Write the VCF header
     *
     * Adds standard GVCF fields to the header
     *
     * @param header a non-null header
     */
    @Override
    public void writeHeader(VCFHeader header) {
        if ( header == null ) throw new IllegalArgumentException("header cannot be null");
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        header.addMetaDataLine(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.MIN_DP_FORMAT_KEY));

        for ( final HomRefBlock partition : GQPartitions ) {
            header.addMetaDataLine(partition.toVCFHeaderLine());
        }

        underlyingWriter.writeHeader(header);
    }

    /**
     * Close this GVCF writer.  Finalizes any pending hom-ref blocks and emits those to the underlyingWriter as well
     */
    @Override
    public void close() {
        close(true);
    }

    /**
     * Horrible work around because there's no clean way to get our VCFWriter closed by the GATK
     *
     * If closeUnderlyingWriter is true, then we'll close the underlying writer, otherwise we'll leave it open
     * so the GATK closes it later
     *
     * @param closeUnderlyingWriter should we leave the underlying writer open or closed?
     */
    public void close(final boolean closeUnderlyingWriter) {
        emitCurrentBlock();
        if ( closeUnderlyingWriter ) underlyingWriter.close();
    }

    /**
     * Add hom-ref site from vc to this gVCF hom-ref state tracking, emitting any pending states if appropriate
     *
     * @param vc a non-null VariantContext
     * @param g a non-null genotype from VariantContext
     * @return a VariantContext to be emitted, or null if non is appropriate
     */
    protected VariantContext addHomRefSite(final VariantContext vc, final Genotype g) {

        if ( nextAvailableStart != -1 ) {
            // don't create blocks while the hom-ref site falls before nextAvailableStart (for deletions)
            if ( vc.getStart() <= nextAvailableStart && vc.getChr().equals(contigOfNextAvailableStart) )
                return null;
            // otherwise, reset to non-relevant
            nextAvailableStart = -1;
            contigOfNextAvailableStart = null;
        }

        final VariantContext result;
        if (genotypeCanBeMergedInCurrentBlock(g)) {
            currentBlock.add(vc.getStart(), g);
            result = null;
        } else {
            result = blockToVCF(currentBlock);
            currentBlock = createNewBlock(vc, g);
        }
        return result;
    }

    private boolean genotypeCanBeMergedInCurrentBlock(final Genotype g) {
        return currentBlock != null && currentBlock.withinBounds(g.getGQ()) && currentBlock.getPloidy() == g.getPloidy()
                && (currentBlock.getMinPLs() == null || !g.hasPL() || (currentBlock.getMinPLs().length == g.getPL().length));
    }

    /**
     * Flush the current hom-ref block, if necessary, to the underlying writer, and reset the currentBlock to null
     */
    private void emitCurrentBlock() {
        if ( currentBlock != null ) {
            // there's actually some work to do
            underlyingWriter.add(blockToVCF(currentBlock));
            currentBlock = null;
        }
    }

    /**
     * Convert a HomRefBlock into a VariantContext
     *
     * @param block the block to convert
     * @return a VariantContext representing the gVCF encoding for this block.
     * It will return {@code null} if input {@code block} is {@code null}, indicating that there
     * is no variant-context to be output into the VCF.
     */
    private VariantContext blockToVCF(final HomRefBlock block) {
        if ( block == null ) return null;

        final VariantContextBuilder vcb = new VariantContextBuilder(block.getStartingVC());
        vcb.attributes(new HashMap<String, Object>(2)); // clear the attributes
        vcb.stop(block.getStop());
        vcb.attribute(VCFConstants.END_KEY, block.getStop());

        // create the single Genotype with GQ and DP annotations
        final GenotypeBuilder gb = new GenotypeBuilder(sampleName, GATKVariantContextUtils.homozygousAlleleList(block.getRef(), block.getPloidy()));
        gb.noAD().noPL().noAttributes(); // clear all attributes

        final int[] minPLs = block.getMinPLs();
        gb.PL(minPLs);
        final int gq = genotypeQualityFromPLs(minPLs);
        gb.GQ(gq);
        gb.DP(block.getMedianDP());
        gb.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, block.getMinDP());

        // This annotation is no longer standard
        //gb.attribute(MIN_GQ_FORMAT_FIELD, block.getMinGQ());

        return vcb.genotypes(gb.make()).make();
    }


    private int genotypeQualityFromPLs(final int[] minPLs) {
        int first = minPLs[0];
        int second = minPLs[1];
        if (first > second) {
            second = first;
            first = minPLs[1];
        }
        for (int i = 3; i < minPLs.length; i++) {
            final int candidate = minPLs[i];
            if (candidate >= second) continue;
            if (candidate <= first) {
                second = first;
                first = candidate;
            } else
                second = candidate;
        }
        return second - first;
    }

    /**
     * Helper function to create a new HomRefBlock from a variant context and current genotype
     *
     * @param vc the VariantContext at the site where want to start the band
     * @param g the genotype of the sample from vc that should be used to initialize the block
     * @return a newly allocated and initialized block containing g already
     */
    private HomRefBlock createNewBlock(final VariantContext vc, final Genotype g) {
        // figure out the GQ limits to use based on the GQ of g
        HomRefBlock partition = null;
        for ( final HomRefBlock maybePartition : GQPartitions ) {
            if ( maybePartition.withinBounds(g.getGQ()) ) {
                partition = maybePartition;
                break;
            }
        }

        if ( partition == null )
            throw new IllegalStateException("GQ " + g + " from " + vc + " didn't fit into any partition");

        // create the block, add g to it, and return it for use
        final HomRefBlock block = new HomRefBlock(vc, partition.getGQLowerBound(), partition.getGQUpperBound(), defaultPloidy);
        block.add(vc.getStart(), g);
        return block;
    }

    /**
     * Add a VariantContext to this writer for emission
     *
     * Requires that the VC have exactly one genotype
     *
     * @param vc a non-null VariantContext
     */
    @Override
    public void add(VariantContext vc) {
        if ( vc == null ) throw new IllegalArgumentException("vc cannot be null");

        if ( sampleName == null )
            sampleName = vc.getGenotype(0).getSampleName();

        if ( ! vc.hasGenotypes() ) {
            throw new IllegalArgumentException("GVCF assumes that the VariantContext has genotypes");
        } else if ( vc.getGenotypes().size() != 1 ) {
            throw new IllegalArgumentException("GVCF assumes that the VariantContext has exactly one genotype but saw " + vc.getGenotypes().size());
        } else {
            if ( currentBlock != null && ! currentBlock.isContiguous(vc) ) {
                // we've made a non-contiguous step (across interval, onto another chr), so finalize
                emitCurrentBlock();
            }

            final Genotype g = vc.getGenotype(0);
            if ( g.isHomRef() && vc.hasAlternateAllele(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) && vc.isBiallelic() ) {
                // create bands
                final VariantContext maybeCompletedBand = addHomRefSite(vc, g);
                if ( maybeCompletedBand != null ) underlyingWriter.add(maybeCompletedBand);
            } else {
                // g is variant, so flush the bands and emit vc
                emitCurrentBlock();
                nextAvailableStart = vc.getEnd();
                contigOfNextAvailableStart = vc.getChr();
                underlyingWriter.add(vc);
            }
        }
    }
}
