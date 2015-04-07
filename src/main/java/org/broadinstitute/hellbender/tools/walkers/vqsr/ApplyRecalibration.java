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
* Copyright 2012-2014 Broad Institute, Inc.
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

package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.engine.walkers.PartitionBy;
import org.broadinstitute.gatk.engine.walkers.PartitionType;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.io.File;
import java.util.*;

/**
 * Apply a score cutoff to filter variants based on a recalibration table
 *
 * <p>
 * This tool performs the second pass in a two-stage process called VQSR; the first pass is performed by the
 * <a href='https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php'>VariantRecalibrator</a> tool.
 * In brief, the first pass consists of creating a Gaussian mixture model by looking at the distribution of annotation
 * values over a high quality subset of the input call set, and then scoring all input variants according to the model.
 * The second pass consists of filtering variants based on score cutoffs identified in the first pass.
 *</p>
 *
 * <p>
 * Using the tranche file and recalibration table generated by the previous step, the ApplyRecalibration tool looks at each variant's VQSLOD value
 * and decides which tranche it falls in. Variants in tranches that fall below the specified truth sensitivity filter level
 * have their FILTER field annotated with the corresponding tranche level. This will result in a call set that is filtered
 * to the desired level but retains the information necessary to increase sensitivity if needed.</p>
 *
 * <p>To be clear, please note that by "filtered", we mean that variants failing the requested tranche cutoff are <b>marked
 * as filtered</b> in the output VCF; they are <b>not discarded</b>.</p>
 *
 * <p>VQSR is probably the hardest part of the Best Practices to get right, so be sure to read the
 * <a href='https://www.broadinstitute.org/gatk/guide/article?id=39'>method documentation</a>,
 * <a href='https://www.broadinstitute.org/gatk/guide/article?id=1259'>parameter recommendations</a> and
 * <a href='https://www.broadinstitute.org/gatk/guide/article?id=2805'>tutorial</a> to really understand what these
 * tools and how to use them for best results on your own data.</p>
 *
 * <h3>Input</h3>
 * <ul>
 * <li>The raw input variants to be filtered.</li>
 * <li>The recalibration table file that was generated by the VariantRecalibrator tool.</li>
 * <li>The tranches file that was generated by the VariantRecalibrator tool.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 * <li>A recalibrated VCF file in which each variant of the requested type is annotated with its VQSLOD and marked as filtered if the score is below the desired quality level.</li>
 * </ul>
 *
 * <h3>Example for filtering SNPs</h3>
 * <pre>
 * java -Xmx3g -jar GenomeAnalysisTK.jar \
 *   -T ApplyRecalibration \
 *   -R reference/human_g1k_v37.fasta \
 *   -input NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf \
 *   --ts_filter_level 99.0 \
 *   -tranchesFile path/to/output.tranches \
 *   -recalFile path/to/output.recal \
 *   -mode SNP \
 *   -o path/to/output.recalibrated.filtered.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 *
 * <ul>
 * <li>The tranche values used in the example above is only a general example. You should determine the level of sensitivity
 * that is appropriate for your specific project. Remember that higher sensitivity (more power to detect variants, yay!) comes
 * at the cost of specificity (more false negatives, boo!). You have to choose at what point you want to set the tradeoff.</li>
 * <li>In order to create the tranche reporting plots (which are only generated for SNPs, not indels!) Rscript needs to be
 * in your environment PATH (this is the scripting version of R, not the interactive version).
 * See <a target="r-project" href="http://www.r-project.org">http://www.r-project.org</a> for more info on how to download and install R.</li>
 * </ul>
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class} )
@PartitionBy(PartitionType.LOCUS)
public class ApplyRecalibration extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {

    public static final String LOW_VQSLOD_FILTER_NAME = "LOW_VQSLOD";
    private final double DEFAULT_VQSLOD_CUTOFF = 0.0;

    /////////////////////////////
    // Inputs
    /////////////////////////////
    /**
     * These calls should be unfiltered and annotated with the error covariates that are intended to use for modeling.
     */
    @Input(fullName="input", shortName = "input", doc="The raw input variants to be recalibrated", required=true)
    public List<RodBinding<VariantContext>> input;
    @Input(fullName="recal_file", shortName="recalFile", doc="The input recal file used by ApplyRecalibration", required=true)
    protected RodBinding<VariantContext> recal;
    @Input(fullName="tranches_file", shortName="tranchesFile", doc="The input tranches file describing where to cut the data", required=false)
    protected File TRANCHES_FILE;

    /////////////////////////////
    // Outputs
    /////////////////////////////
    @Output( doc="The output filtered and recalibrated VCF file in which each variant is annotated with its VQSLOD value")
    private VariantContextWriter vcfWriter = null;

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName="ts_filter_level", shortName="ts_filter_level", doc="The truth sensitivity level at which to start filtering", required=false)
    protected Double TS_FILTER_LEVEL = null;
    @Advanced
    @Argument(fullName="lodCutoff", shortName="lodCutoff", doc="The VQSLOD score below which to start filtering", required=false)
    protected Double VQSLOD_CUTOFF = null;

    /**
     * For this to work properly, the -ignoreFilter argument should also be applied to the VariantRecalibration command.
     */
    @Argument(fullName="ignore_filter", shortName="ignoreFilter", doc="If specified, the recalibration will be applied to variants marked as filtered by the specified filter name in the input VCF file", required=false)
    private String[] IGNORE_INPUT_FILTERS = null;
    @Argument(fullName="ignore_all_filters", shortName="ignoreAllFilters", doc="If specified, the variant recalibrator will ignore all input filters. Useful to rerun the VQSR from a filtered output file.", required=false)
    private boolean IGNORE_ALL_FILTERS = false;
    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't output filtered loci after applying the recalibration", required=false)
    protected boolean EXCLUDE_FILTERED = false;
    @Argument(fullName = "mode", shortName = "mode", doc = "Recalibration mode to employ: 1.) SNP for recalibrating only SNPs (emitting indels untouched in the output VCF); 2.) INDEL for indels; and 3.) BOTH for recalibrating both SNPs and indels simultaneously.", required = false)
    public VariantRecalibratorArgumentCollection.Mode MODE = VariantRecalibratorArgumentCollection.Mode.SNP;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    final private List<Tranche> tranches = new ArrayList<>();
    final private Set<String> inputNames = new HashSet<>();
    final private Set<String> ignoreInputFilterSet = new TreeSet<>();

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        if( TS_FILTER_LEVEL != null ) {
            for ( final Tranche t : Tranche.readTranches(TRANCHES_FILE) ) {
                if ( t.ts >= TS_FILTER_LEVEL ) {
                    tranches.add(t);
                }
                logger.info(String.format("Read tranche " + t));
            }
            Collections.reverse(tranches); // this algorithm wants the tranches ordered from best (lowest truth sensitivity) to worst (highest truth sensitivity)
        }

        for( final RodBinding rod : input ) {
            inputNames.add( rod.getName() );
        }

        if( IGNORE_INPUT_FILTERS != null ) {
            ignoreInputFilterSet.addAll( Arrays.asList(IGNORE_INPUT_FILTERS) );
        }

        // setup the header fields
        final Set<VCFHeaderLine> hInfo = new HashSet<>();
        hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit(), inputNames));
        addVQSRStandardHeaderLines(hInfo);
        final TreeSet<String> samples = new TreeSet<>();
        samples.addAll(SampleUtils.getUniqueSamplesFromRods(getToolkit(), inputNames));

        if( TS_FILTER_LEVEL != null ) {
            // if the user specifies both ts_filter_level and lodCutoff then throw a user error
            if( VQSLOD_CUTOFF != null ) {
                throw new UserException("Arguments --ts_filter_level and --lodCutoff are mutually exclusive. Please only specify one option.");
            }

            if( tranches.size() >= 2 ) {
                for( int iii = 0; iii < tranches.size() - 1; iii++ ) {
                    final Tranche t = tranches.get(iii);
                    hInfo.add(new VCFFilterHeaderLine(t.name, String.format("Truth sensitivity tranche level for " + t.model.toString() + " model at VQS Lod: " + t.minVQSLod + " <= x < " + tranches.get(iii + 1).minVQSLod)));
                }
            }
            if( tranches.size() >= 1 ) {
                hInfo.add(new VCFFilterHeaderLine(tranches.get(0).name + "+", String.format("Truth sensitivity tranche level for " + tranches.get(0).model.toString() + " model at VQS Lod < " + tranches.get(0).minVQSLod)));
            } else {
                throw new UserException("No tranches were found in the file or were above the truth sensitivity filter level " + TS_FILTER_LEVEL);
            }

            logger.info("Keeping all variants in tranche " + tranches.get(tranches.size()-1));
        } else {
            if( VQSLOD_CUTOFF == null ) {
                VQSLOD_CUTOFF = DEFAULT_VQSLOD_CUTOFF;
            }
            hInfo.add(new VCFFilterHeaderLine(LOW_VQSLOD_FILTER_NAME, "VQSLOD < " + VQSLOD_CUTOFF));
            logger.info("Keeping all variants with VQSLOD >= " + VQSLOD_CUTOFF);
        }

        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    public static void addVQSRStandardHeaderLines(final Set<VCFHeaderLine> hInfo) {
        hInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VQS_LOD_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.CULPRIT_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.POSITIVE_LABEL_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.NEGATIVE_LABEL_KEY));
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        if( tracker == null ) { // For some reason RodWalkers get map calls with null trackers
            return 1;
        }

        final List<VariantContext> VCs =  tracker.getValues(input, context.getLocation());
        final List<VariantContext> recals =  tracker.getValues(recal, context.getLocation());

        for( final VariantContext vc : VCs ) {

            if( VariantDataManager.checkVariationClass( vc, MODE ) && (IGNORE_ALL_FILTERS || vc.isNotFiltered() || ignoreInputFilterSet.containsAll(vc.getFilters())) ) {

                final VariantContext recalDatum = getMatchingRecalVC(vc, recals);
                if( recalDatum == null ) {
                    throw new UserException("Encountered input variant which isn't found in the input recal file. Please make sure VariantRecalibrator and ApplyRecalibration were run on the same set of input variants. First seen at: " + vc );
                }

                final String lodString = recalDatum.getAttributeAsString(GATKVCFConstants.VQS_LOD_KEY, null);
                if( lodString == null ) {
                    throw new UserException("Encountered a malformed record in the input recal file. There is no lod for the record at: " + vc );
                }
                final double lod;
                try {
                    lod = Double.valueOf(lodString);
                } catch (NumberFormatException e) {
                    throw new UserException("Encountered a malformed record in the input recal file. The lod is unreadable for the record at: " + vc );
                }

                VariantContextBuilder builder = new VariantContextBuilder(vc);

                // Annotate the new record with its VQSLOD and the worst performing annotation
                builder.attribute(GATKVCFConstants.VQS_LOD_KEY, lod);
                builder.attribute(GATKVCFConstants.CULPRIT_KEY, recalDatum.getAttribute(GATKVCFConstants.CULPRIT_KEY));
                if ( recalDatum.hasAttribute(GATKVCFConstants.POSITIVE_LABEL_KEY))
                    builder.attribute(GATKVCFConstants.POSITIVE_LABEL_KEY, true);
                if ( recalDatum.hasAttribute(GATKVCFConstants.NEGATIVE_LABEL_KEY))
                    builder.attribute(GATKVCFConstants.NEGATIVE_LABEL_KEY, true);

                final String filterString = generateFilterString(lod);

                if( filterString.equals(VCFConstants.PASSES_FILTERS_v4) ) {
                    builder.passFilters();
                } else {
                    builder.filters(filterString);
                }

                final VariantContext outputVC = builder.make();
                if( !EXCLUDE_FILTERED || outputVC.isNotFiltered() ) {
                    vcfWriter.add( outputVC );
                }
            } else { // valid VC but not compatible with this mode, so just emit the variant untouched
                vcfWriter.add( vc );
            }
        }

        return 1; // This value isn't used for anything
    }

    /**
     * Generate the VCF filter string for this record based on the provided lod score
     * @param lod non-null double
     * @return the String to use as the VCF filter field
     */
    protected String generateFilterString( final double lod ) {
        String filterString = null;
        if( TS_FILTER_LEVEL != null ) {
            for( int i = tranches.size() - 1; i >= 0; i-- ) {
                final Tranche tranche = tranches.get(i);
                if( lod >= tranche.minVQSLod ) {
                    if( i == tranches.size() - 1 ) {
                        filterString = VCFConstants.PASSES_FILTERS_v4;
                    } else {
                        filterString = tranche.name;
                    }
                    break;
                }
            }

            if( filterString == null ) {
                filterString = tranches.get(0).name+"+";
            }
        } else {
            filterString = ( lod < VQSLOD_CUTOFF ? LOW_VQSLOD_FILTER_NAME : VCFConstants.PASSES_FILTERS_v4 );
        }

        return filterString;
    }

    private static VariantContext getMatchingRecalVC(final VariantContext target, final List<VariantContext> recalVCs) {
        for( final VariantContext recalVC : recalVCs ) {
            if ( target.getEnd() == recalVC.getEnd() ) {
                return recalVC;
            }
        }

        return null;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    public Integer reduceInit() {
        return 1; // This value isn't used for anything
    }

    public Integer reduce( final Integer mapValue, final Integer reduceSum ) {
        return 1; // This value isn't used for anything
    }

    public Integer treeReduce( final Integer lhs, final Integer rhs ) {
        return 1; // This value isn't used for anything
    }

    public void onTraversalDone( final Integer reduceSum ) {
    }
}

