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

package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.engine.arguments.DbsnpArgumentCollection;
import org.broadinstitute.gatk.engine.arguments.GenotypeCalculationArgumentCollection;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.genotyper.*;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.GeneralPloidyFailOverAFCalculatorProvider;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.genotyper.IndexedSampleList;
import org.broadinstitute.gatk.utils.genotyper.SampleList;
import org.broadinstitute.gatk.utils.genotyper.SampleListUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Perform joint genotyping on gVCF files produced by HaplotypeCaller
 *
 * <p>
 * GenotypeGVCFs merges gVCF records that were produced as part of the Best Practices workflow for variant discovery
 * (see Best Practices documentation for more details) using the '-ERC GVCF' or '-ERC BP_RESOLUTION' mode of the
 * HaplotypeCaller, or result from combining such gVCF files using CombineGVCFs. This tool performs the multi-sample
 * joint aggregation step and merges the records together in a sophisticated manner: at each position of the input
 * gVCFs, this tool will combine all spanning records, produce correct genotype likelihoods, re-genotype the newly
 * merged record, and then re-annotate it.</p>
 *
 * <h3>Input</h3>
 * <p>
 * One or more HaplotypeCaller gVCFs to genotype.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined, genotyped VCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T GenotypeGVCFs \
 *   -R reference.fasta \
 *   --variant sample1.g.vcf \
 *   --variant sample2.g.vcf \
 *   -o output.vcf
 * </pre>
 *
 * <h3>Caveat</h3>
 * <p>Only gVCF files produced by HaplotypeCaller (or CombineGVCFs) can be used as input for this tool. Some other
 * programs produce files that they call gVCFs but those lack some important information (accurate genotype likelihoods
 * for every position) that GenotypeGVCFs requires for its operation.</p>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool is able to handle any ploidy (or mix of ploidies) intelligently; there is no need to specify ploidy
 * for non-diploid organisms.</p>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=-10,stop=10))
@SuppressWarnings("unused")
public class GenotypeGVCFs extends RodWalker<VariantContext, VariantContextWriter> implements AnnotatorCompatible, TreeReducible<VariantContextWriter> {

    /**
     * The gVCF files to merge together
     */
    @Input(fullName="variant", shortName = "V", doc="One or more input gVCF files", required=true)
    public List<RodBindingCollection<VariantContext>> variantCollections;
    final private List<RodBinding<VariantContext>> variants = new ArrayList<>();

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    @Argument(fullName="includeNonVariantSites", shortName="allSites", doc="Include loci found to be non-variant after genotyping", required=false)
    public boolean INCLUDE_NON_VARIANTS = false;

    /**
     * Uniquify all sample names (intended for use with multiple inputs for the same sample)
     */
    @Hidden
    @Advanced
    @Argument(fullName="uniquifySamples", shortName="uniquifySamples", doc="Assume duplicate samples are present and uniquify all names with '.variant' and file number index")
    public boolean uniquifySamples = false;

   @ArgumentCollection
    public GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    /**
     * Which annotations to recompute for the combined output VCF file.
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to recompute.  The single value 'none' removes the default annotations", required=false)
    protected List<String> annotationsToUse = new ArrayList<>(Arrays.asList(new String[]{"InbreedingCoeff", "FisherStrand", "QualByDepth", "ChromosomeCounts", "StrandOddsRatio"}));

    /**
     * Which groups of annotations to add to the output VCF file. The single value 'none' removes the default group. See
     * the VariantAnnotator -list argument to view available groups. Note that this usage is not recommended because
     * it obscures the specific requirements of individual annotations. Any requirements that are not met (e.g. failing
     * to provide a pedigree file for a pedigree-based annotation) may cause the run to fail.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected String[] annotationGroupsToUse = {};


    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    public RodBinding<VariantContext> getDbsnpRodBinding() { return dbsnp.dbsnp; }

    // the genotyping engine
    private UnifiedGenotypingEngine genotypingEngine;
    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    public List<RodBinding<VariantContext>> getCompRodBindings() { return Collections.emptyList(); }
    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }
    public boolean alwaysAppendDbsnpId() { return false; }


    public void initialize() {
        boolean inputsAreTagged = false;

        // collect the actual rod bindings into a list for use later
        for ( final RodBindingCollection<VariantContext> variantCollection : variantCollections ) {
            variants.addAll(variantCollection.getRodBindings());
            if (uniquifySamples) {
                for (final RodBinding<VariantContext> rb : variantCollection.getRodBindings()) {
                    //are inputs passed in with -V:fileTag ?
                    if (!rb.getTags().isEmpty()) inputsAreTagged = true;
                }
            }
        }
        //RodBinding tags are used in sample uniquification
        if (inputsAreTagged)
            logger.warn("Output uniquified VCF may not be suitable for input to CombineSampleData because input VCF(s) contain tags.");

        final GenomeAnalysisEngine toolkit = getToolkit();
        final Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(toolkit, variants);

        final GATKVariantContextUtils.GenotypeMergeType mergeType;
        if(uniquifySamples) {
            mergeType = GATKVariantContextUtils.GenotypeMergeType.UNIQUIFY;
        }
        else
            mergeType = GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE;

        final SampleList samples = new IndexedSampleList(SampleUtils.getSampleList(vcfRods, mergeType));
        // create the genotyping engine
        genotypingEngine = new UnifiedGenotypingEngine(createUAC(), samples, toolkit.getGenomeLocParser(), GeneralPloidyFailOverAFCalculatorProvider.createThreadSafeProvider(toolkit, genotypeArgs, logger),
                toolkit.getArguments().BAQMode);
        // create the annotation engine
        annotationEngine = new VariantAnnotatorEngine(Arrays.asList(annotationGroupsToUse), annotationsToUse, Collections.<String>emptyList(), this, toolkit);

        // take care of the VCF headers
        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions());
        headerLines.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());

        // add headers for annotations added by this tool
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        if ( dbsnp != null && dbsnp.dbsnp.isBound() )
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);

        final Set<String> sampleNameSet = SampleListUtils.asSet(samples);
        final VCFHeader vcfHeader = new VCFHeader(headerLines, sampleNameSet);
        vcfWriter.writeHeader(vcfHeader);

        //now that we have all the VCF headers, initialize the annotations (this is particularly important to turn off RankSumTest dithering in integration tests)
        annotationEngine.invokeAnnotationInitializationMethods(headerLines);

        logger.info("Notice that the -ploidy parameter is ignored in " + getClass().getSimpleName() + " tool as this is automatically determined by the input variant files");
    }

    public VariantContext map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return null;

        final GenomeLoc loc = ref.getLocus();
        final VariantContext combinedVC = ReferenceConfidenceVariantContextMerger.merge(tracker.getPrioritizedValue(variants, loc), loc, INCLUDE_NON_VARIANTS ? ref.getBase() : null, true, uniquifySamples, annotationEngine);
        if ( combinedVC == null )
            return null;
        return regenotypeVC(tracker, ref, combinedVC);
    }

    /**
     * Re-genotype (and re-annotate) a combined genomic VC
     *
     * @param tracker        the ref tracker
     * @param ref            the ref context
     * @param originalVC     the combined genomic VC
     * @return a new VariantContext or null if the site turned monomorphic and we don't want such sites
     */
    protected VariantContext regenotypeVC(final RefMetaDataTracker tracker, final ReferenceContext ref, final VariantContext originalVC) {
        if ( originalVC == null ) throw new IllegalArgumentException("originalVC cannot be null");

        VariantContext rawResult = originalVC;

        // only re-genotype polymorphic sites
        if ( rawResult.isVariant() ) {
            VariantContext regenotypedVC = genotypingEngine.calculateGenotypes(rawResult);
            if ( ! isProperlyPolymorphic(regenotypedVC) ) {
                if (!INCLUDE_NON_VARIANTS)
                    return null;
            }
            else {
                rawResult = addGenotypingAnnotations(rawResult.getAttributes(), regenotypedVC);
            }
        }

        //At this point we should already have DP and AD annotated
        VariantContext result = annotationEngine.finalizeAnnotations(rawResult, originalVC);
        //do trimming after allele-specific annotation reduction or the mapping is difficult
        result = GATKVariantContextUtils.reverseTrimAlleles(result);

        // if it turned monomorphic then we either need to ignore or fix such sites
        boolean createRefGTs = false;
        if ( result.isMonomorphicInSamples() ) {
            if ( !INCLUDE_NON_VARIANTS )
                return null;
            createRefGTs = true;
        }

        // Re-annotate and fix/remove some of the original annotations.
        // Note that the order of these actions matters and is different for polymorphic and monomorphic sites.
        // For polymorphic sites we need to make sure e.g. the SB tag is sent to the annotation engine and then removed later.
        // For monomorphic sites we need to make sure e.g. the hom ref genotypes are created and only then are passed to the annotation engine.
        // We could theoretically make 2 passes to re-create the genotypes, but that gets extremely expensive with large sample sizes.
        if ( createRefGTs ) {
            result = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, true)).make();
            result = annotationEngine.annotateContext(tracker, ref, null, result);
        } else {
            result = annotationEngine.annotateContext(tracker, ref, null, result);
            result = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, false)).make();
        }

        return result;
    }

    /**
     * Determines whether the provided VariantContext has real alternate alleles
     *
     * @param vc  the VariantContext to evaluate
     * @return true if it has proper alternate alleles, false otherwise
     */
    private boolean isProperlyPolymorphic(final VariantContext vc) {
        return ( vc != null &&
                !vc.getAlternateAlleles().isEmpty() &&
                (!vc.isBiallelic() ||
                        (!vc.getAlternateAllele(0).equals(Allele.SPAN_DEL) &&
                                !vc.getAlternateAllele(0).equals(GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED))
                )
        );
    }

    /**
     * Add genotyping-based annotations to the new VC
     *
     * @param originalAttributes the non-null annotations from the original VC
     * @param newVC the new non-null VC
     * @return a non-null VC
     */
    private VariantContext addGenotypingAnnotations(final Map<String, Object> originalAttributes, final VariantContext newVC) {
        // we want to carry forward the attributes from the original VC but make sure to add the MLE-based annotations
        final Map<String, Object> attrs = new HashMap<>(originalAttributes);
        attrs.put(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        attrs.put(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        if (newVC.hasAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY))
            attrs.put(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY, newVC.getAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY));

        return new VariantContextBuilder(newVC).attributes(attrs).make();
    }


    /**
     * Cleans up genotype-level annotations that need to be updated.
     * 1. move MIN_DP to DP if present
     * 2. propagate DP to AD if not present
     * 3. remove SB if present
     * 4. change the PGT value from "0|1" to "1|1" for homozygous variant genotypes
     * 5. move GQ to RGQ if the site is monomorphic
     *
     * @param VC            the VariantContext with the Genotypes to fix
     * @param createRefGTs  if true we will also create proper hom ref genotypes since we assume the site is monomorphic
     * @return a new set of Genotypes
     */
    private List<Genotype> cleanupGenotypeAnnotations(final VariantContext VC, final boolean createRefGTs) {
        final GenotypesContext oldGTs = VC.getGenotypes();
        final List<Genotype> recoveredGs = new ArrayList<>(oldGTs.size());
        for ( final Genotype oldGT : oldGTs ) {
            final Map<String, Object> attrs = new HashMap<>(oldGT.getExtendedAttributes());

            final GenotypeBuilder builder = new GenotypeBuilder(oldGT);
            int depth = oldGT.hasDP() ? oldGT.getDP() : 0;

            // move the MIN_DP to DP
            if ( oldGT.hasExtendedAttribute("MIN_DP") ) {
                depth = Integer.parseInt((String)oldGT.getAnyAttribute("MIN_DP"));
                builder.DP(depth);
                attrs.remove("MIN_DP");
            }

            // move the GQ to RGQ
            if ( createRefGTs && oldGT.hasGQ() ) {
                builder.noGQ();
                attrs.put(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, oldGT.getGQ());
            }

            // remove SB
            attrs.remove("SB");

            // update PGT for hom vars
            if ( oldGT.isHomVar() && oldGT.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) ) {
                attrs.put(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, "1|1");
            }

            // create AD if it's not there
            if ( !oldGT.hasAD() && VC.isVariant() ) {
                final int[] AD = new int[VC.getNAlleles()];
                AD[0] = depth;
                builder.AD(AD);
            }

            if ( createRefGTs ) {
                final int ploidy = oldGT.getPloidy();
                final List<Allele> refAlleles = Collections.nCopies(ploidy,VC.getReference());

                //keep 0 depth samples and 0 GQ samples as no-call
                if (depth > 0 && oldGT.hasGQ() && oldGT.getGQ() > 0) {
                    builder.alleles(refAlleles);
                }

                // also, the PLs are technically no longer usable
                builder.noPL();
            }

            recoveredGs.add(builder.noAttributes().attributes(attrs).make());
        }
        return recoveredGs;
    }

    private void checkRODtags() {

    }

    /**
     * Creates a UnifiedArgumentCollection with appropriate values filled in from the arguments in this walker
     * @return a complete UnifiedArgumentCollection
     */
    private UnifiedArgumentCollection createUAC() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.genotypeArgs = genotypeArgs.clone();
        return uac;
    }

    public VariantContextWriter reduceInit() {
        return vcfWriter;
    }

    public VariantContextWriter reduce(final VariantContext vc, final VariantContextWriter writer) {
        if ( vc != null )
            writer.add(vc);
        return writer;
    }

    @Override
    public VariantContextWriter treeReduce(final VariantContextWriter lhs, final VariantContextWriter rhs) {
        return lhs;
    }

    @Override
    public void onTraversalDone(final VariantContextWriter writer) {}
}
