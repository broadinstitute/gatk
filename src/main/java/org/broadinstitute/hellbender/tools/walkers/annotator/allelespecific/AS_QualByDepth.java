package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.QualByDepth;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

/**
 * Allele-specific call confidence normalized by depth of sample reads supporting the allele
 *
 * <p>This annotation puts the variant confidence QUAL score into perspective by normalizing for the amount of coverage available. Because each read contributes a little to the QUAL score, variants in regions with deep coverage can have artificially inflated QUAL scores, giving the impression that the call is supported by more evidence than it really is. To compensate for this, we normalize the variant confidence by depth, which gives us a more objective picture of how well supported the call is.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The QD is the QUAL score normalized by allele depth (AD) for a variant. For a single sample, the HaplotypeCaller calculates the QD by taking QUAL/AD. For multiple samples, HaplotypeCaller and GenotypeGVCFs calculate the QD by taking QUAL/AD of samples with a non hom-ref genotype call. The reason we leave out the samples with a hom-ref call is to not penalize the QUAL for the other samples with the variant call.</p>
 * <h4>Here is a single-sample example:</h4>
 * <pre>2	37629	.	C	G	1063.77	.	AC=2;AF=1.00;AN=2;DP=31;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=58.50;QD=34.32;SOR=2.376	GT:AD:DP:GQ:PL:QSS	1/1:0,31:31:93:1092,93,0:0,960</pre>
 <p>QUAL/AD = 1063.77/31 = 34.32 = QD</p>
 * <h4>Here is a multi-sample example:</h4>
 * <pre>10	8046	.	C	T	4107.13	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-3.717;DP=1063;FS=1.616;MLEAC=1;MLEAF=0.167;QD=11.54
 GT:AD:DP:GQ:PL:QSS	0/0:369,4:373:99:0,1007,12207:10548,98	    0/0:331,1:332:99:0,967,11125:9576,27	    0/1:192,164:356:99:4138,0,5291:5501,4505</pre>
 * <p>QUAL/AD = 4107.13/356 = 11.54 = QD</p>
 * <p>Note that currently, when HaplotypeCaller is run with `-ERC GVCF`, the QD calculation is invoked before AD itself has been calculated, due to a technical constraint. In that case, HaplotypeCaller uses the number of overlapping reads from the haplotype likelihood calculation in place of AD to calculate QD, which generally yields a very similar number. This does not cause any measurable problems, but can cause some confusion since the number may be slightly different than what you would expect to get if you did the calculation manually. For that reason, this behavior will be modified in an upcoming version.</p>
 *
 * <h3>Caveat</h3>
 * <p>This annotation can only be calculated for sites for which at least one sample was genotyped as carrying a variant allele.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_AS_QualByDepth.php">AS_QualByDepth</a></b> outputs a version of this annotation that includes all alternate alleles in a single calculation.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_Coverage.php">Coverage</a></b> gives the filtered depth of coverage for each sample and the unfiltered depth across all samples.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_DepthPerAlleleBySample.php">DepthPerAlleleBySample</a></b> calculates depth of coverage for each allele per sample (AD).</li>
 * </ul>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Allele-specific call confidence normalized by depth of sample reads supporting the allele (AS_QD)")
public class AS_QualByDepth extends InfoFieldAnnotation implements ReducibleAnnotation, AS_StandardAnnotation {

    @Override
    public List<String> getKeyNames() { return Arrays.asList(GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY); }

    @Override
    public String getRawKeyName() { return GATKVCFConstants.AS_QUAL_KEY; }


    @Override
    public List<VCFInfoHeaderLine> getRawDescriptions() {
        //We only have the finalized key name here because the raw key is internal to GenotypeGVCFs and won't get output in any VCF
        return getDescriptions();
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods ) {
        return Collections.emptyMap();
    }

    /**
     * Note: There is no raw annotation for AS_QualByDepth and thus this method does nothing.
     *       We expect the "AS_QUAL" key to be generated by genotypeGVCFs during genotyping.
     *
     * @param ref the reference context for this annotation
     * @param vc the variant context to annotate
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample
     * @return
     */
    @Override
    public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final ReadLikelihoods<Allele> likelihoods) {
        return Collections.emptyMap();
    }

    /**
     * Note there is no raw annotation data for AS_QualByDepth and thus data cannot be combined
     *
     * @param allelesList   The merged allele list across all variants being combined/merged
     * @param listOfRawData The raw data for all the variants being combined/merged
     * @return
     */
    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    public Map<String, Object> combineRawData(List<Allele> allelesList, List<ReducibleAnnotationData<?>>  listOfRawData) {
        return null;
    }


    /**
     * Uses the "AS_QUAL" key, which must be computed by the genotyping engine in GenotypeGVCFs, to
     * calculate the final AS_QD annotation on the read.
     *
     * @param vc -- contains the final set of alleles, possibly subset by GenotypeGVCFs
     * @param originalVC -- used to get all the alleles for all gVCFs
     * @return
     */
    @Override
    public Map<String, Object> finalizeRawData(VariantContext vc, VariantContext originalVC) {
        //we need to use the AS_QUAL value that was added to the VC by the GenotypingEngine
        if ( !vc.hasAttribute(GATKVCFConstants.AS_QUAL_KEY) ) {
            return null;
        }

        final GenotypesContext genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.isEmpty() ) {
            return null;
        }

        final List<Integer> standardDepth = getAlleleDepths(genotypes);
        if (standardDepth == null) { //all no-calls and homRefs
            return null;
        }

        //Parse the VC's allele-specific qual values
        List<Object> alleleQualObjList = vc.getAttributeAsList(GATKVCFConstants.AS_QUAL_KEY);
        if (alleleQualObjList.size() != vc.getNAlleles() -1) {
            throw new IllegalStateException("Number of AS_QUAL values doesn't match the number of alternate alleles.");
        }
        List<Double> alleleQualList = new ArrayList<>();
        for (final Object obj : alleleQualObjList) {
            alleleQualList.add(Double.parseDouble(obj.toString()));
        }

        // Don't normalize indel length for AS_QD because it will only be called from GenotypeGVCFs, never UG
        List<Double> QDlist = new ArrayList<>();
        double refDepth = (double)standardDepth.get(0);
        for (int i = 0; i < alleleQualList.size(); i++) {
            double AS_QD = -10.0 * alleleQualList.get(i) / ((double)standardDepth.get(i+1) + refDepth); //+1 to skip the reference field of the AD, add ref counts to each to match biallelic case
            // Hack: see note in the fixTooHighQD method below
            AS_QD = QualByDepth.fixTooHighQD(AS_QD);
            QDlist.add(AS_QD);
        }

        final Map<String, Object> map = new HashMap<>();
        map.put(getKeyNames().get(0), AnnotationUtils.encodeValueList(QDlist, "%.2f"));
        return map;
    }

    private List<Integer> getAlleleDepths(final GenotypesContext genotypes) {
        int numAlleles = -1;
        for (final Genotype genotype : genotypes) {
            if (genotype.hasAD()) {
                numAlleles = genotype.getAD().length;
                break;
            }
        }
        if (numAlleles == -1) { //no genotypes have AD
            return null;
        }
        Integer[] alleleDepths = new Integer[numAlleles];
        for (int i = 0; i < alleleDepths.length; i++) {
            alleleDepths[i] = 0;
        }
        for (final Genotype genotype : genotypes) {
            // we care only about genotypes with variant alleles
            if ( !genotype.isHet() && !genotype.isHomVar() ) {
                continue;
            }

            // if we have the AD values for this sample, let's make sure that the variant depth is greater than 1!
            if ( genotype.hasAD() ) {
                final int[] AD = genotype.getAD();
                final int totalADdepth = (int) MathUtils.sum(AD);
                if ( totalADdepth - AD[0] > 1 ) {
                    for (int i = 0; i < AD.length; i++) {
                        alleleDepths[i] += AD[i];
                    }
                }
            }
        }
        return Arrays.asList(alleleDepths);
    }
}
