package org.broadinstitute.hellbender.tools.walkers.linkedreads;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections4.Predicate;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.tools.walkers.variantutils.VariantsToTable;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils.AlleleMapper;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;

/**
 * Evaluate site-level concordance of an single-sample called VCF against a joint-called VCF.
 *
 * <p>This tool evaluates two variant callsets against each other and produces a six-column summary metrics table. The summary:</p>
 *
 * <ul>
 *     <li>stratifies SNP and INDEL calls,</li>
 *     <li>tallies true-positive, false-positive and false-negative calls,</li>
 *     <li>and calculates sensitivity and precision.</li>
 * </ul>
 *
 * <p>The tool assumes all records in the --truth VCF are passing truth variants. For the -eval VCF, the tool uses only unfiltered passing calls.</p>
 *
 * <p>Optionally, the tool can be set to produce VCFs of the following variant records, annotated with each variant's concordance status:</p>
 * <ul>
 *     <li>True positives and false negatives (i.e. all variants in the truth VCF): useful for calculating sensitivity</li>
 *     <li>True positives and false positives (i.e. all variants in the eval VCF): useful for obtaining a training data
 *     set for machine learning classifiers of artifacts</li>
 * </ul>
 *
 * <p>These output VCFs can be passed to {@link VariantsToTable} to produce a TSV file for statistical analysis in R
 * or Python.</p>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 * gatk Concordance \
 *   -R reference.fa \
 *   -eval eval.vcf \
 *   --truth truth.vcf \
 *   --summary summary.tsv
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = SingleSampleVsJointCallConcordance.USAGE_SUMMARY,
        oneLineSummary = SingleSampleVsJointCallConcordance.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class SingleSampleVsJointCallConcordance extends Concordance {
    static final String USAGE_SUMMARY = "Evaluate concordance of an single-sample called VCF against a joint-called VCF";
    static final String USAGE_ONE_LINE_SUMMARY = "Evaluate concordance of an single-sample called VCF against a joint-called VCF";

    @Override
    protected boolean areVariantsAtSameLocusConcordant(final VariantContext truth, final VariantContext eval) {
        final Allele jointRef = truth.getReference();
        final Allele singleSampleRef = eval.getReference();
        if (jointRef.length() < singleSampleRef.length()) return false;

        final AlleleMapper alleleMapping = GATKVariantContextUtils.resolveIncompatibleAlleles(jointRef, eval, new LinkedHashSet<>());
        final boolean jointHasAllAltAlleles = eval.getAlternateAlleles().stream().allMatch(a -> truth.hasAllele(alleleMapping.remap(a)));

        return jointHasAllAltAlleles;
    }

    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return vc -> ! vc.isSymbolicOrSV();
    }

    @Override
    protected Predicate<VariantContext> makeEvalVariantFilter() {
        return vc -> ! vc.isSymbolicOrSV();
    }
}
