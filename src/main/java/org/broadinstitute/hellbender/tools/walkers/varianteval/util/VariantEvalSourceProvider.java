package org.broadinstitute.hellbender.tools.walkers.varianteval.util;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.utils.samples.SampleDB;

import java.util.Collection;
import java.util.List;
import java.util.Set;

// TODO: We need to retain FeatureContext, or some way to allow these to query inputs that are not
// part of the driving variants (which are directly accessible in the apply method since they're included in the
// collection passed to apply)
public interface VariantEvalSourceProvider {

    default Logger getLogger() {return LogManager.getLogger(this.getClass()); }

    SAMSequenceDictionary getProviderReferenceDictionary();

    String getNameForInput(FeatureInput<VariantContext> input);

    boolean getMergeEvals();

    // see IntervalStratification
    boolean getHasIntervalsDefined();

    VCFHeader getHeaderForEvalFeatures(FeatureInput<VariantContext> evals);

    List<VariantContext> getOverlappingKnowns(int start);

    List<Feature> getOverlappingKnownCNVs();

    List<VariantContext> getOverlappingGolds(int start);

    List<Feature> getOverlappingFromIntervals();

    double getMinPhaseQuality();

    int getSamplePloidy();
    double getMendelianViolationQualThreshold();

    String getAllSampleName();
    String getAllFamilyName();

    List<FeatureInput<VariantContext>> getKnowns();

    List<FeatureInput<VariantContext>> getEvals();

    boolean isSubsettingToSpecificSamples();
    Set<String> getSampleNamesForEvaluation();

    Set<String> getFamilyNamesForEvaluation();

    int getNumberOfSamplesForEvaluation();
    Set<String> getSampleNamesForStratification();

    Set<String> getFamilyNamesForStratification();

    List<FeatureInput<VariantContext>> getComps();

    Set<SortableJexlVCMatchExp> getJexlExpressions();

    Set<String> getContigNames();

    boolean ignoreAC0Sites();

    SampleDB getSampleDB();
    long getnProcessedLoci();

    FeatureInput<Feature> getKnownCNVsFile();
}
