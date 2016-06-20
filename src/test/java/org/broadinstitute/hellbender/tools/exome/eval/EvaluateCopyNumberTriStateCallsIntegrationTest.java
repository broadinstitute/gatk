package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetArgumentCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriStateAllele;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.GenotypeCopyNumberTriStateSegments;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Field;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Integration tests for {@link EvaluateCopyNumberTriStateCalls}.
 */
public class EvaluateCopyNumberTriStateCallsIntegrationTest extends CommandLineProgramTest {

    private static final File CALLS_TEST_FILE = new File(CommandLineProgramTest.getTestDataDir(),"exome/eval/eval-calls.vcf.gz");
    private static final File TRUTH_TEST_FILE = new File(CommandLineProgramTest.getTestDataDir(),"exome/eval/eval-truth.vcf.gz");
    private static final File TARGETS_FILE = new File(CommandLineProgramTest.getTestDataDir(), "exome/eval/eval-targets.tsv");

    @Test(dataProvider = "otherOptionsProvider")
    public void testRun(final EvaluationFiltersArgumentCollection filteringOptions) {
        final File vcfOutput = createTempFile("evaluate-vcf-out", ".vcf");
        final File sampleOutput = createTempFile("evaluate-sample-out", ".tab");
        final File caseOutput = createTempFile("evaluate-case-out", ".tab");
        runCommandLine(TRUTH_TEST_FILE, CALLS_TEST_FILE, TARGETS_FILE, vcfOutput, sampleOutput, caseOutput, filteringOptions);
        checkOutputTruthConcordance(TRUTH_TEST_FILE, TARGETS_FILE, vcfOutput, filteringOptions);
        checkOutputCallsWithOverlappingTruthConcordance(TRUTH_TEST_FILE, CALLS_TEST_FILE, TARGETS_FILE, vcfOutput, filteringOptions);
        checkOutputCallsWithoutOverlappingTruthConcordance(TRUTH_TEST_FILE, CALLS_TEST_FILE, TARGETS_FILE, vcfOutput, filteringOptions);
        checkOutputTargetNumbers(TARGETS_FILE, vcfOutput);
        checkCaseOutputConsistency(vcfOutput, caseOutput);
        checkSampleOutputConsistency(vcfOutput, sampleOutput);
    }

    @DataProvider(name = "otherOptionsProvider")
    public Object[][] otherOptionsProvider() {
        final List<Object[]> result = new ArrayList<>();
        final EvaluationFiltersArgumentCollection filterOptions = new EvaluationFiltersArgumentCollection();
        for (final boolean applyMACFilter : Arrays.asList(true, false)) {
            for (final int minCallTargetCount : Arrays.asList(1, 15)) {
                for (final int minCallQuality : Arrays.asList(0, 50)) {
                    for (final double maxCallFreq : Arrays.asList(0.0, 0.10)) {
                        for (final boolean applyMATFilter : Arrays.asList(true, false)) {
                            for (final double maxTruthFreq : Arrays.asList(1.0, 0.15)) {
                                for (final int minTruthLen : Arrays.asList(1, 5)) {
                                    for (final double minTruthQual : Arrays.asList(0.0, 300.0)) {
                                        filterOptions.applyMultiAllelicTruthFilter = applyMATFilter;
                                        filterOptions.applyMultiAllelicCalledFilter = applyMACFilter;
                                        filterOptions.minimumCalledSegmentLength = minCallTargetCount;
                                        filterOptions.minimumCalledSegmentQuality = minCallQuality;
                                        filterOptions.maximumCalledEventFrequency = maxCallFreq;
                                        filterOptions.maximumTruthEventFrequency = maxTruthFreq;
                                        filterOptions.minimumTruthSegmentLength = minTruthLen;
                                        filterOptions.minimumTruthSegmentQuality = minTruthQual;
                                        result.add(new Object[]{filterOptions.clone()});
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    private void checkSampleOutputConsistency(final File vcfOutput, final File sampleOutput) {
        final List<VariantContext> outputVariants = readVCFFile(vcfOutput);
        final List<EvaluationSampleSummaryRecord> summaries;
        try (final EvaluationSampleSummaryReader reader = new EvaluationSampleSummaryReader(sampleOutput)) {
            summaries = reader.toList();
        } catch (final IOException ex) {
            throw new RuntimeException("could not open the sample-output", ex);
        }
        final Set<String> samples = outputVariants.get(0).getSampleNames();
        final Map<String, List<EvaluationSampleSummaryRecord>> summaryBySample = summaries.stream().collect(Collectors.groupingBy(EvaluationSampleSummaryRecord::getSample));
        Assert.assertTrue(summaryBySample.keySet().contains(EvaluateCopyNumberTriStateCalls.DEFAULT_OVERALL_SUMMARY_SAMPLE_NAME));
        Assert.assertEquals(samples.size(), summaryBySample.keySet().size() - 1);
        Assert.assertTrue(summaryBySample.keySet().containsAll(samples));
        summaryBySample.values().forEach(vs -> Assert.assertTrue(vs.size() == 1));
        final int[] totalByClass = new int[EvaluationClass.values().length];
        for (final String sample : samples) {
            final Map<EvaluationClass, List<EvaluationClass>> evs = outputVariants.stream()
                    .filter(v -> v.getFilters().isEmpty() || (v.getFilters().size() == 1 && v.getFilters().contains(VCFConstants.PASSES_FILTERS_v4)))
                    .map(v -> v.getGenotype(sample))
                    .filter(g -> g.getFilters() == null || g.getFilters().isEmpty() || g.getFilters().equals(VCFConstants.PASSES_FILTERS_v4))
                    .map(g -> GATKProtectedVariantContextUtils.getAttributeAsString(g, VariantEvaluationContext.EVALUATION_CLASS_KEY, null))
                    .filter(Objects::nonNull)
                    .map(EvaluationClass::parseString)
                    .collect(Collectors.groupingBy(ev -> ev));
            for (final EvaluationClass clazz : EvaluationClass.values()) {
                evs.putIfAbsent(clazz, Collections.emptyList());
            }
            final EvaluationSampleSummaryRecord outputRecord = summaryBySample.get(sample).get(0);
            int sampleTotal = 0;
            for (final EvaluationClass ev : EvaluationClass.values()) {
                final int classCount = evs.get(ev).size();
                sampleTotal += classCount;
                totalByClass[ev.ordinal()] += classCount;
                Assert.assertEquals(outputRecord.get(ev), evs.get(ev).size());
            }
            Assert.assertEquals(outputRecord.getTotal(), sampleTotal);
        }
        final EvaluationSampleSummaryRecord overallRecord = summaryBySample.get(EvaluateCopyNumberTriStateCalls.DEFAULT_OVERALL_SUMMARY_SAMPLE_NAME).get(0);
        for (final EvaluationClass clazz : EvaluationClass.values()) {
            Assert.assertEquals(overallRecord.get(clazz), totalByClass[clazz.ordinal()]);
        }
    }

    private void checkCaseOutputConsistency(final File vcfOutput, final File caseOutput) {
        // For now we don't test the content of this file as is generated only for courtesy (debugging).
        // If it proves to be useful and some bugs on its content start to craw up perhaps then we should add some tests
        // here.
        Assert.assertTrue(caseOutput.exists());
    }

    private void checkOutputTargetNumbers(final File targetsFile, final File vcfOutput) {
        final List<VariantContext> outputVariants = readVCFFile(vcfOutput);
        final TargetCollection<Target> targets = TargetArgumentCollection.readTargetCollection(targetsFile);
        for (final VariantContext context : outputVariants) {
            final List<Target> overlappingTargets = targets.targets(context);
            Assert.assertEquals(context.getAttributeAsInt(GenotypeCopyNumberTriStateSegments.NUMBER_OF_TARGETS_KEY, -1), overlappingTargets.size());
        }
    }

    private void checkOutputTruthConcordance(final File truthFile, final File targetsFile, final File vcfOutput, final EvaluationFiltersArgumentCollection filteringOptions) {
        final List<VariantContext> truthVariants = readVCFFile(truthFile);
        final List<VariantContext> outputVariants = readVCFFile(vcfOutput);
        final Set<String> outputSamples = outputVariants.get(0).getSampleNames();
        final TargetCollection<Target> targets = TargetArgumentCollection.readTargetCollection(targetsFile);

        for (final VariantContext truth : truthVariants) {
            final List<Target> overlappingTargets = targets.targets(truth);
            final List<VariantContext> overlappingOutput = outputVariants.stream()
                    .filter(vc -> new SimpleInterval(vc).overlaps(truth))
                    .collect(Collectors.toList());

            if (overlappingTargets.isEmpty()) {
                Assert.assertTrue(overlappingOutput.isEmpty());
                continue;
            }
            Assert.assertFalse(overlappingOutput.isEmpty());
            Assert.assertEquals(overlappingOutput.stream().filter(vc -> new SimpleInterval(truth).equals(new SimpleInterval(vc))).count(), 1);
            @SuppressWarnings("all")
            final Optional<VariantContext> prospectiveMatchingOutput = overlappingOutput.stream()
                    .filter(vc -> new SimpleInterval(truth).equals(new SimpleInterval(vc)))
                    .findFirst();
            Assert.assertTrue(prospectiveMatchingOutput.isPresent());
            final VariantContext matchingOutput = prospectiveMatchingOutput.get();
            final int[] truthAC = calculateACFromTruth(truth);
            final long truthAN = MathUtils.sum(truthAC);
            final double[] truthAF = IntStream.of(Arrays.copyOfRange(truthAC, 1, truthAC.length)).mapToDouble(d -> d / (double) truthAN).toArray();
            Assert.assertEquals(matchingOutput.getAttributeAsInt(VariantEvaluationContext.TRUTH_ALLELE_NUMBER_KEY, -1), truthAN);
            assertEquals(GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(matchingOutput, VariantEvaluationContext.TRUTH_ALLELE_FREQUENCY_KEY, () -> new double[2], 0.0),
                    truthAF, 0.001);
            assertOutputVariantFilters(filteringOptions, overlappingTargets, matchingOutput, truthAF);
            for (final String sample : outputSamples) {
                final Genotype outputGenotype = matchingOutput.getGenotype(sample);
                final Genotype truthGenotype = truth.getGenotype(sample);
                final int truthCN = GATKProtectedVariantContextUtils.getAttributeAsInt(truthGenotype, ConvertGSVariantsToSegments.GS_COPY_NUMBER_FORMAT, -1);
                final int truthGT = GATKProtectedVariantContextUtils.getAttributeAsInt(outputGenotype, VariantEvaluationContext.TRUTH_GENOTYPE_KEY, -1);
                final Object truthQualObject = outputGenotype.getAnyAttribute(VariantEvaluationContext.TRUTH_QUALITY_KEY);
                Assert.assertNotNull(truthQualObject, "" + truthGenotype);
                final double truthQual = Double.parseDouble(String.valueOf(truthQualObject));
                if (truthQual < filteringOptions.minimumTruthSegmentQuality) {
                    Assert.assertEquals(outputGenotype.getFilters(),EvaluationFilter.LowQuality.acronym);
                } else {
                    Assert.assertTrue(outputGenotype.getFilters() == null || outputGenotype.getFilters().equals(VCFConstants.PASSES_FILTERS_v4));
                }
                final double[] truthPosteriors = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(truthGenotype, ConvertGSVariantsToSegments.GS_COPY_NUMBER_POSTERIOR,
                        () -> new double[0], Double.NEGATIVE_INFINITY);
                final double truthDelPosterior = MathUtils.log10SumLog10(truthPosteriors, 0, EvaluateCopyNumberTriStateCalls.REFERENCE_COPY_NUMBER_DEFAULT);
                final double truthRefPosterior = truthPosteriors[EvaluateCopyNumberTriStateCalls.REFERENCE_COPY_NUMBER_DEFAULT];
                final double truthDupPosterior = truthPosteriors.length < EvaluateCopyNumberTriStateCalls.REFERENCE_COPY_NUMBER_DEFAULT ? Double.NEGATIVE_INFINITY : MathUtils.log10SumLog10(truthPosteriors, EvaluateCopyNumberTriStateCalls.REFERENCE_COPY_NUMBER_DEFAULT + 1, truthPosteriors.length);

                final CopyNumberTriStateAllele truthAllele = truthGT == -1 ? null : CopyNumberTriStateAllele.ALL_ALLELES.get(truthGT);
                if (truthCN < EvaluateCopyNumberTriStateCalls.REFERENCE_COPY_NUMBER_DEFAULT) {
                    Assert.assertEquals(truthAllele, CopyNumberTriStateAllele.DEL);
                    Assert.assertEquals(truthQual * -.1, MathUtils.log10SumLog10(new double[] {truthRefPosterior, truthDupPosterior}), 0.01);
                } else if (truthCN > EvaluateCopyNumberTriStateCalls.REFERENCE_COPY_NUMBER_DEFAULT) {
                    Assert.assertEquals(truthAllele, CopyNumberTriStateAllele.DUP);
                    Assert.assertEquals(truthQual * -.1, MathUtils.log10SumLog10(new double[] {truthDelPosterior, truthRefPosterior}), 0.01, "" + truthGenotype + " " + outputGenotype);
                } else {
                    Assert.assertEquals(truthAllele, CopyNumberTriStateAllele.REF);
                    Assert.assertEquals(truthQual * -.1, MathUtils.log10SumLog10(new double[] {truthDelPosterior, truthDupPosterior}), 0.01);
                }
                final double outputTruthFraction = GATKProtectedVariantContextUtils.getAttributeAsDouble(outputGenotype, VariantEvaluationContext.TRUTH_COPY_FRACTION_KEY, -1);
                final double inputTruthFraction = GATKProtectedVariantContextUtils.getAttributeAsDouble(truthGenotype, ConvertGSVariantsToSegments.GS_COPY_NUMBER_FRACTION, -1);
                Assert.assertEquals(outputTruthFraction, inputTruthFraction, 0.01);
            }
        }
    }

    private void assertOutputVariantFilters(EvaluationFiltersArgumentCollection filteringOptions, List<Target> overlappingTargets, VariantContext matchingOutput, double[] truthAF) {
        final Set<String> filters = new HashSet<>(matchingOutput.getFilters());
        filters.remove(VCFConstants.PASSES_FILTERS_v4); // don't if PASS is kept as a filter, just making sure that is not.
        if (overlappingTargets.size() < filteringOptions.minimumTruthSegmentLength) {
            Assert.assertTrue(filters.remove(EvaluationFilter.ShortEvent.name()));
        }
        if (truthAF[CopyNumberTriStateAllele.DEL.index() - 1] > 0 && truthAF[CopyNumberTriStateAllele.DUP.index() - 1] > 0 && filteringOptions.applyMultiAllelicTruthFilter) {
            Assert.assertTrue(filters.remove(EvaluationFilter.MultiAllelicTruth.name()));
        }
        if (MathUtils.sum(truthAF) > filteringOptions.maximumTruthEventFrequency) {
            Assert.assertTrue(filters.remove(EvaluationFilter.CommonEvent.name()));
        }
        Assert.assertTrue(filters.isEmpty(), "Remaining filters not explained: " + filters.stream().collect(Collectors.joining(",")));
    }

    private void checkOutputCallsWithoutOverlappingTruthConcordance(final File truthFile, final File callsFile, final File targetsFile, final File vcfOutput, final EvaluationFiltersArgumentCollection filteringOptions) {
        final List<VariantContext> truthVariants = readVCFFile(truthFile);
        final List<VariantContext> outputVariants = readVCFFile(vcfOutput);
        final List<VariantContext> callsVariants = readVCFFile(callsFile);
        final Set<String> outputSamples = outputVariants.get(0).getSampleNames();
        final TargetCollection<Target> targets = TargetArgumentCollection.readTargetCollection(targetsFile);

        for (final VariantContext call : callsVariants) {
            final List<Target> overlappingTargets = targets.targets(call);
            final List<VariantContext> overlappingOutput = outputVariants.stream()
                    .filter(vc -> new SimpleInterval(vc).overlaps(call))
                    .collect(Collectors.toList());
            final List<VariantContext> overlappingTruth = truthVariants.stream()
                    .filter(vc -> new SimpleInterval(vc).overlaps(call))
                    .collect(Collectors.toList());

            if (!overlappingTruth.isEmpty()) {
                continue;
            }

            @SuppressWarnings("all")
            final Optional<VariantContext> matchingOutputOptional = overlappingOutput.stream()
                    .filter(vc -> new SimpleInterval(call).equals(new SimpleInterval(vc)))
                    .findAny();
            final VariantContext matchingOutput = matchingOutputOptional.get();
            final int sampleCallsCount[] = new int[CopyNumberTriStateAllele.ALL_ALLELES.size()];
            for (final String sample : outputSamples) {
                final Genotype outputGenotype = matchingOutput.getGenotype(sample);
                final Genotype callGenotype = call.getGenotype(sample);
                final Allele expectedCall = callGenotype.getAllele(0).isCalled() ? CopyNumberTriStateAllele.valueOf(callGenotype.getAllele(0)) : null;
                final Allele actualCall = outputGenotype.getAllele(0).isCalled() ? CopyNumberTriStateAllele.valueOf(outputGenotype.getAllele(0)) : null;
                Assert.assertEquals(expectedCall, actualCall);
                final boolean expectedDiscovered = GenotypeCopyNumberTriStateSegments.DISCOVERY_TRUE.equals(GATKProtectedVariantContextUtils.getAttributeAsString(callGenotype, GenotypeCopyNumberTriStateSegments.DISCOVERY_KEY, "N"));
                final boolean actualDiscovered = GenotypeCopyNumberTriStateSegments.DISCOVERY_TRUE.equals(GATKProtectedVariantContextUtils.getAttributeAsString(callGenotype, GenotypeCopyNumberTriStateSegments.DISCOVERY_KEY, "N"));
                Assert.assertEquals(actualDiscovered, expectedDiscovered);
                final int[] expectedCounts = new int[CopyNumberTriStateAllele.ALL_ALLELES.size()];
                if (expectedCall.isCalled() && actualDiscovered) {
                    expectedCounts[CopyNumberTriStateAllele.valueOf(expectedCall).index()]++;
                }
                if (outputGenotype.hasExtendedAttribute(VariantEvaluationContext.CALLED_ALLELE_COUNTS_KEY)) {
                    Assert.assertEquals(GATKProtectedVariantContextUtils.getAttributeAsIntArray(outputGenotype, VariantEvaluationContext.CALLED_ALLELE_COUNTS_KEY, () -> new int[CopyNumberTriStateAllele.ALL_ALLELES.size()], 0), expectedCounts);
                }
                if (outputGenotype.hasExtendedAttribute(VariantEvaluationContext.CALLED_SEGMENTS_COUNT_KEY)) {
                    Assert.assertEquals(GATKProtectedVariantContextUtils.getAttributeAsInt(outputGenotype, VariantEvaluationContext.CALLED_SEGMENTS_COUNT_KEY, -1), expectedCall.isCalled() && actualDiscovered ? 1 : 0);
                }
                final String evalClass = GATKProtectedVariantContextUtils.getAttributeAsString(outputGenotype, VariantEvaluationContext.EVALUATION_CLASS_KEY, null);
                Assert.assertEquals(evalClass, expectedCall.isCalled() && actualDiscovered && expectedCall.isNonReference() ? EvaluationClass.UNKNOWN_POSITIVE.acronym : null);
                if (expectedCall.isCalled()) {
                    sampleCallsCount[CopyNumberTriStateAllele.valueOf(expectedCall).index()]++;
                }
                Assert.assertEquals(GATKProtectedVariantContextUtils.getAttributeAsDouble(outputGenotype, VariantEvaluationContext.CALL_QUALITY_KEY, 0.0), callGQ(callGenotype), 0.01);
            }
            final int expectedAN = (int) MathUtils.sum(sampleCallsCount);
            final int observedAN = matchingOutput.getAttributeAsInt(VariantEvaluationContext.CALLS_ALLELE_NUMBER_KEY, -1);
            Assert.assertEquals(observedAN, expectedAN);
            final double[] expectedAF = Arrays.copyOfRange(IntStream.of(sampleCallsCount).mapToDouble(i -> expectedAN > 0 ? i / (double) expectedAN : 0.0)
                    .toArray(), 1, sampleCallsCount.length);
            final double[] observedAF = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(matchingOutput, VariantEvaluationContext.CALLS_ALLELE_FREQUENCY_KEY, () -> new double[matchingOutput.getAlternateAlleles().size()], 0.0);
            Assert.assertNotNull(observedAF);
            assertEquals(observedAF, expectedAF, 0.01);
            Assert.assertEquals(matchingOutput.getAttributeAsInt(VariantEvaluationContext.TRUTH_ALLELE_NUMBER_KEY, -1), 0);
        }
    }

    private void checkOutputCallsWithOverlappingTruthConcordance(final File truthFile, final File callsFile, final File targetsFile, final File vcfOutput, final EvaluationFiltersArgumentCollection filteringOptions) {
        final List<VariantContext> truthVariants = readVCFFile(truthFile);
        final List<VariantContext> outputVariants = readVCFFile(vcfOutput);
        final List<VariantContext> callsVariants = readVCFFile(callsFile);
        final Set<String> outputSamples = outputVariants.get(0).getSampleNames();
        final TargetCollection<Target> targets = TargetArgumentCollection.readTargetCollection(targetsFile);

        for (final VariantContext truth : truthVariants) {
            final List<Target> overlappingTargets = targets.targets(truth);
            final List<VariantContext> overlappingOutput = outputVariants.stream()
                    .filter(vc -> new SimpleInterval(vc).overlaps(truth))
                    .collect(Collectors.toList());
            final List<VariantContext> overlappingCalls = callsVariants.stream()
                    .filter(vc -> new SimpleInterval(vc).overlaps(truth))
                    .collect(Collectors.toList());

            if (overlappingTargets.isEmpty()) {
                Assert.assertTrue(overlappingOutput.isEmpty());
                continue;
            }

            @SuppressWarnings("all")
            final VariantContext matchingOutput = overlappingOutput.stream()
                    .filter(vc -> new SimpleInterval(truth).equals(new SimpleInterval(vc)))
                    .findAny().get();
            final int sampleCallsCount[] = new int[CopyNumberTriStateAllele.ALL_ALLELES.size()];
            for (final String sample : outputSamples) {
                final Genotype outputGenotype = matchingOutput.getGenotype(sample);
                final List<Pair<VariantContext, Genotype>> sampleCalls = overlappingCalls.stream()
                        .map(vc -> new ImmutablePair<>(vc, vc.getGenotype(sample)))
                        .filter(p -> GenotypeCopyNumberTriStateSegments.DISCOVERY_TRUE.equals(p.getRight().getExtendedAttribute(GenotypeCopyNumberTriStateSegments.DISCOVERY_KEY)))
                        .filter(p -> callPassFilters(p.getLeft(), p.getRight(), targets, filteringOptions))
                        .collect(Collectors.toList());
                final int[] expectedCounts = new int[CopyNumberTriStateAllele.ALL_ALLELES.size()];
                sampleCalls.forEach(p -> {
                    expectedCounts[CopyNumberTriStateAllele.valueOf(p.getRight().getAllele(0)).index()]++;
                });
                final int[] actualCounts = GATKProtectedVariantContextUtils.getAttributeAsIntArray(outputGenotype, VariantEvaluationContext.CALLED_ALLELE_COUNTS_KEY, () -> new int[CopyNumberTriStateAllele.ALL_ALLELES.size()], 0);
                Assert.assertEquals(actualCounts, expectedCounts, Arrays.toString(actualCounts) + " " + Arrays.toString(expectedCounts));
                final int expectedTotalCount = (int) MathUtils.sum(expectedCounts);
                final int actualTotalCount = GATKProtectedVariantContextUtils.getAttributeAsInt(outputGenotype, VariantEvaluationContext.CALLED_SEGMENTS_COUNT_KEY, -1);
                Assert.assertEquals(actualTotalCount, expectedTotalCount);
                final int expectedTargetCount = sampleCalls.stream().mapToInt(p -> targets.targetCount(p.getLeft())).sum();
                final int observedTargetCount = GATKProtectedVariantContextUtils.getAttributeAsInt(outputGenotype, VariantEvaluationContext.CALLED_TARGET_COUNT_KEY, -1);
                Assert.assertEquals(observedTargetCount, expectedTargetCount);
                final Allele truthCallAllele = outputTruthAllele(outputGenotype);
                final boolean isMixed = IntStream.of(actualCounts).filter(i -> i > 0).count() > 1;
                final String evalClass = GATKProtectedVariantContextUtils.getAttributeAsString(outputGenotype, VariantEvaluationContext.EVALUATION_CLASS_KEY, null);
                if (sampleCalls.size() > 0 && !isMixed) {
                    final Pair<VariantContext, Genotype> bestCall = sampleCalls.stream()
                            .sorted((p1, p2) -> -Double.compare(callGQ(p1.getRight()), callGQ(p2.getRight())))
                            .findFirst().get();
                    final CopyNumberTriStateAllele expectedCall = CopyNumberTriStateAllele.valueOf(bestCall.getRight().getAllele(0));
                    final CopyNumberTriStateAllele actualCall = CopyNumberTriStateAllele.valueOf(outputGenotype.getAllele(0));
                    Assert.assertEquals(actualCall, expectedCall);
                    sampleCallsCount[expectedCall.index()]++;
                    if (!truthCallAllele.isReference()) {
                        if (truthCallAllele.equals(actualCall)) {
                            Assert.assertEquals(evalClass, EvaluationClass.TRUE_POSITIVE.acronym);
                        } else if (!truthCallAllele.isNoCall()) {
                            Assert.assertEquals(evalClass, EvaluationClass.DISCORDANT_POSITIVE.acronym);
                        } else {
                            Assert.assertNull(evalClass);
                        }
                    } else if (truthCallAllele.isReference()) {
                        Assert.assertEquals(evalClass, EvaluationClass.FALSE_POSITIVE.acronym);
                    }
                } else {
                    Assert.assertEquals(Allele.NO_CALL, outputGenotype.getAllele(0));
                    if (sampleCalls.isEmpty()) {
                        Assert.assertEquals(evalClass, !truthCallAllele.isReference() && truthCallAllele.isCalled() ? EvaluationClass.FALSE_NEGATIVE.acronym : null);
                        Assert.assertEquals(GATKProtectedVariantContextUtils.getAttributeAsDouble(outputGenotype, VariantEvaluationContext.CALL_QUALITY_KEY, -1), 0.0);
                    } else {
                        Assert.assertEquals(evalClass, !truthCallAllele.isReference() && truthCallAllele.isCalled() ? EvaluationClass.MIXED_POSITIVE.acronym : EvaluationClass.FALSE_POSITIVE.acronym);
                        final Pair<VariantContext, Genotype> bestCall = sampleCalls.stream()
                                .sorted((p1, p2) -> -Double.compare(callGQ(p1.getRight()), callGQ(p2.getRight())))
                                .findFirst().get();
                        Assert.assertEquals(GATKProtectedVariantContextUtils.getAttributeAsDouble(outputGenotype, VariantEvaluationContext.CALL_QUALITY_KEY, -1), callGQ(bestCall.getRight()));
                    }
                }
            }
            final int expectedAN = (int) MathUtils.sum(sampleCallsCount);
            final int observedAN = matchingOutput.getAttributeAsInt(VariantEvaluationContext.CALLS_ALLELE_NUMBER_KEY, -1);
            Assert.assertEquals(observedAN, expectedAN);
            final double[] expectedAF = Arrays.copyOfRange(IntStream.of(sampleCallsCount).mapToDouble(i -> expectedAN > 0 ? i / (double) expectedAN : 0.0)
                    .toArray(), 1, sampleCallsCount.length);
            final double[] observedAF = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(matchingOutput, VariantEvaluationContext.CALLS_ALLELE_FREQUENCY_KEY, () -> new double[matchingOutput.getAlternateAlleles().size()], 0.0);
            Assert.assertNotNull(observedAF);
            assertEquals(observedAF, expectedAF, 0.01);
        }
    }

    private Allele outputTruthAllele(final Genotype truthGenotype) {
        final int index = GATKProtectedVariantContextUtils.getAttributeAsInt(truthGenotype, VariantEvaluationContext.TRUTH_GENOTYPE_KEY, -1);
        if (index == -1) {
            return Allele.NO_CALL;
        } else {
            return CopyNumberTriStateAllele.ALL_ALLELES.get(index);
        }
    }

    private boolean callPassFilters(final VariantContext callVariant, final Genotype callGenotype, final TargetCollection<Target> targets, final EvaluationFiltersArgumentCollection filteringOptions) {
        final double[] AF = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(callVariant, VCFConstants.ALLELE_FREQUENCY_KEY, () -> new double[CopyNumberTriStateAllele.ALTERNATIVE_ALLELES.size()], 0.0);
        if (filteringOptions.applyMultiAllelicCalledFilter) {
            if (AF[0] > 0 && AF[1] > 0) {
                return false;
            }
        }
        if (filteringOptions.maximumCalledEventFrequency < MathUtils.sum(AF)) {
            return false;
        }
        if (filteringOptions.minimumCalledSegmentLength > targets.targetCount(callVariant)) {
            return false;
        }
        if (filteringOptions.minimumCalledSegmentQuality > callGQ(callGenotype)) {
            return false;
        }
        return true;
    }

    private double callGQ(final Genotype callGenotype) {
        final int[] pls = callGenotype.getPL();
        final int bestIndex = MathUtils.minElementIndex(pls);
        final int bestPL = pls[bestIndex];
        int secondBestPL = bestIndex == 0 ? pls[1] : pls[0];
        for (int i = 0; i < pls.length; i++) {
            if (i == bestIndex) {
                continue;
            } else if (secondBestPL > pls[i]) {
                secondBestPL = pls[i];
            }
        }
        return secondBestPL - bestPL;
    }


    private void assertEquals(final double[] a, final double[] b, final double epsilon) {
        Assert.assertEquals(a.length, b.length);
        for (int i = 0; i < a.length; i++) {
            Assert.assertEquals(a[i], b[i], epsilon);
        }
    }

    private int[] calculateACFromTruth(final VariantContext vc) {
        final int[] counts = new int[CopyNumberTriStateAllele.ALL_ALLELES.size()];
        for (final Genotype genotype : vc.getGenotypes()) {
            final int cn = GATKProtectedVariantContextUtils.getAttributeAsInt(genotype, ConvertGSVariantsToSegments.GS_COPY_NUMBER_FORMAT, -1);
            if (cn == -1) {
                continue;
            }
            if (cn < 2) {
                counts[CopyNumberTriStateAllele.DEL.index()]++;
            } else if (cn > 2) {
                counts[CopyNumberTriStateAllele.DUP.index()]++;
            } else {
                counts[CopyNumberTriStateAllele.REF.index()]++;
            }
        }
        return counts;
    }

    private List<VariantContext> readVCFFile(final File vcf) {
        final List<VariantContext> result = new ArrayList<>();
        try (final VCFFileReader reader = new VCFFileReader(vcf, false)) {
            reader.forEach(result::add);
        }
        return result;
    }

    private void runCommandLine(final File truth, final File calls, final File targets,
                               final File vcfOutput, final File sampleOutput, final File caseOutput,
                               final EvaluationFiltersArgumentCollection filtersOptions) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + EvaluateCopyNumberTriStateCalls.TRUTH_FILE_SHORT_NAME);
        arguments.add(truth.getAbsolutePath());
        arguments.add("-" + EvaluateCopyNumberTriStateCalls.CALLS_FILE_SHORT_NAME);
        arguments.add(calls.getAbsolutePath());
        arguments.add("-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME);
        arguments.add(targets.getAbsolutePath());
        arguments.add("-includeOverall");
        for (final Field field : filtersOptions.getClass().getFields()) {
            final Argument annotation = field.getAnnotation(Argument.class);
            if (annotation == null) {
                continue;
            }
            final Class<?> type = field.getType();
            try {
                if (type.equals(Boolean.class) || type.equals(Boolean.TYPE)) {
                    if (field.get(filtersOptions).equals(Boolean.TRUE)) {
                        arguments.add("-" + annotation.shortName());
                    }
                } else {
                    arguments.add("-" + annotation.shortName());
                    arguments.add(String.valueOf(field.get(filtersOptions)));
                }
            } catch (final IllegalAccessException e) {
                throw new RuntimeException(e);
            }
        }
        if (vcfOutput != null) {
            arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
            arguments.add(vcfOutput.getAbsolutePath());
        }
        if (sampleOutput != null) {
            arguments.add("-" + EvaluateCopyNumberTriStateCalls.SAMPLE_SUMMARY_OUTPUT_SHORT_NAME);
            arguments.add(sampleOutput.getAbsolutePath());
        }
        if (caseOutput != null) {
            arguments.add("-" + EvaluateCopyNumberTriStateCalls.DETAIL_CALL_OUTPUT_SHORT_NAME);
            arguments.add(caseOutput.getAbsolutePath());
        }
        runCommandLine(arguments);
    }

}
