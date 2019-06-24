package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by valentin on 11/3/16.
 */
public class AssessSTRs extends STRAnalysisWalker {

    @Argument(doc = "output Context table", shortName = "O", fullName = "output")
    public VariantContextWriter output;

    @Argument(doc = "output summary", shortName = "summary", fullName = "summaryOutput", optional = true)
    public File summaryOutput;

    @Argument(doc = "whether missing truth means hom-ref.", shortName = "missingTruthMeansHomeRef", optional = true)
    public boolean missingTruthMeansHomeRef = false;

    private String lastOutputContig;
    private long lastOutputPosition;

    @Argument(doc = "minimum GQ to trust a hom/ref call", optional = true)
    public int minCallQual = 30;

    @Argument(doc = "minimum Qual to trust a truth call", optional = true)
    public int minTruthQual = 30;

    private static final String EVALUATION_CLASS_TAG = "CLASS";
    private static final String TRUTH_ALLELE_TAG = "TruthAlleles";
    private static final String STR_UNIT_TAG = "STRUnit";
    private static final String STR_UNIT_LENGTH_TAG = "STRUnitLength";
    private static final String STR_UNIT_MA_RC_TAG = "STRReferenceUnitRepeatCount";
    private static final String LOW_QUAL_FILTER = "LowQual";

    private final STREvaluationCounter strCounter = new STREvaluationCounter();

    private final STREvaluationCounter nonStrCounter = new STREvaluationCounter();

    public void composeAndSetOutputHeader() {
        final VCFHeader header = getHeaderForVariants();

        header.addMetaDataLine(new VCFInfoHeaderLine(EVALUATION_CLASS_TAG, 1, VCFHeaderLineType.String, "Evaluation class"));
        header.addMetaDataLine(new VCFInfoHeaderLine(STR_UNIT_TAG, 1, VCFHeaderLineType.String, "STR Repeat Unit"));
        header.addMetaDataLine(new VCFInfoHeaderLine(STR_UNIT_MA_RC_TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "STR Reference Repeat Count"));
        header.addMetaDataLine(new VCFInfoHeaderLine(STR_UNIT_LENGTH_TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "STR Repeat Unit Length"));
        header.addMetaDataLine(new VCFInfoHeaderLine(TRUTH_ALLELE_TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "STR Type"));
        header.addMetaDataLine(new VCFFilterHeaderLine(LOW_QUAL_FILTER, "LowQual"));
        output.writeHeader(header);
    }

    @Override
    public void apply(final VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final List<VariantContext> truths = getTruthFeatures(featureContext);
        final List<VariantContext> calls = getCallFeatures(featureContext);
        final VariantContext truth;
        final VariantContext call;
        if (isTruth(variant)) {
            truth = variant;
            call = findBestMatch(truth, calls);
        } else { // isCall(variant) == true.
            call = variant;
            truth = findBestMatch(call, truths);
        }

        if (call != null && truth == variant) {
            return; // in case of matching pair we process when variant == call, to avoid repeats.
        }
        final VariantContext out = compare(call, truth, referenceContext);
        output(out);

    }

    private void output(final VariantContext out) {
        if (out != null) {
            final String contig = out.getContig();
            final long start = out.getStart();
            final long length = out.hasAttribute(STR_UNIT_TAG)
                    ? out.getAttributeAsInt(STR_UNIT_LENGTH_TAG, 0) * out.getAttributeAsInt(STR_UNIT_MA_RC_TAG, 0)
                    : out.getReference().length() - 1;
            if (out.hasAttribute(EVALUATION_CLASS_TAG)) {
                final STREvaluationClass clazz = STREvaluationClass.valueOf(out.getAttributeAsString(EVALUATION_CLASS_TAG, null));
                final STREvaluationCounter counter = out.hasAttribute(STR_UNIT_TAG) ? strCounter : nonStrCounter;
                counter.add(clazz);
            }
            output.add(out);
        }
    }

    public String onTraversalSuccess() {
        logger.info("Number of STR context found: " + strCounter.total());
        output.close();

        try (final PrintWriter summaryWriter = new PrintWriter(new FileWriter(summaryOutput))) {
            summaryWriter.println(String.join("\t", "EVAL_CLASS", "STR", "INDEL"));
            for (final STREvaluationClass evalClass : STREvaluationClass.values()) {
                summaryWriter.println(String.join("\t", evalClass.name(), "" + strCounter.get(evalClass), "" + nonStrCounter.get(evalClass)));
            }
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(summaryOutput, ex);
        }
        return "SUCCESS";
    }

    private VariantContext compare(final VariantContext call, final VariantContext truth, final ReferenceContext ref) {

        final STREvaluationCall rawCall = STREvaluationCall.fromVariantContext(call, ref, false);
        final STREvaluationCall rawTruth = STREvaluationCall.fromVariantContext(truth, ref, missingTruthMeansHomeRef);
        if (rawCall == null || rawTruth == null) {
            logger.debug("Funny situation at " + ref.getInterval());
            return null;
        }
        // make sure they have the same reference allele sequence:
        final Pair<STREvaluationCall, STREvaluationCall> processedPair = STREvaluationCall.harmonize(rawCall, rawTruth);
        final STREvaluationCall processedCall = processedPair.getLeft();
        final STREvaluationCall processedTruth = processedPair.getRight();

        STREvaluationClass evalClass = null;
        VariantContext outputVariantContext = null;
        STRContext strContext;

        if (!processedCall.isCalled) {
            if (processedTruth.isCalled) {
                strContext = composeContext(ref, processedTruth);
                if (strContext != null) {
                    evalClass = processedTruth.isHomRef() ? STREvaluationClass.TRUE_NEGATIVE : STREvaluationClass.FALSE_NEGATIVE;
                    outputVariantContext = truth;
                } else if (processedTruth.containsIndel()) {
                    evalClass = processedTruth.isHomRef() ? STREvaluationClass.TRUE_NEGATIVE : STREvaluationClass.FALSE_NEGATIVE;
                    outputVariantContext = truth;
                } else {
                    return null;
                }
            } else {
                return null;
            }
        } else if (!processedTruth.isCalled) {
                strContext = composeContext(ref, processedCall);
                if (strContext != null) {
                    evalClass = processedCall.isHomRef() ? STREvaluationClass.UNKNOWN_NEGATIVE : STREvaluationClass.UNKNOWN_POSITIVE;
                    outputVariantContext = truth;
                } else if (processedCall.containsIndel()) {
                    evalClass = processedTruth.isHomRef() ? STREvaluationClass.UNKNOWN_NEGATIVE : STREvaluationClass.UNKNOWN_POSITIVE;
                    outputVariantContext = truth;
                } else {
                    return null;
                }
        } else {
            strContext = composeContext(ref, processedTruth);
            if (processedTruth.isHomRef()) {
                if (processedCall.isHomRef()) {
                    if (strContext == null) {
                        return null;
                    }
                    evalClass = STREvaluationClass.TRUE_NEGATIVE;
                    outputVariantContext = call;
                } else {
                    if (strContext == null && !containsIndel(call)) {
                        return null;
                    }
                    evalClass = STREvaluationClass.FALSE_POSITIVE;
                    outputVariantContext = call;
                }
            } else if (processedCall.isHomRef()) {
                if (strContext == null && !containsIndel(truth)) {
                    return null;
                }
                evalClass = STREvaluationClass.FALSE_NEGATIVE;
                outputVariantContext = truth;
            } else if (processedTruth.isHet()) {
                if (strContext == null && !containsIndel(truth)) {
                    return null;
                }
                final VariantContextBuilder vcb = new VariantContextBuilder(call);
                vcb.attribute(TRUTH_ALLELE_TAG, processedTruth.hetAlternative());
                outputVariantContext = vcb.make();
                if (processedCall.isHet()) {
                    evalClass = processedCall.hetAlternative().equals(processedTruth.hetAlternative()) ? STREvaluationClass.TRUE_POSITIVE_HET : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_HET;
                } else if (processedCall.isHomVar()) {
                    evalClass = processedTruth.hetAlternative().equals(processedCall.homAlternative()) ? STREvaluationClass.DISCORDANT_POSITIVE_HET_TO_HOM : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_HET_TO_HOM;
                } else if (processedCall.isCompositeVar()) {
                    evalClass = processedCall.alleles.contains(processedTruth.hetAlternative()) ? STREvaluationClass.DISCORDANT_POSITIVE_HET_TO_COMPOSITE : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_HET_TO_COMPOSITE;
                }
            } else if (processedTruth.isHomVar()) {
                if (strContext == null && !containsIndel(truth)) {
                    return null;
                }
                final VariantContextBuilder vcb = new VariantContextBuilder(call);
                vcb.attribute(TRUTH_ALLELE_TAG, processedTruth.homAlternative());
                outputVariantContext = vcb.make();
                if (processedCall.isHomVar()) {
                    evalClass = processedCall.homAlternative().equals(processedTruth.homAlternative()) ? STREvaluationClass.TRUE_POSITIVE_HOM : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_HOM;
                } else if (processedCall.isHet()) {
                    evalClass = processedCall.hetAlternative().equals(processedTruth.homAlternative()) ? STREvaluationClass.DISCORDANT_POSITIVE_HOM_TO_HET : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_HOM_TO_HET;
                } else if (processedCall.isCompositeVar()) {
                    evalClass = processedCall.alleles.contains(processedTruth.homAlternative()) ? STREvaluationClass.DISCORDANT_POSITIVE_HOM_TO_COMPOSITE : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_HOM_TO_COMPOSITE;
                }
            } else {// if (truth.isCompositeVar()) {
                if (strContext == null && !containsIndel(truth)) {
                    return null;
                }
                final VariantContextBuilder vcb = new VariantContextBuilder(call);
                vcb.attribute(TRUTH_ALLELE_TAG, processedTruth.alternatives().stream().map(Object::toString).collect(Collectors.toList()));
                outputVariantContext = vcb.make();
                if (processedCall.isHomVar()) {
                    evalClass = processedTruth.alleles.contains(processedCall.homAlternative()) ? STREvaluationClass.DISCORDANT_POSITIVE_COMPOSITIVE_TO_HOM : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_COMPOSITE_TO_HOM;
                } else if (processedCall.isHet()) {
                    evalClass = processedTruth.alleles.contains(processedCall.hetAlternative()) ? STREvaluationClass.DISCORDANT_POSITIVE_COMPOSITE_TO_HET : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_COMPOSITE_TO_HOM;
                } else if (processedCall.isCompositeVar()) {
                    evalClass = processedTruth.alleles.stream().allMatch(a -> processedCall.alleles.contains(a)) ? STREvaluationClass.TRUE_POSITIVE_COMPOSITE : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_COMPOSITE;
                }
            }
        }

        return finalizeOutputVariantContext(ref, processedCall, processedTruth, evalClass, outputVariantContext, strContext);
    }

    private VariantContext finalizeOutputVariantContext(ReferenceContext ref, STREvaluationCall processedCall, STREvaluationCall processedTruth, STREvaluationClass evalClass, VariantContext outputVariantContext, STRContext strContext) {
        if (outputVariantContext != null && ref.getBase() != 'N') {
            final VariantContextBuilder vcb = new VariantContextBuilder(outputVariantContext);
            if (outputVariantContext.getStart() != ref.getInterval().getStart()) {
                if (!isNonVariantBlock(outputVariantContext)) {
                    throw new IllegalArgumentException("unexpeced condition");
                }
                vcb.chr(ref.getInterval().getContig());
                vcb.start(ref.getInterval().getStart());
                final Allele reference = Allele.create(ref.getBase(), true);
                vcb.alleles(Collections.singletonList(reference));
                vcb.genotypes(outputVariantContext.getGenotypes().stream()
                        .map(g -> new GenotypeBuilder(g).alleles(Collections.nCopies(g.getAlleles().size(), reference)).make())
                        .collect(Collectors.toList()));
            }
            if (processedTruth.confidence < minTruthQual) {
                vcb.filter(LOW_QUAL_FILTER);
            } else if (processedCall.confidence < minCallQual && evalClass != STREvaluationClass.FALSE_NEGATIVE) { // false negative are not excused by lack of confidence, other eval class
                vcb.filter(LOW_QUAL_FILTER);
            }


            if (evalClass != null) {
                vcb.attribute(EVALUATION_CLASS_TAG, evalClass.name());
            }
            if (strContext != null) {
                vcb.attribute(STR_UNIT_TAG, strContext.getAlleles().getRepeatUnitString());
                vcb.attribute(STR_UNIT_LENGTH_TAG, strContext.getAlleles().getRepeatUnitLength());
                vcb.attribute(STR_UNIT_MA_RC_TAG, strContext.getAlleles().getReferenceRepeatCount());
                strCounter.add(evalClass);
            } else {
                nonStrCounter.add(evalClass);
            }

            return vcb.make();
        } else {
            return null;
        }
    }

    private STRContext composeContext(final ReferenceContext ref, final STREvaluationCall truth) {
        final STRContext result = model.composeContext(ref);
        if (result == null) {
            return null;
        } else if (!result.getAlleles().compatibleAllele(truth.reference)) {
            return null;
        } else {
            return truth.alleles.stream().allMatch(a -> result.getAlleles().compatibleAllele(a)) ? result : null;
        }
    }

    private STRContext composeContext(final ReferenceContext ref, final VariantContext vc, final VariantContext truthVC) {
        STRContext result = model.composeContext(ref, vc);
        if (result == null) {
            result = model.composeContext(ref, truthVC);
        }
        return result;
    }

    private boolean isCalled(final VariantContext vc) {
        final List<Allele> alleles = vc.getGenotypes().get(0).getAlleles();
        if (alleles.size() != 2) {
            return false;
        }
        int normalCount = 0;
        for (final Allele allele : alleles) {
            if (allele.isCalled() && !allele.isSymbolic()) {
                normalCount++;
            }
        }
        return normalCount == 2;
    }

    private boolean containsIndel(final VariantContext truthVC) {
        boolean alleleContainsIndel = false;
        if (truthVC.getAlleles().get(0).length() > 1) {
            alleleContainsIndel = true;
        } else {
            for (final Allele allele : truthVC.getAlternateAlleles()) {
                if (allele.isSymbolic()) {
                    continue;
                } else if (allele.length() != 1) {
                    alleleContainsIndel = true;
                    break;
                }
            }
        }
        return alleleContainsIndel;
    }


    public static boolean isNonVariantBlock(final VariantContext vc) {
        if (vc == null) {
            return false;
        } else if (vc.getAlternateAlleles().size() == 0) {
            return true;
        } else if (vc.getAlternateAlleles().size() == 1) {
            return vc.getAlternateAllele(0).equals(Allele.NON_REF_ALLELE);
        } else {
            return false;
        }
    }






}
