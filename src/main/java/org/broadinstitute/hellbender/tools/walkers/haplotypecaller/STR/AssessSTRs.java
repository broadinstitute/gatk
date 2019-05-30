package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import jdk.nashorn.internal.ir.annotations.Reference;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.GATKVCFUtils;
import org.broadinstitute.hellbender.engine.walkers.*;
import org.broadinstitute.hellbender.utils.collections.Pair;
import org.broadinstitute.hellbender.utils.commandline.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by valentin on 11/3/16.
 */
public class AssessSTRs extends MultiVariantWalker {


    public MultiVariantInputArgumentCollection]

    @Override
    protected MultiVariantInputArgumentCollection getMultiVariantInputArgumentCollection() {
        return new MultiVariantInputArgumentCollection() {
            @Override
            public List<String> getDrivingVariantPaths() {
                if (truth != null) {
                    return Arrays.asList(input.getFeaturePath(), truth.getFeaturePath());
                } else {
                    return Collections.singletonList(input.getFeaturePath());
                }
            }
        };
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Argument(doc = "input VCF", shortName = "V", fullName = "variant")
    public FeatureInput<VariantContext> input;

    @Argument(doc = "truth VCF calls", shortName = "truth", fullName = "truthVariant", optional = true)
    public FeatureInput<VariantContext> truth;

    @ArgumentCollection
    public STRModel model = new STRModel();

    @Argument(doc = "output Context table", shortName = "O", fullName = "output")
    public VariantContextWriter output;

    @Argument(doc = "output summary", shortName = "summary", fullName = "summaryOutput")
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

    @Override
    public void onTraversalStart() {

        final List<VCFHeader> inputHeaders = features.getAllVariantHeaders();

        final List<String> samples = inputHeaders.stream().filter(Objects::nonNull)
                .flatMap(h -> h.getSampleNamesInOrder().stream())
                .sorted()
                .distinct()
                .collect(Collectors.toList());
        final VCFHeader header = new VCFHeader(new HashSet<>(), samples);
        header.addMetaDataLine(new VCFInfoHeaderLine(EVALUATION_CLASS_TAG, 1, VCFHeaderLineType.String, "Evaluation class"));
        header.addMetaDataLine(new VCFInfoHeaderLine(STR_UNIT_TAG, 1, VCFHeaderLineType.String, "STR Repeat Unit"));
        header.addMetaDataLine(new VCFInfoHeaderLine(STR_UNIT_MA_RC_TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "STR Reference Repeat Count"));
        header.addMetaDataLine(new VCFInfoHeaderLine(STR_UNIT_LENGTH_TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "STR Repeat Unit Length"));
        header.addMetaDataLine(new VCFInfoHeaderLine(TRUTH_ALLELE_TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "STR Type"));
        header.addMetaDataLine(new VCFFilterHeaderLine(LOW_QUAL_FILTER, "LowQual"));
        inputHeaders.stream().filter(Objects::nonNull)
                .flatMap(h ->
                    Stream.concat(h.getInfoHeaderLines().stream(),
                            Stream.concat(h.getFormatHeaderLines().stream(), h.getFormatHeaderLines().stream())))
                .forEach(header::addMetaDataLine);
        output.writeHeader(header);
    }

    public VariantContext eval(final VariantContext input, final VariantContext truth, final ReferenceContext ref) {
        final STREvaluationCall rawCall = STREvaluationCall.fromVariantContext(input, ref, false);
        final STREvaluationCall rawTruth = STREvaluationCall.fromVariantContext(truth, ref, missingTruthMeansHomeRef);
        if (rawCall == null || rawTruth == null) {
            logger.debug("Funny situation at " + ref.getInterval();
            return null;
        }
        // make sure they have the same reference allele sequence:
        final Pair<STREvaluationCall, STREvaluationCall> processedCalls = STREvaluationCall.harmonize(rawCall, rawTruth);
        final STREvaluationCall actualCall = processedCalls.getLeft();
        final STREvaluationCall expectedCall = processedCalls.getRight();

        STREvaluationClass evalClass = null;
        VariantContext outputVariantContext = null;
        STRContext strContext;
        final STREvaluationCall bestCall = expectedCall.isCalled
                ? expectedCall
                : (actualCall.isCalled ? actualCall : null);

        if (!bestCall.containsIndel()) {
            return null;
        }

        strContext = composeContext(ref, bestCall);
        if (strContext == null) {
            return null;
        }

        if (!actualCall.isCalled) {
            evalClass = expectedCall.isHomRef() ? STREvaluationClass.TRUE_NEGATIVE : STREvaluationClass.FALSE_NEGATIVE;
            outputVariantContext = truth;
        } else if (!expectedCall.isCalled) {
            evalClass = actualCall.isHomRef() ? STREvaluationClass.UNKNOWN_NEGATIVE : STREvaluationClass.UNKNOWN_POSITIVE;
            outputVariantContext = truth;
        } else {
            if (expectedCall.isHomRef()) {
                if (call.isHomRef()) {
                    if (strContext == null) {
                        return null;
                    }
                    evalClass = STREvaluationClass.TRUE_NEGATIVE;
                    outputVariantContext = vc;
                } else {
                    if (strContext == null && !containsIndel(vc)) {
                        return null;
                    }
                    evalClass = STREvaluationClass.FALSE_POSITIVE;
                    outputVariantContext = vc;
                }
            } else if (call.isHomRef()) {
                if (strContext == null && !containsIndel(truthVC)) {
                    return null;
                }
                evalClass = STREvaluationClass.FALSE_NEGATIVE;
                outputVariantContext = truthVC;
            } else if (truth.isHet()) {
                if (strContext == null && !containsIndel(truthVC)) {
                    return null;
                }
                final VariantContextBuilder vcb = new VariantContextBuilder(vc);
                vcb.attribute(TRUTH_ALLELE_TAG, truth.hetAlternative());
                outputVariantContext = vcb.make();
                if (call.isHet()) {
                    evalClass = call.hetAlternative().equals(truth.hetAlternative()) ? STREvaluationClass.TRUE_POSITIVE_HET : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_HET;
                } else if (call.isHomVar()) {
                    evalClass = truth.hetAlternative().equals(call.homAlternative()) ? STREvaluationClass.DISCORDANT_POSITIVE_HET_TO_HOM : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_HET_TO_HOM;
                } else if (call.isCompositeVar()) {
                    evalClass = call.alleles.contains(truth.hetAlternative()) ? STREvaluationClass.DISCORDANT_POSITIVE_HET_TO_COMPOSITE : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_HET_TO_COMPOSITE;
                }
            } else if (truth.isHomVar()) {
                if (strContext == null && !containsIndel(truthVC)) {
                    return null;
                }
                final VariantContextBuilder vcb = new VariantContextBuilder(vc);
                vcb.attribute(TRUTH_ALLELE_TAG, truth.homAlternative());
                outputVariantContext = vcb.make();
                if (call.isHomVar()) {
                    evalClass = call.homAlternative().equals(truth.homAlternative()) ? STREvaluationClass.TRUE_POSITIVE_HOM : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_HOM;
                } else if (call.isHet()) {
                    evalClass = call.hetAlternative().equals(truth.homAlternative()) ? STREvaluationClass.DISCORDANT_POSITIVE_HOM_TO_HET : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_HOM_TO_HET;
                } else if (call.isCompositeVar()) {
                    evalClass = call.alleles.contains(truth.homAlternative()) ? STREvaluationClass.DISCORDANT_POSITIVE_HOM_TO_COMPOSITE : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_HOM_TO_COMPOSITE;
                }
            } else {// if (truth.isCompositeVar()) {
                if (strContext == null && !containsIndel(truthVC)) {
                    return null;
                }
                final VariantContextBuilder vcb = new VariantContextBuilder(vc);
                vcb.attribute(TRUTH_ALLELE_TAG, truth.alternatives().stream().map(Object::toString).collect(Collectors.toList()));
                outputVariantContext = vcb.make();
                if (call.isHomVar()) {
                    evalClass = truth.alleles.contains(call.homAlternative()) ? STREvaluationClass.DISCORDANT_POSITIVE_COMPOSITIVE_TO_HOM : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_COMPOSITE_TO_HOM;
                } else if (call.isHet()) {
                    evalClass = truth.alleles.contains(call.hetAlternative()) ? STREvaluationClass.DISCORDANT_POSITIVE_COMPOSITE_TO_HET : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_COMPOSITE_TO_HOM;
                } else if (call.isCompositeVar()) {
                    evalClass = truth.alleles.stream().allMatch(a -> call.alleles.contains(a)) ? STREvaluationClass.TRUE_POSITIVE_COMPOSITE : STREvaluationClass.DISCORDANT_POSITIVE_ALLELE_COMPOSITE;
                }
            }
        }

        if (outputVariantContext != null) {
            final VariantContextBuilder vcb = new VariantContextBuilder(outputVariantContext);
            if (outputVariantContext.getStart() != ref.getLocus().getStart()) {
                if (!isNonVariantBlock(outputVariantContext)) {
                    throw new IllegalArgumentException("unexpeced condition");
                }
                vcb.chr(ref.getLocus().getContig());
                vcb.start(ref.getLocus().getStart());
                final Allele reference = Allele.create(ref.getBase(), true);
                vcb.alleles(Collections.singletonList(reference));
                vcb.genotypes(outputVariantContext.getGenotypes().stream()
                        .map(g -> new GenotypeBuilder(g).alleles(Collections.nCopies(g.getAlleles().size(), reference)).make())
                        .collect(Collectors.toList()));
            }
            final boolean filtered;
            if (truth.confidence < minTruthQual) {
                filtered = true;
                vcb.filter(LOW_QUAL_FILTER);
            } else if (call.confidence < minCallQual && evalClass != STREvaluationClass.FALSE_NEGATIVE) { // false negative are not excused by lack of confidence, other eval class
                filtered = true;
                vcb.filter(LOW_QUAL_FILTER);
            } else {
                filtered = false;
            }

            if (vcb.getAlleles().get(0).getBases()[0] == 'N') {
                return null;
            }
            if (evalClass != null) {
                vcb.attribute(EVALUATION_CLASS_TAG, evalClass.name());
            }
            if (strContext != null) {
                vcb.attribute(STR_UNIT_TAG, strContext.getAlleles().getRepeatUnitString());
                vcb.attribute(STR_UNIT_LENGTH_TAG, strContext.getAlleles().getRepeatUnitLength());
                vcb.attribute(STR_UNIT_MA_RC_TAG, strContext.getAlleles().getReferenceRepeatCount());
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

    private VariantContext findTruthVC(final List<VariantContext> truths, final VariantContext vc) {
        if (truths == null || truths.isEmpty()) {
            return null;
        }
        for (final VariantContext truth : truths) {
            if (vc.getStart() == truth.getStart()) {
                return truth;
            }
        }
        for (final VariantContext truth : truths) {
            if (truth.getAlternateAlleles().isEmpty() || (truth.getAlternateAlleles().size() == 1 && truth.getAlternateAlleles().get(0).equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE))) {
                return truth;
            }
        }
        return null;
    }

    private int truthGQ(final VariantContext truth) {
        if (truth == null) {
            return 0;
        } else if (truth.getGenotype(0).hasGQ()) {
            return truth.getGenotype(0).getGQ();
        } else {
            return (int) Math.floor(truth.getPhredScaledQual());
        }
    }

    public static boolean isNonVariantBlock(final VariantContext vc) {
        if (vc == null) {
            return false;
        } else if (vc.getAlternateAlleles().size() == 0) {
            return true;
        } else if (vc.getAlternateAlleles().size() == 1) {
            return vc.getAlternateAllele(0).equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        } else {
            return false;
        }
    }


    @Override
    public Integer reduceInit() {
        return 0;
    }
    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {

    }

    @Override
    public Integer reduce(final VariantContext value, final Integer count) {
        if (value == null) {
            return count;
        }
        final String valueContig = value.getContig();
        final long valueStart = value.getStart();
        if (valueContig.equals(lastOutputContig) && valueStart <= lastOutputPosition) {
            return count;
        }
        lastOutputContig = valueContig;
        final long valueLength = value.hasAttribute(STR_UNIT_TAG)
                ? value.getAttributeAsInt(STR_UNIT_LENGTH_TAG, 0) * value.getAttributeAsInt(STR_UNIT_MA_RC_TAG, 0)
                : value.getReference().length() - 1;
        lastOutputPosition = valueStart + valueLength;
        if (value.hasAttribute(EVALUATION_CLASS_TAG)) {
            final STREvaluationClass clazz = STREvaluationClass.valueOf(value.getAttributeAsString(EVALUATION_CLASS_TAG, null));
            final STREvaluationCounter counter = value.hasAttribute(STR_UNIT_TAG) ? strCounter : nonStrCounter;
            counter.add(clazz);
        }
        output.add(value);
        return count + 1;
    }

    @Override
    public void onTraversalDone(final Integer count) {
        super.onTraversalDone(count);
        logger.info("Number of STR context found: " + count);
        output.close();

        try (final PrintWriter summaryWriter = new PrintWriter(new FileWriter(summaryOutput))) {
            summaryWriter.println(String.join("\t", "EVAL_CLASS", "STR", "INDEL"));
            for (final STREvaluationClass evalClass : STREvaluationClass.values()) {
                summaryWriter.println(String.join("\t", evalClass.name(), "" + strCounter.get(evalClass), "" + nonStrCounter.get(evalClass)));
            }
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(summaryOutput, ex);
        }
    }


}
