package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.PostProcessGermlineCNVCalls;
import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LocatableCopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Helper class for {@link PostProcessGermlineCNVCalls} for single sample postprocessing
 * of {@link GermlineCNVCaller} calls.
 *
 * This class takes in a {@link IntegerCopyNumberStateCollection} a sample name, and a sequence dictionary for
 * validating order of the processed calls. It is used to write a VCF header and genotyping
 * {@link LocatableCopyNumberPosteriorDistribution} records.
 */
public final class GermlineCNVPostProcessingEngine {

    private static final String VARIANT_PREFIX = "CNV";

    //VCF FORMAT header keys
    /**
     * Copy number maximum a posteriori value
     */
    public static final String CN = "CN";
    /**
     * Copy number log posterior (in Phred-scale)
     */
    public static final String CNLP = "CNLP";
    /**
     * Genotype call quality
     */
    public static final String CNQ = "CNQ";

    private final List<Allele> alleles;
    private final IntegerCopyNumberStateCollection integerCopyNumberStateCollection;
    private final String sampleName;
    private final VariantContextWriter outputWriter;
    private final SAMSequenceDictionary samSequenceDictionary;

    /**
     * We store the previously analyzed record for validating that records are ordered according to
     * the sequence dictionary
     */
    private LocatableCopyNumberPosteriorDistribution lastPosteriorRecord = null;

    /**
     * Constructor for {@link PostProcessGermlineCNVCalls} postprocessor
     *
     * @param outputWriter variant context writer
     * @param integerCopyNumberStateCollection collection of copy number states considered by post processor
     * @param sampleName sample name
     * @param samSequenceDictionary sequence dictionary used for interval order validation
     */
    public GermlineCNVPostProcessingEngine(final VariantContextWriter outputWriter,
                                           final IntegerCopyNumberStateCollection integerCopyNumberStateCollection,
                                           final String sampleName,
                                           final SAMSequenceDictionary samSequenceDictionary) {
        this.outputWriter = Utils.nonNull(outputWriter);
        this.integerCopyNumberStateCollection = Utils.nonNull(integerCopyNumberStateCollection);
        Utils.validate(integerCopyNumberStateCollection.size() > 2,
                "There must be at least 3 copy number states present.");
        this.alleles = integerCopyNumberStateCollection.getAlleles();
        this.sampleName = Utils.nonEmpty(sampleName);
        this.samSequenceDictionary = Utils.nonNull(samSequenceDictionary);
    }

    /**
     * Compose the header of the variant context
     *
     * @param vcfDefaultToolHeaderLines default header lines of the VCF generation tool
     */
    public void composeVariantContextHeader(final Set<VCFHeaderLine> vcfDefaultToolHeaderLines) {
        final VCFHeader result = new VCFHeader(Collections.emptySet(), Collections.singletonList(sampleName));

        /* add VCF version */
        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));

        /* add default tool header lines */
        vcfDefaultToolHeaderLines.forEach(line -> result.addMetaDataLine(line));

        /* header lines related to genotype formatting */
        result.addMetaDataLine(new VCFFormatHeaderLine(CN, 1,
                VCFHeaderLineType.Integer, "Copy number maximum a posteriori value"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CNLP, VCFHeaderLineCount.A,
                VCFHeaderLineType.Integer, "Copy number log posterior (in Phred-scale) rounded down"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CNQ, 1,
                VCFHeaderLineType.Integer, "Genotype call quality as the difference between" +
                " the best and second best phred-scaled log posterior scores"));

        /* INFO header lines */
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                VCFHeaderLineType.Integer, "End coordinate of this variant"));
        outputWriter.writeHeader(result);
    }

    /**
     * Write variant context fields given all posterior records from a single gCNV output chunk
     *
     * @param locatableRecordsChunk records from a gCNV output chunk
     */
    public void writeChunkedVariantContext(final List<LocatableCopyNumberPosteriorDistribution> locatableRecordsChunk) {
        for (final LocatableCopyNumberPosteriorDistribution posteriorRecord: locatableRecordsChunk) {
            final VariantContext variantContext = composeVariantContext(posteriorRecord, VARIANT_PREFIX);
            outputWriter.add(variantContext);
            //check that analyzed posterior records are ordered according to the provided sequence dictionary
            if (lastPosteriorRecord != null) {
                if (IntervalUtils.compareLocatables(lastPosteriorRecord, posteriorRecord, samSequenceDictionary) > 0) {
                    throw new UserException.BadInput("Intervals in the chunk directories are not" +
                            " ordered according to the sequence dictionary");
                }
            }
            lastPosteriorRecord = posteriorRecord;
        }
    }

    /**
     * Compose a variant context given a posterior record for an interval
     *
     * @param locatableCopyNumberPosteriorDistribution a posterior record to genotype
     * @param variantPrefix a variant prefix
     * @return composed variant context
     */
    private VariantContext composeVariantContext(
            final LocatableCopyNumberPosteriorDistribution locatableCopyNumberPosteriorDistribution,
            final String variantPrefix) {
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.alleles(alleles);
        variantContextBuilder.chr(locatableCopyNumberPosteriorDistribution.getContig());
        variantContextBuilder.start(locatableCopyNumberPosteriorDistribution.getStart());
        variantContextBuilder.stop(locatableCopyNumberPosteriorDistribution.getEnd());
        variantContextBuilder.id(String.format(variantPrefix + "_%s_%d_%d",
                locatableCopyNumberPosteriorDistribution.getContig(),
                locatableCopyNumberPosteriorDistribution.getStart(),
                locatableCopyNumberPosteriorDistribution.getEnd()));
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final List<Integer> copyNumberPLVector =
                getCopyNumberPLVector(locatableCopyNumberPosteriorDistribution, integerCopyNumberStateCollection);
        final int copyNumberMAP =
                calculateMAPCopyNumberState(locatableCopyNumberPosteriorDistribution, integerCopyNumberStateCollection);
        final int GQ = calculateGenotypeQuality(locatableCopyNumberPosteriorDistribution, integerCopyNumberStateCollection);
        genotypeBuilder.attribute(CN, copyNumberMAP);
        genotypeBuilder.attribute(CNLP, copyNumberPLVector);
        genotypeBuilder.attribute(CNQ, GQ);

        final Genotype genotype = genotypeBuilder.make();

        //Add allele information to the variant context
        variantContextBuilder.attribute(VCFConstants.END_KEY, locatableCopyNumberPosteriorDistribution.getEnd());

        variantContextBuilder.genotypes(genotype);
        return variantContextBuilder.make();
    }

    /**
     * Compute copy number log posterior scores
     *
     * @param locatableCopyNumberPosteriorDistribution copy number posterior locatable record
     * @param integerCopyNumberStateCollection copy number state collection
     * @return a list of copy number log posterior scores
     */
    @VisibleForTesting
    static List<Integer> getCopyNumberPLVector(
            final LocatableCopyNumberPosteriorDistribution locatableCopyNumberPosteriorDistribution,
            final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        final Double largestLogProb = integerCopyNumberStateCollection.getCopyNumberStates().stream()
                .mapToDouble(state -> locatableCopyNumberPosteriorDistribution.getCopyNumberPosteriors().getCopyNumberPosterior(state))
                .max().getAsDouble();
        final List<Integer> unnormalizedPL = integerCopyNumberStateCollection.getCopyNumberStates().stream()
                .mapToDouble(state -> locatableCopyNumberPosteriorDistribution.getCopyNumberPosteriors()
                        .getCopyNumberPosterior(state) - largestLogProb)
                .mapToInt(GermlineCNVPostProcessingEngine::convertLogProbabilityToPhredScore)
                .boxed()
                .collect(Collectors.toList());
        return unnormalizedPL;
    }

    /**
     * Compute MAP copy number state
     *
     * @param locatableCopyNumberPosteriorDistribution copy number posterior locatable record
     * @param integerCopyNumberStateCollection copy number state collection
     * @return MAP copy number state
     */
    @VisibleForTesting
    static int calculateMAPCopyNumberState(
            final LocatableCopyNumberPosteriorDistribution locatableCopyNumberPosteriorDistribution,
            final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        final Optional<IntegerCopyNumberState> copyNumberStateMAP = integerCopyNumberStateCollection
                .getCopyNumberStates().stream()
                .max(Comparator.comparingDouble(state ->
                        locatableCopyNumberPosteriorDistribution.getCopyNumberPosteriors().getCopyNumberPosterior(state)));
        return copyNumberStateMAP.get().getCopyNumber();
    }

    /**
     * Calculate genotype quality score defined as a difference of second smallest and smallest copy number log
     * posterior scores in phred-scale
     *
     * @param locatableCopyNumberPosteriorDistribution copy number posterior locatable record
     * @param integerCopyNumberStateCollection copy number state collection
     * @return genotype quality score
     */
    @VisibleForTesting
    static int calculateGenotypeQuality(final LocatableCopyNumberPosteriorDistribution locatableCopyNumberPosteriorDistribution,
                                        final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        final List<Integer> sortedPosteriors = getCopyNumberPLVector(locatableCopyNumberPosteriorDistribution,
                integerCopyNumberStateCollection).stream().sorted().collect(Collectors.toList());
        //sanity check
        Utils.validate(sortedPosteriors.get(0) == 0, "Something went wrong. Smallest copy" +
                " number posterior score must be 0");
        return sortedPosteriors.get(1);
    }

    private static int convertLogProbabilityToPhredScore(final double posteriorProbInLogSpace) {
        return (int) FastMath.floor(-10.0 * posteriorProbInLogSpace * MathUtils.LOG10_OF_E);
    }
}
