package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;
import org.broadinstitute.hellbender.tools.copynumber.PostprocessGermlineCNVCalls;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LocatableCopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LocatableIntegerCopyNumber;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Helper class for {@link PostprocessGermlineCNVCalls} for single-sample postprocessing of
 * {@link GermlineCNVCaller} calls into genotyped intervals.
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlineCNVIntervalVariantComposer extends GermlineCNVAbstractVariantComposer {

    /* VCF FORMAT header keys */

    /**
     * Copy number maximum a posteriori value
     */
    static final String CN = "CN";

    /**
     * Copy number log posterior (in Phred-scale)
     */
    static final String CNLP = "CNLP";

    /**
     * Genotype call quality
     */
    static final String CNQ = "CNQ";

    private final IntegerCopyNumberStateCollection integerCopyNumberStateCollection;
    private final IntegerCopyNumberState refAutosomalCopyNumberState;
    private final Set<String> allosomalContigSet;

    /**
     * Constructor for {@link PostprocessGermlineCNVCalls} Postprocessor
     *
     * @param outputWriter variant context writer
     * @param integerCopyNumberStateCollection collection of copy-number states considered by post processor
     * @param sampleName sample name
     * @param refAutosomalCopyNumberState ref copy-number state on autosomal contigs
     * @param allosomalContigSet set of allosomal contigs (ref copy-number allele be chosen according to
     *                           given contig baseline copy-number states)
     */
    public GermlineCNVIntervalVariantComposer(final VariantContextWriter outputWriter,
                                              final String sampleName,
                                              final IntegerCopyNumberStateCollection integerCopyNumberStateCollection,
                                              final IntegerCopyNumberState refAutosomalCopyNumberState,
                                              final Set<String> allosomalContigSet) {
        super(outputWriter, sampleName);
        this.integerCopyNumberStateCollection = Utils.nonNull(integerCopyNumberStateCollection);
        Utils.validate(integerCopyNumberStateCollection.size() > 2,
                "There must be at least 3 copy number states present.");
        this.refAutosomalCopyNumberState = Utils.nonNull(refAutosomalCopyNumberState);
        this.allosomalContigSet = Utils.nonNull(allosomalContigSet);
    }

    @Override
    public void composeVariantContextHeader(final Set<VCFHeaderLine> vcfDefaultToolHeaderLines) {
        final VCFHeader result = new VCFHeader(Collections.emptySet(), Collections.singletonList(sampleName));

        /* add VCF version */
        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));

        /* add default tool header lines */
        vcfDefaultToolHeaderLines.forEach(result::addMetaDataLine);

        /* header lines related to genotype formatting */
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1,
                VCFHeaderLineType.Integer, "Genotype"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CN, 1,
                VCFHeaderLineType.Integer, "Copy number maximum a posteriori value"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CNLP, VCFHeaderLineCount.A,
                VCFHeaderLineType.Integer, "Copy number log posterior (in Phred-scale) rounded down"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CNQ, 1,
                VCFHeaderLineType.Integer, "Genotype call quality as the difference between" +
                " the best and second best phred-scaled log posterior scores"));

        /* INFO header lines */
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                VCFHeaderLineType.Integer, "End coordinate of the variant"));
        outputWriter.writeHeader(result);
    }

    /**
     * Write variant context fields given all posterior records from a single gCNV output shard
     *
     * @param copyNumberPosteriorRecordsShard copy-number posterior records from a gCNV output shard
     * @param baselineCopyNumberRecordsShard baseline copy-number records from a gCNV output shard
     */
    public void writeVariantContext(final List<LocatableCopyNumberPosteriorDistribution> copyNumberPosteriorRecordsShard,
                                    final List<LocatableIntegerCopyNumber> baselineCopyNumberRecordsShard) {
        Utils.validateArg(Utils.nonNull(copyNumberPosteriorRecordsShard).size() ==
                Utils.nonNull(baselineCopyNumberRecordsShard).size(), "Copy-number posterior records list has a " +
                "different size than the baseline copy-number records list.");
        final int numRecords = copyNumberPosteriorRecordsShard.size();
        for (int recordIndex = 0; recordIndex < numRecords; recordIndex++) {
            final LocatableCopyNumberPosteriorDistribution copyNumberPosteriorRecord =
                    copyNumberPosteriorRecordsShard.get(recordIndex);
            final LocatableIntegerCopyNumber baselineCopyNumberRecord =
                    baselineCopyNumberRecordsShard.get(recordIndex);
            Utils.validateArg(copyNumberPosteriorRecord.getInterval().equals(baselineCopyNumberRecord.getInterval()),
                    "Copy-number posterior record (%s) and baseline copy-number records (%s) must be given for " +
                            "the same genomic interval.");
            final VariantContext variantContext = composeVariantContext(copyNumberPosteriorRecord,
                    baselineCopyNumberRecord, VARIANT_PREFIX);
            outputWriter.add(variantContext);
        }
    }

    /**
     * Compose a variant context given a posterior record for an interval
     *
     * @param locatableCopyNumberPosteriorDistribution a posterior record to genotype
     * @param variantPrefix a variant prefix
     * @return composed variant context
     */
    @VisibleForTesting
    VariantContext composeVariantContext(
            final LocatableCopyNumberPosteriorDistribution locatableCopyNumberPosteriorDistribution,
            final LocatableIntegerCopyNumber locatableBaselineIntegerCopyNumber,
            final String variantPrefix) {
        final List<Integer> copyNumberPLVector =
                getCopyNumberPLVector(locatableCopyNumberPosteriorDistribution, integerCopyNumberStateCollection);
        final int copyNumberMAP =
                calculateMAPCopyNumberState(locatableCopyNumberPosteriorDistribution, integerCopyNumberStateCollection);
        final int GQ = calculateGenotypeQuality(locatableCopyNumberPosteriorDistribution, integerCopyNumberStateCollection);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.alleles(ALL_ALLELES);
        variantContextBuilder.chr(locatableCopyNumberPosteriorDistribution.getContig());
        variantContextBuilder.start(locatableCopyNumberPosteriorDistribution.getStart());
        variantContextBuilder.stop(locatableCopyNumberPosteriorDistribution.getEnd());
        variantContextBuilder.id(String.format(variantPrefix + "_%s_%d_%d",
                locatableCopyNumberPosteriorDistribution.getContig(),
                locatableCopyNumberPosteriorDistribution.getStart(),
                locatableCopyNumberPosteriorDistribution.getEnd()));

        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final String contig = locatableCopyNumberPosteriorDistribution.getContig();
        final IntegerCopyNumberState refCopyNumber = allosomalContigSet.contains(contig)
                ? locatableBaselineIntegerCopyNumber.getIntegerCopyNumberState()
                : refAutosomalCopyNumberState;

        final Allele allele;
        if (copyNumberMAP > refCopyNumber.getCopyNumber()) {
            allele = DUP_ALLELE;
        } else if (copyNumberMAP < refCopyNumber.getCopyNumber()) {
            allele = DEL_ALLELE;
        } else {
            allele = REF_ALLELE;
        }
        genotypeBuilder.alleles(Collections.singletonList(allele));
        genotypeBuilder.attribute(CN, copyNumberMAP);
        genotypeBuilder.attribute(CNLP, copyNumberPLVector);
        genotypeBuilder.attribute(CNQ, GQ);
        final Genotype genotype = genotypeBuilder.make();

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
                .mapToInt(GermlineCNVIntervalVariantComposer::convertLogProbabilityToPhredScore)
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
