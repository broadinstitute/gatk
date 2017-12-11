package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nullable;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Helper class for single sample gCNV postprocessing
 */
public class GCNVPostProcessor {

    /**
     * VCF miscellaneous header keys
     */
    public static final String SAMPLE_NAME = "sample_name";

    /**
     * VCF FORMAT header keys
     */
    public static final String CN = "CN";
    public static final String CNLP = "CNLP";
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
    private CopyNumberPosteriorLocatableRecord lastPosteriorRecord = null;

    /**
     * gCNV postprocessor constructor
     *
     * @param outputWriter variant context writer
     * @param integerCopyNumberStateCollection collection of copy number states considered by post processor
     * @param sampleName sample name
     * @param samSequenceDictionary sequence dictionary used for interval order validation
     */
    GCNVPostProcessor(final VariantContextWriter outputWriter,
                                final IntegerCopyNumberStateCollection integerCopyNumberStateCollection,
                                final String sampleName,
                                final SAMSequenceDictionary samSequenceDictionary) {
        //TODO pass a list of intervals and add progress meter
        //ProgressMeter progressMeter = new ProgressMeter(1.0);
        //progressMeter.start();
        this.outputWriter = Utils.nonNull(outputWriter);
        this.integerCopyNumberStateCollection = Utils.nonNull(integerCopyNumberStateCollection);
        Utils.validate(integerCopyNumberStateCollection.size() > 2, "There must be at least 3 copy number states present");
        this.alleles = integerCopyNumberStateCollection.getAlleles();
        this.sampleName = sampleName;
        this.samSequenceDictionary = Utils.nonNull(samSequenceDictionary);
    }

    /**
     * Compose the header of the variant context
     *
     * @param commandLine command line of the VCF generation tool
     */
    void composeVariantContextHeader(@Nullable final String commandLine) {
        final VCFHeader result = new VCFHeader(Collections.emptySet(), Arrays.asList(sampleName));

        /* add VCF version */
        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));

        /* add command line */
        if (commandLine != null) {
            result.addMetaDataLine(new VCFHeaderLine("command", commandLine));
        }

        /* add miscellaneous header keys */
        result.addMetaDataLine(new VCFHeaderLine(SAMPLE_NAME, sampleName));

        /* header lines related to genotype formatting */
        result.addMetaDataLine(new VCFFormatHeaderLine(CN, 1,
                VCFHeaderLineType.Integer, "Copy number maximum a posteriori value"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CNLP, VCFHeaderLineCount.A,
                VCFHeaderLineType.Integer, "Copy number log posterior (in Phred-scale) rounded down"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CNQ, 1,
                VCFHeaderLineType.Integer, "Genotype call quality as the difference between the best and second best PL"));

        /* INFO header lines*/
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                VCFHeaderLineType.Integer, "End coordinate of this variant"));
        outputWriter.writeHeader(result);
    }

    /**
     * Write variant context fields given all posterior records from a single gCNV output chunk
     *
     * @param copyNumberPosteriorLocatableRecordsChunk records from a gCNV output chunk
     * @param variantPrefix variant id prefix
     */
    void writeChunkedVariantContext(final List<CopyNumberPosteriorLocatableRecord> copyNumberPosteriorLocatableRecordsChunk,
                                    final String variantPrefix) {
        for (CopyNumberPosteriorLocatableRecord posteriorRecord: copyNumberPosteriorLocatableRecordsChunk) {
            final VariantContext variantContext = composeVariantContext(posteriorRecord, variantPrefix);
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
     * @param copyNumberPosteriorLocatableRecord a posterior record to genotype
     * @param variantPrefix a variant prefix
     * @return composed variant context
     */
    private VariantContext composeVariantContext(final CopyNumberPosteriorLocatableRecord copyNumberPosteriorLocatableRecord,
                                                 final String variantPrefix) {
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.alleles(alleles);
        variantContextBuilder.chr(copyNumberPosteriorLocatableRecord.getContig());
        variantContextBuilder.start(copyNumberPosteriorLocatableRecord.getStart());
        variantContextBuilder.stop(copyNumberPosteriorLocatableRecord.getEnd());
        variantContextBuilder.id(String.format(variantPrefix + "_%s_%d_%d",
                copyNumberPosteriorLocatableRecord.getContig(),
                copyNumberPosteriorLocatableRecord.getStart(),
                copyNumberPosteriorLocatableRecord.getEnd()));
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final List<Integer> copyNumberPLVector = getCopyNumberPLVector(copyNumberPosteriorLocatableRecord, integerCopyNumberStateCollection);
        final int copyNumberMAP = calculateMAPCopyNumberState(copyNumberPosteriorLocatableRecord, integerCopyNumberStateCollection);
        final int GQ = calculateGenotypeQuality(copyNumberPosteriorLocatableRecord, integerCopyNumberStateCollection);
        genotypeBuilder.attribute(CN, copyNumberMAP);
        genotypeBuilder.attribute(CNLP, copyNumberPLVector);
        genotypeBuilder.attribute(CNQ, GQ);

        final Genotype genotype = genotypeBuilder.make();

        //Add allele information to the variant context
        variantContextBuilder.attribute(VCFConstants.END_KEY, copyNumberPosteriorLocatableRecord.getEnd());

        variantContextBuilder.genotypes(genotype);
        return variantContextBuilder.make();
    }

    @VisibleForTesting
    static List<Integer> getCopyNumberPLVector(final CopyNumberPosteriorLocatableRecord copyNumberPosteriorLocatableRecord,
                                                     final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        final List<Integer> unnormalizedPL = integerCopyNumberStateCollection.getCopyNumberStates().stream()
                .mapToDouble(state -> copyNumberPosteriorLocatableRecord.getCopyNumberPosteriors().getCopyNumberPosterior(state))
                .mapToInt(GCNVPostProcessor::convertLogPosteriorToPL).boxed().collect(Collectors.toList());
        final Optional<Integer> smallestPL = unnormalizedPL.stream().min(Integer::compare);
        return unnormalizedPL.stream().mapToInt(p -> p - smallestPL.get()).boxed().collect(Collectors.toList());
    }

    @VisibleForTesting
    static int calculateMAPCopyNumberState(final CopyNumberPosteriorLocatableRecord copyNumberPosteriorLocatableRecord,
                                      final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        final Optional<IntegerCopyNumberState> copyNumberStateMAP = integerCopyNumberStateCollection.getCopyNumberStates()
                .stream().max((a1, a2) -> Double.compare(copyNumberPosteriorLocatableRecord.getCopyNumberPosteriors().getCopyNumberPosterior(a1),
                copyNumberPosteriorLocatableRecord.getCopyNumberPosteriors().getCopyNumberPosterior(a2)));
        return copyNumberStateMAP.get().getCopyNumber();
    }

    @VisibleForTesting
    static int calculateGenotypeQuality(final CopyNumberPosteriorLocatableRecord copyNumberPosteriorLocatableRecord,
                                        final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        final double[] sortedPosteriors = integerCopyNumberStateCollection.getCopyNumberStates().stream()
                .mapToDouble(state -> copyNumberPosteriorLocatableRecord.getCopyNumberPosteriors().getCopyNumberPosterior(state))
                .sorted().toArray();
        final int bestQuality = (int) FastMath.floor(convertLogPosteriorToPL(sortedPosteriors[sortedPosteriors.length - 1]));
        final int secondBestQuality = (int) FastMath.floor(convertLogPosteriorToPL(sortedPosteriors[sortedPosteriors.length - 2]));
        return secondBestQuality - bestQuality;
    }

    private static int convertLogPosteriorToPL(final double posteriorProbInLogSpace) {
        return (int) FastMath.floor(-10.0 * posteriorProbInLogSpace * MathUtils.LOG10_OF_E);
    }
}
