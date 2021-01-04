package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntervalCopyNumberGenotypingData;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class DepthEvidenceAggregator {

    private final List<VCFFileReader> posteriorsReaders;
    private final SAMSequenceDictionary dictionary;
    private final List<String> samples;
    private final List<IntegerCopyNumberState> copyStates;
    private final int numCopyStates;

    private String currentContig;
    private List<IntervalTree<Map<String,double[]>>> currentPosteriorsTreeList;

    public DepthEvidenceAggregator(final List<VCFFileReader> posteriorsReaders,
                                   final List<String> samples,
                                   final SAMSequenceDictionary dictionary) {
        Utils.nonNull(posteriorsReaders);
        Utils.nonNull(samples);
        Utils.nonNull(dictionary);
        this.posteriorsReaders = posteriorsReaders;
        this.samples = samples;
        this.dictionary = dictionary;
        final VariantContext exampleVariant = posteriorsReaders.get(0).iterator().next();
        copyStates = getCopyNumberStates(exampleVariant);
        numCopyStates = copyStates.size();
        validateCopyStates();
    }

    public SVCallRecordDepthPosterior apply(final SVCallRecord call) {
        Utils.nonNull(call);
        if (!call.getContigA().equals(currentContig)) {
            currentContig = call.getContigA();
            currentPosteriorsTreeList = posteriorsReaders.stream().map(this::getCurrentPosteriorsTree).collect(Collectors.toList());
        }
        if (!(call.getType().equals(StructuralVariantType.DEL)  || call.getType().equals(StructuralVariantType.DUP)
                || call.getType().equals(StructuralVariantType.BND))) {
            return null;
        }
        if (!call.getContigA().equals(call.getContigB())) {
            return null;
        }
        return SVCallRecordDepthPosterior.create(copyStates, samples, call, currentPosteriorsTreeList);
    }

    private void validateCopyStates() {
        for (int i = 1; i < posteriorsReaders.size(); i++) {
            final List<IntegerCopyNumberState> otherCopyStates = getCopyNumberStates(posteriorsReaders.get(i).iterator().next());
            if (!copyStates.equals(otherCopyStates)) {
                throw new UserException.BadInput("CNV VCFs do not contain identical copy number states.");
            }
        }
    }

    private List<IntegerCopyNumberState> getCopyNumberStates(final VariantContext variant) {
        final SimpleInterval interval = new SimpleInterval(variant.getContig(), variant.getStart(), variant.getEnd());
        return parseCnvGenotype(variant.getGenotype(0), interval)
                .getCopyNumberPosteriorDistribution().getIntegerCopyNumberStateList();
    }

    private Iterator<VariantContext> queryContig(final String contig, final VCFFileReader reader) {
        final SAMSequenceRecord contigRecord = dictionary.getSequence(contig);
        if (contigRecord == null) {
            throw new UserException.MissingContigInSequenceDictionary(contig, dictionary);
        }
        return reader.query(contig, 1, contigRecord.getSequenceLength());
    }

    private IntervalTree<Map<String,double[]>> getCurrentPosteriorsTree(final VCFFileReader reader) {
        final Iterator<VariantContext> posteriorsIter = queryContig(currentContig, reader);
        final IntervalTree<Map<String,double[]>> tree = new IntervalTree<>();
        while (posteriorsIter.hasNext()) {
            final VariantContext variant = posteriorsIter.next();
            final SimpleInterval interval = new SimpleInterval(variant.getContig(), variant.getStart(), variant.getEnd() - 1);
            final Map<String,double[]> samplePosteriorMap = new HashMap<>(SVUtils.hashMapCapacity(samples.size()));
            for (final Genotype genotype : variant.getGenotypes()) {
                final IntervalCopyNumberGenotypingData data = parseCnvGenotype(genotype, interval);
                final double[] posteriors = new double[numCopyStates];
                int i = 0;
                for (final IntegerCopyNumberState state : copyStates) {
                    posteriors[i] = data.getCopyNumberPosteriorDistribution().getCopyNumberPosterior(state);
                    i++;
                }
                samplePosteriorMap.put(genotype.getSampleName(), posteriors);
            }
            tree.put(interval.getStart(), interval.getEnd(), samplePosteriorMap);
        }
        return tree;
    }

    private IntervalCopyNumberGenotypingData parseCnvGenotype(final Genotype genotype, final SimpleInterval interval) {
        final String posteriorString = (String) genotype.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY);
        Utils.nonNull(posteriorString, "Missing required field " + GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY);
        final String[] posteriorStrings = posteriorString.split(",");

        //Posteriors reported as integer phred-scaled likelihoods and need to be renormalized
        final double[] approximatePosteriors = new double[posteriorStrings.length];
        double total = 0;
        for (int i = 0; i < posteriorStrings.length; i++) {
            final int likelihood = Integer.valueOf(posteriorStrings[i]);
            approximatePosteriors[i] = QualityUtils.qualToErrorProb(likelihood);
            total += approximatePosteriors[i];
        }

        final Map<IntegerCopyNumberState,Double> copyNumberPosteriors = new HashMap<>(SVUtils.hashMapCapacity(posteriorStrings.length));
        for (int i = 0; i < posteriorStrings.length; i++) {
            copyNumberPosteriors.put(new IntegerCopyNumberState(i), FastMath.log(Math.max(approximatePosteriors[i] / total, Double.MIN_VALUE)));
        }

        final String neutralCopyState = (String) genotype.getExtendedAttribute(GATKSVVCFConstants.NEUTRAL_COPY_NUMBER_KEY);
        Utils.nonNull(neutralCopyState, "Missing required field " + GATKSVVCFConstants.NEUTRAL_COPY_NUMBER_KEY);

        return new IntervalCopyNumberGenotypingData(
                interval,
                new CopyNumberPosteriorDistribution(copyNumberPosteriors),
                new IntegerCopyNumberState(Integer.valueOf(neutralCopyState))
        );
    }
}
