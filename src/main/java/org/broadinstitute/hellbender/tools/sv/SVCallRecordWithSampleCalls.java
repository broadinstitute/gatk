package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVSegmentVariantComposer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.SVCallRecordCodec;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.util.*;
import java.util.stream.Collectors;

public class SVCallRecordWithSampleCalls extends SVCallRecordWithEvidence {

    private final Map<String, Integer> samplesWithCopyState;

    private final static List<String> gCNVCallerAttributes = Arrays.asList(
            SVCluster.END_CONTIG_ATTRIBUTE,
            SVCluster.SVLEN_ATTRIBUTE,
            SVCluster.SVTYPE_ATTRIBUTE
    );

    public static SVCallRecordWithSampleCalls create(final VariantContext variant) {
        Utils.nonNull(variant);
        Utils.validate(variant.getAttributes().keySet().containsAll(gCNVCallerAttributes), "Call is missing attributes");
        final String startContig = variant.getContig();
        final int start = variant.getStart();
        final String endContig = variant.getAttributeAsString(SVCluster.END_CONTIG_ATTRIBUTE, "NA");
        final int end = variant.getEnd();
        final StructuralVariantType type = variant.getStructuralVariantType();
        final List<String> algorithms = variant.getAttributeAsStringList(SVCluster.ALGORITHMS_ATTRIBUTE, "NA");
        final String strands = variant.getAttributeAsString(SVCluster.STRANDS_ATTRIBUTE, "0");
        if (strands.length() != 2) {
            throw new IllegalArgumentException("Strands field is not 2 characters long");
        }
        final String startStrandChar = strands.substring(0, 1);
        if (!startStrandChar.equals(SVCallRecordCodec.STRAND_PLUS) && !startStrandChar.equals(SVCallRecordCodec.STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid start strand not found");
        }
        final String endStrandChar = strands.substring(1, 2);
        if (!endStrandChar.equals(SVCallRecordCodec.STRAND_PLUS) && !endStrandChar.equals(SVCallRecordCodec.STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid end strand not found");
        }
        final boolean startStrand = startStrandChar.equals(SVCallRecordCodec.STRAND_PLUS);
        final boolean endStrand = endStrandChar.equals(SVCallRecordCodec.STRAND_PLUS);
        final int length = end - start + 1;
        final Map<String, Integer> samplesWithCopyState = variant.getGenotypes().stream()
                .filter(Genotype::isCalled)
                .collect(Collectors.toMap(Genotype::getSampleName, g -> Integer.parseInt(g.getExtendedAttribute(GermlineCNVSegmentVariantComposer.CN, HomoSapiensConstants.DEFAULT_PLOIDY).toString())));
        return new SVCallRecordWithSampleCalls(startContig, start, startStrand, endContig, end, endStrand, type, length, algorithms, samplesWithCopyState);
    }

    public static SVCallRecordWithSampleCalls createDepthOnlyFromGCNV(final VariantContext variant, final double minQuality) {
        Utils.nonNull(variant);
        final Map<String, Integer> samplesWithCopyState = variant.getGenotypes().stream()
                .filter(Genotype::isCalled)
                .filter(g -> Integer.valueOf((String)g.getExtendedAttribute(GermlineCNVSegmentVariantComposer.QS)) >= minQuality)
                .collect(Collectors.toMap(Genotype::getSampleName, g -> Integer.parseInt(g.getExtendedAttribute(GermlineCNVSegmentVariantComposer.CN, HomoSapiensConstants.DEFAULT_PLOIDY).toString())));
        if (samplesWithCopyState.isEmpty()) return null;
        final List<String> algorithms = Collections.singletonList(SVCluster.DEPTH_ALGORITHM);

        //TODO : use new vcfs to get actual allele
        final int copyNumber = Integer.valueOf((String)variant.getGenotypes().get(0).getExtendedAttribute(GermlineCNVSegmentVariantComposer.CN));
        if (copyNumber == 2) return null;
        final boolean isDel = copyNumber < 2;
        final boolean startStrand = isDel ? true : false;
        final boolean endStrand = isDel ? false : true;
        final StructuralVariantType type = isDel ? StructuralVariantType.DEL : StructuralVariantType.DUP;

        final String startContig = variant.getContig();
        final int start = variant.getStart();
        final int end = variant.getEnd() + 1; // TODO this is a bug with gCNV vcf generation
        final int length = end - start;
        return new SVCallRecordWithSampleCalls(startContig, start, startStrand, startContig, end, endStrand, type, length, algorithms, samplesWithCopyState);
    }

    public SVCallRecordWithSampleCalls(final String startContig,
                                       final int start,
                                       final boolean startStrand,
                                       final String endContig,
                                       final int end,
                                       final boolean endStrand,
                                       final StructuralVariantType type,
                                       final int length,
                                       final List<String> algorithms,
                                       final Map<String, Integer> sampleCopyStates) {
        super(startContig, start, startStrand, endContig, end, endStrand, type, length, algorithms, sampleCopyStates.keySet());
        samplesWithCopyState = sampleCopyStates;
    }
}
