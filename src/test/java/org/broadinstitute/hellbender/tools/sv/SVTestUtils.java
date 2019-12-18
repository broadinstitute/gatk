package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SVTestUtils {
    final static int chr1Length = 249250621;

    final static SAMSequenceDictionary dict = new SAMSequenceDictionary(
            Arrays.asList(new SAMSequenceRecord("chr1", chr1Length),
                    new SAMSequenceRecord("chr10", 135534747),
                    new SAMSequenceRecord("chrX", 155270560)));

    final static int start = 10001;

    final static int length = 10000;

    //separated from end of call1 by defragmenter padding
    final static int start2 = (start + length -1) + (int)Math.round(length * SVDepthOnlyCallDefragmenterTest.defaultDefragmenter.getPaddingFraction());

    final static Genotype sample1 = GenotypeBuilder.create("sample1", Collections.singletonList(Allele.create("<"+ GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL+">", false)));

    final static Genotype sample2 = GenotypeBuilder.create("sample2", Collections.singletonList(Allele.create("<"+GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP+">", false)));

    final static SVCallRecordWithEvidence rightEdgeCall = new SVCallRecordWithEvidence("chr1", chr1Length - 99, true,
            "chr1", chr1Length, true,
            StructuralVariantType.CNV, 100,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence leftEdgeCall = new SVCallRecordWithEvidence("chr1", 1, true,
            "chr1", length, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence nonDepthOnly = new SVCallRecordWithEvidence("chr1", start, true,
            "chr1", start + length -1, true,
            StructuralVariantType.CNV, length,
            Arrays.asList(GATKSVVCFConstants.DEPTH_ALGORITHM, "PE"),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence call1 = new SVCallRecordWithEvidence("chr1", start, true,
            "chr1", start + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence call2 = new SVCallRecordWithEvidence("chr1", start2, true,
            "chr1", start2 + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static Genotype sample3 = GenotypeBuilder.create("sample3", Collections.singletonList(Allele.create("<"+GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL+">", false)));

    final static SVCallRecordWithEvidence sameBoundsSampleMismatch = new SVCallRecordWithEvidence("chr1", start, true,
            "chr1", start + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample3),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static List<Genotype> threeGenotypes = Arrays.asList(sample1, sample2, sample3);

    final static SVCallRecordWithEvidence inversion = new SVCallRecordWithEvidence("chr1", 10000, true, "chr1", 20000, true,
            StructuralVariantType.INV, 10001, Arrays.asList("SR", "PE"), Collections.singletonList(sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    static final SVCallRecordWithEvidence depthOnly = new SVCallRecordWithEvidence("chr1", 10000, true, "chr1", 20000, true,
            StructuralVariantType.CNV, 10001, Collections.singletonList("depth"), Collections.singletonList(sample1),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    static final SVCallRecordWithEvidence depthAndStuff = new SVCallRecordWithEvidence("chr1", 10000, true, "chr1", 20000, true,
                    StructuralVariantType.CNV, 10001, Arrays.asList("depth", "PE"), Collections.singletonList(sample2),
                    Collections.emptyList(), Collections.emptyList(), Collections.emptyList());
}
