package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;

import java.util.*;

public class SVTestUtils {
    final static int chr1Length = 249250621;

    final static SAMSequenceDictionary dict = new SAMSequenceDictionary(
            Arrays.asList(new SAMSequenceRecord("chr1", chr1Length),
                    new SAMSequenceRecord("chr10", 135534747),
                    new SAMSequenceRecord("chrX", 155270560)));

    final private static GenomeLocParser glParser = new GenomeLocParser(SVTestUtils.dict);

    final static int start = 10001;

    final static int length = 10000;

    final static int length1b = (int)Math.round(SVClusterEngine.MIN_RECIPROCAL_OVERLAP_DEPTH * length);

    final static int start1b = start + length - length1b;

    //separated from end of call1 by defragmenter padding (in bin space, according to bins defined by targetIntervals below)
    final static int start2 = start + length*14/9;

    final static int start3 = start + (int)Math.round(Math.floor((1-SVClusterEngine.getMinReciprocalOverlap())*length));

    //make intervals like xxxx----xxxx----xxxx----xxxx----xxxx
    final static List<GenomeLoc> targetIntervals = new ArrayList<>(
            Arrays.asList(glParser.createGenomeLoc("chr1", 1, length/2),  //left edge call
                    glParser.createGenomeLoc("chr1", start, start + length/9),  //start of call1
                    glParser.createGenomeLoc("chr1", start + length*2/9, start + length*3/9),
                    glParser.createGenomeLoc("chr1", start + length*4/9, start + length*5/9),
                    glParser.createGenomeLoc("chr1", start + length*6/9, start + length*7/9),
                    glParser.createGenomeLoc("chr1", start + length*8/9, start + length),  //end of call1
                    glParser.createGenomeLoc("chr1", start + length*10/9, start + length*11/9),  //padding for call1
                    glParser.createGenomeLoc("chr1", start + length*12/9, start + length*13/9),  //padding for call1
                    glParser.createGenomeLoc("chr1", start2, start2 + length/9),
                    glParser.createGenomeLoc("chr1", start2 + length*2/9, start2 + length*3/9),
                    glParser.createGenomeLoc("chr1", start2 + length*4/9, start2 + length*5/9),
                    glParser.createGenomeLoc("chr1", start2 + length*6/9, start2 + length*7/9),
                    glParser.createGenomeLoc("chr1", start2 + length*8/9, start2 + length),
                    glParser.createGenomeLoc("chr1", start2 + length, start2 + length + length/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*2/9, start2 + length+ length*3/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*4/9, start2 + length+ length*5/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*6/9, start2 + length+ length*7/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*8/9, start2 + length+ length),
                    glParser.createGenomeLoc("chr1", chr1Length - 99, chr1Length),
                    glParser.createGenomeLoc("chr10", start, start + length - 1)))
            ;



    final static GenotypeBuilder gb1 = new GenotypeBuilder("sample1", Collections.singletonList(Allele.create("<"+ GATKSVVCFConstants.SYMB_ALT_STRING_DEL +">", false)));
    final static GenotypeBuilder gb2 = new GenotypeBuilder("sample1", Collections.singletonList(Allele.create("<"+ GATKSVVCFConstants.SYMB_ALT_STRING_DEL +">", false)));

    public final static Genotype sample1 = gb1.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1).make();

    final static Genotype sample1_CN0 = gb2.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0).make();

    final static Map<String, Object> sample2AttributeMap = Collections.singletonMap(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 3);

    final static Map<String, Object> sample3AttributeMap = Collections.singletonMap(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0);

    final static Genotype sample2 = GenotypeBuilder.create("sample2",
            Collections.singletonList(Allele.create("<"+GATKSVVCFConstants.SYMB_ALT_STRING_DUP +">", false)), sample2AttributeMap);

    final static SVCallRecordWithEvidence rightEdgeCall = new SVCallRecordWithEvidence("chr1", chr1Length - 99, true,
            "chr1", chr1Length, true,
            StructuralVariantType.CNV, 100,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence leftEdgeCall = new SVCallRecordWithEvidence("chr1", 1, true,
            "chr1", length/2, true,
            StructuralVariantType.CNV, length/2,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence nonDepthOnly = new SVCallRecordWithEvidence("chr1", start, true,
            "chr1", start + length -1, true,
            StructuralVariantType.CNV, length,
            Arrays.asList(GATKSVVCFConstants.DEPTH_ALGORITHM, "PE"),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    public final static SVCallRecordWithEvidence call1 = new SVCallRecordWithEvidence("chr1", start, true,
            "chr1", start + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence call1b = new SVCallRecordWithEvidence("chr1", start1b, true,
            "chr1", start1b + length1b + 1, true,
            StructuralVariantType.CNV, length1b,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence call2 = new SVCallRecordWithEvidence("chr1", start2, true,
            "chr1", start2 + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence call3 = new SVCallRecordWithEvidence("chr1", start2+length, true,
            "chr1", start2 + length + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence call4_chr10 = new SVCallRecordWithEvidence("chr10", start, true,
            "chr10", start + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence call1_CN1 = new SVCallRecordWithEvidence("chr1", start, true,
            "chr1", start + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence call2_CN0 = new SVCallRecordWithEvidence("chr1", start2, true,
            "chr1", start2 + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1_CN0),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static Genotype sample3 = GenotypeBuilder.create("sample3",
            Collections.singletonList(Allele.create("<"+GATKSVVCFConstants.SYMB_ALT_STRING_DEL +">", false)), sample3AttributeMap);

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

    static final SVCallRecordWithEvidence depthAndStuff2 = new SVCallRecordWithEvidence("chr1", 10000 - SVClusterEngine.MAX_BREAKEND_CLUSTERING_WINDOW, true, "chr1", 20000, true,
            StructuralVariantType.CNV, 10001 + SVClusterEngine.MAX_BREAKEND_CLUSTERING_WINDOW, Arrays.asList("depth", "PE"), Collections.singletonList(sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    static final SVCallRecordWithEvidence depthAndStuff3 = new SVCallRecordWithEvidence("chr1", 10000 + SVClusterEngine.MAX_BREAKEND_CLUSTERING_WINDOW, true, "chr1", 20000, true,
            StructuralVariantType.CNV, 10001 - SVClusterEngine.MAX_BREAKEND_CLUSTERING_WINDOW, Arrays.asList("depth", "PE"), Collections.singletonList(sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    final static SVCallRecordWithEvidence overlapsCall1 = new SVCallRecordWithEvidence("chr1", start3, true,
            "chr1", start3 + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Arrays.asList(sample1, sample2),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());
}
