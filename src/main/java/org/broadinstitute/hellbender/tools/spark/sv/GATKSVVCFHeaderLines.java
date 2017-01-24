package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.util.HashMap;
import java.util.Map;

public class GATKSVVCFHeaderLines {

    public static Map<String, VCFHeaderLine> vcfHeaderLines = new HashMap<>();

    public static final String FORMAT_FIELD_SEPARATOR = "%";

    // VCF standard SV header lines
    // todo: add these and the other standard SV info fields from the VCF spec to htsjdk VCFStandardHeaderLines
    public static final String SVTYPE = "SVTYPE";
    public static final String SVLEN = "SVLEN";

    // GATK-SV specific header lines
    public static final String TOTAL_MAPPINGS = "TOTAL_MAPPINGS";
    public static final String HQ_MAPPINGS = "HQ_MAPPINGS";
    public static final String MAPPING_QUALITIES = "MAPPING_QUALITIES";
    public static final String ALIGN_LENGTHS = "ALIGN_LENGTHS";
    public static final String MAX_ALIGN_LENGTH = "MAX_ALIGN_LENGTH";
    public static final String ASSEMBLY_IDS = "ASSEMBLY_IDS";
    public static final String CONTIG_IDS = "CONTIG_IDS";
    public static final String INSERTED_SEQUENCE = "INSERTED_SEQUENCE";
    public static final String INSERTED_SEQUENCE_MAPPINGS = "INSERTED_SEQUENCE_MAPPINGS";
    public static final String HOMOLOGY = "HOMOLOGY";
    public static final String HOMOLOGY_LENGTH = "HOMOLOGY_LENGTH";
    public static final String INV_3_TO_5 = "INV_3_TO_5";
    public static final String INV_5_TO_3 = "INV_5_TO_3";

    static {
        vcfHeaderLines.put(SVTYPE, new VCFInfoHeaderLine(SVTYPE, 1, VCFHeaderLineType.String, "Type of structural variant"));
        vcfHeaderLines.put(SVLEN, new VCFInfoHeaderLine(SVLEN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Difference in length between REF and ALT alleles"));

        vcfHeaderLines.put(TOTAL_MAPPINGS, new VCFInfoHeaderLine(TOTAL_MAPPINGS, 1, VCFHeaderLineType.Integer, "Number of contig alignments that support the variant"));
        vcfHeaderLines.put(HQ_MAPPINGS, new VCFInfoHeaderLine(HQ_MAPPINGS, 1, VCFHeaderLineType.Integer, "Number of high-quality contig alignments that support the variant"));
        vcfHeaderLines.put(MAPPING_QUALITIES, new VCFInfoHeaderLine(MAPPING_QUALITIES, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Mapping qualities of the contig alignments that support the variant"));
        vcfHeaderLines.put(ALIGN_LENGTHS, new VCFInfoHeaderLine(ALIGN_LENGTHS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Minimum lengths of the flanking aligned region from each contig alignment"));
        vcfHeaderLines.put(MAX_ALIGN_LENGTH, new VCFInfoHeaderLine(MAX_ALIGN_LENGTH, 1, VCFHeaderLineType.Integer, "Maximum of the minimum aligned lengths of flanking regions from each contig alignment"));

        // todo: create an alternate assembly file and link to it with breakpoint IDs according to the VCF spec
        vcfHeaderLines.put(ASSEMBLY_IDS, new VCFInfoHeaderLine(ASSEMBLY_IDS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "IDs of the assemblies that produced each contig alignment"));
        vcfHeaderLines.put(CONTIG_IDS, new VCFInfoHeaderLine(CONTIG_IDS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "IDs of the contigs that produced each alignment"));

        vcfHeaderLines.put(INSERTED_SEQUENCE, new VCFInfoHeaderLine(INSERTED_SEQUENCE, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Inserted sequence at the breakpoint"));
        vcfHeaderLines.put(INSERTED_SEQUENCE_MAPPINGS, new VCFInfoHeaderLine(INSERTED_SEQUENCE_MAPPINGS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Alignments of inserted sequence"));

        vcfHeaderLines.put(HOMOLOGY, new VCFInfoHeaderLine(HOMOLOGY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Homologous sequence from contig at the breakpoint"));
        vcfHeaderLines.put(HOMOLOGY_LENGTH, new VCFInfoHeaderLine(HOMOLOGY_LENGTH, 1, VCFHeaderLineType.Integer, "Length of homologous sequence"));
        vcfHeaderLines.put(INV_3_TO_5, new VCFInfoHeaderLine(INV_3_TO_5, 0, VCFHeaderLineType.Flag, "Whether the event represents a 3' to 5' inversion"));
        vcfHeaderLines.put(INV_5_TO_3, new VCFInfoHeaderLine(INV_5_TO_3, 0, VCFHeaderLineType.Flag, "Whether the event represents a 5' to 3' inversion"));
    }

    enum SVTYPES {
        DEL, INS, DUP, INV, CNV, BND
    }
}
