package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.util.HashMap;
import java.util.Map;

public class GATKSVVCFHeaderLines {

    public static Map<String, VCFHeaderLine> vcfHeaderLines = new HashMap<>();

    // VCF standard reserved keys
    // todo: add these and the other standard SV info fields from the VCF spec to htsjdk VCFStandardHeaderLines
    public static final String SVTYPE = "SVTYPE";
    public static final String SVLEN = "SVLEN";

    public static final String VCF_ALT_ALLELE_STRING_DEL = "<DEL>";
    public static final String VCF_ALT_ALLELE_STRING_INV = "<INV>";
    public static final String VCF_ALT_ALLELE_STRING_INS = "<INS>";
    public static final String VCF_ALT_ALLELE_STRING_DUP = "<DUP>";

    // GATK-SV specific header lines
    // applicable to all precise records all the time
    public static final String CONTIG_NAMES = "CTG_NAMES";
    public static final String TOTAL_MAPPINGS = "TOTAL_MAPPINGS";
    public static final String MAPPING_QUALITIES = "MAPPING_QUALITIES";
    public static final String ALIGN_LENGTHS = "ALIGN_LENGTHS";
    public static final String HQ_MAPPINGS = "HQ_MAPPINGS";
    public static final String MAX_ALIGN_LENGTH = "MAX_ALIGN_LENGTH";
    public static final String VARIANT_ID_FIELD_SEPARATOR = "_";

    // applicable to all precise records when available
    public static final String INSERTED_SEQUENCE = "INSERTED_SEQUENCE";
    public static final String INSERTED_SEQUENCE_MAPPINGS = "INSERTED_SEQUENCE_MAPPINGS";
    public static final String HOMOLOGY = "HOMOLOGY";
    public static final String HOMOLOGY_LENGTH = "HOMOLOGY_LENGTH";

    // type specific: inversion
    public static final String INV33 = "INV33";
    public static final String INV55 = "INV55";

    // type specific: tandem duplication
    public static final String DUP_REPEAT_UNIT_REF_SPAN = "DUP_REPEAT_UNIT_REF_SPAN";
    public static final String DUP_SEQ_CIGARS = "DUP_SEQ_CIGARS";
    public static final String DUPLICATION_NUMBERS = "DUP_NUM";
    public static final String DUP_ANNOTATIONS_IMPRECISE = "DUP_ANNOTATIONS_IMPRECISE";
    public static final String TANDUP_CONTRACTION_STRING = "CONTRACTION";
    public static final String TANDUP_EXPANSION_STRING = "EXPANSION";
    public static final String TANDUP_CONTRACTION_ID_START_STRING = "DEL-DUPLICATION-TANDEM-CONTRACTION";
    public static final String TANDUP_EXPANSION_ID_START_STRING = "INS-DUPLICATION-TANDEM-EXPANSION";

    static {
        vcfHeaderLines.put(SVTYPE, new VCFInfoHeaderLine(SVTYPE, 1, VCFHeaderLineType.String, "Type of structural variant"));
        vcfHeaderLines.put(SVLEN, new VCFInfoHeaderLine(SVLEN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Difference in length between REF and ALT alleles"));

        vcfHeaderLines.put(TOTAL_MAPPINGS, new VCFInfoHeaderLine(TOTAL_MAPPINGS, 1, VCFHeaderLineType.Integer, "Number of contig alignments that support the variant"));
        vcfHeaderLines.put(HQ_MAPPINGS, new VCFInfoHeaderLine(HQ_MAPPINGS, 1, VCFHeaderLineType.Integer, "Number of high-quality contig alignments that support the variant"));
        vcfHeaderLines.put(MAPPING_QUALITIES, new VCFInfoHeaderLine(MAPPING_QUALITIES, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Mapping qualities of the contig alignments that support the variant"));
        vcfHeaderLines.put(ALIGN_LENGTHS, new VCFInfoHeaderLine(ALIGN_LENGTHS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Minimum lengths of the flanking aligned region from each contig alignment"));
        vcfHeaderLines.put(MAX_ALIGN_LENGTH, new VCFInfoHeaderLine(MAX_ALIGN_LENGTH, 1, VCFHeaderLineType.Integer, "Maximum of the minimum aligned lengths of flanking regions from each contig alignment"));

        // todo: create an alternate assembly file and link to it with breakpoint IDs according to the VCF spec
        vcfHeaderLines.put(CONTIG_NAMES, new VCFInfoHeaderLine(CONTIG_NAMES, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Name of contigs that evidenced this variant, formatted as \"asm%06d:tig%05d\""));

        vcfHeaderLines.put(INSERTED_SEQUENCE, new VCFInfoHeaderLine(INSERTED_SEQUENCE, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Inserted sequence at the breakpoint"));
        vcfHeaderLines.put(INSERTED_SEQUENCE_MAPPINGS, new VCFInfoHeaderLine(INSERTED_SEQUENCE_MAPPINGS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Alignments of inserted sequence"));

        vcfHeaderLines.put(HOMOLOGY, new VCFInfoHeaderLine(HOMOLOGY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Homologous sequence from contig at the breakpoint"));
        vcfHeaderLines.put(HOMOLOGY_LENGTH, new VCFInfoHeaderLine(HOMOLOGY_LENGTH, 1, VCFHeaderLineType.Integer, "Length of homologous sequence"));
        vcfHeaderLines.put(INV33, new VCFInfoHeaderLine(INV33, 0, VCFHeaderLineType.Flag, "Whether the event represents a 3' to 5' inversion"));
        vcfHeaderLines.put(INV55, new VCFInfoHeaderLine(INV55, 0, VCFHeaderLineType.Flag, "Whether the event represents a 5' to 3' inversion"));

        vcfHeaderLines.put(DUP_REPEAT_UNIT_REF_SPAN, new VCFInfoHeaderLine(DUP_REPEAT_UNIT_REF_SPAN, 1, VCFHeaderLineType.String, "Reference span of the suspected repeated unit in a tandem duplication"));
        vcfHeaderLines.put(DUP_SEQ_CIGARS, new VCFInfoHeaderLine(DUP_SEQ_CIGARS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,
                "CIGARs of the repeated sequence on the locally-assembled contigs when aligned to " + DUP_REPEAT_UNIT_REF_SPAN + " (currently only available for repeats when " + DUP_ANNOTATIONS_IMPRECISE + " is false)"));
        vcfHeaderLines.put(DUPLICATION_NUMBERS, new VCFInfoHeaderLine(DUPLICATION_NUMBERS, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "Number of times the sequence is duplicated on reference and on the alternate alleles"));
        vcfHeaderLines.put(DUP_ANNOTATIONS_IMPRECISE, new VCFInfoHeaderLine(DUP_ANNOTATIONS_IMPRECISE, 0, VCFHeaderLineType.Flag, "Whether the duplication annotations are from an experimental optimization procedure"));
    }


}
