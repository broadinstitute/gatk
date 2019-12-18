package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVSegmentVariantComposer;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

// TODO: 7/24/17 the structure of this file is resembling that of {@link GATKVCFHeaderLines}, should we move this there?
public final class GATKSVVCFHeaderLines {

    public static VCFInfoHeaderLine getInfoLine(final String id) { return infoLines.get(id); }
    public static Set<VCFInfoHeaderLine> getInfoLines() { return new HashSet<>(infoLines.values()); }
    public static VCFFormatHeaderLine getFormatLine(final String id) { return formatLines.get(id); }
    public static Set<VCFFormatHeaderLine> getFormatLines() { return new HashSet<>(formatLines.values()); }
    public static VCFFilterHeaderLine getFilterLine(final String id) { return filterLines.get(id); }
    public static Set<VCFFilterHeaderLine> getFilterLines() { return new HashSet<>(filterLines.values());  }

    private static final Map<String, VCFInfoHeaderLine> infoLines = new LinkedHashMap<>(30);
    private static final Map<String, VCFFormatHeaderLine> formatLines = new LinkedHashMap<>(5);
    private static final Map<String, VCFFilterHeaderLine> filterLines = new LinkedHashMap<>(2);


    private static void addFormatLine(final VCFFormatHeaderLine line) {
        Utils.nonNull(line);
        formatLines.put(line.getID(), line);
    }

    private static void addInfoLine(final VCFInfoHeaderLine line) {
        Utils.nonNull(line);
        infoLines.put(line.getID(), line);
    }

    private static void addFilterLine(final VCFFilterHeaderLine line) {
        Utils.nonNull(line);
        filterLines.put(line.getID(), line);
    }

    // todo htsjdk should have these defined
    public static VCFSimpleHeaderLine getSymbAltAlleleLine(final String id) { return symbAltAlleleLines.get(id); }
    public static Set<VCFSimpleHeaderLine> getSymbAltAlleleLines() { return new HashSet<>(symbAltAlleleLines.values()); }
    private static final Map<String, VCFSimpleHeaderLine> symbAltAlleleLines = new LinkedHashMap<>(10);
    private static void addSymbAltAlleleLine(final VCFSimpleHeaderLine line) {
        Utils.nonNull(line);
        symbAltAlleleLines.put(line.getID(), line);
    }

    static {

        addSymbAltAlleleLine(new VCFSimpleHeaderLine(GATKVCFConstants.SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG,
                GATKSVVCFConstants.SYMB_ALT_STRING_INV, "Inversion of reference sequence"));
        addSymbAltAlleleLine(new VCFSimpleHeaderLine(GATKVCFConstants.SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG,
                GATKSVVCFConstants.SYMB_ALT_STRING_DEL, "Deletion relative to the reference"));
        addSymbAltAlleleLine(new VCFSimpleHeaderLine(GATKVCFConstants.SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG,
                GATKSVVCFConstants.SYMB_ALT_STRING_INS, "Insertion of novel sequence relative to the reference"));
        addSymbAltAlleleLine(new VCFSimpleHeaderLine(GATKVCFConstants.SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG,
                GATKSVVCFConstants.SYMB_ALT_STRING_DUP, "Region of elevated copy number relative to the reference"));
        addSymbAltAlleleLine(new VCFSimpleHeaderLine(GATKVCFConstants.SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG,
                GATKSVVCFConstants.SYMB_ALT_STRING_INVDUP, "Region of elevated copy number relative to the reference, with some copies inverted"));
        addSymbAltAlleleLine(new VCFSimpleHeaderLine(GATKVCFConstants.SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG,
                GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR, "Complex rearrangement of reference sequence"));

        // descriptions on INFO fields that are available for each record
        addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.SVTYPE,
                1,
                VCFHeaderLineType.String,
                "Type of structural variant"));
        
        // optional INFO annotations applicable to all variants
        {
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.SVLEN,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.Integer,
                    "Difference in length between REF and ALT alleles"));
            
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.SEQ_ALT_HAPLOTYPE,
                    VCFHeaderLineCount.A,
                    VCFHeaderLineType.Character, "Alt haplotype sequence, one per alt allele"));

            // TODO: 3/9/18 INSLEN missing (see ticket 4382)
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INSERTED_SEQUENCE,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.String,
                    "Inserted sequence at the breakpoint"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INSERTED_SEQUENCE_LENGTH,
                    VCFHeaderLineCount.A,
                    VCFHeaderLineType.Integer,
                    "Length of inserted sequence (note for duplication records, this does not count the extra copies of the duplicated sequence)"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INSERTED_SEQUENCE_MAPPINGS,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.String,
                    "Alignments of inserted sequence"));

            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HOMOLOGY,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.String,
                    "Sequence of base pair identical micro-homology at event breakpoints"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HOMOLOGY_LENGTH,
                    1,
                    VCFHeaderLineType.Integer,
                    "Length of base pair identical micro-homology at event breakpoints"));

            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.READ_PAIR_SUPPORT,
                    1,
                    VCFHeaderLineType.Integer,
                    "Number of discordant read pairs supporting the variant"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.SPLIT_READ_SUPPORT,
                    1,
                    VCFHeaderLineType.Integer,
                    "Number of split read supplementary mappings supporting the variant"));

            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.EXTERNAL_CNV_CALLS,
                    1,
                    VCFHeaderLineType.String,
                    "Comma-delimited list of external copy number calls that overlap with this variant in format ID:CN:CNQ"));

            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.LINK,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.String,
                    "ID(s) of other record(s) linked to current record"));
        }
        
        // for variants-detected from assembly
        {// todo: create an alternate assembly file and link to it with breakpoint IDs according to the VCF spec
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CONTIG_NAMES,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.String,
                    "Name of local assembly contigs supporting this variant, formatted as \"asm%06d:tig%05d\""));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TOTAL_MAPPINGS,
                    1,
                    VCFHeaderLineType.Integer,
                    "Number of local assembly contigs supporting the variant, i.e. number of entries in " + GATKSVVCFConstants.CONTIG_NAMES));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.MAPPING_QUALITIES,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.Integer,
                    "Mapping qualities of the contig alignments that support the variant"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HQ_MAPPINGS,
                    1,
                    VCFHeaderLineType.Integer,
                    "Number of high-quality contig alignments that support the variant"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.ALIGN_LENGTHS,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.Integer,
                    "Minimum lengths of the flanking aligned region from each contig alignment"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.MAX_ALIGN_LENGTH,
                    1,
                    VCFHeaderLineType.Integer,
                    "Maximum of the values listed in " + GATKSVVCFConstants.ALIGN_LENGTHS));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CTG_GOOD_NONCANONICAL_MAPPING,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.String,
                    "Good mapping of evidence contigs as listed in " + GATKSVVCFConstants.CONTIG_NAMES +
                            " to non-canonical reference chromosomes that could potentially offer better explanation of the assembly contig without the SV record." +
                            " One for each evidence assembly contig, if available; otherwise a \".\"." +
                            " If no evidence contig has such mapping, this annotation is omitted for the record."));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.LINK,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.String,
                    "ID(s) of other variants that are linked to this record, " +
                            "i.e. they constitute a larger more complex variant"));
        }

        // type specific
        {
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BND_MATEID_STR,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.String,
                    "ID(s) for mate(s) of a BND record")); // todo: technically there could be multiple mates, but currently we only have case for 1 mate

            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.IMPRECISE,
                    0,
                    VCFHeaderLineType.Flag,
                    "Imprecise structural variation"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CIPOS,
                    2,
                    VCFHeaderLineType.Integer,
                    "Confidence interval around POS for imprecise variants"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CIEND,
                    2,
                    VCFHeaderLineType.Integer,
                    "Confidence interval around END for imprecise variants"));

            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INV33,
                    0,
                    VCFHeaderLineType.Flag,
                    "Whether the event represents a 3' to 5' breakpoint"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INV55,
                    0,
                    VCFHeaderLineType.Flag,
                    "Whether the event represents a 5' to 3' breakpoint"));

            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN,
                    1,
                    VCFHeaderLineType.String,
                    "Reference span of the suspected repeated unit in a tandem duplication"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DUP_SEQ_CIGARS,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.String,
                    "CIGARs of the repeated sequence on the locally-assembled contigs when aligned to " +
                            GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN +
                            " (currently only available for repeats when " + GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE + " is false)"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DUPLICATION_NUMBERS,
                    VCFHeaderLineCount.R,
                    VCFHeaderLineType.Integer,
                    "Number of times the sequence is duplicated on reference and on the alternate alleles"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE,
                    0, VCFHeaderLineType.Flag,
                    "Whether the duplication annotations are from an experimental optimization procedure"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DUP_IMPRECISE_AFFECTED_RANGE,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.String,
                    "Affected reference range for duplications annotated with the flag " + GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DUP_ORIENTATIONS,
                    VCFHeaderLineCount.A,
                    VCFHeaderLineType.String,
                    "Orientations of the duplicated sequence on alt allele relative to the copy on ref;" +
                            " one group for each alt allele (currently only available for inverted duplication variants)"));

            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DUP_TAN_CONTRACTION_STRING,
                    0,
                    VCFHeaderLineType.Flag,
                    "Tandem repeats contraction compared to reference"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DUP_TAN_EXPANSION_STRING,
                    0,
                    VCFHeaderLineType.Flag,
                    "Tandem repeats expansion compared to reference"));

            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CPX_EVENT_ALT_ARRANGEMENTS,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.String,
                    "For CPX variants only; specifies how reference segments given in " + GATKSVVCFConstants.CPX_SV_REF_SEGMENTS + " are re-arranged"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CPX_SV_REF_SEGMENTS,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.String,
                    "For CPX variants only; segments of reference that are rearranged"));
            addInfoLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CPX_EVENT_KEY,
                    VCFHeaderLineCount.UNBOUNDED,
                    VCFHeaderLineType.String,
                    "ID(s) of " + GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR + "(s) events from which current simple variant record is extracted"));
        }

        // format lines
        {
            addFormatLine(new VCFFormatHeaderLine(GATKSVVCFConstants.COPY_NUMBER_FORMAT,
                    1, // TODO: 7/3/18 Spec 4.3 has this example value at bottom of page 12, but what about multi-allelic sites?
                    VCFHeaderLineType.Integer,
                    "Copy number genotype for imprecise events"));

            addFormatLine(new VCFFormatHeaderLine(GATKSVVCFConstants.COPY_NUMBER_QUALITY_FORMAT,
                    1, // TODO: 7/3/18 Spec 4.3 has this example value at bottom of page 12, but what about multi-allelic sites?
                    VCFHeaderLineType.Float,
                    "Copy number genotype quality for imprecise events"));
        }

        // filter lines
        {
            addFilterLine(new VCFFilterHeaderLine(GATKSVVCFConstants.ASSEMBLY_BASED_VARIANT_MQ_FILTER_KEY,
                    "Assembly evidence based record that whose maximum value specified in " + GATKSVVCFConstants.MAPPING_QUALITIES + " is lower than user specified threshold"));

            addFilterLine(new VCFFilterHeaderLine(GATKSVVCFConstants.ASSEMBLY_BASED_VARIANT_ALN_LENGTH_FILTER_KEY,
                    "Assembly evidence based record that whose " + GATKSVVCFConstants.MAPPING_QUALITIES + " value is lower than user specified threshold"));

            addFilterLine(new VCFFilterHeaderLine(GATKSVVCFConstants.LOW_QS_SCORE_FILTER_KEY,
                    "Depth-only copy number record whose " + GermlineCNVSegmentVariantComposer.QS + " value is lower than user specified threshold"));

            addFilterLine(new VCFFilterHeaderLine(GATKSVVCFConstants.FREQUENCY_FILTER_KEY,
                    "Depth-only copy number record whose " + VCFConstants.ALLELE_FREQUENCY_KEY + " value is higher than user specified threshold"));
        }
    }
}
