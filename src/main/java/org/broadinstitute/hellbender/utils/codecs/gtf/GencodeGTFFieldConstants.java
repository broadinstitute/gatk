package org.broadinstitute.hellbender.utils.codecs.gtf;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * This class contains enums for the fields of the Gencode GTF format (https://www.gencodegenes.org/pages/data_format.html).
 * The intention of this class is to store a list of spec compliant values for each field to reference should they be needed
 * in funcotator. In general however this list should not be considered proscriptive and these enums should generally only
 * be instantiated if necessary to preserve the spec compliance of funcotator against future releases of Gencode.
 */
public class GencodeGTFFieldConstants {

    /**
     * A list of required fields from the GencodeGTF codec. These are the fields that are currently first-class members of the
     * GencodeGTFFeature and are handled specially handled by the GencodeGTFCodec.
     */
    public static final String GENE_ID = "gene_id";
    public static final String TRANSCRIPT_ID = "transcript_id";
    public static final String GENE_TYPE = "gene_type";
    public static final String GENE_BIOTYPE = "gene_biotype";
    public static final String GENE_STATUS = "gene_status";
    public static final String GENE_NAME = "gene_name";
    public static final String TRANSCRIPT_TYPE = "transcript_type";
    public static final String TRANSCRIPT_BIOTYPE = "transcript_biotype";
    public static final String TRANSCRIPT_STATUS = "transcript_status";
    public static final String TRANSCRIPT_NAME = "transcript_name";
    public static final String LEVEL = "level";

    // Mandatory except in Gene and Transcript lines (not enforced by this codec)
    public static final String EXON_NUMBER = "exon_number";
    public static final String EXON_ID = "exon_id";

    public class GencodeOptionalFields {
        /**
         * Spec optional gencode fields:
         */
        public static final String TAG = "tag";
        public static final String CCDSID = "ccdsid";
        public static final String HAVANA_GENE = "havana_gene";
        public static final String HAVANA_TRANSCRIPT = "havana_transcript";
        public static final String PROTEIN_ID = "protein_id";
        public static final String ONT = "ont";
        public static final String TRANSCRIPT_SUPPORT_LEVEL = "transcript_support_level";
        public static final String REMAP_STATUS = "remap_status";
        public static final String REMAP_ORIGINAL_ID = "remap_original_id";
        public static final String REMAP_ORIGINAL_LOCATION = "remap_original_location";
        public static final String REMAP_NUM_MAPPINGS = "remap_num_mappings";
        public static final String REMAP_TARGET_STATUS = "remap_target_status";
        public static final String REMAP_SUBSTITUTED_MISSING_TARGET = "remap_substituted_missing_target";
        public static final String HGNC_ID = "hgnc_id";
        public static final String MGI_ID = "mgi_id";
    }

    /**
     * Biotype / transcript type for the transcript or gene represented in a feature.
     * This is a tag of some biological function associated with a feature.
     *
     * The values here are not exhaustive, but should be used as a reference for already used / known types.
     *
     * For more information, see:
     *     https://www.gencodegenes.org/data_format.html
     *     https://en.wikipedia.org/wiki/General_feature_format
     */
    public enum KnownGeneBiotype {
        // Immunoglobulin (Ig) variable chain and T-cell receptor (TcR) genes imported or annotated according to the IMGT (http://www.imgt.org/)
        IG_C_GENE("IG_C_gene"),
        IG_D_GENE("IG_D_gene"),
        IG_J_GENE("IG_J_gene"),
        IG_LV_GENE("IG_LV_gene"),
        IG_V_GENE("IG_V_gene"),
        TR_C_GENE("TR_C_gene"),
        TR_J_GENE("TR_J_gene"),
        TR_V_GENE("TR_V_gene"),
        TR_D_GENE("TR_D_gene"),

        // Inactivated immunoglobulin gene.
        IG_PSEUDOGENE("IG_pseudogene"),
        IG_C_PSEUDOGENE("IG_C_pseudogene"),
        IG_J_PSEUDOGENE("IG_J_pseudogene"),
        IG_V_PSEUDOGENE("IG_V_pseudogene"),
        TR_V_PSEUDOGENE("TR_V_pseudogene"),
        TR_J_PSEUDOGENE("TR_J_pseudogene"),

        // Non-coding RNA predicted using sequences from Rfam (http://rfam.xfam.org/) and miRBase (http://www.mirbase.org/)
        MT_RRNA("Mt_rRNA"),
        MT_TRNA("Mt_tRNA"),
        MIRNA("miRNA"),
        MISC_RNA("misc_RNA"),
        RRNA("rRNA"),

        SCRNA("scRNA"),
        SNRNA("snRNA"),
        SNORNA("snoRNA"),
        RIBOZYME("ribozyme"),
        SRNA("sRNA"),
        SCARNA("scaRNA"),

        // ENSEMBL-Specific values:
        TRNA("tRNA"),
        TMRNA("tmRNA"),

        // Non-coding RNA predicted to be pseudogene by the Ensembl pipeline
        MT_TRNA_PSEUDOGENE("Mt_tRNA_pseudogene"),
        TRNA_PSEUDOGENE("tRNA_pseudogene"),
        SNORNA_PSEUDOGENE("snoRNA_pseudogene"),
        SNRNA_PSEUDOGENE("snRNA_pseudogene"),
        SCRNA_PSEUDOGENE("scRNA_pseudogene"),
        RRNA_PSEUDOGENE("rRNA_pseudogene"),
        MISC_RNA_PSEUDOGENE("misc_RNA_pseudogene"),
        MIRNA_PSEUDOGENE("miRNA_pseudogene"),

        // To be Experimentally Confirmed. This is used for non-spliced EST clusters that have polyA features. This category has been specifically created for the ENCODE project to highlight regions that could indicate the presence of protein coding genes that require experimental validation, either by 5' RACE or RT-PCR to extend the transcripts, or by confirming expression of the putatively-encoded peptide with specific antibodies.
        TEC("TEC"),

        // If the coding sequence (following the appropriate reference) of a transcript finishes >50bp from a downstream splice site then it is tagged as NMD. If the variant does not cover the full reference coding sequence then it is annotated as NMD if NMD is unavoidable i.e. no matter what the exon structure of the missing portion is the transcript will be subject to NMD.
        NONSENSE_MEDIATED_DECAY("nonsense_mediated_decay"),

        // Transcript that has polyA features (including signal) without a prior stop codon in the CDS, i.e. a non-genomic polyA tail attached directly to the CDS without 3' UTR. These transcripts are subject to degradation.
        NON_STOP_DECAY("non_stop_decay"),

        // Alternatively spliced transcript believed to contain intronic sequence relative to other, coding, variants.
        RETAINED_INTRON("retained_intron"),

        // Contains an open reading frame (ORF).
        PROTEIN_CODING("protein_coding"),

        // Not translated in the reference genome owing to a SNP/DIP but in other individuals/haplotypes/strains the transcript is translated. Replaces the polymorphic_pseudogene transcript biotype.
        PROTEIN_CODING_LOF("protein_coding_LoF"),

        // Transcript that belongs to a protein_coding gene and doesn't contain an ORF. Replaces the processed_transcript transcript biotype in protein_coding genes.
        PROTEIN_CODING_CDS_NOT_DEFINED("protein_coding_CDS_not_defined"),

        // Doesn't contain an ORF.
        PROCESSED_TRANSCRIPT("processed_transcript"),

        // Transcript which is known from the literature to not be protein coding.
        NON_CODING("non_coding"),

        // Transcript believed to be protein coding, but with more than one possible open reading frame.
        AMBIGUOUS_ORF("ambiguous_orf"),

        // Long non-coding transcript in introns of a coding gene that does not overlap any exons.
        SENSE_INTRONIC("sense_intronic"),

        // Long non-coding transcript that contains a coding gene in its intron on the same strand.
        SENSE_OVERLAPPING("sense_overlapping"),

        // Has transcripts that overlap the genomic span (i.e. exon or introns) of a protein-coding locus on the opposite strand.
        ANTISENSE("antisense"),
        ANTISENSE_RNA("antisense_RNA"),

        KNOWN_NCRNA("known_ncrna"),

        // Have homology to proteins but generally suffer from a disrupted coding sequence and an active homologous gene can be found at another locus. Sometimes these entries have an intact coding sequence or an open but truncated ORF, in which case there is other evidence used (for example genomic polyA stretches at the 3' end) to classify them as a pseudogene. Can be further classified as one of the following.
        PSEUDOGENE("pseudogene"),

        // Pseudogene that lack introns and is thought to arise from reverse transcription of mRNA followed by reinsertion of DNA into the genome.
        PROCESSED_PSEUDOGENE("processed_pseudogene"),

        // Pseudogene owing to a SNP/DIP but in other individuals/haplotypes/strains the gene is translated.
        POLYMORPHIC_PSEUDOGENE("polymorphic_pseudogene"),

        // Pseudogene owing to a reverse transcribed and re-inserted sequence.
        RETROTRANSPOSED("retrotransposed"),

        // Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression.
        TRANSCRIBED_PROCESSED_PSEUDOGENE("transcribed_processed_pseudogene"),
        TRANSCRIBED_UNPROCESSED_PSEUDOGENE("transcribed_unprocessed_pseudogene"),
        TRANSCRIBED_UNITARY_PSEUDOGENE("transcribed_unitary_pseudogene"),

        // Pseudogene that has mass spec data suggesting that it is also translated.
        TRANSLATED_PROCESSED_PSEUDOGENE("translated_processed_pseudogene"),
        TRANSLATED_UNPROCESSED_PSEUDOGENE("translated_unprocessed_pseudogene"),

        // A species specific unprocessed pseudogene without a parent gene, as it has an active orthologue in another species.
        UNITARY_PSEUDOGENE("unitary_pseudogene"),

        // Pseudogene that can contain introns since produced by gene duplication.
        UNPROCESSED_PSEUDOGENE("unprocessed_pseudogene"),

        // Used to tag mistakes in the public databases (Ensembl/SwissProt/Trembl)
        ARTIFACT("artifact"),

        // Long, intervening noncoding (linc) RNA that can be found in evolutionarily conserved, intergenic regions.
        LINCRNA("lincRNA"),
        LNCRNA("lncRNA"),

        // Unspliced lncRNA that is several kb in size.
        MACRO_LNCRNA("macro_lncRNA"),

        // Transcript where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR.
        THREE_PRIME_OVERLAPPING_NCRNA("3prime_overlapping_ncRNA"),

        // Otherwise viable coding region omitted from this alternatively spliced transcript because the splice variation affects a region coding for a protein domain.
        DISRUPTED_DOMAIN("disrupted_domain"),

        // Short non coding RNA gene that forms part of the vault ribonucleoprotein complex.
        VAULTRNA("vaultRNA"),

        // A non-coding locus that originates from within the promoter region of a protein-coding gene, with transcription proceeding in the opposite direction on the other strand.
        BIDIRECTIONAL_PROMOTER_LNCRNA("bidirectional_promoter_lncRNA");

        @SuppressWarnings("unchecked")
        private static final Map<String, KnownGeneBiotype> VALUE_MAP =
                Arrays.stream(values()).collect(Collectors.toMap(v -> v.serialized.toLowerCase(), v -> v));

        private final String serialized;

        KnownGeneBiotype(final String serializedValue) {
            serialized = serializedValue;
        }

        @Override
        public String toString() {
            return serialized;
        }

        private static final Map<String, String> SPECIAL_CASE_STRING_VALUE_MAP = createSpecialCaseMap();

        public static KnownGeneBiotype getEnum(final String s) {
            String lowerS = s.toLowerCase();

            // Handle special cases:
            lowerS = SPECIAL_CASE_STRING_VALUE_MAP.getOrDefault(lowerS, lowerS);

            if ( VALUE_MAP.containsKey(lowerS) ){
                return VALUE_MAP.get(lowerS);
            }
            throw new IllegalArgumentException("Unexpected value: " + s);
        }

        /**
         * Create a special case map for alternate field names for known {@link KnownGeneBiotype}s.
         */
        private static Map<String, String> createSpecialCaseMap() {
            final Map<String, String> map = new HashMap<>();

            // From ENSEMBLE GTF files:
            map.put("ncrna", "non_coding");

            return map;
        }

    }

    /**
     * Indication of whether a feature is new, tenatative, or already known.
     *
     * This attribute was removed after release 25.
     *
     * For more information, see:
     *     https://www.gencodegenes.org/data_format.html
     *     https://en.wikipedia.org/wiki/General_feature_format
     */
    public enum GeneTranscriptStatus {
        KNOWN,
        NOVEL,
        PUTATIVE
    }

    /**
     * Status of how a position was annotated / verified:
     *
     *      1 - verified locus
     *      2 - manually annotated locus
     *      3 - automatically annotated locus
     *
     * For more information, see:
     *     https://www.gencodegenes.org/data_format.html
     *     https://en.wikipedia.org/wiki/General_feature_format
     */
    public enum LocusLevel {
        /** Verified locus */
        VERIFIED("1"),

        /** Manually annotated locus */
        MANUALLY_ANNOTATED("2"),

        /** Automatically annotated locus */
        AUTOMATICALLY_ANNOTATED("3");

        @SuppressWarnings("unchecked")
        private static final Map<String, LocusLevel> VALUE_MAP =
                Arrays.stream(values()).collect(Collectors.toMap(v -> v.serialized.toLowerCase(), v -> v));

        private final String serialized;

        LocusLevel(final String serializedValue) {
            serialized = serializedValue;
        }

        @Override
        public String toString() {
            return serialized;
        }

        public static LocusLevel getEnum(final String s) {
            final String lowerS = s.toLowerCase();
            if ( VALUE_MAP.containsKey(lowerS) ){
                return VALUE_MAP.get(lowerS);
            }
            throw new IllegalArgumentException("Unexpected value: " + s);
        }
    }

    /**
     * Additional relevant information appended to a feature.
     *
     * For more information, see:
     *     https://www.gencodegenes.org/data_format.html
     *     https://en.wikipedia.org/wiki/General_feature_format
     *     https://www.gencodegenes.org/pages/tags.html
     */
    public enum FeatureTag {
        /** 3' end extended based on RNA-seq data. */
        THREE_PRIME_NESTED_SUPPORTED_EXTENSION("3_nested_supported_extension"),

        /** 3' end extended based on RNA-seq data. */
        THREE_PRIME_STANDARD_SUPPORTED_EXTENSION("3_standard_supported_extension"),

        /** annotated based on RNA-seq data. */
        FOURFIVEFOUR_RNA_SEQ_SUPPORTED("454_RNA_Seq_supported"),

        /** 5' end extended based on RNA-seq data. */
        FIVE_PRIME_NESTED_SUPPORTED_EXTENSION("5_nested_supported_extension"),

        /** 5' end extended based on RNA-seq data. */
        FIVE_PRIME_STANDARD_SUPPORTED_EXTENSION("5_standard_supported_extension"),

        /** shares an identical CDS but has alternative 5' UTR with respect to a reference variant. */
        ALTERNATIVE_3_UTR("alternative_3_UTR"),

        /** shares an identical CDS but has alternative 3' UTR with respect to a reference variant. */
        ALTERNATIVE_5_UTR("alternative_5_UTR"),

        // --------------------------------------------------------------------------------------------------------
        // Please note that the ordering of the APPRIS_* tags is also used in sorting here.  Do not re-order!
        // --------------------------------------------------------------------------------------------------------
        /** Transcript expected to code for the main functional isoform based on a range of protein features (APPRIS pipeline). */
        APPRIS_PRINCIPAL("appris_principal"),

        /** (This flag corresponds to the older flag "appris_principal") Where the transcript expected to code for the main */
        APPRIS_PRINCIPAL_1("appris_principal_1"),

        /** (This flag corresponds to the older flag "appris_candidate_ccds") Where the APPRIS core modules are unable to choose a */
        APPRIS_PRINCIPAL_2("appris_principal_2"),

        /** Where the APPRIS core modules are unable to choose a clear principal variant and there more than one of the variants */
        APPRIS_PRINCIPAL_3("appris_principal_3"),

        /** (This flag corresponds to the Ensembl 78 flag "appris_candidate_longest_ccds") Where the APPRIS core modules are unable */
        APPRIS_PRINCIPAL_4("appris_principal_4"),

        /** (This flag corresponds to the Ensembl 78 flag "appris_candidate_longest_seq") Where the APPRIS core modules are unable */
        APPRIS_PRINCIPAL_5("appris_principal_5"),

        /** Candidate transcript(s) models that are conserved in at least three tested non-primate species. */
        APPRIS_ALTERNATIVE_1("appris_alternative_1"),

        /** Candidate transcript(s) models that appear to be conserved in fewer than three tested non-primate species. */
        APPRIS_ALTERNATIVE_2("appris_alternative_2"),

        /** where there is no 'appris_principal' variant, the candidate with highest APPRIS score is selected as the primary */
        APPRIS_CANDIDATE_HIGHEST_SCORE("appris_candidate_highest_score"),

        /** the "appris_candidate" transcripts where there are several CCDS, in this case APPRIS labels the longest CCDS. */
        APPRIS_CANDIDATE_LONGEST_CCDS("appris_candidate_longest_ccds"),

        /** the "appris_candidate" transcript that has an unique CCDS. */
        APPRIS_CANDIDATE_CCDS("appris_candidate_ccds"),

        /** where there is no "appris_candidate_ccds" or "appris_candidate_longest_ccds" variant, the longest protein of the */
        APPRIS_CANDIDATE_LONGEST_SEQ("appris_candidate_longest_seq"),

        /** where there is no 'appris_principal' variant, the longest of the 'appris_candidate' variants is selected as the primary */
        APPRIS_CANDIDATE_LONGEST("appris_candidate_longest"),

        /** where there is no single 'appris_principal' variant the main functional isoform will be translated from one of the */
        APPRIS_CANDIDATE("appris_candidate"),

        /** identifies a subset of representative transcripts for each gene; prioritises full-length protein coding transcripts */
        BASIC("basic"),

        /** Transcript contains two confidently annotated CDSs. Support may come from eg proteomic data, cross-species conservation */
        BICISTRONIC("bicistronic"),

        /** Transcript 5' end overlaps ENCODE or Fantom CAGE cluster. */
        CAGE_SUPPORTED_TSS("CAGE_supported_TSS"),

        /** member of the consensus CDS gene set, confirming coding regions between ENSEMBL, UCSC, NCBI and HAVANA. */
        CCDS("CCDS"),

        /** The coding region end could not be confirmed. */
        CDS_END_NF("cds_end_NF"),

        /** The coding region start could not be confirmed. */
        CDS_START_NF("cds_start_NF"),

        /** Transcript QC checked using dotplot to identify features eg splice junctions, end of homology. */
        DOTTER_CONFIRMED("dotter_confirmed"),

        /** an upstream ATG is used where a downstream ATG seems more evolutionary conserved. */
        DOWNSTREAM_ATG("downstream_ATG"),

        /** Transcript was tested and confirmed experimentally. */
        EXP_CONF("exp_conf"),

        /** locus consists of non-overlapping transcript fragments either because of genome assembly issues (i.e., gaps or */
        FRAGMENTED_LOCUS("fragmented_locus"),

        /** Transcript model contains all possible in-frame exons supported by homology, experimental evidence or conservation, but */
        INFERRED_EXON_COMBINATION("inferred_exon_combination"),

        /** Transcript model is not supported by a single piece of transcript evidence. May be supported by multiple fragments of */
        INFERRED_TRANSCRIPT_MODEL("inferred_transcript_model"),

        /** Transcript supported by transcript evidence that, while ampping best-in-genome, shows regions of poor sequence quality. */
        LOW_SEQUENCE_QUALITY("low_sequence_quality"),

        /** the mRNA end could not be confirmed. */
        MRNA_END_NF("mRNA_end_NF"),

        /** the mRNA start could not be confirmed. */
        MRNA_START_NF("mRNA_start_NF"),

        /** the transcript belongs to the MANE Select Plus Clinical data set. The Matched Annotation from NCBI and EMBL-EBI project (MANE) is a collaboration between Ensembl-GENCODE and RefSeq to select a default transcript per human protein coding locus that is representative of biology, well-supported, expressed and conserved. This transcript set matches GRCh38 and is 100% identical between RefSeq and Ensembl-GENCODE for 5' UTR, CDS, splicing and 3' UTR. The Plus Clinical transcripts are chosen to supplement MANE Select when needed for clinical variant reporting.*/
        MANE_PLUS_CLINICAL("MANE_Plus_Clinical"),

        /** the transcript belongs to the MANE Select data set. The Matched Annotation from NCBI and EMBL-EBI project (MANE) is a collaboration between Ensembl-GENCODE and RefSeq to select a default transcript per human protein coding locus that is representative of biology, well-supported, expressed and conserved. This transcript set matches GRCh38 and is 100% identical between RefSeq and Ensembl-GENCODE for 5' UTR, CDS, splicing and 3' UTR. */
        MANE_SELECT("MANE_Select"),

        /** in-frame type of variation where, at the acceptor site, some variants splice after the first AG and others after the */
        NAGNAG_SPLICE_SITE("NAGNAG_splice_site"),

        /** the locus is a host for small non-coding RNAs. */
        NCRNA_HOST("ncRNA_host"),

        /** annotated based on RNA-seq data. */
        NESTED_454_RNA_SEQ_SUPPORTED("nested_454_RNA_Seq_supported"),

        /** the transcript looks like it is subject to NMD but publications, experiments or conservation support the translation of */
        NMD_EXCEPTION("NMD_exception"),

        /** codon if the transcript were longer but cannot currently be annotated as NMD as does not fulfil all criteria - most */
        NMD_LIKELY_IF_EXTENDED("NMD_likely_if_extended"),

        /** the CDS has a non-ATG start and its validity is supported by publication or conservation. */
        NON_ATG_START("non_ATG_start"),

        /** the transcript has a non-canonical splice site conserved in other species. */
        NON_CANONICAL_CONSERVED("non_canonical_conserved"),

        /** the transcript has a non-canonical splice site explained by a genomic sequencing error. */
        NON_CANONICAL_GENOME_SEQUENCE_ERROR("non_canonical_genome_sequence_error"),

        /** the transcript has a non-canonical splice site explained by other reasons. */
        NON_CANONICAL_OTHER("non_canonical_other"),

        /** the transcript has a non-canonical splice site explained by a SNP. */
        NON_CANONICAL_POLYMORPHISM("non_canonical_polymorphism"),

        /** the transcript has a non-canonical splice site that needs experimental confirmation. */
        NON_CANONICAL_TEC("non_canonical_TEC"),

        /** the transcript has a non-canonical splice site explained by a U12 intron (i.e. AT-AC splice site). */
        NON_CANONICAL_U12("non_canonical_U12"),

        /** a splice variant for which supporting evidence has not been submitted to databases, i.e. the model is based on */
        NON_SUBMITTED_EVIDENCE("non_submitted_evidence"),

        /** a transcript is supported by evidence from same species paralogous loci. */
        NOT_BEST_IN_GENOME_EVIDENCE("not_best_in_genome_evidence"),

        /** evidence from other species was used to build model. */
        NOT_ORGANISM_SUPPORTED("not_organism_supported"),

        /** protein-coding locus with no paralogues or orthologs. */
        ORPHAN("orphan"),

        /** exon(s) of the locus overlap exon(s) of a readthrough transcript or a transcript belonging to another locus. */
        OVERLAPPING_LOCUS("overlapping_locus"),

        /** a low confidence upstream ATG existing in other coding variant would lead to NMD in this trancript, that uses the high */
        OVERLAPPING_UORF("overlapping_uORF"),

        /** annotation in the pseudo-autosomal region, which is duplicated between chromosomes X and Y. */
        PAR("PAR"),

        /** member of the pseudogene set predicted by YALE, UCSC and HAVANA. */
        PSEUDO_CONSENS("pseudo_consens"),

        /** a transcript that overlaps two or more independent loci but is considered to belong to a third, separate locus. */
        READTHROUGH_TRANSCRIPT("readthrough_transcript"),

        /** locus overlaps a sequence error or an assembly error in the reference genome that affects its annotation (e.g., 1 or */
        REFERENCE_GENOME_ERROR("reference_genome_error"),

        /** internal intron of CDS portion of transcript is retained. */
        RETAINED_INTRON_CDS("retained_intron_CDS"),

        /** final intron of CDS portion of transcript is retained. */
        RETAINED_INTRON_FINAL("retained_intron_final"),

        /** first intron of CDS portion of transcript is retained. */
        RETAINED_INTRON_FIRST("retained_intron_first"),

        /** protein-coding locus created via retrotransposition. */
        RETROGENE("retrogene"),

        /** Transcript supported by RNAseq data and not supported by mRNA or EST evidence. */
        RNA_SEQ_SUPPORTED_ONLY("RNA_Seq_supported_only"),

        /** Transcript annotated based on mixture of RNA-seq data and EST/mRNA/protein evidence. */
        RNA_SEQ_SUPPORTED_PARTIAL("RNA_Seq_supported_partial"),

        /** Transcript that contains a CDS that has a translation initiation site supported by Ribosomal Profiling data. */
        RP_SUPPORTED_TIS("RP_supported_TIS"),

        /** contains a selenocysteine. */
        SELENO("seleno"),

        /** a processed pseudogene with one or more introns still present. These are likely formed through the retrotransposition */
        SEMI_PROCESSED("semi_processed"),

        /** Transcript contains at least 1 non-canonical splice junction that is associated with a known or novel genome sequence */
        SEQUENCE_ERROR("sequence_error"),

        /** Transcript whose coding sequence contains an internal stop codon that does not cause the translation termination. */
        STOP_CODON_READTHROUGH("stop_codon_readthrough"),

        /** Transcript created or extended using assembled RNA-seq long reads. */
        TAGENE("TAGENE"),

        /** an upstream ATG exists when a downstream ATG is better supported. */
        UPSTREAM_ATG("upstream_ATG"),

        /** a low confidence upstream ATG existing in other coding variant would lead to NMD in this trancript, that uses the high */
        UPSTREAM_UORF("upstream_uORF"),

        /**most representative transcript of the gene. This will be the MANE_Select transcript if there is one, or a transcript chosen by an Ensembl algorithm otherwise. */
        ENSEMBL_CANNONICAL("Ensembl_canonical"),

        /** protein-coding gene that has a readthrough transcript. */
        READTHROUGH_GENE("readthrough_gene"),

        /** annotated on an artifactual duplicate region of the genome assembly. */
        ARTIFACTUAL_DUPLICATION("artifactual_duplication");

        @SuppressWarnings("unchecked")
        private static final Map<String, FeatureTag> VALUE_MAP =
                Arrays.stream(values()).collect(Collectors.toMap(v -> v.serialized.toLowerCase(), v -> v));

        private final String serialized;

        FeatureTag(final String serializedValue) {
            serialized = serializedValue;
        }

        @Override
        public String toString() {
            return serialized;
        }

        public static FeatureTag getEnum(final String s) {
            final String lowerS = s.toLowerCase();
            if ( VALUE_MAP.containsKey(lowerS) ){
                return VALUE_MAP.get(lowerS);
            }
            throw new IllegalArgumentException("Unexpected value: " + s);
        }
    }

    /**
     * Transcript score according to how well mRNA and EST alignments match over its full length.
     *
     * For more information, see:
     *     https://www.gencodegenes.org/data_format.html
     *     https://en.wikipedia.org/wiki/General_feature_format
     */
    public enum TranscriptSupportLevel {
        /** all splice junctions of the transcript are supported by at least one non-suspect mRNA */
        ALL_MRNA_VERIFIED("1"),

        /** the best supporting mRNA is flagged as suspect or the support is from multiple ESTs */
        BEST_MRNA_SUSPECT("2"),

        /** the only support is from a single EST */
        SINGLE_EST_SUPPORT("3"),

        /** the best supporting EST is flagged as suspect */
        BEST_EST_SUSPECT("4"),

        /** no single transcript supports the model structure */
        NO_SINGLE_TRANSCRIPT_SUPPORT("5"),

        /** the transcript was not analyzed */
        NA("NA");

        @SuppressWarnings("unchecked")
        private static final Map<String, TranscriptSupportLevel> VALUE_MAP =
                Arrays.stream(values()).collect(Collectors.toMap(v -> v.serialized.toLowerCase(), v -> v));

        private final String serialized;

        TranscriptSupportLevel(final String serializedValue) {
            serialized = serializedValue;
        }

        @Override
        public String toString() {
            return serialized;
        }

        public static TranscriptSupportLevel getEnum(final String s) {
            final String lowerS = s.toLowerCase();
            if ( VALUE_MAP.containsKey(lowerS) ){
                return VALUE_MAP.get(lowerS);
            }
            throw new IllegalArgumentException("Unexpected value: " + s);
        }
    }

    /**
     * Attribute that indicates the status of the mapping.
     *
     * For more information, see:
     *     https://www.gencodegenes.org/data_format.html
     *     https://en.wikipedia.org/wiki/General_feature_format
     *     http://www.gencodegenes.org/releases/grch37_mapped_releases.html#attrib
     */
    public enum RemapStatus {
        /**
         * Gene or transcript completely mapped to the target genome with all features intact.
         */
        FULL_CONTIG("full_contig"),

        /**
         * Gene or transcript completely mapped to the target genome with insertions in some features. These are usually small insertions.
         */
        FULL_FRAGMENT("full_fragment"),

        /**
         * Gene or transcript partially mapped to the target genome.
         */
        PARTIAL("partial"),

        /**
         * Gene or transcript did not map to the target genome.
         */
        DELETED("deleted"),

        /**
         * The source sequence is not in the assembly alignments. This will occur with alt loci genes if the alignments only contain the primary assembly.
         */
        NO_SEQ_MAP("no_seq_map"),

        /**
         * Transcripts in the gene mapped to multiple locations.
         */
        GENE_CONFLICT("gene_conflict"),

        /**
         * Transcripts caused gene length to change by more than 50%. This is to detect mapping to processed pseudogenes and mapping across tandem gene duplications.
         */
        GENE_SIZE_CHANGE("gene_size_change"),

        /**
         * Gene is from a small, automatic (ENSEMBL source) non-coding RNA. Taken from the target annotation.
         */
        AUTOMATIC_SMALL_NCRNA_GENE("automatic_small_ncrna_gene"),

        /**
         * Gene is from an automatic process (ENSEMBL source). Taken from the target annotation.
         */
        AUTOMATIC_GENE("automatic_gene"),

        /**
         * Pseudogene annotations (excluding polymorphic).
         */
        PSEUDOGENE("pseudogene");

        @SuppressWarnings("unchecked")
        private static final Map<String, RemapStatus> VALUE_MAP =
                Arrays.stream(values()).collect(Collectors.toMap(v -> v.serialized.toLowerCase(), v -> v));

        private final String serialized;

        RemapStatus(final String serializedValue) { serialized = serializedValue; }

        @Override
        public String toString() {
            return serialized;
        }

        public static RemapStatus getEnum(final String s) {
            final String lowerS = s.toLowerCase();
            if ( VALUE_MAP.containsKey(lowerS) ){
                return VALUE_MAP.get(lowerS);
            }
            throw new IllegalArgumentException("Unexpected value: " + s);
        }
    }

    /**
     * Attribute that compares the mapping to the existing target annotations.
     *
     * For more information, see:
     *     https://www.gencodegenes.org/data_format.html
     *     https://en.wikipedia.org/wiki/General_feature_format
     *     http://www.gencodegenes.org/releases/grch37_mapped_releases.html#attrib
     */
    public enum RemapTargetStatus {

        /**
         * Gene or transcript was not in target annotations.
         */
        NEW("new"),

        /**
         * Gene or transcript exists in source and target genome, however source was not mapped.
         */
        LOST("lost"),

        /**
         * Gene or transcript overlaps previous version of annotation on target genome.
         */
        OVERLAP("overlap"),

        /**
         * Gene or transcript exists in target, however source mapping is to a different location. This is often mappings to a gene family members or pseudogenes.
         */
        NONOVERLAP("nonOverlap");

        @SuppressWarnings("unchecked")
        private static final Map<String, RemapTargetStatus> VALUE_MAP =
                Arrays.stream(values()).collect(Collectors.toMap(v -> v.serialized.toLowerCase(), v -> v));

        private final String serialized;

        RemapTargetStatus(final String serializedValue) {
            serialized = serializedValue;
        }

        @Override
        public String toString() {
            return serialized;
        }

        public static RemapTargetStatus getEnum(final String s) {
            final String lowerS = s.toLowerCase();
            if ( VALUE_MAP.containsKey(lowerS) ){
                return VALUE_MAP.get(lowerS);
            }
            throw new IllegalArgumentException("Unexpected value: " + s);
        }
    }
}
