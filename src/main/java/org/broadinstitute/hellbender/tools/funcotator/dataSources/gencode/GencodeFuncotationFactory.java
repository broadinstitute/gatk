package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.*;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.gencode.*;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.collections.Sets;

import java.nio.file.Path;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeature.FeatureTag.*;

/**
 * A factory to create {@link GencodeFuncotation}s.
 * This is a high-level object that interfaces with the internals of {@link Funcotator}.
 * Created by jonn on 8/30/17.
 */
public class GencodeFuncotationFactory extends DataSourceFuncotationFactory {

    //==================================================================================================================
    // Private Static Members:
    /** Standard Logger.  */
    protected static final Logger logger = LogManager.getLogger(GencodeFuncotationFactory.class);

    /**
     * The window around splice sites to mark variants as {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification#SPLICE_SITE}.
     */
    private static final int spliceSiteVariantWindowBases = 2;

    /**
     * Number of bases to the left and right of a variant in which to calculate the GC content.
     */
    private static final int gcContentWindowSizeBases = 200;

    /**
     * The window around a variant to include in the reference context annotation.
     * Also used for context from which to get surrounding codon changes and protein changes.
     */
    final static private int referenceWindow = 10;

    /**
     * List of valid Appris Ranks used for sorting funcotations to get the "best" one.z
     */
    private static final HashSet<GencodeGtfGeneFeature.FeatureTag> apprisRanks = new HashSet<>(
            Arrays.asList(
                APPRIS_PRINCIPAL,
                APPRIS_PRINCIPAL_1,
                APPRIS_PRINCIPAL_2,
                APPRIS_PRINCIPAL_3,
                APPRIS_PRINCIPAL_4,
                APPRIS_PRINCIPAL_5,
                APPRIS_ALTERNATIVE_1,
                APPRIS_ALTERNATIVE_2,
                APPRIS_CANDIDATE_HIGHEST_SCORE,
                APPRIS_CANDIDATE_LONGEST_CCDS,
                APPRIS_CANDIDATE_CCDS,
                APPRIS_CANDIDATE_LONGEST_SEQ,
                APPRIS_CANDIDATE_LONGEST,
                APPRIS_CANDIDATE
            )
    );

    /**
     * The set of {@link GencodeFuncotation.VariantClassification} types that are valid for coding regions.
     */
    private static final Set<GencodeFuncotation.VariantClassification> codingVariantClassifications =
            Sets.newHashSet(Arrays.asList(GencodeFuncotation.VariantClassification.MISSENSE,
                            GencodeFuncotation.VariantClassification.NONSENSE,
                            GencodeFuncotation.VariantClassification.NONSTOP,
                            GencodeFuncotation.VariantClassification.SILENT,
                            GencodeFuncotation.VariantClassification.IN_FRAME_DEL,
                            GencodeFuncotation.VariantClassification.IN_FRAME_INS,
                            GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS,
                            GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL,
                            GencodeFuncotation.VariantClassification.START_CODON_SNP,
                            GencodeFuncotation.VariantClassification.START_CODON_INS,
                            GencodeFuncotation.VariantClassification.START_CODON_DEL));


    //==================================================================================================================
    // Private Members:

    /**
     * ReferenceSequenceFile for the transcript reference file.
     */
    private final ReferenceDataSource transcriptFastaReferenceDataSource;

    /**
     * Map between transcript IDs and the IDs from the FASTA file to look up the transcript.
     * This is necessary because of the way the FASTA file contigs are named.
     */
    private final Map<String, MappedTranscriptIdInfo> transcriptIdMap;

    /**
     * The mode to select the "best" transcript (i.e. the transcript with detailed information) from the list of
     * possible transcripts.
     *
     * For more information on the specifics of the differences go here:
     * https://gatkforums.broadinstitute.org/gatk/discussion/4220/what-is-the-difference-between-tx-mode-best-effect-vs-canonical
     */
    private final FuncotatorArgumentDefinitions.TranscriptSelectionMode transcriptSelectionMode;

    /**
     * {@link List} of Transcript IDs that the user has requested that we annotate.
     * If this list is empty, will default to keeping ALL transcripts.
     * Otherwise, only transcripts on this list will be annotated.
     */
    private final Set<String> userRequestedTranscripts;

    /**
     * The {@link Path} from which we will read the sequences for the coding regions in given transcripts.
     */
    private final Path gencodeTranscriptFastaFile;

    //==================================================================================================================
    // Constructors:

    public GencodeFuncotationFactory(final Path gencodeTranscriptFastaFile, final String version) {
        this(gencodeTranscriptFastaFile, version, FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE, new HashSet<>(), new LinkedHashMap<>());
    }

    public GencodeFuncotationFactory(final Path gencodeTranscriptFastaFile, final String version, final Set<String> userRequestedTranscripts) {
        this(gencodeTranscriptFastaFile, version, FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE, userRequestedTranscripts, new LinkedHashMap<>());
    }

    public GencodeFuncotationFactory(final Path gencodeTranscriptFastaFile, final String version,final FuncotatorArgumentDefinitions.TranscriptSelectionMode transcriptSelectionMode) {
        this(gencodeTranscriptFastaFile, version, transcriptSelectionMode, new HashSet<>(), new LinkedHashMap<>());
    }

    public GencodeFuncotationFactory(final Path gencodeTranscriptFastaFile,
                                     final String version,
                                     final FuncotatorArgumentDefinitions.TranscriptSelectionMode transcriptSelectionMode,
                                     final Set<String> userRequestedTranscripts) {
        this(gencodeTranscriptFastaFile, version, transcriptSelectionMode, userRequestedTranscripts, new LinkedHashMap<>());
    }

    public GencodeFuncotationFactory(final Path gencodeTranscriptFastaFile,
                                     final String version,
                                     final FuncotatorArgumentDefinitions.TranscriptSelectionMode transcriptSelectionMode,
                                     final Set<String> userRequestedTranscripts,
                                     final LinkedHashMap<String, String> annotationOverrides) {

        this.gencodeTranscriptFastaFile = gencodeTranscriptFastaFile;

        transcriptFastaReferenceDataSource = ReferenceDataSource.of(gencodeTranscriptFastaFile);
        transcriptIdMap = createTranscriptIdMap(transcriptFastaReferenceDataSource);

        this.transcriptSelectionMode = transcriptSelectionMode;

        this.version = version;

        // Go through each requested transcript and remove the version numbers from them if they exist:
        this.userRequestedTranscripts = new HashSet<>();
        for ( final String transcript : userRequestedTranscripts ) {
            this.userRequestedTranscripts.add( getTranscriptIdWithoutVersionNumber(transcript) );
        }

        // Initialize overrides / defaults:
        initializeAnnotationOverrides( annotationOverrides );
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public void close() {
        transcriptFastaReferenceDataSource.close();
    }

    @Override
    public String getName() {
        return "Gencode";
    }

    //TODO: Add in version getter/setter:

    @Override
    public LinkedHashSet<String> getSupportedFuncotationFields() {

            return new LinkedHashSet<>(Arrays.asList(
                    getName() + "_" + getVersion() + "_hugoSymbol",
                    getName() + "_" + getVersion() + "_ncbiBuild",
                    getName() + "_" + getVersion() + "_chromosome",
                    getName() + "_" + getVersion() + "_start",
                    getName() + "_" + getVersion() + "_end",
                    getName() + "_" + getVersion() + "_variantClassification",
                    getName() + "_" + getVersion() + "_secondaryVariantClassification",
                    getName() + "_" + getVersion() + "_variantType",
                    getName() + "_" + getVersion() + "_refAllele",
                    getName() + "_" + getVersion() + "_tumorSeqAllele1",
                    getName() + "_" + getVersion() + "_tumorSeqAllele2",
                    getName() + "_" + getVersion() + "_genomeChange",
                    getName() + "_" + getVersion() + "_annotationTranscript",
                    getName() + "_" + getVersion() + "_transcriptStrand",
                    getName() + "_" + getVersion() + "_transcriptExon",
                    getName() + "_" + getVersion() + "_transcriptPos",
                    getName() + "_" + getVersion() + "_cDnaChange",
                    getName() + "_" + getVersion() + "_codonChange",
                    getName() + "_" + getVersion() + "_proteinChange",
                    getName() + "_" + getVersion() + "_gcContent",
                    getName() + "_" + getVersion() + "_referenceContext",
                    getName() + "_" + getVersion() + "_otherTranscripts"
            ));
    }

    @Override
    /**
     * Attempts to treat the given features as {@link GencodeGtfFeature} objects in order to
     * create funcotations for the given variant and reference.
     */
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {
        final List<Funcotation> outputFuncotations = new ArrayList<>();

        // If we have features we need to annotate, go through them and create annotations:
        if ( featureList.size() > 0 ) {
            for ( final Allele altAllele : variant.getAlternateAlleles() ) {
                for ( final Feature feature : featureList ) {

                    // Get the kind of feature we want here:
                    if ( (feature != null) && GencodeGtfGeneFeature.class.isAssignableFrom(feature.getClass()) ) {
                        final List<GencodeFuncotation> gencodeFuncotationList = createFuncotations(variant, altAllele, (GencodeGtfGeneFeature) feature, referenceContext);

                        // Now we have to filter out the output gencodeFuncotations if they are not on the list the user provided:
                        filterAnnotationsByUserTranscripts( gencodeFuncotationList, userRequestedTranscripts );

                        // Add the filtered funcotations here:
                        outputFuncotations.addAll(gencodeFuncotationList);
                    }
                    // TODO: Actually you may want to put another IGR creation here for now...  This may be a more difficult thing if we determine it in here.  There is no way to know if these are IGRs or simply not included in this particular data set.
                }
            }
        }
        else {
            // This is an IGR.  Only bother with it if the User has not asked for a specific transcript (because IGRs
            // by definition have no associated transcript).
            if ( userRequestedTranscripts.size() == 0 ) {
                outputFuncotations.addAll(createIgrFuncotations(variant, referenceContext));
            }
        }

        // Now we set the override values for each annotation:
        for ( final Funcotation funcotation : outputFuncotations ) {
            funcotation.setFieldSerializationOverrideValues( annotationOverrideMap );
        }

        return outputFuncotations;
    }

    @Override
    /**
     * {@inheritDoc}
     * This method should never be called on a {@link GencodeFuncotationFactory}, as that would imply a time-loop.
     */
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList, final List<GencodeFuncotation> gencodeFuncotations) {
        throw new GATKException("This method should never be called on a "+ this.getClass().getName());
    }

    @Override
    public FuncotatorArgumentDefinitions.DataSourceType getType() {
        return FuncotatorArgumentDefinitions.DataSourceType.GENCODE;
    }

    //==================================================================================================================
    // Static Methods:

    /**
     * Filter the given list of {@link GencodeFuncotation} to only contain those funcotations that have transcriptIDs that
     * appear in the given {@code selectedTranscripts}.
     * Ignores transcript version numbers.
     * @param funcotations The {@link List} of {@link GencodeFuncotation} to filter.
     * @param selectedTranscripts The {@link Set} of transcript IDs to keep in the given {@code funcotations}.
     */
    static void filterAnnotationsByUserTranscripts( final List<GencodeFuncotation> funcotations,
                                                    final Set<String> selectedTranscripts ) {
        if ( (selectedTranscripts != null) && (!selectedTranscripts.isEmpty()) ) {
            funcotations.removeIf( f -> !isFuncotationInTranscriptList(f, selectedTranscripts) );
        }
    }

    /**
     * Determines whether the given {@code funcotation} has a transcript ID that is in the given {@code acceptableTranscripts}.
     * Ignores transcript version numbers.
     * @param funcotation The {@link GencodeFuncotation} to check against the set of {@code acceptableTranscripts}.
     * @param acceptableTranscripts The {@link Set} of transcript IDs that are OK to keep.
     * @return {@code true} if funcotation.annotationTranscript is in {@code acceptableTranscripts} (ignoring transcript version); {@code false} otherwise.
     */
    static boolean isFuncotationInTranscriptList( final GencodeFuncotation funcotation,
                                                  final Set<String> acceptableTranscripts ) {
        if ( funcotation.getAnnotationTranscript() != null ) {
            return acceptableTranscripts.contains( getTranscriptIdWithoutVersionNumber(funcotation.getAnnotationTranscript()) );
        }
        else {
            return false;
        }
    }

    /**
     * Removes the transcript ID version number from the given transcript ID (if it exists).
     * @param transcriptId The transcript from which to remove the version number.
     * @return The {@link String} corresponding to the given {@code transcriptId} without a version number.
     */
    static String getTranscriptIdWithoutVersionNumber( final String transcriptId ) {
        return transcriptId.replaceAll("\\.\\d+$", "");
    }

    /**
     * Creates a map of Transcript IDs for use in looking up transcripts from the FASTA dictionary for the GENCODE Transcripts.
     * We include the start and stop codons in the transcripts so we can handle start/stop codon variants.
     * @param fastaReference The {@link ReferenceDataSource} corresponding to the Transcript FASTA file for this GENCODE dataset.
     * @return A {@link Map} of {@link String} -> {@link MappedTranscriptIdInfo} which maps real transcript IDs to the information about that transcript in the transcript FASTA file.
     */
    @VisibleForTesting
    static Map<String, MappedTranscriptIdInfo> createTranscriptIdMap(final ReferenceDataSource fastaReference) {

        final Map<String, MappedTranscriptIdInfo> idMap = new HashMap<>();

        for ( final SAMSequenceRecord sequence : fastaReference.getSequenceDictionary().getSequences() ) {

            final MappedTranscriptIdInfo transcriptInfo = createMappedTranscriptIdInfo( sequence );

            // The names in the file are actually in a list with | between each sequence name.
            // We need to split the names and add them to the dictionary so we can resolve them to the full
            // sequence name as it appears in the file:
            for ( final String transcriptId : sequence.getSequenceName().split("\\|") ) {
                idMap.put(transcriptId, transcriptInfo);
            }
        }

        return idMap;
    }

    /**
     * Creates a {@link MappedTranscriptIdInfo} object based on the given {@link SAMSequenceRecord}.
     * This method is a helper method to get information out of a GENCODE transcript FASTA file easily.
     * This method assumes that {@code sequence} is a sequence from a GENCODE transcript FASTA file.
     * @param sequence The {@link SAMSequenceRecord} from which to create the {@link MappedTranscriptIdInfo}.
     * @return A populated {@link MappedTranscriptIdInfo} object based on the given {@link SAMSequenceRecord}.
     */
    private static MappedTranscriptIdInfo createMappedTranscriptIdInfo( final SAMSequenceRecord sequence ) {

        final MappedTranscriptIdInfo transcriptIdInfo = new MappedTranscriptIdInfo();

        final Pattern utrPattern = Pattern.compile("UTR[35]:(\\d+)-(\\d+)");
        final Pattern cdsPattern = Pattern.compile("CDS:(\\d+)-(\\d+)");

        boolean has3pUtr = false;
        boolean has5pUtr = false;

        // Now let's go through the sequence name and pull out the salient features for each field:
        for (final String field : sequence.getSequenceName().split("\\|")) {
            if ((field.length() > 4) && (field.substring(0, 5).equals("UTR5:"))) {
                final Matcher m = utrPattern.matcher(field);
                m.find();
                transcriptIdInfo.fivePrimeUtrStart = Integer.valueOf(m.group(1));
                transcriptIdInfo.fivePrimeUtrEnd = Integer.valueOf(m.group(2));
                has5pUtr = true;
            } else if ((field.length() > 4) && (field.substring(0, 5).equals("UTR3:"))) {
                final Matcher m = utrPattern.matcher(field);
                m.find();
                transcriptIdInfo.threePrimeUtrStart = Integer.valueOf(m.group(1));
                transcriptIdInfo.threePrimeUtrEnd = Integer.valueOf(m.group(2));
                has3pUtr = true;
            } else if ((field.length() > 3) && (field.substring(0, 4).equals("CDS:"))) {
                final Matcher m = cdsPattern.matcher(field);
                m.find();
                transcriptIdInfo.codingSequenceStart = Integer.valueOf(m.group(1));
                transcriptIdInfo.codingSequenceEnd = Integer.valueOf(m.group(2));
            }
        }

        transcriptIdInfo.mapKey = sequence.getSequenceName();
        transcriptIdInfo.has3pUtr = has3pUtr;
        transcriptIdInfo.has5pUtr = has5pUtr;

        return transcriptIdInfo;
    }

    /**
     * Get the coding sequence from the GENCODE Transcript FASTA file for a given {@code transcriptId}.
     * This will get ONLY the coding sequence for the given {@code transcriptId} and will not include any UTRs.
     * @param transcriptId The ID of the transcript to get from the FASTA file.
     * @param transcriptIdMap A map from transcriptId to MappedTranscriptIdInfo, which tells us how to pull information for the given {@code transcriptId} out of the given {@code transcriptFastaReferenceDataSource}.
     * @param transcriptFastaReferenceDataSource A {@link ReferenceDataSource} for the GENCODE transcript FASTA file.
     * @return The coding sequence for the given {@code transcriptId} as represented in the GENCODE transcript FASTA file.
     */
    private static String getCodingSequenceFromTranscriptFasta( final String transcriptId,
                                                                final Map<String, MappedTranscriptIdInfo> transcriptIdMap,
                                                                final ReferenceDataSource transcriptFastaReferenceDataSource) {

        final MappedTranscriptIdInfo transcriptMapIdAndMetadata = transcriptIdMap.get(transcriptId);

        if ( transcriptMapIdAndMetadata == null ) {
            throw new UserException.BadInput( "Unable to find the given Transcript ID in our transcript list (not in given transcript FASTA file): " + transcriptId );
        }

        final SimpleInterval transcriptInterval = new SimpleInterval(
                transcriptMapIdAndMetadata.mapKey,
                transcriptMapIdAndMetadata.codingSequenceStart,
                transcriptMapIdAndMetadata.codingSequenceEnd
        );

        return transcriptFastaReferenceDataSource.queryAndPrefetch( transcriptInterval ).getBaseString();
    }

    /**
     * Returns whether a variant is in a coding region based on its primary and secondary {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification}.
     * @param varClass Primary {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} of a variant.
     * @param secondaryVarClass Secondary {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} of a variant.
     * @return {@code true} if the corresponding variant is in a coding region; {@code false} otherwise.
     */
    @VisibleForTesting
    static boolean isVariantInCodingRegion(final GencodeFuncotation.VariantClassification varClass,
                                           final GencodeFuncotation.VariantClassification secondaryVarClass ) {

        Utils.nonNull( varClass );

        if ( varClass == GencodeFuncotation.VariantClassification.SPLICE_SITE ) {
            Utils.nonNull( secondaryVarClass );
            return secondaryVarClass != GencodeFuncotation.VariantClassification.INTRON;
        }
        else {
            return codingVariantClassifications.contains(varClass);
        }
    }

    //==================================================================================================================
    // Normal Methods:

    /**
     * Creates a {@link List} of {@link GencodeFuncotation}s based on the given {@link VariantContext}, {@link Allele}, and {@link GencodeGtfGeneFeature}.
     * @param variant The variant to annotate.
     * @param altAllele The allele of the given variant to annotate.
     * @param gtfFeature The GTF feature on which to base annotations.
     * @return A {@link List} of {@link GencodeFuncotation}s for the given variant, allele and gtf feature.
     */
    @VisibleForTesting
    List<GencodeFuncotation> createFuncotations(final VariantContext variant, final Allele altAllele, final GencodeGtfGeneFeature gtfFeature, final ReferenceContext reference) {
        // For each applicable transcript, create an annotation.

        final List<GencodeFuncotation> outputFuncotations = new ArrayList<>();

        // Go through and annotate all our non-best transcripts:
        final List<String> otherTranscriptsCondensedAnnotations = new ArrayList<>();
        for ( final GencodeGtfTranscriptFeature transcript : gtfFeature.getTranscripts() ) {

            // Check if this transcript has the `basic` tag:
            final boolean isBasic = transcript.getOptionalFields().stream()
                    .filter( f -> f.getName().equals("tag") )
                    .filter( f -> f.getValue() instanceof GencodeGtfFeature.FeatureTag )
                    .filter( f -> f.getValue().equals(GencodeGtfFeature.FeatureTag.BASIC) )
                    .count() > 0;

            // Only annotate on the `basic` transcripts:
            if ( isBasic ) {

                // Make sure that we have the transcript in our list.
                // If we don't, warn the user and do not add any funcotations:
                if ( !transcriptIdMap.containsKey(transcript.getTranscriptId()) ) {
                    logger.warn("Coding sequence for given transcript ID (" + transcript.getTranscriptId() + ") is missing in file(" + gencodeTranscriptFastaFile.toUri().toString() + "): Skipping.");
                    continue;
                }

                // Try to create the annotation:
                try {
                    final GencodeFuncotation gencodeFuncotation = createGencodeFuncotationOnTranscript(variant, altAllele, gtfFeature, reference, transcript);

                    // Add it into our transcript:
                    outputFuncotations.add(gencodeFuncotation);

                }
                catch ( final FuncotatorUtils.TranscriptCodingSequenceException ex ) {
                    //TODO: This should never happen, but needs to be here for some transcripts, such as HG19 MUC16 ENST00000599436.1, where the transcript sequence itself is of length not divisible by 3! (3992)
                    //      There may be other erroneous transcripts too.
                    otherTranscriptsCondensedAnnotations.add("ERROR_ON_" + transcript.getTranscriptId());

                    logger.warn("Unable to create a GencodeFuncotation on transcript " + transcript.getTranscriptId() + " for variant: " +
                            variant.getContig() + ":" + variant.getStart() + "-" + variant.getEnd() + "(" + variant.getReference() + " -> " + altAllele + ")"
                    );
                }
            }
        }

        if (outputFuncotations.size() > 0) {
            // Get our "Best Transcript" from our list:
            sortFuncotationsByTranscriptForOutput(outputFuncotations);
            final GencodeFuncotation bestFuncotation = outputFuncotations.remove(0);

            // Now convert the other transcripts into summary strings:
            for ( final GencodeFuncotation funcotation : outputFuncotations ) {
                otherTranscriptsCondensedAnnotations.add(condenseGencodeFuncotation(funcotation));
            }

            // Set our `other transcripts` annotation in our best funcotation:
            bestFuncotation.setOtherTranscripts(otherTranscriptsCondensedAnnotations);

            // Add our best funcotation to the output:
            outputFuncotations.add(bestFuncotation);
        }

        return outputFuncotations;
    }

    /**
     * Create a {@link GencodeFuncotation} for a given variant and transcript.
     * @param variant The {@link VariantContext} to annotate.
     * @param altAllele The alternate {@link Allele} to annotate.
     * @param gtfFeature The corresponding {@link GencodeGtfFeature} from which to create annotations.
     * @param reference The {@link ReferenceContext} for the given {@link VariantContext}.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which the given {@code variant} occurs.
     * @return A {@link GencodeFuncotation}
     */
    private GencodeFuncotation createGencodeFuncotationOnTranscript(final VariantContext variant,
                                                                    final Allele altAllele,
                                                                    final GencodeGtfGeneFeature gtfFeature,
                                                                    final ReferenceContext reference,
                                                                    final GencodeGtfTranscriptFeature transcript) {
        final GencodeFuncotation gencodeFuncotation;

        // We only fully process protein-coding regions.
        // For other gene types, we do trivial processing and label them as either LINCRNA or RNA:
        // TODO: Add more types to Variant Classification to be more explicit and a converter function to go from gene type to variant classification.
        if ( gtfFeature.getGeneType() != GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING ) {

            // Setup the "trivial" fields of the gencodeFuncotation:
            final GencodeFuncotationBuilder gencodeFuncotationBuilder = createGencodeFuncotationBuilderWithTrivialFieldsPopulated(variant, altAllele, gtfFeature, transcript);

            if ( gtfFeature.getGeneType() == GencodeGtfFeature.GeneTranscriptType.LINCRNA ) {
                gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.LINCRNA);
            }
            else {
                gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.RNA);
            }

            // Get GC Content:
            gencodeFuncotationBuilder.setGcContent( calculateGcContent( reference, gcContentWindowSizeBases ) );

            gencodeFuncotation = gencodeFuncotationBuilder.build();
        }
        else {

            // TODO: check for complex INDEL and warn and skip.

            // Find the sub-feature of transcript that contains our variant:
            final GencodeGtfFeature containingSubfeature = getContainingGtfSubfeature(variant, transcript);

            // Determine what kind of region we're in and handle it in it's own way:
            if ( containingSubfeature == null ) {
                // We have an IGR variant
                gencodeFuncotation = createIgrFuncotation(variant.getReference(), altAllele, reference);
            }
            else if ( GencodeGtfExonFeature.class.isAssignableFrom(containingSubfeature.getClass()) ) {
                // We have a coding region variant
                gencodeFuncotation = createCodingRegionFuncotation(variant, altAllele, gtfFeature, reference, transcript, (GencodeGtfExonFeature) containingSubfeature);
            }
            else if ( GencodeGtfUTRFeature.class.isAssignableFrom(containingSubfeature.getClass()) ) {
                // We have a UTR variant
                gencodeFuncotation = createUtrFuncotation(variant, altAllele, reference, gtfFeature, transcript, (GencodeGtfUTRFeature) containingSubfeature);
            }
            else if ( GencodeGtfTranscriptFeature.class.isAssignableFrom(containingSubfeature.getClass()) ) {
                // We have an intron variant
                gencodeFuncotation = createIntronFuncotation(variant, altAllele, reference, gtfFeature, transcript, reference);
            }
            else {
                // Uh-oh!  Problemz.
                throw new GATKException.ShouldNeverReachHereException("Unable to determine type of variant-containing subfeature: " + containingSubfeature.getClass().getName());
            }
        }
        return gencodeFuncotation;
    }

    /**
     * Create a {@link GencodeFuncotation} for a {@code variant} that occurs in a coding region in a {@code transcript}.
     * @param variant The {@link VariantContext} for which to create a {@link GencodeFuncotation}.
     * @param altAllele The {@link Allele} in the given {@code variant} for which to create a {@link GencodeFuncotation}.
     * @param gtfFeature The {@link GencodeGtfGeneFeature} in which the given {@code variant} occurs.
     * @param reference The {@link ReferenceContext} for the current data set.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which the given {@code variant} occurs.
     * @param exon The {@link GencodeGtfExonFeature} in which the given {@code variant} occurs.
     * @return A {@link GencodeFuncotation} containing information about the given {@code variant} given the corresponding {@code exon}.
     */
    private GencodeFuncotation createCodingRegionFuncotation(final VariantContext variant,
                                                             final Allele altAllele,
                                                             final GencodeGtfGeneFeature gtfFeature,
                                                             final ReferenceContext reference,
                                                             final GencodeGtfTranscriptFeature transcript,
                                                             final GencodeGtfExonFeature exon) {

        // Setup the "trivial" fields of the gencodeFuncotation:
        final GencodeFuncotationBuilder gencodeFuncotationBuilder = createGencodeFuncotationBuilderWithTrivialFieldsPopulated(variant, altAllele, gtfFeature, transcript);

        // Get the list of exons by their locations so we can use them to determine our location in the transcript and get
        // the transcript code itself:
        final List<? extends Locatable> exonPositionList = getSortedExonAndStartStopPositions(transcript);

        // Set up our SequenceComparison object so we can calculate some useful fields more easily
        // These fields can all be set without knowing the alternate allele:
        final SequenceComparison sequenceComparison = createSequenceComparison(variant, altAllele, reference, transcript, exonPositionList, transcriptIdMap, transcriptFastaReferenceDataSource);

        final GencodeFuncotation.VariantType variantType = getVariantType(variant.getReference(), altAllele);

        // OK, now that we have our SequenceComparison object set up we can continue the annotation:

        // Set the exon number:
        gencodeFuncotationBuilder.setTranscriptExonNumber( exon.getExonNumber() );

        // Set the Codon and Protein changes:
        final String codonChange = FuncotatorUtils.getCodonChangeString( sequenceComparison );
        final String proteinChange = FuncotatorUtils.getProteinChangeString( sequenceComparison );

        gencodeFuncotationBuilder.setCodonChange(codonChange)
                                 .setProteinChange(proteinChange)
                                 .setGcContent( sequenceComparison.getGcContent() );

        // Set the Variant Classification:
        final GencodeFuncotation.VariantClassification varClass = createVariantClassification( variant, altAllele, variantType, exon, sequenceComparison );
        GencodeFuncotation.VariantClassification secondaryVarClass = null;
        gencodeFuncotationBuilder.setVariantClassification( varClass );
        if ( varClass == GencodeFuncotation.VariantClassification.SPLICE_SITE ) {
            secondaryVarClass = getVariantClassificationForCodingRegion(variant, altAllele, variantType, sequenceComparison );
            gencodeFuncotationBuilder.setSecondaryVariantClassification( secondaryVarClass );
        }

        // Set the Coding DNA change:
        final boolean isInCodingRegion = isVariantInCodingRegion( varClass, secondaryVarClass );

        // Only set cDNA change if we have something in the actual coding region:
        if ( isInCodingRegion ) {
            gencodeFuncotationBuilder.setcDnaChange(FuncotatorUtils.getCodingSequenceChangeString(sequenceComparison));
        }

        // Set the reference context with the bases from the sequence comparison:
        gencodeFuncotationBuilder.setReferenceContext( sequenceComparison.getReferenceBases() );

        return gencodeFuncotationBuilder.build();
    }

    /**
     * Gets a list of locatables representing the start codon, exons, and stop codon containing coding regions within the given {@code transcript}.
     * These exons are sorted by exon-number order.
     * @param transcript A {@link GencodeGtfTranscriptFeature} from which to pull the exons.
     * @return A list of {@link Locatable} objects representing the exons in the given {@code transcript} in the order in which the appear in the expressed protein.
     */
    @VisibleForTesting
    static List<? extends Locatable> getSortedExonAndStartStopPositions(final GencodeGtfTranscriptFeature transcript) {

        // Sort by exon number first:
        transcript.getExons().sort((lhs, rhs) -> lhs.getExonNumber() < rhs.getExonNumber() ? -1 : (lhs.getExonNumber() > rhs.getExonNumber() ) ? 1 : 0 );

        final List<Locatable> exonList = new ArrayList<>(transcript.getExons().size());
        for ( final GencodeGtfExonFeature exon : transcript.getExons() ) {

            // Add in a CDS region:
            if ( exon.getCds() != null ) {

                // If we have a start codon that is not in the CDS for some reason,
                // we need to add it to our list:
                if (exon.getStartCodon() != null) {
                    if ( !exon.getCds().contains(exon.getStartCodon()) ) {
                        exonList.add( exon.getStartCodon() );
                    }
                }

                exonList.add( exon.getCds() );

                // If we have a stop codon that is not in the CDS for some reason,
                // we need to add it to our list:
                if (exon.getStopCodon() != null) {
                    if ( !exon.getCds().contains(exon.getStopCodon()) ) {
                        exonList.add( exon.getStopCodon() );
                    }
                }

            }
            else if (exon.getStartCodon() != null) {
                exonList.add( exon.getStartCodon() );
            }
            else if ( exon.getStopCodon() != null ) {
                exonList.add( exon.getStopCodon() );
            }
        }
        return exonList;
    }

    /**
     * Gets the {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} of the given {@code altAllele} for the given {@code variant}.
     * @param variant The {@link VariantContext} to classify.
     * @param altAllele The {@link Allele} of the given {@code variant} to classify.
     * @param variantType The {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantType} of the given {@code variant}.
     * @param exon The {@link GencodeGtfExonFeature} in which the given {@code variant} occurs.
     * @param sequenceComparison The {@link org.broadinstitute.hellbender.tools.funcotator.SequenceComparison} for the given {@code variant}.
     * @return A {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} based on the given {@code allele}, {@code variant}, {@code exon}, and {@code sequenceComparison}.
     */
    @VisibleForTesting
    static GencodeFuncotation.VariantClassification createVariantClassification(final VariantContext variant,
                                                                                final Allele altAllele,
                                                                                final GencodeFuncotation.VariantType variantType,
                                                                                final GencodeGtfExonFeature exon,
                                                                                final SequenceComparison sequenceComparison ){

        Utils.nonNull(variant);
        Utils.nonNull(altAllele);
        Utils.nonNull(variantType);
        Utils.nonNull(exon);
        Utils.nonNull(sequenceComparison);

        // Get our start position:
        final int startPos = sequenceComparison.getAlleleStart();

        // Determine end position based on whichever allele is longer:
        final int endPos;
        if ( altAllele.length() >= variant.getReference().length()  ) {
            endPos = sequenceComparison.getAlleleStart() + altAllele.length() - 1;
        }
        else {
            endPos = sequenceComparison.getAlleleStart() + variant.getReference().length() - 1;
        }

        // Calculate the number of inserted bases so we can account for them in the splice site calculations:
        final int numInsertedBases = (altAllele.length() > variant.getReference().length()) ? altAllele.length() - variant.getReference().length() : 0;

        GencodeFuncotation.VariantClassification varClass = null;

        boolean hasBeenClassified = false;

        // Check for non-stop first:
        if ( (exon.getStopCodon() != null) && (exon.getStopCodon().overlaps(variant)) ) {

            boolean foundStop = false;

            for (int i = 0; i < sequenceComparison.getAlignedCodingSequenceAlternateAllele().length(); i+=3 ){
                final String codon = sequenceComparison.getAlignedCodingSequenceAlternateAllele().substring(i, i+3);
                if (FuncotatorUtils.getEukaryoticAminoAcidByCodon(codon) == AminoAcid.STOP_CODON) {
                    foundStop = true;
                    break;
                }
            }

            if ( !foundStop ) {
                varClass = GencodeFuncotation.VariantClassification.NONSTOP;
                hasBeenClassified = true;
            }
        }

        // Now check all other cases:
        if ( !hasBeenClassified ) {

            // Check for splice site variants.
            // Here we check to see if a splice site comes anywhere within `spliceSiteVariantWindowBases` of a variant.
            // We add and subtract 1 from the end points because the positons are 1-based & inclusive.
            if ( (((startPos - spliceSiteVariantWindowBases + 1) <= exon.getStart()) && (exon.getStart() <= (spliceSiteVariantWindowBases - numInsertedBases + endPos - 1))) ||
                 (((startPos - spliceSiteVariantWindowBases + 1) <= exon.getEnd()  ) && (exon.getEnd()   <= (spliceSiteVariantWindowBases - numInsertedBases + endPos - 1))) ) {
                varClass = GencodeFuncotation.VariantClassification.SPLICE_SITE;
            }
            else if ((exon.getStartCodon() != null) && (exon.getStartCodon().overlaps(variant))) {
                switch (variantType) {
                    case INS:
                        varClass = GencodeFuncotation.VariantClassification.START_CODON_INS;
                        break;
                    case DEL:
                        varClass = GencodeFuncotation.VariantClassification.START_CODON_DEL;
                        break;
                    default:
                        varClass = GencodeFuncotation.VariantClassification.START_CODON_SNP;
                        break;
                }
            }
            else {
                varClass = getVariantClassificationForCodingRegion(variant, altAllele, variantType, sequenceComparison);
            }
        }

        return varClass;
    }

    /**
     * Get the {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} for a given {@code variant}/{@code allele} in a coding region.
     * @param variant The {@link VariantContext} to classify.
     * @param altAllele The {@link Allele} of the given {@code variant} to classify.
     * @param variantType The {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantType} of the given {@code variant}.
     * @param sequenceComparison The {@link org.broadinstitute.hellbender.tools.funcotator.SequenceComparison} for the given {@code variant}.
     * @return A {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} based on the given {@code allele}, {@code variant}, {@code exon}, and {@code sequenceComparison}.
     */
    private static GencodeFuncotation.VariantClassification getVariantClassificationForCodingRegion(final VariantContext variant,
                                                                                             final Allele altAllele,
                                                                                             final GencodeFuncotation.VariantType variantType,
                                                                                             final SequenceComparison sequenceComparison) {
        final GencodeFuncotation.VariantClassification varClass;

        if (variantType == GencodeFuncotation.VariantType.INS) {
            if ( GATKVariantContextUtils.isFrameshift(variant.getReference(), altAllele)) {
                varClass = GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS;
            }
            else {
                varClass = GencodeFuncotation.VariantClassification.IN_FRAME_INS;
            }
        }
        else if (variantType == GencodeFuncotation.VariantType.DEL) {
            if (GATKVariantContextUtils.isFrameshift(variant.getReference(), altAllele)) {
                varClass = GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL;
            }
            else {
                varClass = GencodeFuncotation.VariantClassification.IN_FRAME_DEL;
            }
        }
        else {
            // This is a SNP/MNP
            // We just check to see what the protein change is to check for MISSENSE/NONSENSE/SILENT:
            varClass = getVarClassFromEqualLengthCodingRegions( sequenceComparison );
        }

        return varClass;
    }

    /**
     * Gets the {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} for a {@code variant} where the reference and alternate
     * alleles are the same length and the variant appears in a coding region.
     * This essentially compares the amino acid sequences of both alleles and returns a value based on the differences between them.
     * @return The {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} corresponding to the given variant / reference allele / alternate allele.
     */
    private static GencodeFuncotation.VariantClassification getVarClassFromEqualLengthCodingRegions(final SequenceComparison sequenceComparison) {

        GencodeFuncotation.VariantClassification varClass = GencodeFuncotation.VariantClassification.SILENT;

        boolean foundStop = false;

        for ( int i = 0; i < sequenceComparison.getAlternateAminoAcidSequence().length(); ++i ) {
            final char altAminoAcid = sequenceComparison.getAlternateAminoAcidSequence().charAt(i);

            if ( FuncotatorUtils.getAminoAcidByLetter(altAminoAcid) == AminoAcid.STOP_CODON ) {
                foundStop = true;
                break;
            }
            else if ( altAminoAcid != sequenceComparison.getReferenceAminoAcidSequence().charAt(i) ) {
                varClass = GencodeFuncotation.VariantClassification.MISSENSE;
            }
        }

        if ( foundStop ) {
            varClass = GencodeFuncotation.VariantClassification.NONSENSE;
        }

        return varClass;
    }

    /**
     * Create a {@link GencodeFuncotation} for a {@code variant} that occurs in an untranslated region in a given {@code transcript}.
     * @param variant The {@link VariantContext} for which to create a {@link GencodeFuncotation}.
     * @param altAllele The {@link Allele} in the given {@code variant} for which to create a {@link GencodeFuncotation}.
     * @param reference The {@link ReferenceContext} for the current data set.
     * @param gtfFeature The {@link GencodeGtfGeneFeature} in which the given {@code variant} occurs.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which the given {@code variant} occurs.
     * @param utr The {@link GencodeGtfUTRFeature} in which the given {@code variant} occurs.
     * @return A {@link GencodeFuncotation} containing information about the given {@code variant} given the corresponding {@code utr}.
     */
    private GencodeFuncotation createUtrFuncotation(final VariantContext variant,
                                                    final Allele altAllele,
                                                    final ReferenceContext reference,
                                                    final GencodeGtfGeneFeature gtfFeature,
                                                    final GencodeGtfTranscriptFeature transcript,
                                                    final GencodeGtfUTRFeature utr) {

        // Setup the "trivial" fields of the gencodeFuncotation:
        final GencodeFuncotationBuilder gencodeFuncotationBuilder = createGencodeFuncotationBuilderWithTrivialFieldsPopulated(variant, altAllele, gtfFeature, transcript);

        // Find which exon this UTR is in:
        for ( final GencodeGtfExonFeature exon : transcript.getExons() ) {
            if ( exon.contains( utr ) ) {
                gencodeFuncotationBuilder.setTranscriptExonNumber( exon.getExonNumber() );
            }
        }

        // Set GC Content:
        gencodeFuncotationBuilder.setGcContent( calculateGcContent( reference, gcContentWindowSizeBases ) );

        // Determine the strand for the variant:
        final Strand strand = Strand.decode( transcript.getGenomicStrand().toString() );
        FuncotatorUtils.assertValidStrand(strand);

        // Get the strand-corrected alleles from the inputs.
        // Also get the reference sequence for the variant region.
        // (spanning the entire length of both the reference and the variant, regardless of which is longer).
        final Allele strandCorrectedRefAllele;
        final Allele strandCorrectedAltAllele;

        if ( strand == Strand.POSITIVE ) {
            strandCorrectedRefAllele = variant.getReference();
            strandCorrectedAltAllele = altAllele;
        }
        else {
            strandCorrectedRefAllele = Allele.create(ReadUtils.getBasesReverseComplement( variant.getReference().getBases() ), true);
            strandCorrectedAltAllele = Allele.create(ReadUtils.getBasesReverseComplement( altAllele.getBases() ), false);
        }

        final String referenceBases = FuncotatorUtils.getBasesInWindowAroundReferenceAllele(strandCorrectedRefAllele, strandCorrectedAltAllele, strand, referenceWindow, reference);

        // Set our reference sequence in the Gencode Funcotation Builder:
        gencodeFuncotationBuilder.setReferenceContext( referenceBases );

        // Set whether it's the 5' or 3' UTR:
        if ( is5PrimeUtr(utr, transcript) ) {

            // We're 5' to the coding region.
            // This means we need to check for de novo starts.

            // Get our coding sequence for this region:
            final List<Locatable> activeRegions = Collections.singletonList(utr);

            final String referenceCodingSequence;
            if ( transcriptFastaReferenceDataSource != null ) {
                referenceCodingSequence = getCodingSequenceFromTranscriptFasta( transcript.getTranscriptId(), transcriptIdMap, transcriptFastaReferenceDataSource);
            }
            else {
                referenceCodingSequence = FuncotatorUtils.getCodingSequence(reference, activeRegions, strand);
            }

            final int codingStartPos = FuncotatorUtils.getStartPositionInTranscript(variant, activeRegions, strand);

            //Check for de novo starts:
            if ( FuncotatorUtils.getEukaryoticAminoAcidByCodon(referenceCodingSequence.substring(codingStartPos, codingStartPos+3) )
                    == AminoAcid.METHIONINE ) {

                // We know we have a new start.
                // Is it in frame or out of frame?
                if ( FuncotatorUtils.isInFrameWithEndOfRegion(codingStartPos, referenceCodingSequence.length()) ) {
                    gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME);
                }
                else {
                    gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME);
                }
            }
            else {
                gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR);
            }
        }
        else {
            gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.THREE_PRIME_UTR);
        }

        return gencodeFuncotationBuilder.build();
    }

    /**
     * Create a {@link GencodeFuncotation} for a {@code variant} that occurs in an intron in the given {@code transcript}.
     * @param variant The {@link VariantContext} for which to create a {@link GencodeFuncotation}.
     * @param altAllele The {@link Allele} in the given {@code variant} for which to create a {@link GencodeFuncotation}.
     * @param reference The {@link ReferenceContext} for the given {@code variant}.
     * @param gtfFeature The {@link GencodeGtfGeneFeature} in which the given {@code variant} occurs.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which the given {@code variant} occurs.
     * @param referenceContext The {@link ReferenceContext} in which the given variant appears.
     * @return A {@link GencodeFuncotation} containing information about the given {@code variant} given the corresponding {@code transcript}.
     */
    private static GencodeFuncotation createIntronFuncotation(final VariantContext variant,
                                                              final Allele altAllele,
                                                              final ReferenceContext reference,
                                                              final GencodeGtfGeneFeature gtfFeature,
                                                              final GencodeGtfTranscriptFeature transcript,
                                                              final ReferenceContext referenceContext) {

        // Setup the "trivial" fields of the gencodeFuncotation:
        final GencodeFuncotationBuilder gencodeFuncotationBuilder = createGencodeFuncotationBuilderWithTrivialFieldsPopulated(variant, altAllele, gtfFeature, transcript);

        // Determine the strand for the variant:
        final Strand strand = Strand.decode( transcript.getGenomicStrand().toString() );
        FuncotatorUtils.assertValidStrand(strand);

        // Get the strand-corrected alleles from the inputs.
        // Also get the reference sequence for the variant region.
        // (spanning the entire length of both the reference and the variant, regardless of which is longer).
        final Allele strandCorrectedRefAllele;
        final Allele strandCorrectedAltAllele;

        if ( strand == Strand.POSITIVE ) {
            strandCorrectedRefAllele = variant.getReference();
            strandCorrectedAltAllele = altAllele;
        }
        else {
            strandCorrectedRefAllele = Allele.create(ReadUtils.getBasesReverseComplement( variant.getReference().getBases() ), true);
            strandCorrectedAltAllele = Allele.create(ReadUtils.getBasesReverseComplement( altAllele.getBases() ), false);
        }

        final String referenceBases = FuncotatorUtils.getBasesInWindowAroundReferenceAllele(strandCorrectedRefAllele, strandCorrectedAltAllele, strand, referenceWindow, referenceContext);

        // Set our reference sequence in the Gencode Funcotation Builder:
        gencodeFuncotationBuilder.setReferenceContext( referenceBases );

        // Set as default INTRON variant classification:
        gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.INTRON);

        // Set GC Content:
        gencodeFuncotationBuilder.setGcContent( calculateGcContent( reference, gcContentWindowSizeBases ) );

        // Need to check if we're within the window for splice site variants:
        final GencodeGtfExonFeature spliceSiteExon = getExonWithinSpliceSiteWindow(variant, transcript, spliceSiteVariantWindowBases);
        if ( spliceSiteExon != null ) {
            // Set the variant classification:
            gencodeFuncotationBuilder.setVariantClassification(GencodeFuncotation.VariantClassification.SPLICE_SITE)
                                     .setSecondaryVariantClassification(GencodeFuncotation.VariantClassification.INTRON);

            // In deletions we have added a base to the front because of VCF requirements, thus we add an
            // offset of 1 to account for that:
            // (TODO: come to think of it this is really bad, because we're tying our parsing / computations to a data format).
            int offsetIndelAdjustment = 0;
            if ( GATKVariantContextUtils.isDeletion(variant.getReference(), altAllele) ) {
                offsetIndelAdjustment = 1;
            }

            gencodeFuncotationBuilder.setCodonChange(
                    FuncotatorUtils.createSpliceSiteCodonChange(variant.getStart(), spliceSiteExon.getExonNumber(), spliceSiteExon.getStart(), spliceSiteExon.getEnd(), strand, offsetIndelAdjustment)
            );
        }

        return gencodeFuncotationBuilder.build();
    }

    /**
     * Gets the {@link GencodeGtfExonFeature} that is within {@code spliceSiteVariantWindowBases} bases of the given {@code variant}.
     * @param variant The {@link VariantContext} to check for position within {@code spliceSiteVariantWindowBases} bases of each {@link GencodeGtfExonFeature} in {@code transcript}.
     * @param transcript The {@link GencodeGtfTranscriptFeature} containing the given {@code variant}.
     * @return The {@link GencodeGtfExonFeature} that is within {@code spliceSiteVariantWindowBases} bases of the given {@code variant}; {@code null} if no such {@link GencodeGtfExonFeature} exists in the given {@code transcript}.
     */
    private static GencodeGtfExonFeature getExonWithinSpliceSiteWindow( final VariantContext variant,
                                                                        final GencodeGtfTranscriptFeature transcript,
                                                                        final int spliceSiteVariantWindowBases ) {
        GencodeGtfExonFeature spliceSiteExon = null;

        for ( final GencodeGtfExonFeature exon : transcript.getExons() ) {
            if ((Math.abs(exon.getStart() - variant.getStart()) <= spliceSiteVariantWindowBases) ||
                    (Math.abs(exon.getEnd() - variant.getStart()) <= spliceSiteVariantWindowBases)) {
                spliceSiteExon = exon;
                break;
            }
        }

        return spliceSiteExon;
    }

    /**
     * Get the subfeature contained in {@code transcript} that contains the given {@code variant}.
     * The returned subfeature will be of type {@link GencodeGtfFeature} with concrete type based on the type of region
     * in which the variant is found:
     *      Found in coding region -> {@link GencodeGtfExonFeature}
     *      Found in UTR ->{@link GencodeGtfUTRFeature}
     *      Found in intron ->{@link GencodeGtfTranscriptFeature}
     *      Not Found in transcript ->{@code null}
     * @param variant A {@link VariantContext} of which to determine the containing subfeature.
     * @param transcript A {@link GencodeGtfTranscriptFeature} in which to find the subfeature containing the given {@code variant}.
     * @return The {@link GencodeGtfFeature} corresponding to the subfeature of {@code transcript} in which the given {@code variant} was found.
     */
    private static GencodeGtfFeature getContainingGtfSubfeature(final VariantContext variant, final GencodeGtfTranscriptFeature transcript) {

        boolean determinedRegionAlready = false;
        GencodeGtfFeature subFeature = null;

        if ( transcript.contains(variant) ) {
            if ( transcript.getUtrs().size() > 0 ) {
                for ( final GencodeGtfUTRFeature utr : transcript.getUtrs() ) {
                    if ( utr.overlaps(variant) ) {
                        subFeature = utr;
                        determinedRegionAlready = true;
                    }
                }
            }

            // Even though we may have an overlapping UTR already, we may be able to find a spot in the transcript
            // where this overlaps something more meaningful.
            // For example, see HG19 - chr19:8959608
            for (final GencodeGtfExonFeature exon : transcript.getExons()) {
                if ((exon.getCds() != null) && (exon.getCds().overlaps(variant))) {
                    subFeature = exon;
                    determinedRegionAlready = true;
                }
                else if ((exon.getStartCodon() != null) && (exon.getStartCodon().overlaps(variant))) {
                    subFeature = exon;
                    determinedRegionAlready = true;
                }
                else if ((exon.getStopCodon() != null) && (exon.getStopCodon().overlaps(variant))) {
                    subFeature = exon;
                    determinedRegionAlready = true;
                }
            }

            if ( !determinedRegionAlready ) {
                subFeature = transcript;
            }
        }

        return subFeature;
    }

    /**
     * Creates a {@link org.broadinstitute.hellbender.tools.funcotator.SequenceComparison} object with the fields populated.
     * @param variant The {@link VariantContext} for the current variant.
     * @param alternateAllele The current alternate {@link Allele} for the variant.
     * @param reference The {@link ReferenceContext} for the current sample set.
     * @param transcript The {@link GencodeGtfTranscriptFeature} for the current gene feature / alt allele.
     * @param exonPositionList A {@link List} of {@link htsjdk.samtools.util.Locatable} objects representing exon positions in the transcript.
     * @return A populated {@link org.broadinstitute.hellbender.tools.funcotator.SequenceComparison} object.
     */
    @VisibleForTesting
    static SequenceComparison createSequenceComparison(final VariantContext variant,
                                                       final Allele alternateAllele,
                                                       final ReferenceContext reference,
                                                       final GencodeGtfTranscriptFeature transcript,
                                                       final List<? extends htsjdk.samtools.util.Locatable> exonPositionList,
                                                       final Map<String, MappedTranscriptIdInfo> transcriptIdMap,
                                                       final ReferenceDataSource transcriptFastaReferenceDataSource) {

        final SequenceComparison sequenceComparison = new SequenceComparison();

        // Get the contig:
        sequenceComparison.setContig(variant.getContig());

        // Get the strand:
        final Strand strand = Strand.decode( transcript.getGenomicStrand().toString() );
        sequenceComparison.setStrand(strand);

        // Get the alleles from the inputs
        // Also get the reference sequence for the variant region
        // (spanning the entire length of both the reference and the variant, regardless of which is longer).
        final Allele refAllele;
        final Allele altAllele;

        // TODO: Make this a parameter:
        final int referenceWindow = 10;
        final String referenceBases;

        final SimpleInterval currentReferenceWindow = reference.getWindow();

        if ( strand == Strand.POSITIVE ) {
            refAllele = variant.getReference();
            altAllele = alternateAllele;

            // Calculate our window to include any extra bases but also have the right referenceWindow:
            final int endWindow = refAllele.length() >= altAllele.length() ? referenceWindow + refAllele.length() - 1: referenceWindow + altAllele.length() - 1;

            // Get the reference sequence:
            referenceBases = new String(reference.getBases(new SimpleInterval(currentReferenceWindow.getContig(), currentReferenceWindow.getStart() - referenceWindow, currentReferenceWindow.getEnd() + endWindow)));
        }
        else {
            refAllele = Allele.create(ReadUtils.getBasesReverseComplement( variant.getReference().getBases() ), true);
            altAllele = Allele.create(ReadUtils.getBasesReverseComplement( alternateAllele.getBases() ), false);

            // Calculate our window to include any extra bases but also have the right referenceWindow:
            final int endWindow = refAllele.length() >= altAllele.length() ? referenceWindow + refAllele.length() - 1: referenceWindow + altAllele.length() - 1;

            // Get the reference sequence:
            referenceBases = ReadUtils.getBasesReverseComplement(reference.getBases(new SimpleInterval(currentReferenceWindow.getContig(), currentReferenceWindow.getStart() - referenceWindow, currentReferenceWindow.getEnd() + endWindow)));
        }

        // Set our reference sequence in the SequenceComparison:
        sequenceComparison.setReferenceWindow( referenceWindow );
        sequenceComparison.setReferenceBases( referenceBases );

        // Set our GC content:
        sequenceComparison.setGcContent( calculateGcContent( reference, gcContentWindowSizeBases ) );

        // Get the coding sequence for the transcript:
        final String transcriptSequence;
        // NOTE: This can't be null because of the Funcotator input args.
        transcriptSequence = getCodingSequenceFromTranscriptFasta( transcript.getTranscriptId(), transcriptIdMap, transcriptFastaReferenceDataSource );

        // Get the transcript sequence as described by the given exonPositionList:
        sequenceComparison.setTranscriptCodingSequence(new ReferenceSequence(transcript.getTranscriptId(),transcript.getStart(),transcriptSequence.getBytes()));

        // Get the ref allele:
        sequenceComparison.setReferenceAllele(refAllele.getBaseString());

        // Get the allele genomic start position:
        sequenceComparison.setAlleleStart(variant.getStart());

        // Get the allele transcript start position (including UTRs):
        sequenceComparison.setTranscriptAlleleStart(
                FuncotatorUtils.getTranscriptAlleleStartPosition( variant.getStart(), transcript.getStart(), transcript.getEnd(), sequenceComparison.getStrand() )
        );

        // Get the coding region start position (in the above computed transcript coding region):
        sequenceComparison.setCodingSequenceAlleleStart(
                FuncotatorUtils.getStartPositionInTranscript( variant, exonPositionList, strand )
        );

        // Get the overlapping exon start / stop as an interval from the given variant:
        //TODO: See the overlap detector for this:
        sequenceComparison.setExonPosition(
                FuncotatorUtils.getOverlappingExonPositions( refAllele, altAllele, variant.getContig(), variant.getStart(), variant.getEnd(), strand, exonPositionList )
        );

        // Get the in-frame start position of the codon containing the given variant:
        sequenceComparison.setAlignedCodingSequenceAlleleStart(FuncotatorUtils.getAlignedPosition(sequenceComparison.getCodingSequenceAlleleStart()));

        // Get the in-frame stop position of the codon containing the given variant:
        sequenceComparison.setAlignedReferenceAlleleStop(
                FuncotatorUtils.getAlignedEndPosition(
                        // Subtract 1 because of the 1-based/inclusive nature of genetic coordinates:
                        sequenceComparison.getCodingSequenceAlleleStart() + refAllele.length() - 1
                )
        );

        // Get the in-frame/codon-aligned region containing the reference allele:
        sequenceComparison.setAlignedReferenceAllele(
                FuncotatorUtils.getAlignedRefAllele(
                        referenceBases,
                        referenceWindow,
                        refAllele,
                        sequenceComparison.getCodingSequenceAlleleStart(),
                        sequenceComparison.getAlignedCodingSequenceAlleleStart())
        );

        // Get the in-frame/codon-aligned CODING region containing the reference allele:
        // NOTE: We are calling this with Strand.POSITIVE because we have already reverse complemented the reference sequence.
        sequenceComparison.setAlignedCodingSequenceReferenceAllele(
                FuncotatorUtils.getAlignedCodingSequenceAllele(
                        sequenceComparison.getTranscriptCodingSequence().getBaseString(),
                        sequenceComparison.getAlignedCodingSequenceAlleleStart(),
                        sequenceComparison.getAlignedReferenceAlleleStop(),
                        refAllele,
                        sequenceComparison.getCodingSequenceAlleleStart(),
                        Strand.POSITIVE )
        );

        // Get the amino acid sequence of the reference allele:
        sequenceComparison.setReferenceAminoAcidSequence(
                FuncotatorUtils.createAminoAcidSequence( sequenceComparison.getAlignedCodingSequenceReferenceAllele() )
        );

        // Get the starting protein position of this variant.
        sequenceComparison.setProteinChangeStartPosition(
                FuncotatorUtils.getProteinChangePosition( sequenceComparison.getAlignedCodingSequenceAlleleStart() )
        );

        // Set our alternate allele:
        sequenceComparison.setAlternateAllele( altAllele.getBaseString() );

        // Set our stop position:
        sequenceComparison.setAlignedAlternateAlleleStop(
                FuncotatorUtils.getAlignedEndPosition(
                        sequenceComparison.getCodingSequenceAlleleStart() + altAllele.length() - 1
                )
        );

        // Get the aligned alternate allele:
        final int alignedRefAlleleStartPos = sequenceComparison.getCodingSequenceAlleleStart() - sequenceComparison.getAlignedCodingSequenceAlleleStart() + 1;
        sequenceComparison.setAlignedAlternateAllele(
                FuncotatorUtils.getAlternateSequence(
                        sequenceComparison.getAlignedReferenceAllele(),
                        alignedRefAlleleStartPos,
                        refAllele,
                        altAllele )
        );

        // Get the aligned coding sequence alternate allele:
        sequenceComparison.setAlignedCodingSequenceAlternateAllele(
                FuncotatorUtils.getAlternateSequence(
                        sequenceComparison.getAlignedCodingSequenceReferenceAllele(),
                        alignedRefAlleleStartPos,
                        refAllele,
                        altAllele )
        );

        // Set our alternate amino acid sequence:
        sequenceComparison.setAlternateAminoAcidSequence(
                FuncotatorUtils.createAminoAcidSequence(sequenceComparison.getAlignedCodingSequenceAlternateAllele())
        );

        // Set our protein end position:
        sequenceComparison.setProteinChangeEndPosition(
                FuncotatorUtils.getProteinChangeEndPosition(sequenceComparison.getProteinChangeStartPosition(), sequenceComparison.getAlignedCodingSequenceAlternateAllele().length())
        );

        return sequenceComparison;
    }

    /**
     * Calculates the fraction of Guanine and Cytosine bases in a window of a given size around a variant.
     * Note: Since Guanine and Cytosine are complementary bases, strandedness makes no difference.
     * @param referenceContext The {@link ReferenceContext} for a variant.  Assumed to already be centered on the variant of interest.  Must not be {@code null}.
     * @param windowSize The number of bases to the left and right of the given {@code variant} to calculate the GC Content.  Must be >=1.
     * @return The fraction of Guanine and Cytosine bases / total bases in a window of size {@code windowSize} around a variant.
     */
    public static double calculateGcContent( final ReferenceContext referenceContext,
                                             final int windowSize) {

        Utils.nonNull( referenceContext );
        ParamUtils.isPositive( windowSize, "Window size must be >= 1." );

        // Create a placeholder for the bases:
        final byte[] bases;

        // Get the bases:
        bases = referenceContext.getBases(windowSize, windowSize);

        // Get the gcCount:
        long gcCount = 0;
        for ( final byte base : bases ) {
            if ( BaseUtils.basesAreEqual(base, BaseUtils.Base.G.base) || BaseUtils.basesAreEqual(base, BaseUtils.Base.C.base) ) {
                ++gcCount;
            }
        }

        // Calculate the ratio itself:
        return ((double)gcCount) / ((double) bases.length);
    }

    /**
     * Creates a {@link GencodeFuncotationBuilder} with some of the fields populated.
     * @param variant The {@link VariantContext} for the current variant.
     * @param altAllele The alternate {@link Allele} we are currently annotating.
     * @param gtfFeature The current {@link GencodeGtfGeneFeature} read from the input feature file.
     * @param transcript The current {@link GencodeGtfTranscriptFeature} containing our {@code alternateAllele}.
     * @return A trivially populated {@link GencodeFuncotationBuilder} object.
     */
     private static GencodeFuncotationBuilder createGencodeFuncotationBuilderWithTrivialFieldsPopulated(final VariantContext variant,
                                                                                                        final Allele altAllele,
                                                                                                        final GencodeGtfGeneFeature gtfFeature,
                                                                                                        final GencodeGtfTranscriptFeature transcript) {

         final GencodeFuncotationBuilder gencodeFuncotationBuilder = new GencodeFuncotationBuilder();
         final Strand strand = Strand.decode(transcript.getGenomicStrand().toString());

         gencodeFuncotationBuilder.setRefAlleleAndStrand(variant.getReference(), strand)
                 .setHugoSymbol(gtfFeature.getGeneName())
                 .setNcbiBuild(gtfFeature.getUcscGenomeVersion())
                 .setChromosome(gtfFeature.getChromosomeName())
                 .setStart(variant.getStart());

         // The end position is inclusive, so we need to make sure we don't double-count the start position (so we subtract 1):
         gencodeFuncotationBuilder.setEnd(variant.getStart() + altAllele.length() - 1)
                 .setVariantType(getVariantType(variant.getReference(), altAllele))
                 .setTumorSeqAllele1(altAllele.getBaseString())
                 .setTumorSeqAllele2(altAllele.getBaseString())
                 .setGenomeChange(getGenomeChangeString(variant, altAllele, gtfFeature))
                 .setAnnotationTranscript(transcript.getTranscriptId())
                 .setTranscriptPos(
                         FuncotatorUtils.getTranscriptAlleleStartPosition(variant.getStart(), transcript.getStart(), transcript.getEnd(), strand)
                 )
                 .setOtherTranscripts(
                    gtfFeature.getTranscripts().stream().map(GencodeGtfTranscriptFeature::getTranscriptId).collect(Collectors.toList())
                 );

         // Check for the optional non-serialized values for sorting:
         // NOTE: This is kind of a kludge:
         gencodeFuncotationBuilder.setLocusLevel( Integer.valueOf(gtfFeature.getLocusLevel().toString()) );

        // Check for existence of Appris Rank and set it:
         gencodeFuncotationBuilder.setApprisRank( getApprisRank( gtfFeature ) );

         // Get the length of the transcript:
         // NOTE: We add 1 because of genomic cordinates:
         gencodeFuncotationBuilder.setTranscriptLength( transcript.getEnd() - transcript.getStart() + 1);return gencodeFuncotationBuilder;
    }

    /**
     * Determines if the given UTR is 3' or 5' of the given transcript.
     * Assumes the UTR is part of the given transcript.
     * @param utr The {@link GencodeGtfUTRFeature} to check for relative location in the given {@link GencodeGtfTranscriptFeature}.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which to check for the given {@code utr}.
     * @return {@code true} if the given {@code utr} is 5' for the given {@code transcript}; {@code false} otherwise.
     */
    private static boolean is5PrimeUtr(final GencodeGtfUTRFeature utr, final GencodeGtfTranscriptFeature transcript) {
        boolean isBefore = true;
        if ( transcript.getGenomicStrand() == GencodeGtfFeature.GenomicStrand.FORWARD ) {
            for ( final GencodeGtfExonFeature exon : transcript.getExons() ) {
                if ( exon.getStart() < utr.getStart()) {
                    isBefore = false;
                    break;
                }
            }
        }
        else {
            for ( final GencodeGtfExonFeature exon : transcript.getExons() ) {
                if ( exon.getStart() > utr.getStart()) {
                    isBefore = false;
                    break;
                }
            }
        }

        return isBefore;
    }

    /**
     * Creates a string representing the genome change given the variant, allele, and gene feature for this variant.
     * NOTE: The genome change will only reflect positions and bases within an exon.
     *       Positions beyond the ends of exons will be changed to the last position in the exon.
     *       Bases beyond the ends of exons will be truncated from the resulting string.
     * @param variant {@link VariantContext} of which to create the change.
     * @param altAllele {@link Allele} representing the alternate allele for this variant.
     * @param gtfFeature {@link GencodeGtfGeneFeature} corresponding to this variant.
     * @return A short {@link String} representation of the genomic change for the given variant, allele, and feature.
     */
    private static String getGenomeChangeString(final VariantContext variant,
                                                final Allele altAllele,
                                                final GencodeGtfGeneFeature gtfFeature) {

        // Check for insertion:
        if ( variant.getReference().length() < altAllele.length() ) {
            final String cleanAltAlleleString = FuncotatorUtils.getNonOverlappingAltAlleleBaseString( variant.getReference(), altAllele, false);

            return "g." + gtfFeature.getChromosomeName() +
                    ":" + variant.getStart() + "_" + (variant.getStart() + 1) + "ins" +
                    cleanAltAlleleString;
        }
        // Check for deletion:
        else if ( variant.getReference().length() > altAllele.length() ) {

            final String cleanAltAlleleString = FuncotatorUtils.getNonOverlappingAltAlleleBaseString(variant.getReference(), altAllele, true);

            final int startPos = variant.getStart() + 1;
            final int endPos = variant.getStart() + variant.getReference().length() - 1;

            if ( startPos == endPos ) {
                return "g." + gtfFeature.getChromosomeName() +
                        ":" + startPos + "del" + cleanAltAlleleString;
            }
            else {
                return "g." + gtfFeature.getChromosomeName() +
                        ":" + startPos + "_" + endPos +
                        "del" + cleanAltAlleleString;
            }
        }
        // Check for SNP:
        else if ( variant.getReference().length() == 1 ) {
            return "g." + gtfFeature.getChromosomeName() +
                    ":" + variant.getStart() +
                    variant.getReference().getBaseString() + ">" + altAllele.getBaseString();
        }
        else {
            return "g." + gtfFeature.getChromosomeName() +
                    ":" + variant.getStart() + "_" + ( variant.getStart() + variant.getReference().length() - 1) +
                    variant.getReference().getBaseString() + ">" + altAllele.getBaseString();
        }
    }

    /**
     * Sort the given list of funcotations such that the list becomes sorted in "best"->"worst" order by each funcotation's
     * transcript.
     * @param funcotationList The {@link List} of {@link GencodeFuncotation} to sort.
     */
    private void sortFuncotationsByTranscriptForOutput( final List<GencodeFuncotation> funcotationList ) {
        //TODO: Make this sort go from "worst" -> "best" so we can just pop the last element off and save some time.

        if ( transcriptSelectionMode == FuncotatorArgumentDefinitions.TranscriptSelectionMode.BEST_EFFECT ) {
            funcotationList.sort(new BestEffectGencodeFuncotationComparator(userRequestedTranscripts));
        }
        else {
            funcotationList.sort(new CannonicalGencodeFuncotationComparator(userRequestedTranscripts));
        }
    }

    /**
     * Creates a {@link List} of {@link GencodeFuncotation}s based on the given {@link VariantContext} with type
     * {@link GencodeFuncotation.VariantClassification#IGR}.
     * @param variant The {@link VariantContext} for which to create {@link Funcotation}s.
     * @param reference The {@link ReferenceContext} against which to compare the given {@link VariantContext}.
     * @return A list of IGR annotations for the given variant.
     */
    private static List<GencodeFuncotation> createIgrFuncotations(final VariantContext variant, final ReferenceContext reference) {
        // for each allele, create an annotation.

        final List<GencodeFuncotation> gencodeFuncotations = new ArrayList<>();

        for ( final Allele allele : variant.getAlternateAlleles() ) {
            gencodeFuncotations.add( createIgrFuncotation(variant.getReference(), allele, reference) );
        }

        return gencodeFuncotations;
    }

    /**
     * Condenses a given {@link GencodeFuncotation} into a string for the `other transcripts` annotation.
     * @param funcotation The {@link GencodeFuncotation} to condense.
     * @return A {@link String} representing the given {@link GencodeFuncotation}.
     */
    private static String condenseGencodeFuncotation( final GencodeFuncotation funcotation ) {
        Utils.nonNull( funcotation );

        final StringBuilder condensedFuncotationStringBuilder = new StringBuilder();

        if ( !funcotation.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR) ) {
            condensedFuncotationStringBuilder.append(funcotation.getHugoSymbol());
            condensedFuncotationStringBuilder.append("_");
            condensedFuncotationStringBuilder.append(funcotation.getAnnotationTranscript());
            condensedFuncotationStringBuilder.append("_");
            condensedFuncotationStringBuilder.append(funcotation.getVariantClassification());

            if ( !(funcotation.getVariantClassification().equals(GencodeFuncotation.VariantClassification.INTRON ) ||
                    ((funcotation.getSecondaryVariantClassification() != null) && funcotation.getSecondaryVariantClassification().equals(GencodeFuncotation.VariantClassification.INTRON))) ) {
                condensedFuncotationStringBuilder.append("_");
                condensedFuncotationStringBuilder.append(funcotation.getProteinChange());
            }
        }
        else {
            //TODO: This is known issue #3849:
            condensedFuncotationStringBuilder.append("IGR_ANNOTATON");
        }
        return condensedFuncotationStringBuilder.toString();
    }

    /**
     * Creates a {@link GencodeFuncotation}s based on the given {@link Allele} with type
     * {@link GencodeFuncotation.VariantClassification#IGR}.
     * Reports reference bases as if they are on the {@link Strand#POSITIVE} strand.
     * @param refAllele The reference {@link Allele} to use for this {@link GencodeFuncotation}.
     * @param altAllele The alternate {@link Allele} to use for this {@link GencodeFuncotation}.
     * @param reference The {@link ReferenceContext} in which the given {@link Allele}s appear.
     * @return An IGR funcotation for the given allele.
     */
    private static GencodeFuncotation createIgrFuncotation(final Allele refAllele,
                                                           final Allele altAllele,
                                                           final ReferenceContext reference){

        final GencodeFuncotationBuilder funcotationBuilder = new GencodeFuncotationBuilder();

        // Get GC Content:
        funcotationBuilder.setGcContent( calculateGcContent( reference, gcContentWindowSizeBases ) );

        funcotationBuilder.setVariantClassification( GencodeFuncotation.VariantClassification.IGR );
        funcotationBuilder.setRefAlleleAndStrand( refAllele, Strand.POSITIVE );
        funcotationBuilder.setTumorSeqAllele1( altAllele.getBaseString() );
        funcotationBuilder.setTumorSeqAllele2( altAllele.getBaseString() );

        final String referenceBases = FuncotatorUtils.getBasesInWindowAroundReferenceAllele(refAllele, altAllele, Strand.POSITIVE, referenceWindow, reference);

        // Set our reference context in the the FuncotatonBuilder:
        funcotationBuilder.setReferenceContext( referenceBases );

        return funcotationBuilder.build();
    }

    /**
     * Determines the variant type based on the given reference allele and alternate allele.
     * @param refAllele The reference {@link Allele} for this variant.
     * @param altAllele The alternate {@link Allele} for this variant.
     * @return A {@link GencodeFuncotation.VariantType} representing the variation type between the given reference and alternate {@link Allele}.
     */
    private static GencodeFuncotation.VariantType getVariantType( final Allele refAllele, final Allele altAllele ) {

        if ( altAllele.length() > refAllele.length() ) {
            return GencodeFuncotation.VariantType.INS;
        }
        else if (altAllele.length() < refAllele.length()) {
            return GencodeFuncotation.VariantType.DEL;
        }
        else {
            // We know they are the same length, now we just need to check one of them:
            switch (refAllele.length()) {
                case 1:  return GencodeFuncotation.VariantType.SNP;
                case 2:  return GencodeFuncotation.VariantType.DNP;
                case 3:  return GencodeFuncotation.VariantType.TNP;
                default: return GencodeFuncotation.VariantType.ONP;
            }
        }
    }

    /**
     * Get the Appris Rank from the given {@link GencodeGtfGeneFeature}.
     * @param gtfFeature The {@link GencodeGtfGeneFeature} from which to get the Appris Rank.
     * @return The highest Appris Rank found in the given {@code gtfFeature}; if no Appris Rank exists, {@code null}.
     */
    private static GencodeGtfFeature.FeatureTag getApprisRank( final GencodeGtfGeneFeature gtfFeature ) {

        // Get our appris tag(s) if it/they exist(s):
        final List<GencodeGtfFeature.FeatureTag> gtfApprisTags = gtfFeature.getOptionalFields().stream()
                .filter( f -> f.getName().equals("tag") )
                .filter( f -> f.getValue() instanceof GencodeGtfFeature.FeatureTag )
                .filter( f -> apprisRanks.contains( f.getValue() ) )
                .map( f -> (GencodeGtfFeature.FeatureTag)f.getValue() ).collect(Collectors.toList());

        if ( gtfApprisTags.isEmpty() ) {
            return null;
        }
        else if ( gtfApprisTags.size() == 1 ) {
            return gtfApprisTags.get(0);
        }
        else {
            // This case should never happen, but just in case we take the highest Appris Rank:
            gtfApprisTags.sort( Comparator.naturalOrder() );
            return gtfApprisTags.get(0);
        }
    }

    //==================================================================================================================
    // Helper Data Types:

    /**
     * Comparator class for Best Effect order.
     * Complex enough that a Lambda would be utter madness.
     */
    static class BestEffectGencodeFuncotationComparator implements Comparator<GencodeFuncotation> {

        final Set<String> userRequestedTranscripts;

        public BestEffectGencodeFuncotationComparator( final Set<String> userRequestedTranscripts ) {
            this.userRequestedTranscripts = userRequestedTranscripts;
        }

        @Override
        public int compare( final GencodeFuncotation a, final GencodeFuncotation b ) {
            // 1)
            // Choose the transcript that is on the custom list specified by the user:
            if ( isFuncotationInTranscriptList(a, userRequestedTranscripts) && (!isFuncotationInTranscriptList(b, userRequestedTranscripts)) ) {
                return -1;
            }
            else if ( (!isFuncotationInTranscriptList(a, userRequestedTranscripts)) && isFuncotationInTranscriptList(b, userRequestedTranscripts) ) {
                return 1;
            }

            // 1.5)
            // Check to see if one is an IGR.  IGR's have only a subset of the information in them, so it's easier to
            // order them if they're IGRs:
            else if ( (b.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR)) &&
                    (!a.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR)) ) {
                return -1;
            }
            else if ( (a.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR)) &&
                    (!b.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR)) ) {
                return 1;
            }

            // 2)
            // Check highest variant classification:
            else if ( a.getVariantClassification().getSeverity() < b.getVariantClassification().getSeverity() ) {
                return -1;
            }
            else if ( a.getVariantClassification().getSeverity() > b.getVariantClassification().getSeverity() ) {
                return 1;
            }

            // 3)
            // Check locus/curation levels:
            if ( (a.getLocusLevel() != null) && (b.getLocusLevel() == null) ) {
                return -1;
            }
            else if ( (a.getLocusLevel() == null ) && (b.getLocusLevel() != null) ) {
                return 1;
            }
            else if ( (a.getLocusLevel() != null) && (b.getLocusLevel() != null) && (!a.getLocusLevel().equals(b.getLocusLevel())) ) {
                return a.getLocusLevel().compareTo( b.getLocusLevel() );
            }

            // 4)
            // Check the appris annotation:
            else if ( (a.getApprisRank() != null) && (b.getApprisRank() == null) ) {
                return -1;
            }
            else if ( (a.getApprisRank() == null ) && (b.getApprisRank() != null) ) {
                return 1;
            }
            else if ( (a.getApprisRank() != null) && (b.getApprisRank() != null) && (!a.getApprisRank().equals(b.getApprisRank())) ) {
                return a.getApprisRank().compareTo( b.getApprisRank() );
            }

            // 5)
            // Check transcript sequence length:
            else if ( (a.getTranscriptLength() != null) && (b.getTranscriptLength() == null) ) {
                return -1;
            }
            else if ( (a.getTranscriptLength() == null ) && (b.getTranscriptLength() != null) ) {
                return 1;
            }
            else if ( (a.getTranscriptLength() != null) && (b.getTranscriptLength() != null) && (!a.getTranscriptLength().equals(b.getTranscriptLength())) ) {
                return b.getTranscriptLength().compareTo( a.getTranscriptLength() );
            }

            // 6)
            // Default to ABC order by transcript name:
            else if ( (a.getAnnotationTranscript() != null) && (b.getAnnotationTranscript() == null) ) {
                return -1;
            }
            else if ( (a.getAnnotationTranscript() == null ) && (b.getAnnotationTranscript() != null) ) {
                return 1;
            }
            // Need a default case in case all the comparison criteria are the same:
            else if ( (a.getAnnotationTranscript() == null ) && (b.getAnnotationTranscript() == null) ) {
                return -1;
            }
            else {
                return a.getAnnotationTranscript().compareTo(b.getAnnotationTranscript());
            }
        }
    }

    /**
     * Comparator class for Cannonical order.
     * Complex enough that a Lambda would be utter madness.
     */
    static class CannonicalGencodeFuncotationComparator implements Comparator<GencodeFuncotation> {

        final Set<String> userRequestedTranscripts;

        public CannonicalGencodeFuncotationComparator( final Set<String> userRequestedTranscripts ) {
            this.userRequestedTranscripts = userRequestedTranscripts;
        }

        @Override
        public int compare( final GencodeFuncotation a, final GencodeFuncotation b ) {

            // 1)
            // Choose the transcript that is on the custom list specified by the user:
            if ( isFuncotationInTranscriptList(a, userRequestedTranscripts) && (!isFuncotationInTranscriptList(b, userRequestedTranscripts)) ) {
                return -1;
            }
            else if ( (!isFuncotationInTranscriptList(a, userRequestedTranscripts)) && isFuncotationInTranscriptList(b, userRequestedTranscripts) ) {
                return 1;
            }

            // 2)
            // Check locus/curation levels:
            if ( (a.getLocusLevel() != null) && (b.getLocusLevel() == null) ) {
                return -1;
            }
            else if ( (a.getLocusLevel() == null ) && (b.getLocusLevel() != null) ) {
                return 1;
            }
            else if ( (a.getLocusLevel() != null) && (b.getLocusLevel() != null) && (!a.getLocusLevel().equals(b.getLocusLevel())) ) {
                return a.getLocusLevel().compareTo( b.getLocusLevel() );
            }

            // 2.5)
            // Check to see if one is an IGR.  IGR's have only a subset of the information in them, so it's easier to
            // order them if they're IGRs:
            if ( (b.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR)) &&
                    (!a.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR)) ) {
                return -1;
            }
            else if ( (a.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR)) &&
                    (!b.getVariantClassification().equals(GencodeFuncotation.VariantClassification.IGR)) ) {
                return 1;
            }

            // 3)
            // Check highest variant classification:
            else if ( a.getVariantClassification().getSeverity() < b.getVariantClassification().getSeverity() ) {
                return -1;
            }
            else if ( a.getVariantClassification().getSeverity() > b.getVariantClassification().getSeverity() ) {
                return 1;
            }

            // 4)
            // Check the appris annotation:
            else if ( (a.getApprisRank() != null) && (b.getApprisRank() == null) ) {
                return -1;
            }
            else if ( (a.getApprisRank() == null ) && (b.getApprisRank() != null) ) {
                return 1;
            }
            else if ( (a.getApprisRank() != null) && (b.getApprisRank() != null) && (!a.getApprisRank().equals(b.getApprisRank())) ) {
                return a.getApprisRank().compareTo( b.getApprisRank() );
            }

            // 5)
            // Check transcript sequence length:
            else if ( (a.getTranscriptLength() != null) && (b.getTranscriptLength() == null) ) {
                return -1;
            }
            else if ( (a.getTranscriptLength() == null ) && (b.getTranscriptLength() != null) ) {
                return 1;
            }
            else if ( (a.getTranscriptLength() != null) && (b.getTranscriptLength() != null) && (!a.getTranscriptLength().equals(b.getTranscriptLength())) ) {
                return b.getTranscriptLength().compareTo( a.getTranscriptLength() );
            }

            // 6)
            // Default to ABC order by transcript name:
            else if ( (a.getAnnotationTranscript() != null) && (b.getAnnotationTranscript() == null) ) {
                return -1;
            }
            else if ( (a.getAnnotationTranscript() == null ) && (b.getAnnotationTranscript() != null) ) {
                return 1;
            }
            // Need a default case in case all the comparison criteria are the same:
            else if ( (a.getAnnotationTranscript() == null ) && (b.getAnnotationTranscript() == null) ) {
                return -1;
            }
            else {
                return a.getAnnotationTranscript().compareTo(b.getAnnotationTranscript());
            }
        }
    }

    /**
     * A simple data object class to hold information about the transcripts in the
     * GENCODE transcript FASTA file.
     */
    @VisibleForTesting
    static class MappedTranscriptIdInfo {
        /**
         * The key in the GENCODE transcript FASTA file to use to get the coding sequence associated with this Transcript.
         */
        String mapKey;

        /**
         * The start position (1-based, inclusive) of the coding sequence in this transcript.
         */
        int codingSequenceStart;

        /**
         * The start position (1-based, inclusive) of the coding sequence in this transcript.
         */
        int codingSequenceEnd;

        /**
         * Whether or not the transcript has a 3' UTR.
         */
        boolean has3pUtr;

        /**
         * The start position (1-based, inclusive) of the 3' UTR in this transcript.
         */
        int threePrimeUtrStart;

        /**
         * The end position (1-based, inclusive) of the 3' UTR in this transcript.
         */
        int threePrimeUtrEnd;

        /**
         * Whether or not the transcript has a 3' UTR.
         */
        boolean has5pUtr;

        /**
         * The start position (1-based, inclusive) of the 5' UTR in this transcript.
         */
        int fivePrimeUtrStart;

        /**
         * The end position (1-based, inclusive) of the 5' UTR in this transcript.
         */
        int fivePrimeUtrEnd;
    }

}
