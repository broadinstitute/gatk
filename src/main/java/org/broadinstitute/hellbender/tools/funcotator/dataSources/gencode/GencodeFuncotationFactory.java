package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.GENCODE.*;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * A factory to create {@link GencodeFuncotation}s.
 * Created by jonn on 8/30/17.
 */
public class GencodeFuncotationFactory extends DataSourceFuncotationFactory {

    //==================================================================================================================

    /**
     * The window around splice sites to mark variants as {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification#SPLICE_SITE}.
     */
    final static private int spliceSiteVariantWindowBases = 2;

    //==================================================================================================================

    /**
     * ReferenceSequenceFile for the transcript reference file.
     */
    private final ReferenceDataSource transcriptFastaReferenceDataSource;

    /**
     * Map between transcript IDs and the IDs from the FASTA file to look up the transcript.
     * This is necessary because of the way the FASTA file contigs are named.
     */
    private final Map<String, MappedTranscriptIdInfo> transcriptIdMap;

    //==================================================================================================================

    public GencodeFuncotationFactory(final File gencodeTranscriptFastaFile) {
        transcriptFastaReferenceDataSource = ReferenceDataSource.of(gencodeTranscriptFastaFile);

        transcriptIdMap = createTranscriptIdMap(transcriptFastaReferenceDataSource);
    }

    //==================================================================================================================

    @Override
    public void close() {
        transcriptFastaReferenceDataSource.close();
    }

    @Override
    public List<String> getSupportedFuncotationFields() {
        return GencodeFuncotation.getSerializedFieldNames();
    }

    @Override
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {
        final List<Funcotation> funcotations = new ArrayList<>();

        // If we have features we need to annotate, go through them and create annotations:
        if ( featureList.size() > 0 ) {
            for ( final Allele altAllele : variant.getAlternateAlleles() ) {
                for ( final Feature feature : featureList ) {

                    // Get the kind of feature we want here:
                    if ( GencodeGtfGeneFeature.class.isAssignableFrom(feature.getClass()) ) {
                        funcotations.addAll(createFuncotations(variant, altAllele, (GencodeGtfGeneFeature) feature, referenceContext));
                    }

                    // NOTE: If we don't have any funcotations for this feature, it's OK.
                    //       However, this means that some other DataSourceFuncotationFactory must be producing a
                    //       funcotation for this variant.
                    //       For it is decreed that all variants must have a funcotation, even if that funcotation be
                    //       empty.
                    // TODO: Actually you may want to put another IGR creation here for now...  This may be a more difficult thing if we determine it in here.  There is no way to know if these are IGRs or simply not included in this particular data set.
                }
            }
        }
        else {
            // This is an IGR.
            funcotations.addAll( createIgrFuncotations(variant, referenceContext) );
        }

        return funcotations;
    }

    //==================================================================================================================

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
            return  (varClass == GencodeFuncotation.VariantClassification.MISSENSE)        ||
                    (varClass == GencodeFuncotation.VariantClassification.NONSENSE)        ||
                    (varClass == GencodeFuncotation.VariantClassification.NONSTOP)         ||
                    (varClass == GencodeFuncotation.VariantClassification.SILENT)          ||
                    (varClass == GencodeFuncotation.VariantClassification.IN_FRAME_DEL)    ||
                    (varClass == GencodeFuncotation.VariantClassification.IN_FRAME_INS)    ||
                    (varClass == GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS) ||
                    (varClass == GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL) ||
                    (varClass == GencodeFuncotation.VariantClassification.START_CODON_SNP) ||
                    (varClass == GencodeFuncotation.VariantClassification.START_CODON_INS) ||
                    (varClass == GencodeFuncotation.VariantClassification.START_CODON_DEL);
        }
    }

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

        final List<GencodeFuncotation> gencodeFuncotations = new ArrayList<>();

        // TODO: instead of getting the best one here, we do them all, then in the renderer we condense them into 1 annotation based on worst effect.
        // Get our "best" transcript:
        final int bestTranscriptIndex = getBestTranscriptIndex(gtfFeature, variant);
        if ( bestTranscriptIndex == -1 ) {
            throw new GATKException("Could not get a good transcript for the given feature: " + gtfFeature.toString());
        }
        final GencodeGtfTranscriptFeature transcript = gtfFeature.getTranscripts().remove(bestTranscriptIndex);

        final GencodeGtfFeature.GeneTranscriptType geneType = gtfFeature.getGeneType();

        // We only fully process protein-coding regions.
        // For other gene types, we do trivial processing and label them as either LINCRNA or RNA:
        // TODO: Add more types to Variant Classification to be more explicit and a converter function to go from gene type to variant classification.
        if ( geneType != GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING) {

            // Setup the "trivial" fields of the gencodeFuncotation:
            final GencodeFuncotation gencodeFuncotation = createGencodeFuncotationWithTrivialFieldsPopulated(variant, altAllele, gtfFeature, transcript);

            if ( geneType == GencodeGtfFeature.GeneTranscriptType.LINCRNA) {
                gencodeFuncotation.setVariantClassification(GencodeFuncotation.VariantClassification.LINCRNA);
            }
            else {
                gencodeFuncotation.setVariantClassification(GencodeFuncotation.VariantClassification.RNA);
            }
            gencodeFuncotations.add(gencodeFuncotation);
        }
        else {
            // Find the sub-feature of transcript that contains our variant:
            final GencodeGtfFeature containingSubfeature = getContainingGtfSubfeature(variant, transcript);

            // Determine what kind of region we're in and handle it in it's own way:
            if (containingSubfeature == null) {

                // We have an IGR variant
                gencodeFuncotations.add( createIgrFuncotation(altAllele) );

            } else if (GencodeGtfExonFeature.class.isAssignableFrom(containingSubfeature.getClass())) {

                // We have a coding region variant
                gencodeFuncotations.add( createCodingRegionFuncotation(variant, altAllele, gtfFeature, reference, transcript, (GencodeGtfExonFeature) containingSubfeature) );

            } else if (GencodeGtfUTRFeature.class.isAssignableFrom(containingSubfeature.getClass())) {

                // We have a UTR variant
                gencodeFuncotations.add( createUtrFuncotation(variant, altAllele, reference, gtfFeature, transcript, (GencodeGtfUTRFeature) containingSubfeature) );

            } else if (GencodeGtfTranscriptFeature.class.isAssignableFrom(containingSubfeature.getClass())) {

                // We have an intron variant
                gencodeFuncotations.add( createIntronFuncotation(variant, altAllele, gtfFeature, transcript) );

            } else {

                // Uh-oh!  Problemz.
                throw new GATKException("Unable to determine type of variant-containing subfeature: " + containingSubfeature.getClass().getName());
            }
        }

        return gencodeFuncotations;
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
        final GencodeFuncotation gencodeFuncotation = createGencodeFuncotationWithTrivialFieldsPopulated(variant, altAllele, gtfFeature, transcript);

        // Set the exon number:
        gencodeFuncotation.setTranscriptExon( exon.getExonNumber() );

        // Get the list of exons by their locations so we can use them to determine our location in the transcript and get
        // the transcript code itself:
        final List<? extends Locatable> exonPositionList = getSortedExonAndStartStopPositions(transcript);

        // Set up our SequenceComparison object so we can calculate some useful fields more easily
        // These fields can all be set without knowing the alternate allele:
        final FuncotatorUtils.SequenceComparison sequenceComparison = createSequenceComparison(variant, altAllele, reference, transcript, exonPositionList, transcriptIdMap, transcriptFastaReferenceDataSource);

        // OK, now that we have our SequenceComparison object set up we can continue the annotation:

        // Set the Codon and Protein changes:
        gencodeFuncotation.setCodonChange(FuncotatorUtils.getCodonChangeString(sequenceComparison));
        gencodeFuncotation.setProteinChange(FuncotatorUtils.getProteinChangeString(sequenceComparison));

        // Set the Variant Classification:
        final GencodeFuncotation.VariantClassification varClass = getVariantClassification( variant, altAllele, gencodeFuncotation.getVariantType(), exon, sequenceComparison );
        gencodeFuncotation.setVariantClassification( varClass );
        if ( varClass == GencodeFuncotation.VariantClassification.SPLICE_SITE ) {
            gencodeFuncotation.setSecondaryVariantClassification(
                    getVariantClassificationForCodingRegion(variant, altAllele, gencodeFuncotation.getVariantType(), sequenceComparison )
            );
        }

        // Set the Coding DNA change:
        final boolean isInCodingRegion = isVariantInCodingRegion( varClass, gencodeFuncotation.getSecondaryVariantClassification() );

        // Only set cDNA change if we have something in the actual coding region:
        if ( isInCodingRegion ) {
            gencodeFuncotation.setcDnaChange(FuncotatorUtils.getCodingSequenceChangeString(sequenceComparison));
        }

        return gencodeFuncotation;
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
     * @param sequenceComparison The {@link org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils.SequenceComparison} for the given {@code variant}.
     * @return A {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} based on the given {@code allele}, {@code variant}, {@code exon}, and {@code sequenceComparison}.
     */
    @VisibleForTesting
    static GencodeFuncotation.VariantClassification getVariantClassification(final VariantContext variant,
                                                                      final Allele altAllele,
                                                                      final GencodeFuncotation.VariantType variantType,
                                                                      final GencodeGtfExonFeature exon,
                                                                      final FuncotatorUtils.SequenceComparison sequenceComparison ){

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
     * @param sequenceComparison The {@link org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils.SequenceComparison} for the given {@code variant}.
     * @return A {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} based on the given {@code allele}, {@code variant}, {@code exon}, and {@code sequenceComparison}.
     */
    private static GencodeFuncotation.VariantClassification getVariantClassificationForCodingRegion(final VariantContext variant,
                                                                                             final Allele altAllele,
                                                                                             final GencodeFuncotation.VariantType variantType,
                                                                                             final FuncotatorUtils.SequenceComparison sequenceComparison) {
        final GencodeFuncotation.VariantClassification varClass;

        if (variantType == GencodeFuncotation.VariantType.INS) {
            if (FuncotatorUtils.isFrameshift(variant.getReference(), altAllele)) {
                varClass = GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS;
            }
            else {
                varClass = GencodeFuncotation.VariantClassification.IN_FRAME_INS;
            }
        }
        else if (variantType == GencodeFuncotation.VariantType.DEL) {
            if (FuncotatorUtils.isFrameshift(variant.getReference(), altAllele)) {
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
    private static GencodeFuncotation.VariantClassification getVarClassFromEqualLengthCodingRegions(final FuncotatorUtils.SequenceComparison sequenceComparison) {

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
        final GencodeFuncotation gencodeFuncotation = createGencodeFuncotationWithTrivialFieldsPopulated(variant, altAllele, gtfFeature, transcript);

        // Find which exon this UTR is in:
        for ( final GencodeGtfExonFeature exon : transcript.getExons() ) {
            if ( exon.contains( utr ) ) {
                gencodeFuncotation.setTranscriptExon( exon.getExonNumber() );
            }
        }

        // Set whether it's the 5' or 3' UTR:
        if ( is5PrimeUtr(utr, transcript) ) {

            // We're 5' to the coding region.
            // This means we need to check for de novo starts.

            // Get our coding sequence for this region:
            final List<Locatable> activeRegions = Collections.singletonList(utr);
            final Strand strand = Strand.toStrand( transcript.getGenomicStrand().toString() );

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
                    gencodeFuncotation.setVariantClassification(GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME);
                }
                else {
                    gencodeFuncotation.setVariantClassification(GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME);
                }
            }
            else {
                gencodeFuncotation.setVariantClassification(GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR);
            }
        }
        else {
            gencodeFuncotation.setVariantClassification(GencodeFuncotation.VariantClassification.THREE_PRIME_UTR);
        }

        return gencodeFuncotation;
    }

    /**
     * Create a {@link GencodeFuncotation} for a {@code variant} that occurs in an intron in the given {@code transcript}.
     * @param variant The {@link VariantContext} for which to create a {@link GencodeFuncotation}.
     * @param altAllele The {@link Allele} in the given {@code variant} for which to create a {@link GencodeFuncotation}.
     * @param gtfFeature The {@link GencodeGtfGeneFeature} in which the given {@code variant} occurs.
     * @param transcript The {@link GencodeGtfTranscriptFeature} in which the given {@code variant} occurs.
     * @return A {@link GencodeFuncotation} containing information about the given {@code variant} given the corresponding {@code transcript}.
     */
    private static GencodeFuncotation createIntronFuncotation(final VariantContext variant,
                                                       final Allele altAllele,
                                                       final GencodeGtfGeneFeature gtfFeature,
                                                       final GencodeGtfTranscriptFeature transcript) {

        // Setup the "trivial" fields of the gencodeFuncotation:
        final GencodeFuncotation gencodeFuncotation = createGencodeFuncotationWithTrivialFieldsPopulated(variant, altAllele, gtfFeature, transcript);

        // Set as default INTRON variant classification:
        gencodeFuncotation.setVariantClassification(GencodeFuncotation.VariantClassification.INTRON);

        // Determine the strand for the variant:
        final Strand strand = Strand.toStrand( transcript.getGenomicStrand().toString() );
        if ( strand == Strand.NONE ) {
            throw new GATKException("Unable to handle NONE strand.");
        }

        // Need to check if we're within the window for splice site variants:
        for ( final GencodeGtfExonFeature exon : transcript.getExons() ) {
            if (( Math.abs( exon.getStart() - variant.getStart() ) <= spliceSiteVariantWindowBases ) ||
                ( Math.abs( exon.getEnd() - variant.getStart() ) <= spliceSiteVariantWindowBases )) {

                // Set the variant classification:
                gencodeFuncotation.setVariantClassification(GencodeFuncotation.VariantClassification.SPLICE_SITE);
                gencodeFuncotation.setSecondaryVariantClassification(GencodeFuncotation.VariantClassification.INTRON);

                // In deletions we have added a base to the front because of VCF requirements, thus we add an
                // offset of 1 to account for that:
                // (TODO: come to think of it this is really bad, because we're tying our parsing / computations to a data format).
                int offsetIndelAdjustment = 0;
                if ( FuncotatorUtils.isDeletion(variant.getReference(), altAllele) ) {
                    offsetIndelAdjustment = 1;
                }

                gencodeFuncotation.setCodonChange(
                        FuncotatorUtils.createSpliceSiteCodonChange(variant.getStart(), exon.getExonNumber(), exon.getStart(), exon.getEnd(), strand, offsetIndelAdjustment)
                );
            }
        }

        return gencodeFuncotation;
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
     * Creates a {@link org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils.SequenceComparison} object with the fields populated.
     * @param variant The {@link VariantContext} for the current variant.
     * @param alternateAllele The current alternate {@link Allele} for the variant.
     * @param reference The {@link ReferenceContext} for the current sample set.
     * @param transcript The {@link GencodeGtfTranscriptFeature} for the current gene feature / alt allele.
     * @param exonPositionList A {@link List} of {@link htsjdk.samtools.util.Locatable} objects representing exon positions in the transcript.
     * @return A populated {@link org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils.SequenceComparison} object.
     */
    @VisibleForTesting
    static FuncotatorUtils.SequenceComparison createSequenceComparison(final VariantContext variant,
                                                                final Allele alternateAllele,
                                                                final ReferenceContext reference,
                                                                final GencodeGtfTranscriptFeature transcript,
                                                                final List<? extends htsjdk.samtools.util.Locatable> exonPositionList,
                                                                final Map<String, MappedTranscriptIdInfo> transcriptIdMap,
                                                                final ReferenceDataSource transcriptFastaReferenceDataSource) {

        final FuncotatorUtils.SequenceComparison sequenceComparison = new FuncotatorUtils.SequenceComparison();

        // Get the contig:
        sequenceComparison.setContig(variant.getContig());

        // Get the strand:
        final Strand strand = Strand.toStrand( transcript.getGenomicStrand().toString() );
        sequenceComparison.setStrand(strand);

        // Get the alleles from the inputs
        // Also get the reference sequence for the variant region
        // (spanning the entire length of both the reference and the variant, regardless of which is longer).
        final Allele refAllele;
        final Allele altAllele;

        // TODO: Make this a parameter:
        final int referenceWindow = 10;
        final String referenceBases;

        if ( strand == Strand.POSITIVE ) {
            refAllele = variant.getReference();
            altAllele = alternateAllele;

            // Calculate our window to include any extra bases but also have the right referenceWindow:
            final int endWindow = refAllele.length() >= altAllele.length() ? referenceWindow + refAllele.length() - 1: referenceWindow + altAllele.length() - 1;

            // Set our reference window:
            reference.setWindow(referenceWindow, endWindow);

            // Get the reference sequence:
            referenceBases = new String(reference.getBases());
        }
        else {
            refAllele = Allele.create(ReadUtils.getBasesReverseComplement( variant.getReference().getBases() ), true);
            altAllele = Allele.create(ReadUtils.getBasesReverseComplement( alternateAllele.getBases() ), false);

            // Calculate our window to include any extra bases but also have the right referenceWindow:
            final int endWindow = refAllele.length() >= altAllele.length() ? referenceWindow + refAllele.length() - 1: referenceWindow + altAllele.length() - 1;

            // Set our reference window:
            reference.setWindow(referenceWindow, endWindow);

            // Get the reference sequence:
            referenceBases = ReadUtils.getBasesReverseComplement(reference.getBases());
        }

        // Set our reference sequence in the SequenceComparison:
        sequenceComparison.setReferenceBases( referenceBases );
        sequenceComparison.setReferenceWindow( referenceWindow );

        // Get the coding sequence for the transcript:
        final String referenceCodingSequence;
        if ( transcriptFastaReferenceDataSource != null ) {
            referenceCodingSequence = getCodingSequenceFromTranscriptFasta( transcript.getTranscriptId(), transcriptIdMap, transcriptFastaReferenceDataSource );
        }
        else {
            referenceCodingSequence = FuncotatorUtils.getCodingSequence( reference, exonPositionList, strand );
        }

        // Get the reference sequence in the coding region as described by the given exonPositionList:
        sequenceComparison.setReferenceCodingSequence(new ReferenceSequence(transcript.getTranscriptId(),transcript.getStart(),referenceCodingSequence.getBytes()));

        // Get the ref allele:
        sequenceComparison.setReferenceAllele(refAllele.getBaseString());

        // Get the allele genomic start position:
        sequenceComparison.setAlleleStart(variant.getStart());

        // Get the allele transcript start position:
        sequenceComparison.setTranscriptAlleleStart(
                FuncotatorUtils.getTranscriptAlleleStartPosition( variant.getStart(), transcript.getStart(), transcript.getEnd(), sequenceComparison.getStrand() )
        );

        // Get the coding region start position (in the above computed reference coding region):
        sequenceComparison.setCodingSequenceAlleleStart(
                FuncotatorUtils.getStartPositionInTranscript( variant, exonPositionList, strand )
        );

        // Get the overlapping exon start / stop as an interval from the given variant:
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
                        sequenceComparison.getReferenceCodingSequence().getBaseString(),
                        sequenceComparison.getAlignedCodingSequenceAlleleStart(),
                        sequenceComparison.getAlignedReferenceAlleleStop(),
                        refAllele,
                        sequenceComparison.getCodingSequenceAlleleStart(),
                        Strand.POSITIVE )
        );

        // TODO: Check from here down for where you should use the coding sequence VS the raw reference alleles:

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
     * Creates a Gencode Funcotation with all trivial fields populated.
     * @param variant The {@link VariantContext} for the current variant.
     * @param altAllele The alternate {@link Allele} we are currently annotating.
     * @param gtfFeature The current {@link GencodeGtfGeneFeature} read from the input feature file.
     * @param transcript The current {@link GencodeGtfTranscriptFeature} containing our {@code alternateAllele}.
     * @return A trivially populated {@link GencodeFuncotation} object.
     */
     private static GencodeFuncotation createGencodeFuncotationWithTrivialFieldsPopulated(final VariantContext variant,
                                                                                  final Allele altAllele,
                                                                                  final GencodeGtfGeneFeature gtfFeature,
                                                                                  final GencodeGtfTranscriptFeature transcript) {
        final GencodeFuncotation gencodeFuncotation = new GencodeFuncotation();

        final Strand strand = Strand.toStrand( transcript.getGenomicStrand().toString() );

         if ( strand == Strand.POSITIVE ) {
             gencodeFuncotation.setRefAllele(variant.getReference().getBaseString());
             gencodeFuncotation.setTranscriptStrand("+");
         }
         else {
             gencodeFuncotation.setRefAllele(
                     ReadUtils.getBasesReverseComplement( variant.getReference().getBases() )
             );
             gencodeFuncotation.setTranscriptStrand("-");
         }

        gencodeFuncotation.setHugoSymbol( gtfFeature.getGeneName() );
        gencodeFuncotation.setNcbiBuild( gtfFeature.getUcscGenomeVersion() );
        gencodeFuncotation.setChromosome( gtfFeature.getChromosomeName() );

        gencodeFuncotation.setStart(variant.getStart());

         // The end position is inclusive, so we need to make sure we don't double-count the start position (so we subtract 1):
         gencodeFuncotation.setEnd(variant.getStart() + altAllele.length() - 1);

         gencodeFuncotation.setVariantType( getVariantType(variant.getReference(), altAllele) );

        gencodeFuncotation.setTumorSeqAllele1( altAllele.getBaseString() );
        gencodeFuncotation.setTumorSeqAllele2( altAllele.getBaseString() );

        gencodeFuncotation.setGenomeChange(getGenomeChangeString(variant, altAllele, gtfFeature));

        gencodeFuncotation.setAnnotationTranscript( transcript.getTranscriptId() );

        gencodeFuncotation.setTranscriptPos(
             FuncotatorUtils.getTranscriptAlleleStartPosition( variant.getStart(), transcript.getStart(), transcript.getEnd(), strand )
        );

        gencodeFuncotation.setOtherTranscripts(
                gtfFeature.getTranscripts().stream().map(GencodeGtfTranscriptFeature::getTranscriptId).collect(Collectors.toList())
        );

        return gencodeFuncotation;
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

//            // We handle the case where we have a coding region slightly differently than when we have only non-coding
//            // deleted bases:
//            if ( exon != null ) {
//                startPos = (variant.getStart() < exon.getStart()) ? exon.getStart() + 1 : variant.getStart() + 1;
//                endPos = (variant.getEnd() > exon.getEnd()) ? exon.getEnd() : (variant.getStart() + variant.getReference().length() - 1);
//
//                if (cleanAltAlleleString.length() > (endPos - startPos + 1)) {
//                    cleanAltAlleleString = cleanAltAlleleString.substring(0, endPos - startPos + 1);
//                }
//            }

            return "g." + gtfFeature.getChromosomeName() +
                    ":" + startPos + "_" + endPos +
                    "del" + cleanAltAlleleString;
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
     * Return the index of the "best" transcript in this gene.
     * @param geneFeature A {@link GencodeGtfGeneFeature} from which to get the index of the "best" transcript.
     * @param variant The {@link VariantContext} for which we want to get the best index.
     * @return The index of the "best" {@link GencodeGtfTranscriptFeature} in the given {@link GencodeGtfGeneFeature}.  Returns -1 if no transcript is present.
     */
    private static int getBestTranscriptIndex(final GencodeGtfGeneFeature geneFeature, final VariantContext variant) {
        if ( geneFeature.getTranscripts().size() == 0 ) {
            return -1;
        }

        for ( int i = 0 ; i < geneFeature.getTranscripts().size() ; ++i ) {
            if ( geneFeature.getTranscripts().get(i).getGenomicPosition().overlaps(variant) ) {
                return i;
            }
        }

        // Oops... we didn't find anything.
        return -1;
    }

    /**
     * Creates a {@link List} of {@link GencodeFuncotation}s based on the given {@link VariantContext} with type
     * {@link GencodeFuncotation.VariantClassification#IGR}.
     * @param variant The variant to annotate.
     * @param reference The reference against which to compare the given variant.
     * @return A list of IGR annotations for the given variant.
     */
    private static List<GencodeFuncotation> createIgrFuncotations(final VariantContext variant, final ReferenceContext reference) {
        // for each allele, create an annotation.

        // TODO: NEED TO FIX THIS LOGIC TO INCLUDE MORE INFO!

        final List<GencodeFuncotation> gencodeFuncotations = new ArrayList<>();

        for ( final Allele allele : variant.getAlternateAlleles() ) {
            gencodeFuncotations.add( createIgrFuncotation(allele) );
        }

        return gencodeFuncotations;
    }

    /**
     * Creates a {@link GencodeFuncotation}s based on the given {@link Allele} with type
     * {@link GencodeFuncotation.VariantClassification#IGR}.
     * @param altAllele The alternate allele to use for this funcotation.
     * @return An IGR funcotation for the given allele.
     */
    private static GencodeFuncotation createIgrFuncotation(final Allele altAllele){
        final GencodeFuncotation gencodeFuncotation = new GencodeFuncotation();

        gencodeFuncotation.setVariantClassification( GencodeFuncotation.VariantClassification.IGR );
        gencodeFuncotation.setTumorSeqAllele1( altAllele.getBaseString() );
        gencodeFuncotation.setTumorSeqAllele2( altAllele.getBaseString() );

        return gencodeFuncotation;
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

    //==================================================================================================================

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
