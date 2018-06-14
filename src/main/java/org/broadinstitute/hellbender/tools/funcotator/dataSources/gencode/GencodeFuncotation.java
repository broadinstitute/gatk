package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeature;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfGeneFeature;

import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A class to represent a Functional Annotation.  Each instance represents the annotations on a single transcript.
 * Created by jonn on 8/22/17.
 *
 * TODO: This will likely need to be renamed TranscriptFuncotation (or will have to implement an interface for TranscriptFuncotation) in order to handle non-gencode transcript sources.
 */
public class GencodeFuncotation implements Funcotation {

    //==================================================================================================================

    //------------------------------------------------------------
    // Fields for serialization:

    // NOTE: The order here matters.
    //       The way the column headers get generated for VCFs is dependent on the order here.
    //       In turn, the order of the fields corresponds to this order.

    private String                  hugoSymbol;                         // TRIVIAL (i.e. by the time we match to a transcript, we have this info regardless to where in the transcript the variant lies)
    private String                  ncbiBuild;                          // TRIVIAL
    private String                  chromosome;                         // TRIVIAL
    private int                     start;                              // TRIVIAL
    private int                     end;                                // TRIVIAL
    private VariantClassification   variantClassification;              //          CDS / UTR / INTRON / IGR
    private VariantClassification   secondaryVariantClassification;     //          CDS / INTRON
    private VariantType             variantType;                        // TRIVIAL
    private String                  refAllele;                          // TRIVIAL
    private String                  tumorSeqAllele2;                    // TRIVIAL

    private String                  genomeChange;                       // TRIVIAL
    private String                  annotationTranscript;               // TRIVIAL
    private String                  transcriptStrand;                   // TRIVIAL
    private Integer                 transcriptExon;                     //           CDS / UTRs
    private Integer                 transcriptPos;                      // TRIVIAL
    private String                  cDnaChange;                         //           CDS
    private String                  codonChange;                        //           CDS
    private String                  proteinChange;                      //           CDS
    private Double                  gcContent;

    private String                  referenceContext;                   // Already calculated.

    private List<String>            otherTranscripts;                   // TRIVIAL

    private String                  dataSourceName;

    //------------------------------------------------------------
    // Non-serialized fields:

    // These are included because they help determine the transcript selection
    private Integer                              locusLevel;
    private GencodeGtfGeneFeature.FeatureTag     apprisRank;
    private Integer                              transcriptLength;
    private String                               version;
    private GencodeGtfFeature.GeneTranscriptType geneTranscriptType;

    //------------------------------------------------------------
    // Fields for overriding serialized values:
    private String hugoSymbolSerializedOverride                     = null;
    private String ncbiBuildSerializedOverride                      = null;
    private String chromosomeSerializedOverride                     = null;
    private String startSerializedOverride                          = null;
    private String endSerializedOverride                            = null;
    private String variantClassificationSerializedOverride          = null;
    private String secondaryVariantClassificationSerializedOverride = null;
    private String variantTypeSerializedOverride                    = null;
    private String refAlleleSerializedOverride                      = null;
    private String tumorSeqAllele1SerializedOverride                = null;
    private String tumorSeqAllele2SerializedOverride                = null;

    private String genomeChangeSerializedOverride         = null;
    private String annotationTranscriptSerializedOverride = null;
    private String transcriptStrandSerializedOverride     = null;
    private String transcriptExonSerializedOverride       = null;
    private String transcriptPosSerializedOverride        = null;
    private String cDnaChangeSerializedOverride           = null;
    private String codonChangeSerializedOverride          = null;
    private String proteinChangeSerializedOverride        = null;
    private String gcContentSerializedOverride            = null;
    private String referenceContextSerializedOverride     = null;
    private String otherTranscriptsSerializedOverride     = null;

    private FuncotationMetadata metadata;
    //==================================================================================================================

    /**
     * Basic constructor for a {@link GencodeFuncotation}.
     *
     * Not public, since calling code should use the GencodeFuncotationBuilder.
     */
    @VisibleForTesting
    GencodeFuncotation() {}

    /**
     * Copy constructor for a {@link GencodeFuncotation}.
     */
    @VisibleForTesting
    GencodeFuncotation(final GencodeFuncotation that) {
        this.hugoSymbol = that.hugoSymbol;                         
        this.ncbiBuild = that.ncbiBuild;
        this.chromosome = that.chromosome;
        this.start = that.start;
        this.end = that.end;
        this.variantClassification = that.variantClassification;
        this.secondaryVariantClassification = that.secondaryVariantClassification;
        this.variantType = that.variantType;
        this.refAllele = that.refAllele;
        this.tumorSeqAllele2 = that.tumorSeqAllele2;
        this.genomeChange = that.genomeChange;
        this.annotationTranscript = that.annotationTranscript;
        this.transcriptStrand = that.transcriptStrand;
        this.transcriptExon = that.transcriptExon;
        this.transcriptPos = that.transcriptPos;
        this.cDnaChange = that.cDnaChange;
        this.codonChange = that.codonChange;
        this.proteinChange = that.proteinChange;
        this.gcContent = that.gcContent;
        this.referenceContext = that.referenceContext;
        this.otherTranscripts = that.otherTranscripts;
        this.dataSourceName = that.dataSourceName;
        this.locusLevel = that.locusLevel;
        this.apprisRank = that.apprisRank;
        this.transcriptLength = that.transcriptLength;
        this.version = that.version;
        this.geneTranscriptType = that.geneTranscriptType;
        this.hugoSymbolSerializedOverride = that.hugoSymbolSerializedOverride;
        this.ncbiBuildSerializedOverride = that.ncbiBuildSerializedOverride;
        this.chromosomeSerializedOverride = that.chromosomeSerializedOverride;
        this.startSerializedOverride = that.startSerializedOverride;
        this.endSerializedOverride = that.endSerializedOverride;
        this.variantClassificationSerializedOverride = that.variantClassificationSerializedOverride;
        this.secondaryVariantClassificationSerializedOverride = that.secondaryVariantClassificationSerializedOverride;
        this.variantTypeSerializedOverride = that.variantTypeSerializedOverride;
        this.refAlleleSerializedOverride = that.refAlleleSerializedOverride;
        this.tumorSeqAllele1SerializedOverride = that.tumorSeqAllele1SerializedOverride;
        this.tumorSeqAllele2SerializedOverride = that.tumorSeqAllele2SerializedOverride;
        this.genomeChangeSerializedOverride = that.genomeChangeSerializedOverride;
        this.annotationTranscriptSerializedOverride = that.annotationTranscriptSerializedOverride;
        this.transcriptStrandSerializedOverride = that.transcriptStrandSerializedOverride;
        this.transcriptExonSerializedOverride = that.transcriptExonSerializedOverride;
        this.transcriptPosSerializedOverride = that.transcriptPosSerializedOverride;
        this.cDnaChangeSerializedOverride = that.cDnaChangeSerializedOverride;
        this.codonChangeSerializedOverride = that.codonChangeSerializedOverride;
        this.proteinChangeSerializedOverride = that.proteinChangeSerializedOverride;
        this.gcContentSerializedOverride = that.gcContentSerializedOverride;
        this.referenceContextSerializedOverride = that.referenceContextSerializedOverride;
        this.otherTranscriptsSerializedOverride = that.otherTranscriptsSerializedOverride;
        this.metadata = that.metadata;
    }

    //==================================================================================================================

    /**
     * @return A new builder for a {@link GencodeFuncotation}.
     */
    public static GencodeFuncotationBuilder getBuilder() {
        return new GencodeFuncotationBuilder();
    }

    //==================================================================================================================

    @Override
    public Allele getAltAllele() {
        return Allele.create(tumorSeqAllele2.getBytes(), false);
    }

    @Override
    public String serializeToVcfString() {
        // Alias for the FIELD_DELIMITER so we can have nicer looking code:
        final String DELIMITER = VcfOutputRenderer.FIELD_DELIMITER;
        //TODO See issue https://github.com/broadinstitute/gatk/issues/4797

        // After the manual string, we check to see if we have an override first and if not we get the set field value:
        final List<String> funcotations = Arrays.asList((hugoSymbolSerializedOverride != null ? hugoSymbolSerializedOverride : (hugoSymbol != null ? hugoSymbol : "")),
                (ncbiBuildSerializedOverride != null ? ncbiBuildSerializedOverride : (ncbiBuild != null ? ncbiBuild : "")),
                (chromosomeSerializedOverride != null ? chromosomeSerializedOverride : (chromosome != null ? chromosome : "")),
                (startSerializedOverride != null ? startSerializedOverride : String.valueOf(start)),
                (endSerializedOverride != null ? endSerializedOverride : String.valueOf(end)),
                (variantClassificationSerializedOverride != null ? variantClassificationSerializedOverride : (variantClassification != null ? variantClassification.toString() : "")),
                (secondaryVariantClassificationSerializedOverride != null ? secondaryVariantClassificationSerializedOverride : (secondaryVariantClassification != null ? secondaryVariantClassification.toString() : "")),
                (variantTypeSerializedOverride != null ? variantTypeSerializedOverride : (variantType != null ? variantType.toString() : "")),
                (refAlleleSerializedOverride != null ? refAlleleSerializedOverride : (refAllele != null ? refAllele : "")),
                // NOTE: Ref allele gets serialized as the tumorSeqAllele1 as well, but we have to account for the override:
                (tumorSeqAllele1SerializedOverride != null ? tumorSeqAllele1SerializedOverride : (refAllele != null ? refAllele : "")),
                (tumorSeqAllele2SerializedOverride != null ? tumorSeqAllele2SerializedOverride : (tumorSeqAllele2 != null ? tumorSeqAllele2 : "")),
                (genomeChangeSerializedOverride != null ? genomeChangeSerializedOverride : (genomeChange != null ? genomeChange : "")),
                (annotationTranscriptSerializedOverride != null ? annotationTranscriptSerializedOverride : (annotationTranscript != null ? annotationTranscript : "")),
                (transcriptStrandSerializedOverride != null ? transcriptStrandSerializedOverride : (transcriptStrand != null ? transcriptStrand : "")),
                (transcriptExonSerializedOverride != null ? transcriptExonSerializedOverride : (transcriptExon != null ? transcriptExon.toString() : "")),
                (transcriptPosSerializedOverride != null ? transcriptPosSerializedOverride : (transcriptPos != null ? transcriptPos.toString() : "")),
                (cDnaChangeSerializedOverride != null ? cDnaChangeSerializedOverride : (cDnaChange != null ? cDnaChange : "")),
                (codonChangeSerializedOverride != null ? codonChangeSerializedOverride : (codonChange != null ? codonChange : "")),
                (proteinChangeSerializedOverride != null ? proteinChangeSerializedOverride : (proteinChange != null ? proteinChange : "")),
                (gcContentSerializedOverride != null ? gcContentSerializedOverride : (gcContent != null ? gcContent.toString() : "")),
                (referenceContextSerializedOverride != null ? referenceContextSerializedOverride : (referenceContext != null ? referenceContext : "")),
                (otherTranscriptsSerializedOverride != null ? otherTranscriptsSerializedOverride : (otherTranscripts != null ? otherTranscripts.stream().map(Object::toString).collect(Collectors.joining(VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER)) : ""))
            );

        return funcotations.stream().map(FuncotatorUtils::sanitizeFuncotationForVcf).collect(Collectors.joining(DELIMITER));
    }

    @Override
    public void setFieldSerializationOverrideValue( final String fieldName, final String overrideValue ) {

        // Cut off the "Gencode" and version number at the start of the string:
        final String shortFieldName = fieldName.replaceAll("^" + getDataSourceName()+ "_" + version + "_", "");

        switch (shortFieldName) {
            case "hugoSymbol":                     hugoSymbolSerializedOverride = overrideValue;                     break;
            case "ncbiBuild":                      ncbiBuildSerializedOverride = overrideValue;                      break;
            case "chromosome":                     chromosomeSerializedOverride = overrideValue;                     break;
            case "start":                          startSerializedOverride = overrideValue;                          break;
            case "end":                            endSerializedOverride = overrideValue;                            break;
            case "variantClassification":          variantClassificationSerializedOverride = overrideValue;          break;
            case "secondaryVariantClassification": secondaryVariantClassificationSerializedOverride = overrideValue; break;
            case "variantType":                    variantTypeSerializedOverride = overrideValue;                    break;
            case "refAllele":                      refAlleleSerializedOverride = overrideValue;                      break;
            case "tumorSeqAllele1":                tumorSeqAllele1SerializedOverride = overrideValue;                break;
            case "tumorSeqAllele2":                tumorSeqAllele2SerializedOverride = overrideValue;                break;
            case "genomeChange":                   genomeChangeSerializedOverride = overrideValue;                   break;
            case "annotationTranscript":           annotationTranscriptSerializedOverride = overrideValue;           break;
            case "transcriptStrand":               transcriptStrandSerializedOverride = overrideValue;               break;
            case "transcriptExon":                 transcriptExonSerializedOverride = overrideValue;                 break;
            case "transcriptPos":                  transcriptPosSerializedOverride = overrideValue;                  break;
            case "cDnaChange":                     cDnaChangeSerializedOverride = overrideValue;                     break;
            case "codonChange":                    codonChangeSerializedOverride = overrideValue;                    break;
            case "proteinChange":                  proteinChangeSerializedOverride = overrideValue;                  break;
            case "gcContent":                      gcContentSerializedOverride = overrideValue;                      break;
            case "referenceContext":               referenceContextSerializedOverride = overrideValue;               break;
            case "otherTranscripts":               otherTranscriptsSerializedOverride = overrideValue;               break;
            default: throw new UserException("Attempted to override invalid field in this GencodeFuncotation: " + fieldName + " (value was: " + overrideValue + ")");
        }
    }

    @Override
    public String getDataSourceName() {
        return dataSourceName;
    }

    @Override
    public LinkedHashSet<String> getFieldNames() {
        return new LinkedHashSet<>(
                Arrays.asList(
                        getDataSourceName() + "_" + version + "_hugoSymbol",
                        getDataSourceName() + "_" + version + "_ncbiBuild",
                        getDataSourceName() + "_" + version + "_chromosome",
                        getDataSourceName() + "_" + version + "_start",
                        getDataSourceName() + "_" + version + "_end",
                        getDataSourceName() + "_" + version + "_variantClassification",
                        getDataSourceName() + "_" + version + "_secondaryVariantClassification",
                        getDataSourceName() + "_" + version + "_variantType",
                        getDataSourceName() + "_" + version + "_refAllele",
                        getDataSourceName() + "_" + version + "_tumorSeqAllele1",
                        getDataSourceName() + "_" + version + "_tumorSeqAllele2",
                        getDataSourceName() + "_" + version + "_genomeChange",
                        getDataSourceName() + "_" + version + "_annotationTranscript",
                        getDataSourceName() + "_" + version + "_transcriptStrand",
                        getDataSourceName() + "_" + version + "_transcriptExon",
                        getDataSourceName() + "_" + version + "_transcriptPos",
                        getDataSourceName() + "_" + version + "_cDnaChange",
                        getDataSourceName() + "_" + version + "_codonChange",
                        getDataSourceName() + "_" + version + "_proteinChange",
                        getDataSourceName() + "_" + version + "_gcContent",
                        getDataSourceName() + "_" + version + "_referenceContext",
                        getDataSourceName() + "_" + version + "_otherTranscripts"
                )
        );
    }

    @Override
    public String getField(final String fieldName) {

        // Allow a user to specify the name of the field, or the fully-qualified name of the field
        // with GencodeFuncotationFactory.DATA_SOURCE_NAME + "_" + version + "_" at the start.
        final String altFieldName = getDataSourceName() + "_" + version + "_" + fieldName;
        final LinkedHashSet<String> fieldNames = getFieldNames();

        if ( fieldNames.contains(fieldName) || fieldNames.contains(altFieldName) ) {
            switch(fieldName.replace(getDataSourceName() + "_" + version + "_", "")) {
                case "hugoSymbol":
                    return (hugoSymbolSerializedOverride != null ? hugoSymbolSerializedOverride : (hugoSymbol != null ? hugoSymbol : ""));
                case "ncbiBuild":
                    return (ncbiBuildSerializedOverride != null ? ncbiBuildSerializedOverride : (ncbiBuild != null ? ncbiBuild : ""));
                case "chromosome":
                    return (chromosomeSerializedOverride != null ? chromosomeSerializedOverride : (chromosome != null ? chromosome : ""));
                case "start":
                    return (startSerializedOverride != null ? startSerializedOverride : String.valueOf(start));
                case "end":
                    return (endSerializedOverride != null ? endSerializedOverride : String.valueOf(end));
                case "variantClassification":
                    return (variantClassificationSerializedOverride != null ? variantClassificationSerializedOverride : (variantClassification != null ? String.valueOf(variantClassification) : ""));
                case "secondaryVariantClassification":
                    return (secondaryVariantClassificationSerializedOverride != null ? secondaryVariantClassificationSerializedOverride : (secondaryVariantClassification != null ? String.valueOf(secondaryVariantClassification) : ""));
                case "variantType":
                    return (variantTypeSerializedOverride != null ? variantTypeSerializedOverride : (variantType != null ? String.valueOf(variantType) : ""));
                case "refAllele":
                    return (refAlleleSerializedOverride != null ? refAlleleSerializedOverride : (refAllele != null ? refAllele : ""));
                case "tumorSeqAllele1":
                    return (tumorSeqAllele1SerializedOverride != null ? tumorSeqAllele1SerializedOverride : (refAllele != null ? refAllele : ""));
                case "tumorSeqAllele2":
                    return (tumorSeqAllele2SerializedOverride != null ? tumorSeqAllele2SerializedOverride : (tumorSeqAllele2 != null ? tumorSeqAllele2 : ""));
                case "genomeChange":
                    return (genomeChangeSerializedOverride != null ? genomeChangeSerializedOverride : (genomeChange != null ? genomeChange : ""));
                case "annotationTranscript":
                    return (annotationTranscriptSerializedOverride != null ? annotationTranscriptSerializedOverride : (annotationTranscript != null ? annotationTranscript : ""));
                case "transcriptStrand":
                    return (transcriptStrandSerializedOverride != null ? transcriptStrandSerializedOverride : (transcriptStrand != null ? transcriptStrand : ""));
                case "transcriptExon":
                    return (transcriptExonSerializedOverride != null ? transcriptExonSerializedOverride : (transcriptExon != null ? String.valueOf(transcriptExon) : ""));
                case "transcriptPos":
                    return (transcriptPosSerializedOverride != null ? transcriptPosSerializedOverride : (transcriptPos != null ? String.valueOf(transcriptPos) : ""));
                case "cDnaChange":
                    return (cDnaChangeSerializedOverride != null ? cDnaChangeSerializedOverride : (cDnaChange != null ? cDnaChange : ""));
                case "codonChange":
                    return (codonChangeSerializedOverride != null ? codonChangeSerializedOverride : (codonChange != null ? codonChange : ""));
                case "proteinChange":
                    return (proteinChangeSerializedOverride != null ? proteinChangeSerializedOverride : (proteinChange != null ? proteinChange : ""));
                case "gcContent":
                    return (gcContentSerializedOverride != null ? gcContentSerializedOverride : (gcContent != null ? String.valueOf(gcContent) : ""));
                case "referenceContext":
                    return (referenceContextSerializedOverride != null ? referenceContextSerializedOverride : (referenceContext != null ? referenceContext : ""));
                case "otherTranscripts":
                    return (otherTranscriptsSerializedOverride != null ? otherTranscriptsSerializedOverride : (otherTranscripts != null ? otherTranscripts.stream().map(Object::toString).collect(Collectors.joining(VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER)) : ""));
            }
        }

        throw new GATKException(this.getClass().getSimpleName() + ": Does not contain field: " + fieldName);
    }

    @Override
    public boolean hasField(final String fieldName) {
        final LinkedHashSet<String> fieldNames = getFieldNames();
        final String altFieldName = getDataSourceName() + "_" + version + "_" + fieldName;
        return ( fieldNames.contains(fieldName) || fieldNames.contains(altFieldName) );
    }

    @Override
    public FuncotationMetadata getMetadata() {
        return this.metadata;
    }

    //==================================================================================================================


    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        GencodeFuncotation that = (GencodeFuncotation) o;

        if (start != that.start) return false;
        if (end != that.end) return false;
        if (hugoSymbol != null ? !hugoSymbol.equals(that.hugoSymbol) : that.hugoSymbol != null) return false;
        if (ncbiBuild != null ? !ncbiBuild.equals(that.ncbiBuild) : that.ncbiBuild != null) return false;
        if (chromosome != null ? !chromosome.equals(that.chromosome) : that.chromosome != null) return false;
        if (variantClassification != that.variantClassification) return false;
        if (secondaryVariantClassification != that.secondaryVariantClassification) return false;
        if (variantType != that.variantType) return false;
        if (refAllele != null ? !refAllele.equals(that.refAllele) : that.refAllele != null) return false;
        if (tumorSeqAllele2 != null ? !tumorSeqAllele2.equals(that.tumorSeqAllele2) : that.tumorSeqAllele2 != null)
            return false;
        if (genomeChange != null ? !genomeChange.equals(that.genomeChange) : that.genomeChange != null) return false;
        if (annotationTranscript != null ? !annotationTranscript.equals(that.annotationTranscript) : that.annotationTranscript != null)
            return false;
        if (transcriptStrand != null ? !transcriptStrand.equals(that.transcriptStrand) : that.transcriptStrand != null)
            return false;
        if (transcriptExon != null ? !transcriptExon.equals(that.transcriptExon) : that.transcriptExon != null)
            return false;
        if (transcriptPos != null ? !transcriptPos.equals(that.transcriptPos) : that.transcriptPos != null)
            return false;
        if (cDnaChange != null ? !cDnaChange.equals(that.cDnaChange) : that.cDnaChange != null) return false;
        if (codonChange != null ? !codonChange.equals(that.codonChange) : that.codonChange != null) return false;
        if (proteinChange != null ? !proteinChange.equals(that.proteinChange) : that.proteinChange != null)
            return false;
        if (gcContent != null ? !gcContent.equals(that.gcContent) : that.gcContent != null) return false;
        if (referenceContext != null ? !referenceContext.equals(that.referenceContext) : that.referenceContext != null)
            return false;
        if (otherTranscripts != null ? !otherTranscripts.equals(that.otherTranscripts) : that.otherTranscripts != null)
            return false;
        if (dataSourceName != null ? !dataSourceName.equals(that.dataSourceName) : that.dataSourceName != null)
            return false;
        if (locusLevel != null ? !locusLevel.equals(that.locusLevel) : that.locusLevel != null) return false;
        if (apprisRank != that.apprisRank) return false;
        if (transcriptLength != null ? !transcriptLength.equals(that.transcriptLength) : that.transcriptLength != null)
            return false;
        if (version != null ? !version.equals(that.version) : that.version != null) return false;
        if (geneTranscriptType != that.geneTranscriptType) return false;
        if (hugoSymbolSerializedOverride != null ? !hugoSymbolSerializedOverride.equals(that.hugoSymbolSerializedOverride) : that.hugoSymbolSerializedOverride != null)
            return false;
        if (ncbiBuildSerializedOverride != null ? !ncbiBuildSerializedOverride.equals(that.ncbiBuildSerializedOverride) : that.ncbiBuildSerializedOverride != null)
            return false;
        if (chromosomeSerializedOverride != null ? !chromosomeSerializedOverride.equals(that.chromosomeSerializedOverride) : that.chromosomeSerializedOverride != null)
            return false;
        if (startSerializedOverride != null ? !startSerializedOverride.equals(that.startSerializedOverride) : that.startSerializedOverride != null)
            return false;
        if (endSerializedOverride != null ? !endSerializedOverride.equals(that.endSerializedOverride) : that.endSerializedOverride != null)
            return false;
        if (variantClassificationSerializedOverride != null ? !variantClassificationSerializedOverride.equals(that.variantClassificationSerializedOverride) : that.variantClassificationSerializedOverride != null)
            return false;
        if (secondaryVariantClassificationSerializedOverride != null ? !secondaryVariantClassificationSerializedOverride.equals(that.secondaryVariantClassificationSerializedOverride) : that.secondaryVariantClassificationSerializedOverride != null)
            return false;
        if (variantTypeSerializedOverride != null ? !variantTypeSerializedOverride.equals(that.variantTypeSerializedOverride) : that.variantTypeSerializedOverride != null)
            return false;
        if (refAlleleSerializedOverride != null ? !refAlleleSerializedOverride.equals(that.refAlleleSerializedOverride) : that.refAlleleSerializedOverride != null)
            return false;
        if (tumorSeqAllele1SerializedOverride != null ? !tumorSeqAllele1SerializedOverride.equals(that.tumorSeqAllele1SerializedOverride) : that.tumorSeqAllele1SerializedOverride != null)
            return false;
        if (tumorSeqAllele2SerializedOverride != null ? !tumorSeqAllele2SerializedOverride.equals(that.tumorSeqAllele2SerializedOverride) : that.tumorSeqAllele2SerializedOverride != null)
            return false;
        if (genomeChangeSerializedOverride != null ? !genomeChangeSerializedOverride.equals(that.genomeChangeSerializedOverride) : that.genomeChangeSerializedOverride != null)
            return false;
        if (annotationTranscriptSerializedOverride != null ? !annotationTranscriptSerializedOverride.equals(that.annotationTranscriptSerializedOverride) : that.annotationTranscriptSerializedOverride != null)
            return false;
        if (transcriptStrandSerializedOverride != null ? !transcriptStrandSerializedOverride.equals(that.transcriptStrandSerializedOverride) : that.transcriptStrandSerializedOverride != null)
            return false;
        if (transcriptExonSerializedOverride != null ? !transcriptExonSerializedOverride.equals(that.transcriptExonSerializedOverride) : that.transcriptExonSerializedOverride != null)
            return false;
        if (transcriptPosSerializedOverride != null ? !transcriptPosSerializedOverride.equals(that.transcriptPosSerializedOverride) : that.transcriptPosSerializedOverride != null)
            return false;
        if (cDnaChangeSerializedOverride != null ? !cDnaChangeSerializedOverride.equals(that.cDnaChangeSerializedOverride) : that.cDnaChangeSerializedOverride != null)
            return false;
        if (codonChangeSerializedOverride != null ? !codonChangeSerializedOverride.equals(that.codonChangeSerializedOverride) : that.codonChangeSerializedOverride != null)
            return false;
        if (proteinChangeSerializedOverride != null ? !proteinChangeSerializedOverride.equals(that.proteinChangeSerializedOverride) : that.proteinChangeSerializedOverride != null)
            return false;
        if (gcContentSerializedOverride != null ? !gcContentSerializedOverride.equals(that.gcContentSerializedOverride) : that.gcContentSerializedOverride != null)
            return false;
        if (referenceContextSerializedOverride != null ? !referenceContextSerializedOverride.equals(that.referenceContextSerializedOverride) : that.referenceContextSerializedOverride != null)
            return false;
        if (otherTranscriptsSerializedOverride != null ? !otherTranscriptsSerializedOverride.equals(that.otherTranscriptsSerializedOverride) : that.otherTranscriptsSerializedOverride != null)
            return false;
        return metadata != null ? metadata.equals(that.metadata) : that.metadata == null;
    }

    @Override
    public int hashCode() {
        int result = hugoSymbol != null ? hugoSymbol.hashCode() : 0;
        result = 31 * result + (ncbiBuild != null ? ncbiBuild.hashCode() : 0);
        result = 31 * result + (chromosome != null ? chromosome.hashCode() : 0);
        result = 31 * result + start;
        result = 31 * result + end;
        result = 31 * result + (variantClassification != null ? variantClassification.hashCode() : 0);
        result = 31 * result + (secondaryVariantClassification != null ? secondaryVariantClassification.hashCode() : 0);
        result = 31 * result + (variantType != null ? variantType.hashCode() : 0);
        result = 31 * result + (refAllele != null ? refAllele.hashCode() : 0);
        result = 31 * result + (tumorSeqAllele2 != null ? tumorSeqAllele2.hashCode() : 0);
        result = 31 * result + (genomeChange != null ? genomeChange.hashCode() : 0);
        result = 31 * result + (annotationTranscript != null ? annotationTranscript.hashCode() : 0);
        result = 31 * result + (transcriptStrand != null ? transcriptStrand.hashCode() : 0);
        result = 31 * result + (transcriptExon != null ? transcriptExon.hashCode() : 0);
        result = 31 * result + (transcriptPos != null ? transcriptPos.hashCode() : 0);
        result = 31 * result + (cDnaChange != null ? cDnaChange.hashCode() : 0);
        result = 31 * result + (codonChange != null ? codonChange.hashCode() : 0);
        result = 31 * result + (proteinChange != null ? proteinChange.hashCode() : 0);
        result = 31 * result + (gcContent != null ? gcContent.hashCode() : 0);
        result = 31 * result + (referenceContext != null ? referenceContext.hashCode() : 0);
        result = 31 * result + (otherTranscripts != null ? otherTranscripts.hashCode() : 0);
        result = 31 * result + (dataSourceName != null ? dataSourceName.hashCode() : 0);
        result = 31 * result + (locusLevel != null ? locusLevel.hashCode() : 0);
        result = 31 * result + (apprisRank != null ? apprisRank.hashCode() : 0);
        result = 31 * result + (transcriptLength != null ? transcriptLength.hashCode() : 0);
        result = 31 * result + (version != null ? version.hashCode() : 0);
        result = 31 * result + (geneTranscriptType != null ? geneTranscriptType.hashCode() : 0);
        result = 31 * result + (hugoSymbolSerializedOverride != null ? hugoSymbolSerializedOverride.hashCode() : 0);
        result = 31 * result + (ncbiBuildSerializedOverride != null ? ncbiBuildSerializedOverride.hashCode() : 0);
        result = 31 * result + (chromosomeSerializedOverride != null ? chromosomeSerializedOverride.hashCode() : 0);
        result = 31 * result + (startSerializedOverride != null ? startSerializedOverride.hashCode() : 0);
        result = 31 * result + (endSerializedOverride != null ? endSerializedOverride.hashCode() : 0);
        result = 31 * result + (variantClassificationSerializedOverride != null ? variantClassificationSerializedOverride.hashCode() : 0);
        result = 31 * result + (secondaryVariantClassificationSerializedOverride != null ? secondaryVariantClassificationSerializedOverride.hashCode() : 0);
        result = 31 * result + (variantTypeSerializedOverride != null ? variantTypeSerializedOverride.hashCode() : 0);
        result = 31 * result + (refAlleleSerializedOverride != null ? refAlleleSerializedOverride.hashCode() : 0);
        result = 31 * result + (tumorSeqAllele1SerializedOverride != null ? tumorSeqAllele1SerializedOverride.hashCode() : 0);
        result = 31 * result + (tumorSeqAllele2SerializedOverride != null ? tumorSeqAllele2SerializedOverride.hashCode() : 0);
        result = 31 * result + (genomeChangeSerializedOverride != null ? genomeChangeSerializedOverride.hashCode() : 0);
        result = 31 * result + (annotationTranscriptSerializedOverride != null ? annotationTranscriptSerializedOverride.hashCode() : 0);
        result = 31 * result + (transcriptStrandSerializedOverride != null ? transcriptStrandSerializedOverride.hashCode() : 0);
        result = 31 * result + (transcriptExonSerializedOverride != null ? transcriptExonSerializedOverride.hashCode() : 0);
        result = 31 * result + (transcriptPosSerializedOverride != null ? transcriptPosSerializedOverride.hashCode() : 0);
        result = 31 * result + (cDnaChangeSerializedOverride != null ? cDnaChangeSerializedOverride.hashCode() : 0);
        result = 31 * result + (codonChangeSerializedOverride != null ? codonChangeSerializedOverride.hashCode() : 0);
        result = 31 * result + (proteinChangeSerializedOverride != null ? proteinChangeSerializedOverride.hashCode() : 0);
        result = 31 * result + (gcContentSerializedOverride != null ? gcContentSerializedOverride.hashCode() : 0);
        result = 31 * result + (referenceContextSerializedOverride != null ? referenceContextSerializedOverride.hashCode() : 0);
        result = 31 * result + (otherTranscriptsSerializedOverride != null ? otherTranscriptsSerializedOverride.hashCode() : 0);
        result = 31 * result + (metadata != null ? metadata.hashCode() : 0);
        return result;
    }

    @Override
    public String toString() {
        return "GencodeFuncotation{" +
                "hugoSymbol='" + hugoSymbol + '\'' +
                ", ncbiBuild='" + ncbiBuild + '\'' +
                ", chromosome='" + chromosome + '\'' +
                ", start=" + start +
                ", end=" + end +
                ", variantClassification=" + variantClassification +
                ", secondaryVariantClassification=" + secondaryVariantClassification +
                ", variantType=" + variantType +
                ", refAllele='" + refAllele + '\'' +
                ", tumorSeqAllele2='" + tumorSeqAllele2 + '\'' +
                ", genomeChange='" + genomeChange + '\'' +
                ", annotationTranscript='" + annotationTranscript + '\'' +
                ", transcriptStrand='" + transcriptStrand + '\'' +
                ", transcriptExon=" + transcriptExon +
                ", transcriptPos=" + transcriptPos +
                ", cDnaChange='" + cDnaChange + '\'' +
                ", codonChange='" + codonChange + '\'' +
                ", proteinChange='" + proteinChange + '\'' +
                ", gcContent=" + gcContent +
                ", referenceContext='" + referenceContext + '\'' +
                ", otherTranscripts=" + otherTranscripts +
                ", dataSourceName=" + dataSourceName +
                ", locusLevel=" + locusLevel +
                ", apprisRank=" + apprisRank +
                ", transcriptLength=" + transcriptLength +
                ", version='" + version + '\'' +
                ", geneTranscriptType=" + geneTranscriptType +
                ", hugoSymbolSerializedOverride='" + hugoSymbolSerializedOverride + '\'' +
                ", ncbiBuildSerializedOverride='" + ncbiBuildSerializedOverride + '\'' +
                ", chromosomeSerializedOverride='" + chromosomeSerializedOverride + '\'' +
                ", startSerializedOverride='" + startSerializedOverride + '\'' +
                ", endSerializedOverride='" + endSerializedOverride + '\'' +
                ", variantClassificationSerializedOverride='" + variantClassificationSerializedOverride + '\'' +
                ", secondaryVariantClassificationSerializedOverride='" + secondaryVariantClassificationSerializedOverride + '\'' +
                ", variantTypeSerializedOverride='" + variantTypeSerializedOverride + '\'' +
                ", refAlleleSerializedOverride='" + refAlleleSerializedOverride + '\'' +
                ", tumorSeqAllele1SerializedOverride='" + tumorSeqAllele1SerializedOverride + '\'' +
                ", tumorSeqAllele2SerializedOverride='" + tumorSeqAllele2SerializedOverride + '\'' +
                ", genomeChangeSerializedOverride='" + genomeChangeSerializedOverride + '\'' +
                ", annotationTranscriptSerializedOverride='" + annotationTranscriptSerializedOverride + '\'' +
                ", transcriptStrandSerializedOverride='" + transcriptStrandSerializedOverride + '\'' +
                ", transcriptExonSerializedOverride='" + transcriptExonSerializedOverride + '\'' +
                ", transcriptPosSerializedOverride='" + transcriptPosSerializedOverride + '\'' +
                ", cDnaChangeSerializedOverride='" + cDnaChangeSerializedOverride + '\'' +
                ", codonChangeSerializedOverride='" + codonChangeSerializedOverride + '\'' +
                ", proteinChangeSerializedOverride='" + proteinChangeSerializedOverride + '\'' +
                ", gcContentSerializedOverride='" + gcContentSerializedOverride + '\'' +
                ", referenceContextSerializedOverride='" + referenceContextSerializedOverride + '\'' +
                ", otherTranscriptsSerializedOverride='" + otherTranscriptsSerializedOverride + '\'' +
                '}';
    }

    //==================================================================================================================

    public String getHugoSymbol() {
        return hugoSymbol;
    }

    public void setHugoSymbol(final String hugoSymbol) {
        this.hugoSymbol = hugoSymbol;
    }

    public String getNcbiBuild() {
        return ncbiBuild;
    }

    public void setNcbiBuild(final String ncbiBuild) {
        this.ncbiBuild = ncbiBuild;
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(final String chromosome) {
        this.chromosome = chromosome;
    }

    public int getStart() {
        return start;
    }

    public void setStart(final int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(final int end) {
        this.end = end;
    }

    public VariantClassification getVariantClassification() {
        return variantClassification;
    }

    public void setVariantClassification(final VariantClassification variantClassification) {
        this.variantClassification = variantClassification;
    }

    public VariantClassification getSecondaryVariantClassification() {
        return secondaryVariantClassification;
    }

    public void setSecondaryVariantClassification(final VariantClassification secondaryVariantClassification) {
        this.secondaryVariantClassification = secondaryVariantClassification;
    }

    public VariantType getVariantType() {
        return variantType;
    }

    public void setVariantType(final VariantType variantType) {
        this.variantType = variantType;
    }

    public String getRefAllele() {
        return refAllele;
    }

    public void setRefAllele(final String refAllele) {
        this.refAllele = refAllele;
    }

    public String getTumorSeqAllele2() {
        return tumorSeqAllele2;
    }

    public void setTumorSeqAllele2(final String tumorSeqAllele2) {
        this.tumorSeqAllele2 = tumorSeqAllele2;
    }

    public String getGenomeChange() {
        return genomeChange;
    }

    public void setGenomeChange(final String genomeChange) {
        this.genomeChange = genomeChange;
    }

    public String getAnnotationTranscript() {
        return annotationTranscript;
    }

    public void setAnnotationTranscript(final String annotationTranscript) {
        this.annotationTranscript = annotationTranscript;
    }

    public String getTranscriptStrand() {
        return transcriptStrand;
    }

    public void setTranscriptStrand(final String transcriptStrand) {
        this.transcriptStrand = transcriptStrand;
    }

    public Integer getTranscriptExonNumber() {
        return transcriptExon;
    }

    public void setTranscriptExonNumber(final Integer transcriptExonNumber) {
        this.transcriptExon = transcriptExonNumber;
    }

    public Integer getTranscriptPos() {
        return transcriptPos;
    }

    public void setTranscriptPos(final Integer transcriptPos) {
        this.transcriptPos = transcriptPos;
    }

    public String getcDnaChange() {
        return cDnaChange;
    }

    public void setcDnaChange(final String cDnaChange) {
        this.cDnaChange = cDnaChange;
    }

    public String getCodonChange() {
        return codonChange;
    }

    public void setCodonChange(final String codonChange) {
        this.codonChange = codonChange;
    }

    public String getProteinChange() {
        return proteinChange;
    }

    public void setProteinChange(final String proteinChange) {
        this.proteinChange = proteinChange;
    }

    public Double getGcContent() {
        return gcContent;
    }

    public void setGcContent(final Double gcContent) {
        this.gcContent = gcContent;
    }

    public List<String> getOtherTranscripts() {
        return otherTranscripts;
    }

    public void setOtherTranscripts(final List<String> otherTranscripts) {
        this.otherTranscripts = otherTranscripts;
    }

    public Integer getLocusLevel() {
        return locusLevel;
    }

    public void setLocusLevel(final Integer locusLevel) {
        this.locusLevel = locusLevel;
    }

    public GencodeGtfGeneFeature.FeatureTag getApprisRank() {
        return apprisRank;
    }

    public void setApprisRank(final GencodeGtfGeneFeature.FeatureTag apprisRank) {
        this.apprisRank = apprisRank;
    }

    public Integer getTranscriptLength() {
        return transcriptLength;
    }

    public void setTranscriptLength(final Integer transcriptLength) {
        this.transcriptLength = transcriptLength;
    }
    public String getReferenceContext() {
        return referenceContext;
    }

    public void setReferenceContext(final String referenceContext) {
        this.referenceContext = referenceContext;
    }

    public void setVersion(final String version) {
        this.version = version;
    }

    public GencodeGtfFeature.GeneTranscriptType getGeneTranscriptType() {
        return geneTranscriptType;
    }

    public void setGeneTranscriptType(final GencodeGtfFeature.GeneTranscriptType geneTranscriptType) {
        this.geneTranscriptType = geneTranscriptType;
    }

    public void setDataSourceName(final String dataSourceName) {
        this.dataSourceName = dataSourceName;
    }

    public void setMetadata(final FuncotationMetadata metadata) {
        this.metadata = metadata;
    }
//==================================================================================================================

    /**
     * The names of the fields in this GencodeFuncotation.
     * Used to access the field map.
     */
    private enum FieldName {
        hugoSymbol,
        ncbiBuild,
        chromosome,
        start,
        end,
        variantClassification,
        secondaryVariantClassification,
        variantType,
        refAllele,
        tumorSeqAllele1,
        tumorSeqAllele2,
        genomeChange,
        annotationTranscript,
        transcriptStrand,
        transcriptExon,
        transcriptPos,
        cDnaChange,
        codonChange,
        proteinChange,
        gcContent,
        referenceContext,
        otherTranscripts;
    }

    //==================================================================================================================

    public enum VariantType {
        INS("INS"),
        DEL("DEL"),
        SNP("SNP"),
        DNP("DNP"),
        TNP("TNP"),
        ONP("ONP"),
        MNP("MNP"),
        NA("NA");

        final private String serialized;

        VariantType(final String serialized) { this.serialized = serialized; }

        @Override
        public String toString() {
            return serialized;
        }
    }

    /**
     * Represents the type and severity of a variant.
     * The lower the {@link VariantClassification#relativeSeverity}, the more severe the variant.
     * Descriptions taken from:
     *     https://gatkforums.broadinstitute.org/gatk/discussion/8815/oncotator-variant-classification-and-secondary-variant-classification
     */
    public enum VariantClassification {

        /** Variant classification could not be determined. */
        COULD_NOT_DETERMINE("",99),

        // Only for Introns:
        /** Variant lies between exons within the bounds of the chosen transcript. */
        INTRON("INTRON", 10),

        // Only for UTRs:
        /** Variant is on the 5'UTR for the chosen transcript. */
        FIVE_PRIME_UTR("FIVE_PRIME_UTR", 6),
        /** Variant is on the 3'UTR for the chosen transcript */
        THREE_PRIME_UTR("THREE_PRIME_UTR", 6),

        // Only for IGRs:
        /** Intergenic region. Does not overlap any transcript. */
        IGR("IGR", 20),
        /** The variant is upstream of the chosen transcript (within 3kb) */
        FIVE_PRIME_FLANK("FIVE_PRIME_FLANK", 15),

        // These can be in Coding regions or Introns (only SPLICE_SITE):
        /** The point mutation alters the protein structure by one amino acid. */
        MISSENSE("MISSENSE", 1),
        /** A premature stop codon is created by the variant. */
        NONSENSE("NONSENSE", 0),
        /** Variant removes stop codon. */
        NONSTOP("NONSTOP", 0),
        /** Variant is in coding region of the chosen transcript, but protein structure is identical (i.e. a synonymous mutation). */
        SILENT("SILENT", 5),
        /** The variant is within a configurable number of bases (default=2) of a splice site. See the secondary classification to determine if it lies on the exon or intron side. */
        SPLICE_SITE("SPLICE_SITE", 4),
        /** Deletion that keeps the sequence in frame (i.e. deletion of a length evenly divisible by 3). */
        IN_FRAME_DEL("IN_FRAME_DEL", 1),
        /** Insertion that keeps the sequence in frame (i.e. insertion of a length evenly divisible by 3). */
        IN_FRAME_INS("IN_FRAME_INS", 1),
        /** Insertion that moves the coding sequence out of frame (i.e. insertion of a length not evenly divisible by 3). */
        FRAME_SHIFT_INS("FRAME_SHIFT_INS", 2),
        /** Deletion that moves the sequence out of frame (i.e. deletion of a length not evenly divisible by 3). */
        FRAME_SHIFT_DEL("FRAME_SHIFT_DEL", 2),
        /** Point mutation that overlaps the start codon. */
        START_CODON_SNP("START_CODON_SNP", 3),
        /** Insertion that overlaps the start codon. */
        START_CODON_INS("START_CODON_INS", 3),
        /** Deletion that overlaps the start codon. */
        START_CODON_DEL("START_CODON_DEL", 3),

        // These can only be in 5' UTRs:
        /** New start codon is created by the given variant using the chosen transcript. However, it is in frame relative to the coded protein. */
        DE_NOVO_START_IN_FRAME("DE_NOVO_START_IN_FRAME", 1),
        /** New start codon is created by the given variant using the chosen transcript. However, it is out of frame relative to the coded protein. */
        DE_NOVO_START_OUT_FRAME("DE_NOVO_START_OUT_FRAME", 0),

        // These are special / catch-all cases:
        /** Variant lies on one of the RNA transcripts. */
        RNA("RNA", 4),
        /** Variant lies on one of the lincRNAs. */
        LINCRNA("LINCRNA", 4);

        /**
         * The relative severity of each {@link VariantClassification}.
         * Lower numbers are considered more severe.
         * Higher numbers are considered less severe.
         */
        final private int relativeSeverity;

        /** The serialized version of this {@link VariantClassification} */
        final private String serialized;

        VariantClassification(final String serialized, final int sev) {
            this.serialized = serialized;
            relativeSeverity = sev;
        }

        /**
         * @return The {@link VariantClassification#relativeSeverity} of {@code this} {@link VariantClassification}.
         */
        public int getSeverity() { return relativeSeverity; }

        @Override
        public String toString() {
            return serialized;
        }
    }
}
