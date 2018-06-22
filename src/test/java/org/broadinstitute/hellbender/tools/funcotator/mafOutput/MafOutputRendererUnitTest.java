package org.broadinstitute.hellbender.tools.funcotator.mafOutput;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.collections.MapUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.*;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;

/**
 * Unit test class for the {@link MafOutputRenderer}.
 * Created by jonn on 1/22/18.
 */
public class MafOutputRendererUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    private final String OTHER_TRANSCRIPTS_RAW_NAMES =
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString()            + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString()            + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString()            + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString()         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString()         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.MISSENSE.toString()                + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.NONSENSE.toString()                + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.SILENT.toString()                  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.SPLICE_SITE.toString()             + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString()  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString() + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.START_CODON_SNP.toString()         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.START_CODON_INS.toString()         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.START_CODON_DEL.toString()         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.NONSTOP.toString()                 + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString()          + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString()         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString()        + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.INTRON.toString()                  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.LINCRNA.toString();

    private final String OTHER_TRANSCRIPTS_MAF_READY_NAMES =
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString())     + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString())     + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString())     + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString())  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString())  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.MISSENSE.toString())         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSENSE.toString())         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SILENT.toString())           + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SPLICE_SITE.toString())      + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString()                                                    + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString()                                                   + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.START_CODON_SNP.toString()                                                           + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.START_CODON_INS.toString()                                                           + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.START_CODON_DEL.toString())  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSTOP.toString())          + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString())   + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString())  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString()) + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.INTRON.toString())           + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.LINCRNA.toString());

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    private List<String> createFieldValuesFromNameList(final String prefix, final List<String> baseFieldList, final int fieldSize) {
        final List<String> outList = new ArrayList<>(baseFieldList.size());

        for ( int i = 0; i < baseFieldList.size() ; ++i ) {
            final String formatString = "%s%0" +
                    ((fieldSize - prefix.length()) > 0 ? fieldSize - prefix.length() : "") +
                    "d";
            outList.add(String.format(formatString, prefix, i+1));
        }

        return outList;
    }

    private MafOutputRenderer createMafOutputRenderer(final File outputFile) {

        final Map<Path, Properties> configData =
                DataSourceUtils.getAndValidateDataSourcesFromPaths(
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        Collections.singletonList(FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER)
                );

        final List<DataSourceFuncotationFactory> dataSourceFuncotationFactories = DataSourceUtils.createDataSourceFuncotationFactoriesForDataSources(
                configData,
                new LinkedHashMap<>(),
                TranscriptSelectionMode.BEST_EFFECT,
                new HashSet<>()
        );

        // Sort the datasources to ensure the same order every time:
        dataSourceFuncotationFactories.sort( Comparator.comparing(DataSourceFuncotationFactory::getName) );

        return new MafOutputRenderer(
                IOUtils.getPath(outputFile.toURI().toString()),
                dataSourceFuncotationFactories,
                new VCFHeader(),
                new LinkedHashMap<>(),
                new LinkedHashMap<>(),
                new HashSet<>()
        );
    }

    private MafOutputRenderer createMafOutputRenderer() {
        return createMafOutputRenderer(getSafeNonExistentFile("TestMafOutputFile"));
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Object[][] provideForTestReplaceFuncotationValuesWithMafCompliantValues() {
        return new Object[][] {
                // Empty maps:
                {
                    new HashMap<>(),
                    new LinkedHashMap<>()
                },
                // Singleton map that doesn't contain a replaceable element:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        {"Unreplaceable_Key", "Unreplaceable_Value"},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"Unreplaceable_Key", "Unreplaceable_Value"},
                                }
                        )
                },
                // Map with multiple elements, none of which are replaced:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        {"Unreplaceable_Key1", "Unreplaceable_Value1"},
                                        {"Unreplaceable_Key2", "Unreplaceable_Value2"},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"Unreplaceable_Key1", "Unreplaceable_Value1"},
                                        {"Unreplaceable_Key2", "Unreplaceable_Value2"},
                                }
                        )
                },
                // Singleton map containing a replaced value:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Other_Transcripts, OTHER_TRANSCRIPTS_RAW_NAMES },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Other_Transcripts, OTHER_TRANSCRIPTS_MAF_READY_NAMES.replaceAll(VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER, MafOutputRendererConstants.OTHER_TRANSCRIPT_DELIMITER)},
                                }
                        )
                },
                // Map with multiple elements, some of which are replaced:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        {"Unreplaceable_Key1", "Unreplaceable_Value1"},
                                        { MafOutputRendererConstants.FieldName_Other_Transcripts, OTHER_TRANSCRIPTS_RAW_NAMES },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"Unreplaceable_Key1", "Unreplaceable_Value1"},
                                        { MafOutputRendererConstants.FieldName_Other_Transcripts, OTHER_TRANSCRIPTS_MAF_READY_NAMES.replaceAll(VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER, MafOutputRendererConstants.OTHER_TRANSCRIPT_DELIMITER)},
                                }
                        )
                },
                // Map with multiple elements, all of which are replaced:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Other_Transcripts, OTHER_TRANSCRIPTS_RAW_NAMES },
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.INS.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "A" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "A" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "AT" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "100" },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Other_Transcripts, OTHER_TRANSCRIPTS_MAF_READY_NAMES.replaceAll(VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER, MafOutputRendererConstants.OTHER_TRANSCRIPT_DELIMITER)},
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.INS.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "-" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "-" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "T" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "101" },
                                }
                        )
                }
        };
    }

    @DataProvider
    private Object[][] provideForAdjustIndelAlleleInformationForMafOutput() {
        return new Object[][] {
                // Empty maps:
                {
                        new LinkedHashMap<>(),
                        new LinkedHashMap<>()
                },
                // Singleton map that doesn't contain a replaceable element:
                {
                        MapUtils.putAll(new LinkedHashMap<String, Object>(),
                                new Object[][] {
                                        {"Unreplaceable_Key", "Unreplaceable_Value"},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"Unreplaceable_Key", "Unreplaceable_Value"},
                                }
                        )
                },
                // Singleton map that contains a variant type that is not replaced:
                {
                        MapUtils.putAll(new LinkedHashMap<String, Object>(),
                                new Object[][] {
                                        {MafOutputRendererConstants.FieldName_Variant_Type, GencodeFuncotation.VariantType.SNP.toString()},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {MafOutputRendererConstants.FieldName_Variant_Type, GencodeFuncotation.VariantType.SNP.toString()},
                                }
                        )
                },
                // Map with multiple elements, none of which are replaced:
                {
                        MapUtils.putAll(new LinkedHashMap<String, Object>(),
                                new Object[][] {
                                        {"Unreplaceable_Key1", "Unreplaceable_Value1"},
                                        {MafOutputRendererConstants.FieldName_Variant_Type, GencodeFuncotation.VariantType.SNP.toString()},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"Unreplaceable_Key1", "Unreplaceable_Value1"},
                                        {MafOutputRendererConstants.FieldName_Variant_Type, GencodeFuncotation.VariantType.SNP.toString()},
                                }
                        )
                },
                // Map contains an INSERTION of size 1 with associated fields:
                {
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.INS.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "A" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "A" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "AT" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "100" },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.INS.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "-" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "-" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "T" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "101" },
                                }
                        )
                },
                // Map contains a DELETION of size 1 with associated fields:
                {
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.DEL.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "AA" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "AA" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "A" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "100" },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.DEL.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "A" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "A" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "-" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "101" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "101" },
                                }
                        )
                },
                // Map contains an INSERTION of size 5 with associated fields:
                {
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.INS.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "AAAAA" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "AAAAA" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "AAAAATTTTT" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "104" },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.INS.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "-" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "-" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "TTTTT" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "101" },
                                }
                        )
                },
                // Map contains a DELETION of size 5 with associated fields:
                {
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.DEL.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "ATTTTT" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "ATTTTT" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "A" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "100" },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.DEL.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "TTTTT" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "TTTTT" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "-" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "101" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "105" },
                                }
                        )
                },
        };
    }

    @DataProvider
    private Object[][] provideForMafTransform() {
        return new Object[][] {
                // Empty Strings:
                {
                    "", "", ""
                },

                // Field has no replacements
                {
                        "SPECIAL_NO_REPLACEY_FIELD", "GARBAGE", "GARBAGE"
                },

                // Field containing replaceable string from another field:
                {
                    MafOutputRendererConstants.FieldName_Variant_Classification, MafOutputRendererConstants.FieldValue_Gencode_Chromosome_Mito, MafOutputRendererConstants.FieldValue_Gencode_Chromosome_Mito
                },

                // Variant classification tests:
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString(),
                        "In_Frame_Del"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString(),
                        "In_Frame_Ins"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString(),
                        "Frame_Shift_Ins"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString(),
                        "Frame_Shift_Del"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.MISSENSE.toString(),
                        "Missense_Mutation"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.NONSENSE.toString(),
                        "Nonsense_Mutation"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.SILENT.toString(),
                        "Silent"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.SPLICE_SITE.toString(),
                        "Splice_Site"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.START_CODON_DEL.toString(),
                        "Translation_Start_Site"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.NONSTOP.toString(),
                        "Nonstop_Mutation"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString(),
                        "5'UTR"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString(),
                        "3'UTR"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString(),
                        "5'Flank"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.INTRON.toString(),
                        "Intron"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.LINCRNA.toString(),
                        "RNA"
                },

                // Chromosome tests:
                {
                    MafOutputRendererConstants.FieldName_Chromosome,
                    MafOutputRendererConstants.FieldValue_Gencode_Chromosome_Mito,
                    MafOutputRendererConstants.FieldValue_Chromosome_Mito
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chr7",
                        "7"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "ChR7",
                        "7"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chracteristic",
                        "acteristic"
                },

                // Other transcripts:
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.MISSENSE.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.MISSENSE.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.NONSENSE.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSENSE.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.SILENT.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SILENT.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.SPLICE_SITE.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SPLICE_SITE.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.START_CODON_DEL.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.START_CODON_DEL.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.NONSTOP.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSTOP.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.INTRON.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.INTRON.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.LINCRNA.toString(),
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.LINCRNA.toString())
                },

                // Fields that we do not replace:
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString(),
                        GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString()
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString(),
                        GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString()
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.START_CODON_SNP.toString(),
                        GencodeFuncotation.VariantClassification.START_CODON_SNP.toString()
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.START_CODON_INS.toString(),
                        GencodeFuncotation.VariantClassification.START_CODON_INS.toString()
                },

                // Multiple other transcripts:
                {
                    MafOutputRendererConstants.FieldName_Other_Transcripts,
                        OTHER_TRANSCRIPTS_RAW_NAMES,
                        OTHER_TRANSCRIPTS_MAF_READY_NAMES
                }
        };
    }

    @DataProvider
    private Object[][] provideForWrite() {

        final VariantContext variant1 = new VariantContextBuilder().chr("chr1").start(25).stop(25).alleles("A", "T").make();
        final VariantContext variant2 = new VariantContextBuilder().chr("chrX").start(1759).stop(1759).alleles("G", "C").make();
        final VariantContext variant3 = new VariantContextBuilder().chr("chr3").start(89247).stop(89248).alleles("AA", "GG").make();
        final VariantContext variant4 = new VariantContextBuilder().chr("chr6").start(888996).stop(888999).alleles("TCCC", "T").make();
        final VariantContext variant5 = new VariantContextBuilder().chr("chr9").start(1067897).stop(1067899).alleles("AGG", "A", "AG", "AT").make();

        final List<VariantContext> variantList = Arrays.asList(variant1, variant2, variant3, variant4, variant5);

        final List<String> baseFieldNameList = new ArrayList<>(createMafOutputRenderer( getSafeNonExistentFile("GARBAGE") ).getDefaultMap().keySet());

        final int fieldSize = 10;
        
        // NOTE: The data field names must match data sources that are checked in for this to work in an expected way:
        final List<List<Funcotation>> funcotationList = Arrays.asList(
                Collections.singletonList(
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("A", baseFieldNameList, fieldSize),
                                Allele.create("T"),
                                GencodeFuncotationFactory.DEFAULT_NAME, null
                        )
                ),
                Collections.singletonList(
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("B", baseFieldNameList, fieldSize),
                                Allele.create("C"),
                                GencodeFuncotationFactory.DEFAULT_NAME, null
                        )
                ),
                Collections.singletonList(
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("C", baseFieldNameList, fieldSize),
                                Allele.create("GG"),
                                GencodeFuncotationFactory.DEFAULT_NAME, null
                        )
                ),
                Collections.singletonList(
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("D", baseFieldNameList, fieldSize),
                                Allele.create("T"),
                                "TestDataSource4", null
                        )
                ),
                Arrays.asList(
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("E", baseFieldNameList, fieldSize),
                                Allele.create("A"),
                                "TestDataSource5", null
                        ),
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("F", baseFieldNameList, fieldSize),
                                Allele.create("AG"),
                                "TestDataSource5", null
                        ),
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("G", baseFieldNameList, fieldSize),
                                Allele.create("AT"),
                                "TestDataSource5", null
                        )
                )
        );

        return new Object[][] {
                {
                    variantList,
                    funcotationList,
                    new File(FuncotatorTestConstants.FUNCOTATOR_TEST_DIR + File.separator + "ExampleMafFileForTests.maf")
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestReplaceFuncotationValuesWithMafCompliantValues")
    public void testReplaceFuncotationValuesWithMafCompliantValues(final Map<String, Object> outputMap, final LinkedHashMap<String, String> expected ) {
        final LinkedHashMap<String, String> compliantMap = createMafOutputRenderer().replaceFuncotationValuesWithMafCompliantValues(outputMap);
        Assert.assertEquals(compliantMap, expected);
    }

    @Test(dataProvider = "provideForMafTransform")
    public void testMafTransform(final String key, final String value, final String expectedValue) {
        final String transformedValue = createMafOutputRenderer().mafTransform(key, value);
        Assert.assertEquals(transformedValue, expectedValue);
    }

    @Test(dataProvider = "provideForAdjustIndelAlleleInformationForMafOutput")
    public void testAdjustIndelAlleleInformationForMafOutput(final LinkedHashMap<String, String> outputMap, final LinkedHashMap<String, String> expected) {
        createMafOutputRenderer().adjustIndelAlleleInformationForMafOutput(outputMap);
        Assert.assertEquals(outputMap, expected);
    }

    @Test(dataProvider = "provideForWrite")
    public void testWrite(final List<VariantContext> variants, final List<List<Funcotation>> funcotations, final File expectedFile) {

        final File outFile = getSafeNonExistentFile("TestMafOutputFile");
        try ( final MafOutputRenderer mafOutputRenderer = createMafOutputRenderer( outFile ) ) {
            for ( int i = 0 ; i < variants.size(); ++i ) {
                final FuncotationMap funcotationMap = FuncotationMap.createNoTranscriptInfo(funcotations.get(i));
                mafOutputRenderer.write(variants.get(i), funcotationMap);
            }
        }

        try {
            IntegrationTestSpec.assertEqualTextFiles(outFile, expectedFile, "#");
        }
        catch (final IOException ex) {
            throw new GATKException("ERROR comparing text files: " + outFile.toURI().toString() + " and " + expectedFile.toURI().toString(), ex);
        }
    }
}
