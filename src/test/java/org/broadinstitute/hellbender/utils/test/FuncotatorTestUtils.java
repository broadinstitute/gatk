package org.broadinstitute.hellbender.utils.test;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationBuilder;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class FuncotatorTestUtils {
    private FuncotatorTestUtils() {}


    /**
     * Since funcotation factories need an instance of {@link FeatureContext} to funcotate, this convenience method can
     *  create a new instance for test methods.
     *
     * @param funcotationFactories {@link List} of {@link DataSourceFuncotationFactory} that should be used to generate the
     *                                         {@link FeatureContext}.  Never {@code null}, but empty list is acceptable.
     * @param dummyToolInstanceName A name to use for the "tool".  Any string will work here.  Never {@code null}.
     * @param interval genomic interval for the result.  Typically, this would be the interval of the variant.  Never {@link null}.
     * @param featureQueryLookahead When querying FeatureDataSources, cache this many extra bases of context beyond
     *                              the end of query intervals in anticipation of future queries. Must be >= 0.  If uncertain, use zero.
     * @param cloudPrefetchBuffer See {@link FeatureManager#FeatureManager(CommandLineProgram, int, int, int, Path)}  If uncertain, use zero.
     * @param cloudIndexPrefetchBuffer See {@link FeatureManager#FeatureManager(CommandLineProgram, int, int, int, Path)}  If uncertain, use zero.
     * @param reference See {@link FeatureManager#FeatureManager(CommandLineProgram, int, int, int, Path)}  If uncertain, use {@code null}.
     * @return a {@link FeatureContext} ready for querying the funcotation factories on the given interval.  Never {@code null}.
     */
    @VisibleForTesting
    public static FeatureContext createFeatureContext(final List<DataSourceFuncotationFactory> funcotationFactories, final String dummyToolInstanceName,
                                                      final SimpleInterval interval, final int featureQueryLookahead, final int cloudPrefetchBuffer,
                                                      final int cloudIndexPrefetchBuffer, final Path reference) {
        Utils.nonNull(funcotationFactories);
        Utils.nonNull(dummyToolInstanceName);
        Utils.nonNull(interval);

        final Map<FeatureInput<? extends Feature>, Class<? extends Feature>> featureInputsWithType =
                funcotationFactories.stream()
                        .collect(Collectors.toMap(ff -> ff.getMainSourceFileAsFeatureInput(), ff -> ff.getAnnotationFeatureClass()));

        return FeatureContext.createFeatureContextForTesting(featureInputsWithType, dummyToolInstanceName, interval,
                featureQueryLookahead, cloudPrefetchBuffer, cloudIndexPrefetchBuffer, reference);
    }

    /**
     * Create a dummy gencode funcotation in the most manual way possible.  This method does not check for consistency
     *  or reasonable values.
     *
     * @param hugoSymbol
     * @param ncbiBuild
     * @param chromosome
     * @param start
     * @param end
     * @param variantClassification
     * @param secondaryVariantClassification
     * @param variantType
     * @param refAllele
     * @param tumorSeqAllele2
     * @param genomeChange
     * @param annotationTranscript
     * @param transcriptStrand
     * @param transcriptExon
     * @param transcriptPos
     * @param cDnaChange
     * @param codonChange
     * @param proteinChange
     * @param gcContent
     * @param referenceContext
     * @param otherTranscripts
     * @param version Use "19" for hg19/b37 and "28" for hg38.
     * @return
     */
    public static GencodeFuncotation createGencodeFuncotation(final String hugoSymbol, final String ncbiBuild,
                                                                             final String chromosome, final int start, final int end,
                                                                             final GencodeFuncotation.VariantClassification variantClassification,
                                                                             final GencodeFuncotation.VariantClassification secondaryVariantClassification,
                                                                             final GencodeFuncotation.VariantType variantType,
                                                                             final String refAllele,
                                                                             final String tumorSeqAllele2, final String genomeChange,
                                                                             final String annotationTranscript, final Strand transcriptStrand,
                                                                             final Integer transcriptExon, final Integer transcriptPos,
                                                                             final String cDnaChange, final String codonChange,
                                                                             final String proteinChange, final Double gcContent,
                                                                             final String referenceContext,
                                                                             final List<String> otherTranscripts, final String version) {

        final GencodeFuncotationBuilder funcotationBuilder = new GencodeFuncotationBuilder();

        funcotationBuilder.setVersion(version);
        funcotationBuilder.setDataSourceName(GencodeFuncotationFactory.DEFAULT_NAME);

        funcotationBuilder.setHugoSymbol( hugoSymbol );
        funcotationBuilder.setNcbiBuild( ncbiBuild );
        funcotationBuilder.setChromosome( chromosome );
        funcotationBuilder.setStart( start );
        funcotationBuilder.setEnd( end );
        funcotationBuilder.setVariantClassification( variantClassification );
        funcotationBuilder.setSecondaryVariantClassification(secondaryVariantClassification);
        funcotationBuilder.setVariantType( variantType );
        funcotationBuilder.setRefAllele(Allele.create(refAllele));
        funcotationBuilder.setTumorSeqAllele2( tumorSeqAllele2 );

        funcotationBuilder.setGenomeChange( genomeChange );
        funcotationBuilder.setAnnotationTranscript( annotationTranscript );
        funcotationBuilder.setStrand( transcriptStrand );
        funcotationBuilder.setTranscriptExonNumber( transcriptExon );
        funcotationBuilder.setTranscriptPos( transcriptPos );
        funcotationBuilder.setcDnaChange( cDnaChange );
        funcotationBuilder.setCodonChange( codonChange );
        funcotationBuilder.setProteinChange( proteinChange );
        funcotationBuilder.setGcContent( gcContent );
        funcotationBuilder.setReferenceContext( referenceContext );
        funcotationBuilder.setOtherTranscripts( otherTranscripts );

        return funcotationBuilder.build();
    }

    /**
     *  Create a variant context with the following fields.  Note that the genotype will be empty, as will
     *  INFO annotations.
     *
     * @param reference E.g. {@link FuncotatorReferenceTestUtils#retrieveHg19Chr3Ref}
     * @param contig Never {@code null}
     * @param start Must be positive.
     * @param end Must be positive.
     * @param refAlleleString Valid string for an allele.  Do not include "*".  Never {@code null}
     * @param altAlleleString Valid string for an allele.  Do not include "*".  Never {@code null}
     * @return a simple, biallelic variant context
     */
    public static VariantContext createSimpleVariantContext(final String reference, final String contig,
                                                            final int start,
                                                            final int end,
                                                            final String refAlleleString,
                                                            final String altAlleleString) {
        Utils.nonNull(contig);
        Utils.nonNull(refAlleleString);
        Utils.nonNull(altAlleleString);
        ParamUtils.isPositive(start, "Invalid start position: " + start);
        ParamUtils.isPositive(end, "Invalid end position: " + end);

        final Allele refAllele = Allele.create(refAlleleString, true);
        final Allele altAllele = Allele.create(altAlleleString);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                reference,
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );

        return variantContextBuilder.make();
    }
}
