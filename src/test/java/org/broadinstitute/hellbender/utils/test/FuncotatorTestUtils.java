package org.broadinstitute.hellbender.utils.test;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationBuilder;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeature;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeatureBaseData;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfTranscriptFeature;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
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
     *  This method is meant for use in test cases where the developer can provide reasonable values.  Please use
     *   caution as it is easy to specify parameters that would result in a funcotation that would never be seen in the wild.
     *
     *  {@code null} is not recommended, but this method will not throw an exception.
     *
     * @param hugoSymbol gene name, such as PIK3CA or TEST_GENE
     * @param ncbiBuild Associated build.  Typically, "hg38" or "hg19"
     * @param chromosome Contig name.  Please note that this method will not enforce the value here.
     * @param start Must be greater than zero.  No other checks provided.  For example, this method will not check whether a chromosome is of an appropriate length for this value.
     * @param end Must be greater than zero.  No other checks provided.  For example, this method will not check whether a chromosome is of an appropriate length for this value.
     * @param variantClassification Variant classification.
     * @param secondaryVariantClassification Typically, only used for {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification#SPLICE_SITE} and other variant classifications that need more information to disambiguate.  For example, a splice site can be in an exon or an intron.  So if the Splice site is in an intron, this field would be {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification#INTRON}
     * @param variantType Variant type to use.  Again, this method will not check consistency.  For example, this method will not throw an exception if you have {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantType#SNP} with an indel-specific variant classification.
     * @param refAllele reference allele
     * @param tumorSeqAllele2 alternate allele
     * @param genomeChange Genome change string.  This method will not check validity.
     * @param annotationTranscript Transcript ID to use.  If you wish to fake a no transcript gencode funcotation, use {@link org.broadinstitute.hellbender.tools.funcotator.FuncotationMap#NO_TRANSCRIPT_AVAILABLE_KEY}
     * @param transcriptStrand This method will not check validity against the annotation transcript.
     * @param transcriptExon 1-based exon number.  In the Funcotator main code, this field takes into account coding direction, so exon 1 will be at the highest genomic position on a negative strand.
     *                       This method will not check validity of the specified value.  Must be >= 1.
     * @param transcriptStartPos 1-based position of the start of the variant in the transcript genome sequence.  In the Funcotator main code, this field takes into account the coding direction, so it is possible to have
     *                      a higher transcript position that is a lower position in genomic coordinate space.
     *                      This method will not check validity against the annotation transcript and transcript exon. Must be >= 1.
     * @param transcriptEndPos 1-based position of the end of the variant in the transcript genome sequence.  In the Funcotator main code, this field takes into account the coding direction, so it is possible to have
     *                      a higher transcript position that is a lower position in genomic coordinate space.
     *                      This method will not check validity against the annotation transcript and transcript exon. Must be >= 1.
     * @param cDnaChange This method will not check validity of the specified value.
     * @param codonChange This method will not check validity of the specified value.
     * @param proteinChange This method will not check validity of the specified value.
     * @param gcContent This method will not check validity.  Must be >= 0.0 and =< 1.0
     * @param referenceContext This method will not check validity of the specified value.  This method will not even check that you are specifying valid bases.
     * @param otherTranscripts This method will not check validity of the specified value.
     * @param version Gencode version.  Use "19" for hg19/b37 and "28" for hg38.  This method will not check validity,
     *                but much of the rendering code will fail if this field is not valid.
     * @return a gencode funcotation representing the above parameters.  Never {@code null}
     */
    public static GencodeFuncotation createGencodeFuncotation(final String hugoSymbol, final String ncbiBuild,
                                                                 final String chromosome, final int start, final int end,
                                                                 final GencodeFuncotation.VariantClassification variantClassification,
                                                                 final GencodeFuncotation.VariantClassification secondaryVariantClassification,
                                                                 final GencodeFuncotation.VariantType variantType,
                                                                 final String refAllele,
                                                                 final String tumorSeqAllele2, final String genomeChange,
                                                                 final String annotationTranscript, final Strand transcriptStrand,
                                                                 final Integer transcriptExon,
                                                                  final Integer transcriptStartPos,
                                                                  final Integer transcriptEndPos,
                                                                 final String cDnaChange, final String codonChange,
                                                                 final String proteinChange, final Double gcContent,
                                                                 final String referenceContext,
                                                                 final List<String> otherTranscripts, final String version) {


        ParamUtils.isPositive(start, "Start position is 1-based and must be greater that zero.");
        ParamUtils.isPositive(end, "End position is 1-based and must be greater that zero.");
        ParamUtils.isPositive(transcriptExon, "Transcript exon is a 1-based index.");
        ParamUtils.isPositive(transcriptStartPos, "Transcript start position is a 1-based index.");
        ParamUtils.isPositive(transcriptEndPos, "Transcript end position is a 1-based index.");
        ParamUtils.inRange(gcContent, 0.0, 1.0, "GC Content must be between 0.0 and 1.0.");

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
        funcotationBuilder.setTranscriptStartPos( transcriptStartPos );
        funcotationBuilder.setTranscriptEndPos( transcriptEndPos );
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
        return createSimpleVariantContext(reference, contig, start, end, refAlleleString, altAlleleString, null);
    }

    /**
     *  Create a variant context with the following fields.  Note that the genotype will be empty, as will
     *  INFO annotations.
     *
     * @param reference E.g. {@link FuncotatorReferenceTestUtils#retrieveHg19Chr3Ref}
     * @param contig Never {@code null}
     * @param start Must be positive.
     * @param end Must be positive.
     * @param refAlleleString Valid {@link String} for an allele.  Do not include "*".  Never {@code null}
     * @param altAlleleString Valid {@link String} for an allele.  Do not include "*".  Never {@code null}
     * @param filterString Valid {@link String} for filters (See VCF 4.2 spec).  May be {@code null}.
     * @return a simple, biallelic variant context
     */
    public static VariantContext createSimpleVariantContext(final String reference, final String contig,
                                                            final int start,
                                                            final int end,
                                                            final String refAlleleString,
                                                            final String altAlleleString,
                                                            final String filterString) {
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

        if ( filterString != null ) {
            variantContextBuilder.filter(filterString);
        }
        else {
            variantContextBuilder.passFilters();
        }

        final VariantContext vc = variantContextBuilder.make();

        return vc;
    }

    /**
     * Create a {@link ReferenceContext} object from a sequence of bases and a location.
     * @param sequence {@link String} of bases from which to create the {@link ReferenceContext} data.
     * @param contigName {@link String} for the contig name for the {@link ReferenceContext} location.
     * @param refStartPos Start postiion (1-based, inclusive) for the location of the {@link ReferenceContext}.
     * @param refEndPos End postiion (1-based, inclusive) for the location of the {@link ReferenceContext}.
     * @return A {@link ReferenceContext} containing the given {@code sequence} at the given location.
     */
    public static ReferenceContext createReferenceContextFromBasesAndLocation(final String sequence,
                                                                              final String contigName,
                                                                              final int refStartPos,
                                                                              final int refEndPos) {
        // Create an in-memory ReferenceContext:
        final SimpleInterval wholeReferenceInterval = new SimpleInterval( contigName, 1, sequence.length() );

        return new ReferenceContext(
                new ReferenceMemorySource(
                        new ReferenceBases(sequence.getBytes(), wholeReferenceInterval),
                        new SAMSequenceDictionary(
                                Collections.singletonList(
                                        new SAMSequenceRecord(contigName, sequence.length())
                                )
                        )
                ),
                new SimpleInterval( contigName, refStartPos, refEndPos )
        );
    }

    /**
     * Creates an artifical GencodeGtfTranscriptFeature for testing with dummy values for all fields except
     * for the contig, start, stop, and strand.
     *
     * @param contig Contig that should be assigned to the new GencodeGtfTranscriptFeature
     * @param start Start position that should be assigned to the new GencodeGtfTranscriptFeature
     * @param stop Stop position that should be assigned to the new GencodeGtfTranscriptFeature
     * @param strand Strand that should be assigned to the new GencodeGtfTranscriptFeature
     * @return A new GencodeGtfTranscriptFeature with the specified contig, start, stop, and strand, and dummy
     *         values for all other fields
     */
    public static GencodeGtfTranscriptFeature createArtificialGencodeGtfTranscriptFeatureForTesting( final String contig, final int start, final int stop, final Strand strand) {
        return (GencodeGtfTranscriptFeature)GencodeGtfTranscriptFeature.create(
                new GencodeGtfFeatureBaseData(
                        2,
                        contig,
                        GencodeGtfFeature.AnnotationSource.ENSEMBL,
                        GencodeGtfFeature.FeatureType.TRANSCRIPT,
                        start,
                        stop,
                        strand,
                        GencodeGtfFeature.GenomicPhase.DOT,
                        "FakeGeneID",
                        "FakeTranscriptID",
                        GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                        null,
                        "FakeGeneName",
                        GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                        null,
                        "FakeTranscriptName",
                        -1,
                        null,
                        GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                        Collections.emptyList(),
                        null)
        );
    }
}
