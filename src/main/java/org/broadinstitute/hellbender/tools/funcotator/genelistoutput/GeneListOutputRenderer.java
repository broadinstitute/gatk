package org.broadinstitute.hellbender.tools.funcotator.genelistoutput;


import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import net.greypanther.natsort.SimpleNaturalComparator;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotationMap;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;
import org.broadinstitute.hellbender.tools.funcotator.OutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.simpletsvoutput.SimpleTsvOutputRenderer;
import org.broadinstitute.hellbender.utils.Utils;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/** This class can
 * - only work on segments.
 * - requires a "genes" funcotation field on each segment.  Additionally, this field must be a comma-delimited list of
 * gene names.  Also, requires, start/end exons, and start/end genes.  For a total of five fields.
 * - requires that the input FuncotationMap in {@link GeneListOutputRenderer#write(VariantContext, FuncotationMap)} only
 * have one transcript ID and that is {@link FuncotationMap#NO_TRANSCRIPT_AVAILABLE_KEY}.
 * - each segment can only have one alternate allele across all of the funcotations.
 *
 * If the above conditions are violated, even for a single funcotationmap/segment, this output renderer will throw an
 *  exception.
 *
 * If a gene is not covered by the segments (i.e. if a gene does not appear in any segment genes funcotation field), the
 *  gene will not be inclulded in the output.
 *
 *  Genes will be sorted alphabetically.
 *
 *  Currently, this output renderer only supports gencode genes.
 *
 */
public class GeneListOutputRenderer extends OutputRenderer {

    @VisibleForTesting
    final public static String CONFIG_RESOURCE = "org/broadinstitute/hellbender/tools/funcotator/gene_list_output.config";
    final private static String GENE_FIELD_NAME = "gene";
    final private static String EXON_FIELD_NAME = "exon";
    final private static String GENE_LIST_OUTPUT_RENDERER_DUMMY_NAME = "GENE_LIST_OUTPUT_RENDERER";
    final private static String GENCODE_GENES_FIELD_PATTERN_STRING = "Gencode_[0-9]+" + GencodeFuncotationFactory.GENES_SUFFIX;
    final private static String GENCODE_START_GENE_FIELD_PATTERN_STRING = "Gencode_[0-9]+" + GencodeFuncotationFactory.START_GENE_SUFFIX;
    final private static String GENCODE_START_EXON_FIELD_PATTERN_STRING = "Gencode_[0-9]+" + GencodeFuncotationFactory.START_EXON_SUFFIX;
    final private static String GENCODE_END_GENE_FIELD_PATTERN_STRING = "Gencode_[0-9]+" + GencodeFuncotationFactory.END_GENE_SUFFIX;
    final private static String GENCODE_END_EXON_FIELD_PATTERN_STRING = "Gencode_[0-9]+" + GencodeFuncotationFactory.END_EXON_SUFFIX;
    private static final String ERR_STR_IS_THIS_SEG_FILE = "Is the input a file of segment variant contexts?";

    // Validation error messages.  Warning:  These are used in regular expressions for testing.
    @VisibleForTesting
    static final String MISSING_GENCODE_FIELD_ERROR_MSG = "Could not find exactly one gencode field match in the funcotation map";
    @VisibleForTesting
    static final String ONE_ALTERNATE_ALLELE_ERR_MSG = "  Only one alternate allele per variant context is accepted.";
    @VisibleForTesting
    static final String INVALID_VARIANT_CONTEXT_ERR_MSG = "  Variant context does not represent a copy number segment";
    @VisibleForTesting
    static final String INCORRECT_NUM_TRANSCRIPT_ERR_MSG = "  Need exactly one transcript ID";
    @VisibleForTesting
    static final String INVALID_TRANSCRIPT_ID_ERR_MSG = "  Invalid transcript ID seen";

    final private SortedMap<Pair<String,String>, Pair<VariantContext, FuncotationMap>> geneExonToVariantFuncotationMap;
    final private SimpleTsvOutputRenderer simpleTsvOutputRenderer;
    final private Pattern gencodeGenesFieldPattern;
    final private Pattern gencodeStartExonFieldPattern;
    final private Pattern gencodeStartGeneFieldPattern;
    final private Pattern gencodeEndGeneFieldPattern;
    final private Pattern gencodeEndExonFieldPattern;

    final private int minNumBasesForValidSegment;

    public GeneListOutputRenderer(final Path outputFilePath,
                                  final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                                  final LinkedHashMap<String, String> unaccountedForOverrideAnnotations,
                                  final Set<String> excludedOutputFields, final String toolVersion){
        this(outputFilePath, unaccountedForDefaultAnnotations, unaccountedForOverrideAnnotations,
                excludedOutputFields, toolVersion, FuncotatorUtils.DEFAULT_MIN_NUM_BASES_FOR_VALID_SEGMENT);
    }

    public GeneListOutputRenderer(final Path outputFilePath,
                                  final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                                  final LinkedHashMap<String, String> unaccountedForOverrideAnnotations,
                                  final Set<String> excludedOutputFields, final String toolVersion,
                                  final int minNumBasesForValidSegment){
        super(toolVersion);
        Utils.nonNull(unaccountedForDefaultAnnotations);
        Utils.nonNull(unaccountedForOverrideAnnotations);
        Utils.nonNull(excludedOutputFields);
        Utils.nonNull(toolVersion);
        Utils.nonNull(outputFilePath);
        IOUtil.assertFileIsWritable(outputFilePath.toFile());

        geneExonToVariantFuncotationMap = new TreeMap<>(new PairedSimpleNaturalComparator<>());

        // When initializing the simpleoutput renderer, make sure that extra funcotation fields not explicitly in the
        //  configuration are not written.
        simpleTsvOutputRenderer = SimpleTsvOutputRenderer.createFromResource(outputFilePath,
            unaccountedForDefaultAnnotations, unaccountedForOverrideAnnotations, excludedOutputFields,
            Paths.get(CONFIG_RESOURCE), toolVersion, false);

        // Cache patterns to determine which field
        gencodeGenesFieldPattern = Pattern.compile(GENCODE_GENES_FIELD_PATTERN_STRING);
        gencodeStartGeneFieldPattern = Pattern.compile(GENCODE_START_GENE_FIELD_PATTERN_STRING);
        gencodeStartExonFieldPattern = Pattern.compile(GENCODE_START_EXON_FIELD_PATTERN_STRING);
        gencodeEndGeneFieldPattern = Pattern.compile(GENCODE_END_GENE_FIELD_PATTERN_STRING);
        gencodeEndExonFieldPattern = Pattern.compile(GENCODE_END_EXON_FIELD_PATTERN_STRING);

        this.minNumBasesForValidSegment = minNumBasesForValidSegment;
    }

    /**
     * Since we have to wait until all segments are processed, the output renderer will do nothing until closed.
     */
    @Override
    public void close() {

        // Flush the map.
        flushAllGenes();

        simpleTsvOutputRenderer.close();
    }

    /** Since this output renderer needs to see all of the input before it can render properly, it will not actually
     * write anything during this call.  The writing is done all at once during {@link GeneListOutputRenderer#close()}
     *
     * @param variant {@link VariantContext} to write to the file.  Must be a segment.
     * @param txToFuncotationMap {@link FuncotationMap} to add to the given {@code variant} on output.
     */
    @Override
    public void write(final VariantContext variant, final FuncotationMap txToFuncotationMap) {

        // Validate that we can work with the incoming data.
        validateAbleToWrite(variant, txToFuncotationMap);

        // For each gene in the genes field, add an entry in the map.
        // HACK: This is more difficult since we do not actually know what the genes et al field names actually are.
        final String genesFieldName = determineGencodeFieldName(txToFuncotationMap, FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY, gencodeGenesFieldPattern);
        final List<String> genes = Arrays.asList(StringUtils.split(txToFuncotationMap.getFieldValue(
                FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY,
                genesFieldName,
                txToFuncotationMap.getAlleles(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY).iterator().next())
                ,","));
        final String startGene = txToFuncotationMap.getFieldValue(
                FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY,
                determineGencodeFieldName(txToFuncotationMap, FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY, gencodeStartGeneFieldPattern),
                txToFuncotationMap.getAlleles(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY).iterator().next());
        final String startExon = txToFuncotationMap.getFieldValue(
                FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY,
                determineGencodeFieldName(txToFuncotationMap, FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY, gencodeStartExonFieldPattern),
                txToFuncotationMap.getAlleles(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY).iterator().next());
        final String endGene = txToFuncotationMap.getFieldValue(
                FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY,
                determineGencodeFieldName(txToFuncotationMap, FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY, gencodeEndGeneFieldPattern),
                txToFuncotationMap.getAlleles(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY).iterator().next());
        final String endExon = txToFuncotationMap.getFieldValue(
                FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY,
                determineGencodeFieldName(txToFuncotationMap, FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY, gencodeEndExonFieldPattern),
                txToFuncotationMap.getAlleles(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY).iterator().next());

        // For the next two steps make sure to copy the funcotation maps, since we may have to alter them downstream.
        // Add start gene and end gene to the map (with exon string).  Ignore the empty gene case, which means a segment
        //  breakpoint overlapped an IGR.
        Stream.of(Pair.of(startGene, startExon), Pair.of(endGene, endExon))
                .filter(p -> !StringUtils.isEmpty(p.getLeft()))
                .forEach(p -> geneExonToVariantFuncotationMap.put(p, Pair.of(variant, FuncotationMap.create(txToFuncotationMap))));

        // Add each gene to the internal map, making sure not to re-add the start and end gene.  Also, do not include
        //  blank genes
        genes.stream().filter(g -> !g.equals(startGene) && !g.equals(endGene))
                .filter(g -> !StringUtils.isEmpty(g))
                .forEach(g -> geneExonToVariantFuncotationMap.put(Pair.of(g, ""), Pair.of(variant, FuncotationMap.create(txToFuncotationMap))));
    }

    private void flushAllGenes() {

        // Create a dummy funcotation that has the gene and exon fields and add it to the funcotation map.
        for (final Pair<String, String> geneExon : geneExonToVariantFuncotationMap.keySet()) {
            final VariantContext segmentVariantContext = geneExonToVariantFuncotationMap.get(geneExon).getLeft();

            final FuncotationMap funcotationMap = geneExonToVariantFuncotationMap.get(geneExon).getRight();

            final Allele altAllele = funcotationMap.get(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY).get(0).getAltAllele();
            final Funcotation newFuncotation = createGeneExonFuncotation(geneExon.getLeft(), geneExon.getRight(), altAllele);

            funcotationMap.add(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY, newFuncotation);
            simpleTsvOutputRenderer.write(segmentVariantContext, funcotationMap);
        }
    }

    private Funcotation createGeneExonFuncotation(final String gene, final String exon, final Allele allele) {
        return TableFuncotation.create(
                Arrays.asList(GENE_FIELD_NAME, EXON_FIELD_NAME), Arrays.asList(gene, exon),
                allele, GENE_LIST_OUTPUT_RENDERER_DUMMY_NAME, null);
    }

    private static String determineGencodeFieldName(final FuncotationMap funcotationMap, final String txId, final Pattern pattern) {
        final Set<String> allPossibleFields = funcotationMap.getFieldNames(txId);
        final List<String> matches = allPossibleFields.stream().map(pattern::matcher).filter(Matcher::matches)
                .map(Matcher::group).collect(Collectors.toList());

        if (matches.size() != 1) {
            throw new UserException.BadInput(MISSING_GENCODE_FIELD_ERROR_MSG + ": " + matches.stream().collect(Collectors.joining(", ")));
        }

        return matches.get(0);
    }

    @VisibleForTesting
    SortedMap<Pair<String, String>, Pair<VariantContext, FuncotationMap>> getGeneExonToVariantFuncotationMap() {
        return geneExonToVariantFuncotationMap;
    }

    /**
     * Rules:
     *
     * - only work on segments.
     * - requires a "genes" funcotation field on each segment.  Additionally, this field must be a comma-delimited list of
     * gene names.  Also, requires, start/end exons, and start/end genes.  For a total of five fields.  This check is
     * not done in this method.  It will get checked in {@link GeneListOutputRenderer#write(VariantContext, FuncotationMap)}
     * - requires that the input FuncotationMap in {@link GeneListOutputRenderer#write(VariantContext, FuncotationMap)} only
     * have one transcript ID and that is {@link FuncotationMap#NO_TRANSCRIPT_AVAILABLE_KEY}.
     * - each segment can only have one alternate allele across all of the funcotations.
     *
     * Will convert illegal argument exceptions into user errors, since this means that the user probably gave the wrong
     *  type of input file.
     */
    @VisibleForTesting
    void validateAbleToWrite(final VariantContext variant, final FuncotationMap txToFuncotationMap) {
        try {
            Utils.validateArg(FuncotatorUtils.isSegmentVariantContext(variant, minNumBasesForValidSegment), ERR_STR_IS_THIS_SEG_FILE + INVALID_VARIANT_CONTEXT_ERR_MSG + ": " + variant);
            Utils.validateArg(txToFuncotationMap.getTranscriptList().size() == 1, ERR_STR_IS_THIS_SEG_FILE + INCORRECT_NUM_TRANSCRIPT_ERR_MSG + ": " + StringUtils.join(txToFuncotationMap.getTranscriptList(), ","));
            Utils.validateArg(txToFuncotationMap.getTranscriptList().get(0).equals(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY), ERR_STR_IS_THIS_SEG_FILE + INVALID_TRANSCRIPT_ID_ERR_MSG + "  (must be no transcript available dummy ID): " + txToFuncotationMap.getTranscriptList().get(0));
            Utils.validateArg(txToFuncotationMap.getAlleles(FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY).size() == 1, ERR_STR_IS_THIS_SEG_FILE + ONE_ALTERNATE_ALLELE_ERR_MSG);
        } catch (final IllegalArgumentException iae) {
            throw new UserException.BadInput(iae.getMessage());
        }
    }

    private final class PairedSimpleNaturalComparator<T extends CharSequence,U extends CharSequence>
            implements Comparator<Pair<T, U>> {

        private final Comparator<T> leftComparator = SimpleNaturalComparator.getInstance();
        private final Comparator<U> rightComparator = SimpleNaturalComparator.getInstance();

        PairedSimpleNaturalComparator() {}

        @Override
        public int compare(final Pair<T, U> o1, final Pair<T, U> o2) {
            final int compareLeft = leftComparator.compare(o1.getLeft(), o2.getLeft());
            if (compareLeft != 0) {
                return compareLeft;
            } else {
                return rightComparator.compare(o1.getRight(), o2.getRight());
            }
        }
    }
}
