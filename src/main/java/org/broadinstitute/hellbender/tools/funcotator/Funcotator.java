package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.mafOutput.MafOutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.transformers.VariantTransformer;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Funcotator (FUNCtional annOTATOR) analyzes given variants for their function (as retrieved from a set of data sources) and produces the analysis in a specified output file.
 *
 * <p>
 *     This tool is a functional annotation tool that allows a user to add annotations to called variants based on a set of data sources, each with its own matching criteria.
 *     Data sources are expected to be in folders that are specified as input arguments.  While multiple data source folders can be specified, <b>no two data sources can have the same name</b>.
 * </p>
 * <h3>Data Source Folders</h3>
 * <p>
 *     In each main data source folder, there should be sub-directories for each individual data source, with further sub-directories for a specific reference (i.e. <i>hg19</i> or <i>hg38</i>).
 *     In the reference-specific data source directory, there is a configuration file detailing information about the data source and how to match it to a variant.  This configuration file is required.
 * </p>
 * <p>
 *     An example of a data source directory is the following:
 *
 *     <pre>
 *         dataSourcesFolder/
 *              Data_Source_1/
 *                  hg19
 *                      data_source_1.config
 *                      data_source_1.data.file.one
 *                      data_source_1.data.file.two
 *                      data_source_1.data.file.three
 *                      ...
 *                   hg38
 *                      data_source_1.config
 *                      data_source_1.data.file.one
 *                      data_source_1.data.file.two
 *                      data_source_1.data.file.three
 *                      ...
 *              Data_Source_2/
 *                  hg19
 *                      data_source_2.config
 *                      data_source_2.data.file.one
 *                      data_source_2.data.file.two
 *                      data_source_2.data.file.three
 *                      ...
 *                   hg38
 *                      data_source_2.config
 *                      data_source_2.data.file.one
 *                      data_source_2.data.file.two
 *                      data_source_2.data.file.three
 *                      ...
 *               ...
 *     </pre>
 *     <b>A gzip of data source files is provided here: <a href="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/">ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/</a>.</b>
 * </p>
 * <h4>User-Defined Data Sources</h4>
 * <p>
 *     Users can define their own data sources by creating a new correctly-formatted data source sub-directory in the main data sources folder.  In this sub-directory, the user must create an additional folder for the reference for which the data source is valid.  If the data source is valid for multiple references, then multiple reference folders should be created.
 *     Inside each reference folder, the user should place the file(s) containing the data for the data source.  Additionally the user <b>must</b> create a configuration file containing metadata about the data source.
 * </p>
 * <p>
 *     There are several formats allowed for data sources, however the two most useful are arbitrarily separated value (XSV) files, such as comma-separated value (CSV), tab-separated value (TSV).  These files contain a table of data that can be matched to a variant by <i>gene name</i>, <i>transcript ID</i>, or <i>genome position</i>.
 *     In the case of <i>gene name</i> and <i>transcript ID</i>, one column must contain the <i>gene name</i> or <i>transcript ID</i> for each row's data.
 *     <ul>
 *         <li>For <i>gene name</i>, when a variant is annotated with a gene name that <i>exactly matches</i> an entry in the gene name column for a row, that row's other fields will be added as annotations to the variant.</li>
 *         <li>For <i>transcript ID</i>, when a variant is annotated with a transcript ID that <i>exactly matches</i> an entry in the transcript ID column for a row, that row's other fields will be added as annotations to the variant.</li>
 *         <li>For <i>genome position</i>, one column must contain the contig ID, another column must contain the start position (1-based, inclusive), and a column must contain the stop position (1-based, inclusive).  The start and stop columns may be the same column.  When a variant is annotated with a genome position that <i>overlaps</i> an entry in the three genome position columns for a row, that row's other fields will be added as annotations to the variant.</li>
 *     </ul>
 * </p>
 * <h4>Configuration File Format</h4>
 * <p>
 *     The configuration file is a standard Java properties-style configuration file with key-value pairs.  This file name <b>must end in .config</b>.
 * </p>
 * <p>
 *     The following is an example of a genome position XSV configuration file (for the ORegAnno data source):
 *     <pre>
 *         name = Oreganno
 *         version = 20160119
 *         src_file = oreganno.tsv
 *         origin_location = http://www.oreganno.org/dump/ORegAnno_Combined_2016.01.19.tsv
 *         preprocessing_script = getOreganno.py
 *
 *         # Supported types:
 *         # simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID
 *         # locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location
 *         # gencode      -- Custom datasource class for GENCODE
 *         # cosmic       -- Custom datasource class for COSMIC
 *         type = locatableXSV
 *
 *         # Required field for GENCODE files.
 *         # Path to the FASTA file from which to load the sequences for GENCODE transcripts:
 *         gencode_fasta_path =
 *
 *         # Required field for simpleXSV files.
 *         # Valid values:
 *         #     GENE_NAME
 *         #     TRANSCRIPT_ID
 *         xsv_key =
 *
 *         # Required field for simpleXSV files.
 *         # The 0-based index of the column containing the key on which to match
 *         xsv_key_column =
 *
 *         # Required field for simpleXSV AND locatableXSV files.
 *         # The delimiter by which to split the XSV file into columns.
 *         xsv_delimiter = \t
 *
 *         # Required field for simpleXSV files.
 *         # Whether to permissively match the number of columns in the header and data rows
 *         # Valid values:
 *         #     true
 *         #     false
 *         xsv_permissive_cols = true
 *
 *         # Required field for locatableXSV files.
 *         # The 0-based index of the column containing the contig for each row
 *         contig_column = 1
 *
 *         # Required field for locatableXSV files.
 *         # The 0-based index of the column containing the start position for each row
 *         start_column = 2
 *
 *         # Required field for locatableXSV files.
 *         # The 0-based index of the column containing the end position for each row
 *         end_column = 3
 *     </pre>
 * </p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>A reference genome sequence.</li>
 *     <li>The version of the reference genome sequence being used (either <i>hg19</i> or <i>hg38</i>).</li>
 *     <li>A VCF of variant calls to annotate.</li>
 *     <li>The path to a folder of data sources formatted for use by Funcotator.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A VCF containing all variants from the input file with added annotation columns corresponding to annotations from each data source that matched a given variant according to that data source's matching criteria.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   ./gatk Funcotator \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -O output.vcf \
 *   --data-sources-path dataSourcesFolder/ \
 *   --ref-version hg19
 * </pre>
 *
 * <h3>Notes</h3>
 * <ul>
 *     <li>This is a beta tool, and as such may generate errors or warnings.</li>
 *     <li>This tool is the GATK analog of <a href="http://portals.broadinstitute.org/oncotator/">Oncotator</a>.</li>
 * </ul>
 *
 * <h3>Known Issues</h3>
 * <p>A complete list of known open issues can be found on <a href="https://github.com/broadinstitute/gatk/issues?q=is%3Aopen+is%3Aissue+label%3AFuncotator">the GATK github entry for funcotator here.</a></p>
 *
 * <h4>Notable Issues as of 2018 Jan 3</h4>
 * <ul>
 *     <li>Only supports VCF for inputs and outputs (<a href="https://github.com/broadinstitute/gatk/issues/3922">Issue 3922</a>).</li>
 *     <li>Only supports a single GENCODE data source (<a href="https://github.com/broadinstitute/gatk/issues/3956">Issue 3956</a>).</li>
 *     <li>Non-GENCODE annotations on multiallelic variants are only rendered properly for the last allele (<a href="https://github.com/broadinstitute/gatk/issues/3896">Issue 3896</a>).</li>
 *     <li>The "other transcripts" annotation is missing for IGR variants (<a href="https://github.com/broadinstitute/gatk/issues/3849">Issue 3849</a>).</li>
 *     <li>Only the "Format" field of input VCF files is preserved in the output VCF file (<a href="https://github.com/broadinstitute/gatk/issues/3895">Issue 3895</a>).</li>
 *     <li>Codon Change and Protein Change for indels spanning splice sites are not properly rendered (<a href="https://github.com/broadinstitute/gatk/issues/3749">Issue 3749</a>).</li>
 *     <li>Does not support Structural Variants (<a href="https://github.com/broadinstitute/gatk/issues/4083">Issue 4083</a>).</li>
 * </ul>
 *
 */
@CommandLineProgramProperties(
        summary = "Create functional annotations on given variants cross-referenced by a given set of data sources.\n" +
                "A GATK functional annotation tool (similar functionality to Oncotator).",
        oneLineSummary = "Functional Annotator",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class Funcotator extends VariantWalker {
    private static final Logger logger = LogManager.getLogger(Funcotator.class);

    /**
     * The current version of {@link Funcotator}.
     */
    public static final String VERSION = "0.0.4";

    //==================================================================================================================
    // Arguments:

    @ArgumentCollection
    private FuncotatorArgumentCollection funcotatorArgs = new FuncotatorArgumentCollection();

    //==================================================================================================================

    private OutputRenderer outputRenderer;

    private FuncotatorEngine funcotatorEngine;

    /**
     * Whether the input variant contigs must be converted to hg19.
     * This is only the case when the input reference is b37 AND when
     * the reference version is hg19 (i.e. {@link FuncotatorArgumentCollection#referenceVersion} == {@link FuncotatorArgumentDefinitions#HG19_REFERENCE_VERSION_STRING}).
     */
    private boolean mustConvertInputContigsToHg19 = false;


    //==================================================================================================================

    /**
     * @return The {@link Funcotator}-specific arguments used to instantiate this {@link Funcotator} instance.
     */
    public FuncotatorArgumentCollection getArguments() {
        return funcotatorArgs;
    }

    @Override
    protected String getVersion() {
        return super.getVersion() + "-" + VERSION;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {

        if (seqValidationArguments.performSequenceDictionaryValidation()) {
            // Ensure that the reference dictionary is a superset of the variant dictionary:
            checkReferenceDictionaryIsSupersetOfVariantDictionary();
        }

        // Next set up our transcript list:
        final Set<String> finalUserTranscriptIdSet = processTranscriptList(funcotatorArgs.userTranscriptIdSet);

        final LinkedHashMap<String, String> annotationDefaultsMap = splitAnnotationArgsIntoMap(funcotatorArgs.annotationDefaults);
        final LinkedHashMap<String, String> annotationOverridesMap = splitAnnotationArgsIntoMap(funcotatorArgs.annotationOverrides);

        // Initialize all of our data sources:
        // Sort data sources to make them process in the same order each time:
        funcotatorArgs.dataSourceDirectories.sort(Comparator.naturalOrder());
        final Map<Path, Properties> configData = DataSourceUtils.getAndValidateDataSourcesFromPaths(funcotatorArgs.referenceVersion, funcotatorArgs.dataSourceDirectories);

        // Create the data sources from the input:
        // This will also create and register the FeatureInputs (created by the Data Sources)
        // with the GATK Engine, so we do not have to plumb them in after the fact.
        final List<DataSourceFuncotationFactory> dataSourceFuncotationFactories = DataSourceUtils.createDataSourceFuncotationFactoriesForDataSources(
                configData,
                annotationOverridesMap,
                funcotatorArgs.transcriptSelectionMode,
                finalUserTranscriptIdSet,
                this,
                funcotatorArgs.lookaheadFeatureCachingInBp
        );

        // Create our engine to do our work and drive this Funcotation train!
        funcotatorEngine = new FuncotatorEngine(
                VcfFuncotationMetadata.create(
                    new ArrayList<>(getHeaderForVariants().getInfoHeaderLines())
                ),
                dataSourceFuncotationFactories
        );

        // Determine which annotations are accounted for (by the funcotation factories) and which are not.
        final LinkedHashMap<String, String> unaccountedForDefaultAnnotations = getUnaccountedForAnnotations( funcotatorEngine.getFuncotationFactories(), annotationDefaultsMap );
        final LinkedHashMap<String, String> unaccountedForOverrideAnnotations = getUnaccountedForAnnotations( funcotatorEngine.getFuncotationFactories(), annotationOverridesMap );

        // Set up our output renderer:
        switch (funcotatorArgs.outputFormatType) {
            case MAF:
                outputRenderer = new MafOutputRenderer(funcotatorArgs.outputFile.toPath(),
                        funcotatorEngine.getFuncotationFactories(),
                        getHeaderForVariants(),
                        unaccountedForDefaultAnnotations,
                        unaccountedForOverrideAnnotations,
                        getDefaultToolVCFHeaderLines().stream().map(Object::toString).collect(Collectors.toCollection(LinkedHashSet::new)),
                        funcotatorArgs.referenceVersion);
                break;
            case VCF:
                outputRenderer = new VcfOutputRenderer(createVCFWriter(funcotatorArgs.outputFile),
                        funcotatorEngine.getFuncotationFactories(),
                        getHeaderForVariants(),
                        unaccountedForDefaultAnnotations,
                        unaccountedForOverrideAnnotations,
                        getDefaultToolVCFHeaderLines());
                break;
            default:
                throw new GATKException("Unsupported output format type specified: " + funcotatorArgs.outputFormatType.toString());
        }
        logger.info("Creating a " + outputFormatType + " file for output: " + outputFile.toURI());

        // Check for reference version (in)compatibility:
        determineReferenceAndDatasourceCompatibility();
    }

    /**
     * Checks to see that the given reference's sequence dictionary is a
     * superset of the given variant file's dictionary.
     *
     * This is a more strict check than the one found in {@link GATKTool#validateSequenceDictionaries()}.
     */
    private void checkReferenceDictionaryIsSupersetOfVariantDictionary() {

        final SAMSequenceDictionary referenceDictionary = getReferenceDictionary();
        final SAMSequenceDictionary variantDictionary = getSequenceDictionaryForDrivingVariants();

        if ( referenceDictionary == null ) {
            throw new UserException.BadInput("Reference fasta sequence dictionary is null!");
        }

        if ( variantDictionary == null ) {
            throw new UserException.BadInput("Funcotator by default requires that the variant input have a sequence dictionary in its header. To disable this safety check, use argument --" + StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME);
        }

        SequenceDictionaryUtils.validateDictionaries(
                "Reference", getReferenceDictionary(),
                "Driving Variants", getSequenceDictionaryForDrivingVariants(),
                true,
                false
                );
    }

    private void determineReferenceAndDatasourceCompatibility() {
        if ( funcotatorArgs.forceB37ToHg19ContigNameConversion ||
                ( funcotatorArgs.referenceVersion.equals(FuncotatorArgumentDefinitions.HG19_REFERENCE_VERSION_STRING) &&
                  FuncotatorUtils.isSequenceDictionaryUsingB37Reference(getSequenceDictionaryForDrivingVariants()) )) {

            // NOTE AND WARNING:
            // hg19 is from ucsc. b37 is from the genome reference consortium.
            // ucsc decided the grc version had some bad data in it, so they blocked out some of the bases, aka "masked" them
            // so the lengths of the contigs are the same, the bases are just _slightly_ different.
            // ALSO, the contig naming convention is different between hg19 and hg38:
            //      hg19 uses contigs of the form "chr1"
            //      b37 uses contigs of the form  "1"
            // This naming convention difference causes a LOT of issues and was a bad idea.

            logger.warn("WARNING: You are using B37 as a reference.  " +
                    "Funcotator will convert your variants to GRCh37, and this will be fine in the vast majority of cases.  " +
                    "There MAY be some errors (e.g. in the Y chromosome, but possibly in other places as well) due to changes between the two references.");

            mustConvertInputContigsToHg19 = true;
        }
    }

    private VariantContext getCorrectVariantContextForReference(final VariantContext variant) {
        if ( mustConvertInputContigsToHg19 ) {
            final VariantContextBuilder vcb = new VariantContextBuilder(variant);
            vcb.chr(FuncotatorUtils.convertB37ContigToHg19Contig(variant.getContig()));
            return vcb.make();
        }
        else {
            return variant;
        }
    }

    @Override
    protected VariantFilter makeVariantFilter() {
        return  variant -> {
            // Ignore variants that have been filtered if the user requests it:
            if ( funcotatorArgs.removeFilteredVariants && variant.isFiltered() ) {
                return false;
            }
            return true;
        };
    }

    @Override
    public VariantTransformer makePostVariantFilterTransformer(){
        return variantContext -> getCorrectVariantContextForReference(variantContext);
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        final ReferenceContext correctReferenceContext;

        // Check to see if we need to revert the ReferenceContext's interval to the original variant interval
        // (This would only happen in the case where we were given b37 variants with hg19 data sources):
        if ( mustConvertInputContigsToHg19 ) {

            // Convert our contig back to B37 here so it matches the variant:
            final SimpleInterval interval = new SimpleInterval(
                    FuncotatorUtils.convertHG19ContigToB37Contig(variant.getContig()), variant.getStart(), variant.getEnd()
            );

            correctReferenceContext = new ReferenceContext(referenceContext, interval);
        }
        else {
            correctReferenceContext = referenceContext;
        }

        // Place the variant on our queue to be funcotated:
        enqueueAndHandleVariant(variant, correctReferenceContext, featureContext);
    }

    @Override
    public Object onTraversalSuccess() {
        return true;
    }

    @Override
    public void closeTool() {
        if ( funcotatorEngine != null) {
            funcotatorEngine.close();
        }

        if ( outputRenderer != null ) {
            outputRenderer.close();
        }
    }

    //==================================================================================================================

    /**
     * Creates a {@link LinkedHashMap} of annotations in the given {@code annotationMap} that do not occur in the given {@code dataSourceFactories}.
     * @param dataSourceFactories {@link List} of {@link DataSourceFuncotationFactory} to check for whether each annotation in the {@code annotationMap} is handled.
     * @param annotationMap {@link Map} (of ANNOTATION_NAME : ANNOTATION_VALUE) to check
     * @return A {@link LinkedHashMap} of annotations in the given {@code annotationMap} that do not occur in the given {@code dataSourceFactories}.
     */
    private LinkedHashMap<String, String> getUnaccountedForAnnotations( final List<DataSourceFuncotationFactory> dataSourceFactories,
                                                                        final Map<String, String> annotationMap ) {
        final LinkedHashMap<String, String> outAnnotations = new LinkedHashMap<>();

        // Check each field in each factory:
        for ( final String field : annotationMap.keySet() ) {
            boolean accountedFor = false;
            for ( final DataSourceFuncotationFactory funcotationFactory : dataSourceFactories ) {

                if ( funcotationFactory.getSupportedFuncotationFields().contains(field) ) {
                    accountedFor = true;
                    break;
                }
            }
            if ( !accountedFor ) {
                outAnnotations.put(field, annotationMap.get(field));
            }
        }

        return outAnnotations;
    }

    /**
     * Creates an annotation on the given {@code variant} or enqueues it to be processed during a later call to this method.
     * @param variant {@link VariantContext} to annotate.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.
     * @param featureContext {@link FeatureContext} corresponding to the given {@code variant}.
     */
    private void enqueueAndHandleVariant(final VariantContext variant, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        final FuncotationMap funcotationMap = funcotatorEngine.createFuncotationMapForVariant(variant, referenceContext, featureContext);

        // At this point there is only one transcript ID in the funcotation map if canonical or best effect are selected
        outputRenderer.write(variant, funcotationMap);
    }

    /**
     * Split each element of the given {@link List} into a key and value.
     * Assumes each element of the given {@link List} is formatted as follows:
     *     KEY:VALUE
     * @param annotationArgs {@link List} of strings formatted KEY:VALUE to turn into a {@link Map}.
     * @return A {@link LinkedHashMap} of KEY:VALUE pairs corresponding to entries in the given list.
     */
    private LinkedHashMap<String, String> splitAnnotationArgsIntoMap( final List<String> annotationArgs ) {

        final LinkedHashMap<String, String> annotationMap = new LinkedHashMap<>();

        for ( final String s : annotationArgs ) {
            final List<String> keyVal = ParsingUtils.split(s, FuncotatorArgumentDefinitions.MAP_NAME_VALUE_DELIMITER);
            if ( keyVal.size() != 2) {
                throw new UserException.BadInput( "Argument annotation incorrectly formatted: " + s );
            }

            annotationMap.put( keyVal.get(0), keyVal.get(1) );
        }

        return annotationMap;
    }

    /**
     * Processes the given {@link Set} into a list of transcript IDs.
     * This is necessary because the command-line input argument is overloaded to be either a file containing transcript
     * IDs (1 per line) OR as a list of transcript IDs.
     * @param rawTranscriptSet {@link Set} of {@link String}s from which to create a list of Transcript IDs.  If of size 1, will try to open as a file.
     * @return A {@link Set} of {@link String} contianing Transcript IDs in which the user is interested.
     */
    private Set<String> processTranscriptList(final Set<String> rawTranscriptSet) {
        if ( rawTranscriptSet.size() == 1 ) {
            final String filePathString = rawTranscriptSet.iterator().next();
            try ( final BufferedReader bufferedReader = Files.newBufferedReader(IOUtils.getPath(filePathString)) ) {
                logger.info("Opened transcript file: " + filePathString);

                // Create a place to put our output:
                final Set<String> transcriptIdSet = new HashSet<>();

                String line = bufferedReader.readLine();
                while ( line != null ) {
                    logger.info("    Adding transcript ID to transcript set: " + line);
                    transcriptIdSet.add(line);
                    line = bufferedReader.readLine();
                }
                logger.info("Transcript parsing complete.");

                return transcriptIdSet;
            }
            catch ( final IOException ex ) {
                logger.warn("Could not open transcript selection list as a file.  Using it as a singleton list of transcript IDs: [" + filePathString + "]");
                return rawTranscriptSet;
            }
        }
        else {
            return rawTranscriptSet;
        }
    }
}
