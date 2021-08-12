package org.broadinstitute.hellbender.tools.funcotator;


import org.apache.http.HttpEntity;
import org.apache.http.HttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClientBuilder;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;

import java.nio.file.Files;
import java.nio.file.Path;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.time.LocalDate;
import java.time.Month;


/**
 * Class to copy a file using {@link CloseableHttpClient}.
 * Operates using paths.
 *
 * Created by Hailey on 8/5/21.
 */
public class FuncotatorDataSourceBundlerHttpClient {

    //==================================================================================================================
    // Standard logger:
    private final static Logger logger = LogManager.getLogger(FuncotatorDataSourceBundlerHttpClient.class);

    //==================================================================================================================
    // Public Static Members:
    public static final String ENSEMBL_TEMPLATE_CONFIG      = "src/test/resources/large/funcotator/funcotator_dataSources/ensembl.template.config";
    public static final String ENSEMBL_CONFIG_NAME          = "ensembl.config";
    public static final String MANIFEST_FILE_NAME           = "MANIFEST.txt";
    public static final String TEMPLATE_CONFIG_FILE_NAME    = "template.config";
    public static final String README_FILE_NAME             = "README.txt";

    //==================================================================================================================
    // Private Static Members:
    private static final int BUFFER_SIZE_BYTES    = 1024 * 1024;

    //==================================================================================================================
    // Private Members:

    // Data variables:
    protected String dsOrganism;
    protected String fileName;
    protected String fastaFileName;
    protected String dsURL;
    protected Path dsPath;
    protected Path dsUnzipPath;
    protected String dsFastaURL;
    protected Path dsFastaPath;
    protected Path dsFastaUnzipPath;
    protected Path indexFilePath;
    protected String baseURL;
    protected String speciesName;
    protected String baseFastaURL;
    protected Path outputDestination;
    protected Path outputUnzippedDest;
    protected Path outputFastaDest;
    protected Path outputFastaUnzipDest;
    protected Path outputIndexDest;
    protected Path configFilePath;
    protected Path metadataFilePath;
    protected boolean overwriteOutputFile;
    protected boolean extractAfterDownload;
    protected String dsGtfReadMeURL;
    protected String dsFastaReadMeURL;
    protected Path dsGtfReadMePath;
    protected Path dsFastaReadMePath;

    // Copy buffer:
    public static byte copyBuffer[] = new byte[BUFFER_SIZE_BYTES];

    //==================================================================================================================
    // Constructors:

    /**
     * {@link FuncotatorDataSourceBundlerHttpClient}
     * This internal constructor is to be used by the class itself.
     * @param dsOrganism The {@link String} representing the chosen organism.
     * @param speciesName The {@link String} representing the chosen division.
     * @param baseURL The {@link String} representing the base url for the chosen organism.
     */
    protected FuncotatorDataSourceBundlerHttpClient(final String dsOrganism, final String speciesName, final String baseURL, final String baseFastaURL) {
        this.dsOrganism             = dsOrganism;
        this.speciesName            = speciesName;
        this.baseURL                = baseURL;
        this.baseFastaURL           = baseFastaURL;
        this.fileName               = FuncotatorDataSourceBundlerUtils.buildMapGetFileName(this.dsOrganism, this.speciesName, false);
        this.fastaFileName          = FuncotatorDataSourceBundlerUtils.buildMapGetFileName(this.dsOrganism, this.speciesName, true);
        this.dsURL                  = setURL(this.baseURL, this.speciesName, this.fileName);
        this.dsPath                 = setPath(this.speciesName, this.fileName);
        this.dsUnzipPath            = setUnzipPath(this.speciesName, this.fileName);
        this.dsFastaURL             = setFastaUrl(this.baseFastaURL, this.speciesName, this.fastaFileName);
        this.dsFastaPath            = setFastaPath(this.speciesName, this.fastaFileName);
        this.dsFastaUnzipPath       = setFastaUnzipPath(this.speciesName, this.fastaFileName);
        this.indexFilePath          = setIndexFilePath(this.speciesName, this.fileName);
        this.outputDestination      = this.dsPath.toAbsolutePath();
        this.outputUnzippedDest     = this.dsUnzipPath.toAbsolutePath();
        this.outputFastaDest        = this.dsFastaPath.toAbsolutePath();
        this.outputFastaUnzipDest   = this.dsFastaUnzipPath.toAbsolutePath();
        this.outputIndexDest        = this.indexFilePath.toAbsolutePath();
        this.configFilePath         = setConfigFilePath(this.speciesName);
        this.metadataFilePath       = setMetadataFilePath(this.speciesName);
        this.dsGtfReadMeURL         = setReadMeURL(this.baseURL, this.speciesName);
        this.dsFastaReadMeURL       = setReadMeURL(this.baseFastaURL, this.speciesName);
        this.dsGtfReadMePath        = setGtfReadMePath(this.speciesName, this.fileName);
        this.dsFastaReadMePath      = setFastaReadMePath(this.speciesName, this.fileName);
    }

    //==================================================================================================================
    // Static Methods:

    /**
     * Create an {@link FuncotatorDataSourceBundlerHttpClient}.
     * @param dsOrganism The {@link String} representing the chosen organism to bundle data sources for.
     * @param speciesName The {@link String} representing the specific division of the organism to bundle data sources for.
     * @param baseURL The {@link String} base url for the specific organism which was chosen.
     * @param baseFastaURL The {@link String} base url for the fasta file for the specific organism which was chosen.
     * @return An {@link FuncotatorDataSourceBundlerHttpClient} initialized to copy the file for the species name to the user's computer.
     */
    public static FuncotatorDataSourceBundlerHttpClient create(final String dsOrganism, final String speciesName, final String baseURL, final String baseFastaURL) {
        return new FuncotatorDataSourceBundlerHttpClient(dsOrganism, speciesName, baseURL, baseFastaURL);
    }

    /**
     * Download the data source given by {@code dsURL} to {@code outputDestination}.
     * @param dsURL The {@link String} representing the url at which the file to download is found.
     * @param outputDestination The {@link Path} representing the output path for the copied file to be put.
     */
    public static void downloadDataSources(final String dsURL, final Path outputDestination) {

        // Creating CloseableHttpClient object to access the webpage and retrieve the file:
        CloseableHttpClient client = HttpClientBuilder.create().build();

        // Creating an HttpGet object to send the request to the server:
        HttpGet request = new HttpGet(dsURL);

        try {
            // Using an HttpResponse class object to catch the response from the server
            HttpResponse response = client.execute(request);
            // The data sent by the server is obtained in this getEntity() function:
            HttpEntity entity = response.getEntity();

            // Extracting the data from the entity object:
            try( final InputStream inputStream = entity.getContent();
                 final OutputStream outputStream = Files.newOutputStream(outputDestination) ) {

                // Perform the copy:
                while (true) {

                    // Read from our input:
                    final int bytesRead = inputStream.read(copyBuffer);
                    if (bytesRead == -1) {
                        break;
                    }

                    // Write to our output:
                    outputStream.write(copyBuffer, 0, bytesRead);
                }
            }
            catch (final IOException ex) {
                throw new UserException("Could not copy file: " + dsURL + " -> " + outputDestination.toUri().toString(), ex);
            }
        }
        catch (final IOException ex) {
            throw new UserException("Could not obtain data from "+ dsURL, ex);
        }
    }

    /**
     * Build an index file for the gtf data source file.
     * @param gtfFilePath The {@link Path} representing the path to the gtf data source file we have downloaded.
     * @param idxFilePath The {@link Path} representing the path where we want our indexed file to be.
     */
    public static void buildIndexFile(Path gtfFilePath, Path idxFilePath) {
        IndexFeatureFile indexer = new IndexFeatureFile();
        indexer.indexGTF(gtfFilePath.toAbsolutePath(), idxFilePath.toAbsolutePath());
    }

    /**
     * Build an index file for the fasta data source file.
     * @param dsFastaUnzipPath The {@link Path} representing the path to the unzipped fasta file of the data source we have downloaded.
     */
    public static void buildFastaIndexFile(Path dsFastaUnzipPath) {
        FastaReferenceWriterBuilder fastaIndexer = new FastaReferenceWriterBuilder();
        fastaIndexer.setFastaFile(dsFastaUnzipPath);
        fastaIndexer.setMakeFaiOutput(true);
        try {
            fastaIndexer.build();
        } catch (IOException e) {
            throw new IllegalArgumentException("Error: ", e);
        }
    }

    /**
     * Build a config file for the data source we have downloaded.
     * @param bundler The {@link FuncotatorDataSourceBundlerHttpClient} which holds all of the variables we will need.
     */
    public static void buildConfigFile(FuncotatorDataSourceBundlerHttpClient bundler) {
        try {
            FileWriter writer = new FileWriter(bundler.configFilePath.toAbsolutePath().toString());
            BufferedWriter buffer = new BufferedWriter(writer);
            buffer.write(
                    "name = Ensembl\n" +
                            "version = 104\n" +
                            "src_file = " + bundler.fileName + "\n" +
                            "origin_location = " + bundler.dsURL + " \n" +
                            "preprocessing_script = FuncotatorDataSourceBundler \n" +
                            "\n" +
                            "# Supported types:\n" +
                            "# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID\n" +
                            "# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location\n" +
                            "# gencode      -- Custom datasource class for GENCODE\n" +
                            "#\tcosmic       -- Custom datasource class for COSMIC\n" +
                            "type = gencode\n" +
                            "\n" +
                            "# Required field for GENCODE files.\n" +
                            "# Path to the FASTA file from which to load the sequences for GENCODE transcripts:\n" +
                            "gencode_fasta_path = " + bundler.dsFastaURL + "\n" +
                            "\n" +
                            "# Required field for simpleXSV files.\n" +
                            "# Valid values:\n" +
                            "#     GENE_NAME\n" +
                            "#     TRANSCRIPT_ID\n" +
                            "xsv_key = GENE_NAME\n" +
                            "\n" +
                            "# Required field for simpleXSV files.\n" +
                            "# The 0-based index of the column containing the key on which to match\n" +
                            "xsv_key_column = 0\n" +
                            "\n" +
                            "# Required field for simpleXSV AND locatableXSV files.\n" +
                            "# The delimiter by which to split the XSV file into columns.\n" +
                            "xsv_delimiter = ,\n" +
                            "\n" +
                            "# Required field for simpleXSV files.\n" +
                            "# Whether to permissively match the number of columns in the header and data rows\n" +
                            "# Valid values:\n" +
                            "#     true\n" +
                            "#     false\n" +
                            "xsv_permissive_cols = true\n" +
                            "\n" +
                            "# Required field for locatableXSV files.\n" +
                            "# The 0-based index of the column containing the contig for each row\n" +
                            "contig_column =\n" +
                            "\n" +
                            "# Required field for locatableXSV files.\n" +
                            "# The 0-based index of the column containing the start position for each row\n" +
                            "start_column =\n" +
                            "\n" +
                            "# Required field for locatableXSV files.\n" +
                            "# The 0-based index of the column containing the end position for each row\n" +
                            "end_column =\n"
            );
            buffer.close();
        } catch (IOException e) {
            throw new UserException("Error. Unable to build file: " + bundler.configFilePath);
        }
    }

    /**
     * Build the template config file in the correct folder.
     * @param bundler The {@link FuncotatorDataSourceBundlerHttpClient} which holds all of the variables we will need.
     */
    public static void buildTemplateConfigFile(FuncotatorDataSourceBundlerHttpClient bundler) {
        try {
            FileWriter writer = new FileWriter(bundler.metadataFilePath.toAbsolutePath() + "/" + TEMPLATE_CONFIG_FILE_NAME);
            BufferedWriter buffer = new BufferedWriter(writer);
            buffer.write(
                    "name = Achilles\n" +
                            "version = 110303\n" +
                            "src_file = achilles_lineage_results.import.txt\n" +
                            "origin_location = UNKNOWN\n" +
                            "preprocessing_script =\n" +
                            "\n" +
                            "# Supported types:\n" +
                            "# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID\n" +
                            "# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location\n" +
                            "# gencode      -- Custom datasource class for GENCODE\n" +
                            "#\tcosmic       -- Custom datasource class for COSMIC\n" +
                            "type = simpleXSV\n" +
                            "\n" +
                            "# Required field for GENCODE files.\n" +
                            "# Path to the FASTA file from which to load the sequences for GENCODE transcripts:\n" +
                            "gencode_fasta_path =\n" +
                            "\n" +
                            "# Required field for simpleXSV files.\n" +
                            "# Valid values:\n" +
                            "#     GENE_NAME\n" +
                            "#     TRANSCRIPT_ID\n" +
                            "xsv_key = GENE_NAME\n" +
                            "\n" +
                            "# Required field for simpleXSV files.\n" +
                            "# The 0-based index of the column containing the key on which to match\n" +
                            "xsv_key_column = 0\n" +
                            "\n" +
                            "# Required field for simpleXSV AND locatableXSV files.\n" +
                            "# The delimiter by which to split the XSV file into columns.\n" +
                            "xsv_delimiter = ,\n" +
                            "\n" +
                            "# Required field for simpleXSV files.\n" +
                            "# Whether to permissively match the number of columns in the header and data rows\n" +
                            "# Valid values:\n" +
                            "#     true\n" +
                            "#     false\n" +
                            "xsv_permissive_cols = true\n" +
                            "\n" +
                            "# Required field for locatableXSV files.\n" +
                            "# The 0-based index of the column containing the contig for each row\n" +
                            "contig_column =\n" +
                            "\n" +
                            "# Required field for locatableXSV files.\n" +
                            "# The 0-based index of the column containing the start position for each row\n" +
                            "start_column =\n" +
                            "\n" +
                            "# Required field for locatableXSV files.\n" +
                            "# The 0-based index of the column containing the end position for each row\n" +
                            "end_column ="
            );
        } catch (IOException e) {
            throw new UserException("Error. Unable to make template config file in location: " + bundler.metadataFilePath + "/" + TEMPLATE_CONFIG_FILE_NAME);
        }
    }

    /**
     * Build a ReadMe file in the correct folder.
     * @param bundler The {@link FuncotatorDataSourceBundlerHttpClient} which holds all of the variables we will need.
     */
    public static void buildReadMeFile(FuncotatorDataSourceBundlerHttpClient bundler) {
        try {
            FileWriter writer = new FileWriter(bundler.metadataFilePath.toAbsolutePath().toString() + "/" + README_FILE_NAME);
            BufferedWriter buffer = new BufferedWriter(writer);
            buffer.write(
                    "################################################################################\n" +
                            "# Funcotator Data Sources Bundler Package README\n" +
                            "################################################################################\n" +
                            "\n" +
                            "+---------------------------------------------+ \n" +
                            "| Data Source Version Information             |\n" +
                            "+---------------------------------------------+ \n" +
                            "\n" +
                            "Version:          0.0." + bundler.getDate() + "\n" +
                            "Use Case:         species name\n" +
                            "Source:           ./gatk -bundler.dsOrganism -species-name bundler.speciesName \n" +
                            "Alternate Source: ./gatk -bundler.dsOrganism -species-name bundler.speciesName \n" +
                            "\n" +
                            "################################################################################\n" +
                            "\n" +
                            "+---------------------------------------------+ \n" +
                            "| README                                      | \n" +
                            "+---------------------------------------------+ \n" +
                            "\n" +
                            "This is a collection of data sources to be used in conjunction with Funcotator\n" +
                            "to annotate data samples for a variety of species. \n" +
                            "\n" +
                            "This folder is a top-level Data Sources Folder for The Broad Institute's \n" +
                            "Funcotator Data Source Bundler tool.  When running Funcotator, pass the path to this directory in\n" +
                            "as a command-line argument:\n" +
                            "\n" +
                            "   ./gatk Funcotator --data-sources-path PATH/TO/THIS/FOLDER ...\n" +
                            "\n" +
                            "For more information on Funcotator, see the GATK development github site:\n" +
                            "\n" +
                            "\thttps://github.com/broadinstitute/gatk\n" +
                            "\n" +
                            "################################################################################\n" +
                            "\n" +
                            "+---------------------------------------------+ \n" +
                            "| Data Sources                                |\n" +
                            "+---------------------------------------------+ \n" +
                            "\n" +
                            "Using this Data Sources Folder will enable the following data sources:\n" +
                            "\n" +
                            " acmg_lof\n" +
                            "--------------------\n" +
                            " Custom list of genes found with (or without) LoF mechanisms in earlier studies from the Laboratory for Molecular Medicine (LMM).  \n" +
                            " This datasource provides information whether a given gene was found in that list and whether there were any LoF variants found by the LMM.\n" +
                            "\n" +
                            " acmg_rec\n" +
                            "--------------------\n" +
                            " The ACMG59 recommendation list for clinical intervention.  See https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/ for more information.\n" +
                            "\n" +
                            " clinvar\n" +
                            "--------------------\n" +
                            " ClinVar is a freely accessible, public archive of reports of the relationships among human variations and phenotypes, with supporting evidence.  \n" +
                            " See https://www.ncbi.nlm.nih.gov/clinvar/intro/\n" +
                            "\n" +
                            " exac\n" +
                            "--------------------\n" +
                            " Variants from the Exome Aggregation Consortium.  See http://exac.broadinstitute.org/\n" +
                            "\n" +
                            " gencode\n" +
                            "--------------------\n" +
                            " The GENCODE project produces high quality reference gene annotation and experimental validation for human and mouse genomes.\n" +
                            "\n" +
                            " lmm_known\n" +
                            "--------------------\n" +
                            " Variants flagged as possibly having pathogenic effects found by the Laboratory for Molecular Medicine (LMM).\n" +
                            "\n" +
                            "  ensembl\n" +
                            "--------------------\n" +
                            "  The ENSEMBL Project produces high quality reference gene annotation and experimental validation for over 50,000 genomes. \n"
            );
            buffer.close();
        } catch (IOException e) {
            throw new UserException("Error. Unable to make ReadMe file in location: " + bundler.metadataFilePath + "/" + README_FILE_NAME);
        }
    }

    /**
     * Build a manifest file in the correct folder.
     * @param bundler The {@link FuncotatorDataSourceBundlerHttpClient} which holds all of the variables we will need.
     */
    public static void buildManifestFile(FuncotatorDataSourceBundlerHttpClient bundler) {
        try {
            FileWriter writer = new FileWriter(bundler.metadataFilePath.toAbsolutePath().toString() + "/" + MANIFEST_FILE_NAME);
            BufferedWriter buffer = new BufferedWriter(writer);
            buffer.write(
                    "Version:          0.0." + bundler.getDate() + "\n" +
                            "Use Case:         " + bundler.speciesName + "\n" +
                            "Source:           ./gatk FuncotatorDataSourceBundler -" + bundler.dsOrganism + "-species-name " + bundler.speciesName + "\n" +
                            "Alternate Source: ./gatk FuncotatorDataSourceBundler -" + bundler.dsOrganism + "-species-name " + bundler.speciesName + "\n"
            );
            buffer.close();
        } catch (IOException e) {
            throw new UserException("Error. Unable to make manifest file in location: " + bundler.metadataFilePath.toString() + "/" + MANIFEST_FILE_NAME);
        }
    }

    //==================================================================================================================
    // Getters / Setters:

    /**
     * @return A copy of the {@link String} used as the data source URL for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public String getDSUrl() {
        return this.dsURL;
    }

    /**
     * @return A copy of the {@link Path} used as the output destination path for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getOutputDestination() {
        return this.outputDestination;
    }

    /**
     * @return A copy of the {@link Path} used as the output destination for the unzipped gtf file for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getOutputUnzippedDest() {
        return this.outputUnzippedDest;
    }

    /**
     * @return A copy of the {@link String} used as the fasta data source URL for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public String getFastaURL() {
        return this.dsFastaURL;
    }

    /**
     * @return A copy of the {@link Path} used as the output fasta desination path for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getFastaOutputDestination() {
        return this.outputFastaDest;
    }

    /**
     * @return A copy of the {@link Path} used as the data source path for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getDSPath() {
        return this.dsPath;
    }

    /**
     * @return A copy of the {@link Path} used as the unzipped data source path for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getDSUnzipPath() {
        return this.dsUnzipPath;
    }

    /**
     * @return A copy of the {@link Path} used as the unzipped fasta data source path for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getFastaUnzipPath() {
        return this.dsFastaUnzipPath;
    }

    /**
     * @return A copy of the {@link Path} used as the index file path for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getIndexPath() {
        return this.indexFilePath;
    }

    /**
     * @return A copy of the {@link Path} which is the path to the config file for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getConfigPath() {
        return this.configFilePath;
    }

    /**
     * @return A copy of the {@link String} which is the url for the gtf ReadMe for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public String getGtfReadMeURL() {
        return this.dsGtfReadMeURL;
    }

    /**
     * @return A copy of the {@link String} which is the url for the fasta ReadMe for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public String getFastaReadMeURL() {
        return this.dsFastaReadMeURL;
    }

    /**
     * @return A copy of the {@link Path} which is the path to the gtf ReadMe file for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getGtfReadMePath() {
        return this.dsGtfReadMePath;
    }

    /**
     * @return A copy of the {@link Path} which is the path to the fasta ReadMe file for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public Path getFastaReadMePath() {
        return this.dsFastaReadMePath;
    }

    /**
     * @return A copy of the {@link String} which is the date for this {@link FuncotatorDataSourceBundlerHttpClient}.
     */
    public String getDate() {
        final LocalDate date = LocalDate.of(2021, Month.AUGUST, 10);
        return String.format("%d%02d%02d", date.getYear(), date.getMonthValue(), date.getDayOfMonth());
    }

    /**
     * Constructs the url for the data source file.
     * @param baseURL The {@link String} representing the base url where the data source will be found
     * @param speciesName The {@link String} representing the chosen species to gather data sources for.
     * @param fileName The {@link String} representing the file name for the data source we need to download.
     * @return The url constructed using the base URL, the species name, and the file name.
     */
    public String setURL(String baseURL, String speciesName, String fileName) {
        return baseURL + speciesName + "/" + fileName + "." + DataSourceUtils.GTF_GZ_EXTENSION;
    }

    /**
     * Constructs the url for the gtf ReadMe for the data source file.
     * @param baseURL The {@link String} representing the base url where the data source will be found.
     * @param speciesName The {@link String} representing the chosen species to gather data sources for.
     * @return The url constructed using the base URL, the species name, and the ReadMe extension.
     */
    public String setReadMeURL(String baseURL, String speciesName) {
        return baseURL + speciesName + "/" + DataSourceUtils.README_EXTENSION;
    }

    /**
     * Constructs the path where the data source to be copied is found.
     * @param speciesName The {@link String} representing the chosen species to gather data sources for.
     * @param fileName The {@link String} representing the file name for the data source we need to download.
     * @return A path to the data source.
     */
    public Path setPath(String speciesName, String fileName) {
        return IOUtils.getPath(speciesName + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + fileName + "." + DataSourceUtils.GTF_GZ_EXTENSION);
    }

    /**
     * Constructs the path where the gtf ReadMe for the data source is found.
     * @param speciesName The {@link String} representing the chosen species to gather data sources for.
     * @param fileName The {@link String} representing the file name for the data source we need to download.
     * @return A path to the data source.
     */
    public Path setGtfReadMePath(String speciesName, String fileName) {
        return IOUtils.getPath(speciesName + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + DataSourceUtils.GTF_README_EXTENSION);
    }

    /**
     * Constructs the path for where the unzipped gtf file should go.
     * @param speciesName The {@link String} representing the chosen species to gather data sources for.
     * @param fileName The {@link String} representing the file name for the data source we need to download.
     * @return A path to the unzipped gtf file.
     */
    public Path setUnzipPath(String speciesName, String fileName) {
        return IOUtils.getPath(speciesName + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + fileName + DataSourceUtils.GTF_UNZIPPED_EXTENSION);
    }

    /**
     * Constructs the url for the fasta file of the data source.
     * @param baseUrl The {@link String} representing the base URL for the fasta file.
     * @param speciesName The {@link String} reprsenting the chosen species to gather data soruces for.
     * @param fileName The {@link String} representing the file name for the data source we need to download.
     * @return The url at which the fasta file can be found.
     */
    public String setFastaUrl(String baseUrl, String speciesName, String fileName) {
        return baseUrl + speciesName + "/" + DataSourceUtils.CDNA_EXTENSION + fileName + "." + DataSourceUtils.FASTA_GZ_EXTENSION;
    }

    /**
     * Constructs the path where the fasta file for the data source to be copied is found.
     * @param speciesName The {@link String} representing the chosen species to gather data sources for.
     * @param fileName The {@link String} representing the file name for the data source we need to download.
     * @return A path to the fasta data source.
     */
    public Path setFastaPath(String speciesName, String fileName) {
        return IOUtils.getPath(speciesName + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + fileName + "." + DataSourceUtils.FASTA_GZ_EXTENSION);
    }

    /**
     * Constructs the path where the fasta ReadMe file for the data source to be copied is found.
     * @param speciesName The {@link String} representing the chosen species to gather data sources for.
     * @param fileName The {@link String} representing the file name for the data source we need to download.
     * @return A path to the fasta data source.
     */
    public Path setFastaReadMePath(String speciesName, String fileName) {
        return IOUtils.getPath(speciesName + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + DataSourceUtils.FASTA_README_EXTENSION);
    }

    /**
     * Constructs the path for where the unzipped fasta file should go.
     * @param speciesName The {@link String} representing the chosen species to gather data sources for.
     * @param fileName The {@link String} representing the file name for the data source we need to download.
     * @return A path to the unzipped fasta file.
     */
    public Path setFastaUnzipPath(String speciesName, String fileName) {
        return IOUtils.getPath(speciesName + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + fileName + DataSourceUtils.FASTA_UNZIPPED_EXTENSION);
    }

    /**
     * Constructs the path for where the index file should go.
     * @param speciesName The {@link String} representing the chosen species to gather data sources for.
     * @param fileName The {@link String} representing the file name for the data source we need to download.
     * @return A path to the index file for this data source.
     */
    public Path setIndexFilePath(String speciesName, String fileName) {
        return IOUtils.getPath(speciesName + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + fileName + "." + DataSourceUtils.GTF_GZ_EXTENSION + DataSourceUtils.IDX_EXTENSION);
    }

    /**
     * Constructs the path where the config file should go.
     * @param speciesName The {@link String} representing the chosen species to gather data sources for.
     * @return The path where the config file will be found.
     */
    public Path setConfigFilePath(String speciesName) {
        return IOUtils.getPath(speciesName + "/" + DataSourceUtils.ENSEMBL_EXTENSION + "/" + speciesName + "/" + ENSEMBL_CONFIG_NAME);
    }

    /**
     * Constructs the path where the manifest file should go.
     * @param speciesName The {@link String} representing the chosen species to gather data sources for.
     * @return The path where the manifest file will be found.
     */
    public Path setMetadataFilePath(String speciesName) {
        return IOUtils.getPath(speciesName + "/");
    }

}
