package org.broadinstitute.hellbender.tools.reference;

import htsjdk.samtools.reference.FastaReferenceWriter;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.alignment.MummerExecutor;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * Display reference comparison as a tab-delimited table and summarize reference differences.
 *
 * <p>This tool generates a MD5-keyed table comparing specified references and does an analysis to summarize the differences
 * between the references provided. Comparisons are made against a "special" reference, specified with the -R
 * argument. Subsequent references to be compared may be specified using the -references-to-compare argument.
 * The table can be directed to a file or standard output using provided command-line arguments.
 * A supplementary table keyed by sequence name can be displayed using the -display-sequences-by-name argument;
 * to display only sequence names for which the references are not consistent, run with the -display-only-differing-sequences
 * argument as well.</p>
 *
 * <h3>Input</h3>
 * <p>
 * The references and MD5 calculation mode.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A TSV file or output stream displaying the MD5-keyed table comparison, and a table analysis displayed to standard output.
 * <pre>
 * MD5	Length	hg19mini.fasta	hg19mini_1renamed.fasta
 * 8c0c38e352d8f3309eabe4845456f274	16000	1	chr1
 * 5f8388fe3fb34aa38375ae6cf5e45b89	16000	2	2
 * 94de808a3a2203dbb02434a47bd8184f	16000	3	3
 * 7d397ee919e379328d8f52c57a54c778	16000	4	4
 *
 * REFERENCE PAIR: hg19mini.fasta, hg19mini_1renamed.fasta
 * Status:
 * 	DIFFER_IN_SEQUENCE_NAMES
 * </pre>
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk CompareReferences \
 *   -R reference1.fasta \
 *   -refcomp reference2.fasta
 *   -O output.fasta \
 *   -md5-calculation-mode USE_DICT
 * </pre>
 *
 */

@CommandLineProgramProperties(
        summary = "Compare multiple references and output a tab-delimited table detailing what the differences are and a summarized analysis of each pair of references.",
        oneLineSummary = "Display reference comparison as a tab-delimited table and summarize reference differences.",
        programGroup = ReferenceProgramGroup.class
)

@DocumentedFeature
@ExperimentalFeature
public class CompareReferences extends GATKTool {



    @Argument(fullName = "references-to-compare", shortName = "refcomp", doc = "Reference sequence file(s) to compare.")
    private List<GATKPath> references;

    /**
     * Output file will be written here.
     *
     * Note: If no output file provided, table will print to standard output.
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "If specified, output reference sequence table in TSV format to this file. Otherwise print it to stdout.", optional = true)
    private GATKPath output;

    @Argument(fullName = "md5-calculation-mode", shortName = "md5-calculation-mode", doc = "MD5CalculationMode indicating method of MD5 calculation.", optional = true)
    private MD5CalculationMode md5CalculationMode = MD5CalculationMode.RECALCULATE_IF_MISSING;

    @Argument(fullName = "display-sequences-by-name", doc = "If provided, the table by sequence name will be printed.", optional = true)
    private boolean displaySequencesByName = false;

    @Argument(fullName = "display-only-differing-sequences", doc = "If provided, only display sequence names that differ in their actual sequence.", optional = true)
    private boolean onlyDisplayDifferingSequences = false;

    @Argument(fullName = "base-comparison", doc = "If provided, any mismatching, same-length sequences will be aligned for a base-comparison.", optional = true)
    private BaseComparisonMode baseComparisonMode = BaseComparisonMode.NO_BASE_COMPARISON;

    @Argument(fullName = "base-comparison-output", doc = "Output directory for base comparison outputs. Required for running base-comparison in FULL_ALIGNMENT or FIND_SNPS_ONLY mode.", optional = true)
    private GATKPath baseComparisonOutputDirectory;




    public enum MD5CalculationMode {
        // use only MD5s found in dictionary; if MD5 missing, crashes
        USE_DICT,
        // use any MD5s found in dictionary and recalculate any missing MD5s
        RECALCULATE_IF_MISSING,
        // recalculate all MD5s, regardless of presence in dictionary
        ALWAYS_RECALCULATE;
    }

    public enum BaseComparisonMode{
        // no base comparison
        NO_BASE_COMPARISON,
        // run the mummer pipeline to generate a snps file containing SNPs and INDELs for any mismatching sequences of same name
        FULL_ALIGNMENT,
        // do a base-by-base comparison of any mistmatching sequences of same name to output a table containing each base mismatch
        FIND_SNPS_ONLY
    }

    private Map<GATKPath, ReferenceDataSource> referenceSources;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        // add data source for -R reference
        referenceSources = new LinkedHashMap<>();
        referenceSources.put(referenceArguments.getReferenceSpecifier(), directlyAccessEngineReferenceDataSource());

        // add data sources for remaining references
        for(GATKPath path : references){
            referenceSources.put(path, ReferenceDataSource.of(path.toPath()));
        }
    }

    @Override
    public void traverse(){
        ReferenceSequenceTable table = new ReferenceSequenceTable(referenceSources, md5CalculationMode);
        logger.info("Building reference sequence table.");
        table.build();
        logger.info("Finished building table.");

        writeTable(table);

        if(displaySequencesByName){
            writeTableBySequenceName(table);
        }

        logger.info("Analyzing table.");
        List<ReferencePair> referencePairs = table.compareAllReferences();
        logger.info("Finished analyzing table.");
        System.out.println("*********************************************************");
        for(ReferencePair pair : referencePairs){
            System.out.println(pair);
        }

        if(referencePairs.size() != 1 && baseComparisonMode != BaseComparisonMode.NO_BASE_COMPARISON){
            throw new UserException.BadInput("");
        }
        if(baseComparisonMode != BaseComparisonMode.NO_BASE_COMPARISON){
            if(baseComparisonOutputDirectory == null) {
                throw new UserException.CouldNotCreateOutputFile(baseComparisonOutputDirectory, "Output directory not provided but required in base comparison mode.");
            }
            if(!Files.exists(baseComparisonOutputDirectory.toPath())){
                throw new UserException.CouldNotCreateOutputFile(baseComparisonOutputDirectory, "Output directory non-existent.");
            }
        }
        ReferencePair refPair = referencePairs.get(0);

        switch (baseComparisonMode) {
            case FULL_ALIGNMENT:
                runFullAlignment(refPair, table);
                break;
            case FIND_SNPS_ONLY:
                runFindSNPS(refPair, table);
                break;
        }
    }


    /**
     * Given a table, write table to output.
     *
     * Note: if no output file specified, displays to standard output.
     *
     * @param table
     */
    private void writeTable(ReferenceSequenceTable table){
        TableColumnCollection columns = new TableColumnCollection(table.getColumnNames());
        try(CompareReferences.CompareReferencesOutputTableWriter writer = output == null
            ? new CompareReferences.CompareReferencesOutputTableWriter(new OutputStreamWriter(System.out), columns)
            : new CompareReferences.CompareReferencesOutputTableWriter(output.toPath(), columns)
        ){
            writer.writeHeaderIfApplies();
            for(ReferenceSequenceTable.TableRow row : table){
                writer.writeRecord(row);
            }
        }
        catch(IOException exception){
             throw (output == null)
             ? new UserException.CouldNotCreateOutputFile("System.out", "Failed to write output table.", exception)
             : new UserException.CouldNotCreateOutputFile(output, "Failed to write output table.", exception);
        }
    }

    /**
     * Given a table, write table by sequence name to standard output
     *
     * @param table
     */
    public void writeTableBySequenceName(ReferenceSequenceTable table){
        System.out.print("*********************************************************\n");
        System.out.print("Name \tMD5 \tReference\n");
        for(String sequenceName : table.getAllSequenceNames()){
            Set<ReferenceSequenceTable.TableRow> rows = table.queryBySequenceName(sequenceName);
            if(onlyDisplayDifferingSequences) {
                if(rows.size() > 1) {
                    displayBySequenceName(rows, sequenceName);
                }
            } else {
                displayBySequenceName(rows, sequenceName);
            }
        }
    }

    private void displayBySequenceName(Set<ReferenceSequenceTable.TableRow> rows, String sequenceName){
        System.out.print(sequenceName);
        for(ReferenceSequenceTable.TableRow row : rows) {
            ReferenceSequenceTable.TableEntry[] entries = row.getEntries();
            System.out.print("\n\t" + row.getMd5() + "\t");

            for(int i = 2; i < entries.length; i++) {
                if(entries[i].getColumnValue().equals(sequenceName)) {
                    System.out.print(entries[i].getColumnName() + "\t");
                }
            }
        }
        System.out.println();
    }

    private void runFullAlignment(ReferencePair refPair, ReferenceSequenceTable table){
        List<File> snpsFiles = new ArrayList<>();
        if(refPair.getAnalysis().contains(ReferencePair.Status.DIFFER_IN_SEQUENCE)) {
            // find the mismatch sequence
            for (String sequenceName : table.getAllSequenceNames()) {
                Set<ReferenceSequenceTable.TableRow> rows = table.queryBySequenceName(sequenceName);
                if (rows.size() == 2) {
                    // generate fasta files for mismatching sequences
                    GATKPath ref1Path = refPair.getRef1();
                    GATKPath ref2Path = refPair.getRef2();

                    String sequenceInRef1Name = ref1Path.toPath().getFileName() + "." + sequenceName;
                    String sequenceInRef2Name = ref2Path.toPath().getFileName() + "." + sequenceName;

                    File ref1TempFastaOutput = IOUtils.createTempFile(sequenceInRef1Name, ".fasta");
                    File ref2TempFastaOutput = IOUtils.createTempFile(sequenceInRef2Name, ".fasta");

                    GATKPath ref1Fasta = generateFastaForSequence(ReferenceDataSource.of(ref1Path.toPath(), true), sequenceName, new GATKPath(ref1TempFastaOutput.toString()));
                    GATKPath ref2Fasta = generateFastaForSequence(ReferenceDataSource.of(ref2Path.toPath(), true), sequenceName, new GATKPath(ref2TempFastaOutput.toString()));

                    // pass fastas into mummer, get back a vcf and add to list
                    MummerExecutor executor = new MummerExecutor();
                    logger.info("Running mummer alignment on sequence " + sequenceName);
                    File tempSnpsDirectory = IOUtils.createTempDir("tempsnps");
                    File mummerOutput = executor.executeMummer(ref1Fasta.toPath().toFile(), ref2Fasta.toPath().toFile(), tempSnpsDirectory, sequenceName);
                    logger.info("Finished running mummer alignment on sequence " + sequenceName);
                    snpsFiles.add(mummerOutput);
                }
            }
            // merge individual snps files
            File snps = new File(baseComparisonOutputDirectory.toPath().toString(), String.format("%s_%s.snps", refPair.getRef1AsString(), refPair.getRef2AsString()));
            try (PrintWriter writer = new PrintWriter(snps)) {
                for (File file : snpsFiles) {
                    try (XReadLines reader = new XReadLines(file)) {
                            for (String line : reader) {
                                writer.write(line + "\n");
                            }
                    } catch (IOException e) {
                        throw new UserException("Error merging show-snps outputs.", e);
                    }
                }
            }
            catch(IOException e){
                throw new UserException("Error writing show-snps output file.");
            }
            File vcf = new File(baseComparisonOutputDirectory.toPath().toString(), String.format("%s_%s.vcf", refPair.getRef1AsString(), refPair.getRef2AsString()));
            try(XReadLines reader = new XReadLines(snps);
                VariantContextWriter writer = createVCFWriter(vcf);
                ReferenceDataSource source = ReferenceDataSource.of(refPair.getRef1().toPath())){
                VCFHeader header = new VCFHeader();
                header.setSequenceDictionary(getReferenceDictionary());
                writer.writeHeader(header);
                MummerIndel currentIndel = null;
                int previousPos = -1;
                    for (String line : reader) {
                        String[] fields = line.split("\\t", -1);
                        String contig = fields[10];
                        int pos = Integer.valueOf(fields[0]);
                        String ref = fields[1];
                        String alt = fields[2];

                        // insertion
                        if (ref.equals(".") && !alt.equals(".")) {
                            if (pos == previousPos && currentIndel.isInsertion) {
                                currentIndel.alt += alt;
                            } else {
                                if (currentIndel != null) {
                                    writer.add(currentIndel.getAsVCFRecord());
                                }
                                String refPadding = new String(source.queryAndPrefetch(new SimpleInterval(contig, pos-1, pos-1)).getBases());
                                currentIndel = new MummerIndel(contig, refPadding, refPadding + alt, pos-1, true, false);
                            }
                        }

                        // deletion
                        else if(!ref.equals(".") && alt.equals(".")){
                            if(pos == previousPos+1 && currentIndel.isDeletion){
                               currentIndel.ref += ref;
                            }
                            else {
                                if(currentIndel != null) {
                                    writer.add(currentIndel.getAsVCFRecord());
                                }
                                String refPadding = new String(source.queryAndPrefetch(new SimpleInterval(contig, pos-1, pos-1)).getBases());
                                currentIndel = new MummerIndel(contig, refPadding + ref, refPadding, pos-1, false, true);
                            }
                        }
                        // snp
                        else {
                            if(currentIndel != null){
                                writer.add(currentIndel.getAsVCFRecord());
                                currentIndel = null;
                            }
                            VariantContextBuilder vcfBuilder = new VariantContextBuilder();
                            VariantContext record = vcfBuilder.chr(contig).start(pos).stop(pos).alleles(ref, alt).make();
                            writer.add(record);
                        }
                        previousPos = pos;
                    }
                    if(currentIndel != null){
                        writer.add(currentIndel.getAsVCFRecord());
                    }
            }
            catch(IOException e){
                throw new UserException("");
            }
            catch(NumberFormatException e){
                throw new UserException("");
            }
        }
        else{
            logger.info("No mismatching sequences found.");
        }
    }

    private void runFindSNPS(ReferencePair refPair, ReferenceSequenceTable table) {
        if(refPair.getAnalysis().contains(ReferencePair.Status.DIFFER_IN_SEQUENCE)) {
            TableColumnCollection columns = new TableColumnCollection(Arrays.asList("Sequence Name", "Position", refPair.getRef1AsString(), refPair.getRef2AsString()));
            File snpsOutput = new File(baseComparisonOutputDirectory.toPath().toString(), String.format("%s_%s_snps.tsv", refPair.getRef1AsString(), refPair.getRef2AsString()));

            try(FindSNPsOnlyTableWriter writer = new FindSNPsOnlyTableWriter(snpsOutput.toPath(), columns); final ReferenceDataSource source1 = ReferenceDataSource.of(refPair.getRef1().toPath(), true);
                final ReferenceDataSource source2 = ReferenceDataSource.of(refPair.getRef2().toPath(), true)){
                writer.writeHeaderIfApplies();

                // find the mismatch sequence
                for (String sequenceName : table.getAllSequenceNames()) {
                    Set<ReferenceSequenceTable.TableRow> rows = table.queryBySequenceName(sequenceName);
                    if (rows.size() == 2) {
                        // if the lengths of the 2 sequences aren't equal, error - can't compare different sequence lengths, probably indel need alignment
                        ReferenceSequenceTable.TableRow[] rowArray = rows.toArray(new ReferenceSequenceTable.TableRow[0]);
                        if (rowArray[0].getLength() != rowArray[1].getLength()) {
                            logger.warn("Sequence lengths are not equal and can't be compared. Consider running in FULL_ALIGNMENT mode.");
                            continue;
                        }

                        int sequenceLength = rowArray[0].getLength();
                        Iterator<Byte> ref1BaseIterator = source1.query(new SimpleInterval(sequenceName, 1, sequenceLength));
                        Iterator<Byte> ref2BaseIterator = source2.query(new SimpleInterval(sequenceName, 1, sequenceLength));
                        int position = 0;

                        // find SNPS and generate SNPRecords
                        while (ref1BaseIterator.hasNext()) {
                            position++;
                            Byte ref1Allele = ref1BaseIterator.next();
                            Byte ref2Allele = ref2BaseIterator.next();
                            if (!ref1Allele.equals(ref2Allele)) {
                                SNPRecord record = new SNPRecord(sequenceName, position, new String(new byte[]{ref1Allele}), new String(new byte[]{ref2Allele}), refPair.getRef1AsString(), refPair.getRef2AsString());
                                writer.writeRecord(record);
                            }
                        }
                    }
                }
            }
            catch(IOException exception){
                throw new UserException.CouldNotCreateOutputFile(baseComparisonOutputDirectory + "/" + snpsOutput, "Failed to write output table.", exception);
            }
        }
        else{
            logger.info("No mismatching sequences found.");
        }
    }

    /**
     * Method to generate a fasta file for an individual sequence in a reference
     *
     * @param source ReferenceDataSource for the reference
     * @param sequenceName target sequence name as a String
     * @param output GATKPath for file output
     * @return output location as a GATKPath
     */
    public static GATKPath generateFastaForSequence(ReferenceDataSource source, String sequenceName, GATKPath output){
        try {
            int sequenceLength = source.getSequenceDictionary().getSequence(sequenceName).getSequenceLength();
            FastaReferenceWriter writer = new FastaReferenceWriterBuilder()
                    .setFastaFile(output.toPath())
                    .setBasesPerLine(80)
                    .build();
            ReferenceSequence seq = source.queryAndPrefetch(new SimpleInterval(sequenceName, 1, sequenceLength));
            writer.addSequence(seq);
            writer.close();
            return output;
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Couldn't create " + output + ", encountered exception: " + e.getMessage(), e);
        }
    }

    @Override
    public void closeTool() {
        for(Map.Entry<GATKPath, ReferenceDataSource> entry : referenceSources.entrySet()){
            if(!entry.getKey().equals(referenceArguments.getReferenceSpecifier())){
                entry.getValue().close();
            }
        }
    }

    /**
     * TableWriter to format and write the table output.
     */
    public static class CompareReferencesOutputTableWriter extends TableWriter<ReferenceSequenceTable.TableRow> {

        public CompareReferencesOutputTableWriter(final Path table, TableColumnCollection columns) throws IOException {
            super(table, columns);
        }

        public CompareReferencesOutputTableWriter(final Writer writer, TableColumnCollection columns) throws IOException {
            super(writer, columns);
        }

        @Override
        protected void composeLine(final ReferenceSequenceTable.TableRow record, final DataLine dataLine) {
            List<String> columnNames = record.getColumnNames();
            ReferenceSequenceTable.TableEntry[] entries = record.getEntries();
            for(int i = 0; i < columnNames.size(); i++){
                dataLine.set(entries[i].getColumnName(), entries[i].getColumnValue());
            }
        }
    }

    /**
     * Minimal class representing a single SNP in a pair of references.
     * Stores the sequence, position of the SNP, and the allele in each reference.
     */
    private static class SNPRecord{
        String sequence;
        int position;
        String ref1Allele;
        String ref2Allele;
        String ref1;
        String ref2;

        public SNPRecord(String sequence, int position, String ref1Allele, String ref2Allele, String ref1, String ref2){
            this.sequence = sequence;
            this.position = position;
            this.ref1Allele = ref1Allele;
            this.ref2Allele = ref2Allele;
            this.ref1 = ref1;
            this.ref2 = ref2;
        }
    }

    /**
     * TableWriter to format and write SNP table output.
     */
    public static class FindSNPsOnlyTableWriter extends TableWriter<SNPRecord> {

        public FindSNPsOnlyTableWriter(final Path table, TableColumnCollection columns) throws IOException {
            super(table, columns);
        }

        @Override
        protected void composeLine(final SNPRecord record, final DataLine dataLine) {
            dataLine.set("Sequence Name", record.sequence)
                    .set("Position", record.position)
                    .set(record.ref1, record.ref1Allele)
                    .set(record.ref2, record.ref2Allele);
        }
    }

    private static class MummerIndel{
        String ref;
        String alt;
        String chr;
        int pos;
        boolean isInsertion;
        boolean isDeletion;

        MummerIndel(String chr, String ref, String alt, int pos, boolean isInsertion, boolean isDeletion){
            this.ref = ref;
            this.alt = alt;
            this.chr = chr;
            this.pos = pos;
            this.isInsertion = isInsertion;
            this.isDeletion = isDeletion;
        }

        public VariantContext getAsVCFRecord(){
            // builder, add fields
            VariantContextBuilder vcfBuilder = new VariantContextBuilder();
            int stopPos = isInsertion ? pos : pos + ref.length()-1;
            VariantContext record = vcfBuilder.chr(chr).start(pos).stop(stopPos).alleles(ref, alt).make();
            return record;
        }


    }
}
