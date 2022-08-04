package org.broadinstitute.hellbender.tools.reference;

import htsjdk.samtools.reference.FastaReferenceWriter;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.FilterFuncotations;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
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

    @Argument(fullName = "base-comparison-output", doc = "", optional = true)
    private GATKPath baseComparisonOutputDirectory;

    @Argument(fullName = "base-comparison", doc = "If provided, any mismatching sequences will be aligned for a base-comparison.", optional = true)
    private BaseComparisonMode baseComparisonMode = BaseComparisonMode.NO_BASE_COMPARISON;


    public enum MD5CalculationMode {
        // use only MD5s found in dictionary; if MD5 missing, crashes
        USE_DICT,
        // use any MD5s found in dictionary and recalculate any missing MD5s
        RECALCULATE_IF_MISSING,
        // recalculate all MD5s, regardless of presence in dictionary
        ALWAYS_RECALCULATE;
    }

    public enum BaseComparisonMode{
        NO_BASE_COMPARISON,
        FULL_ALIGNMENT,
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
        if(baseComparisonMode != BaseComparisonMode.NO_BASE_COMPARISON && baseComparisonOutputDirectory == null){
            throw new UserException.CouldNotCreateOutputFile(baseComparisonOutputDirectory, "Output directory non existent.");
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
            // only do alignment if exactly 2 inputs (1 ReferencePair)



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

    public void runFullAlignment(ReferencePair refPair, ReferenceSequenceTable table){

        // create the fastas (see below)
        // mummer executor -- sep class, pass in directory containing mummer executables
        // execute mummer on each pair of fastas --> list of VCFs
        // merge into 1 VCF, output to final location

        // create FeatureDataSource<VariantContext>
        // VariantContextWriter (writes final VCF) -- call createVCFWriter(mummeroutput.vcf)
        // create blank vcf header --> setSequenceDictionary(-R sequence dictionary) --> writer.writeHeader(header)
        // for each VCF
            // create FDS
            // loop over for each in FDS to get each variant record
            // writer.add(variant context)
        // close OR try w/r

    }


    public void runFindSNPS(ReferencePair refPair, ReferenceSequenceTable table) {
        List<SNPRecord> records = new ArrayList<>();
        if(refPair.getAnalysis().contains(ReferencePair.Status.DIFFER_IN_SEQUENCE)) {
            // reference data sources for both references -- initialized to not drop bases
            final ReferenceDataSource source1 = ReferenceDataSource.of(refPair.getRef1().toPath(), true);
            final ReferenceDataSource source2 = ReferenceDataSource.of(refPair.getRef2().toPath(), true);

            // find the mismatch sequence
            for (String sequenceName : table.getAllSequenceNames()) {
                Set<ReferenceSequenceTable.TableRow> rows = table.queryBySequenceName(sequenceName);
                if (rows.size() == 2) {
                    // if the lengths of the 2 sequences aren't equal, error - can't compare different sequence lengths, probably indel need alignment
                    ReferenceSequenceTable.TableRow[] rowArray = rows.toArray(new ReferenceSequenceTable.TableRow[0]);
                    if(rowArray[0].getLength() != rowArray[1].getLength()) {
                        throw new UserException("Sequence lengths are not equal and can't be compared. Consider running in FULL_ALIGNMENT mode.");
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
                        if(!ref1Allele.equals(ref2Allele)) {
                            SNPRecord record = new SNPRecord(sequenceName, position, new String(new byte[]{ref1Allele}), new String(new byte[]{ref2Allele}), refPair.getRef1AsString(), refPair.getRef2AsString());
                            records.add(record);
                        }
                    }
                }
                // write SNP tsv file
                TableColumnCollection columns = new TableColumnCollection(Arrays.asList("Sequence Name", "Position", refPair.getRef1AsString(), refPair.getRef2AsString()));
                try(FindSNPsOnlyTableWriter writer = new FindSNPsOnlyTableWriter(baseComparisonOutputDirectory.toPath(), columns)){
                    writer.writeHeaderIfApplies();
                    for(SNPRecord record : records){
                        writer.writeRecord(record);
                    }
                }
                catch(IOException exception){
                    throw new UserException.CouldNotCreateOutputFile(baseComparisonOutputDirectory, "Failed to write output table.", exception);
                }
            }
            source1.close();
            source2.close();
        }
    }

    public static void generateFastaForSequence(ReferenceDataSource source, String sequenceName, GATKPath output){
        try {
            int sequenceLength = source.getSequenceDictionary().getSequence(sequenceName).getSequenceLength();
            FastaReferenceWriter writer = new FastaReferenceWriterBuilder()
                    .setFastaFile(output.toPath())
                    .setBasesPerLine(80)
                    .build();
            ReferenceSequence seq = source.queryAndPrefetch(new SimpleInterval(sequenceName, 1, sequenceLength));
            writer.addSequence(seq);
            writer.close();
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Couldn't create " + output + ", encountered exception: " + e.getMessage(), e);
        }
    }

    public void doBaseComparison(ReferencePair refPair, ReferenceSequenceTable table){
        // if mismatching MD5s found between the two inputs
        if(refPair.getAnalysis().contains(ReferencePair.Status.DIFFER_IN_SEQUENCE)){
            // find the mismatch sequence
            for(String sequenceName : table.getAllSequenceNames()) {
                Set<ReferenceSequenceTable.TableRow> rows = table.queryBySequenceName(sequenceName);
                if (rows.size() == 2) {
                    // generate fasta files for mismatching sequences
                    GATKPath ref1Path = refPair.getRef1();
                    GATKPath ref2Path = refPair.getRef2();

                    String sequenceInRef1Name = ref1Path.toPath().getFileName() + "." + sequenceName;
                    String sequenceInRef2Name = ref2Path.toPath().getFileName() + "." + sequenceName;

                        /*File ref1SequenceOutput = IOUtils.createTempFileInDirectory(sequenceInRef1Name, ".fasta", sequenceOutputDirectory);
                        File ref2SequenceOutput = IOUtils.createTempFileInDirectory(sequenceInRef2Name, ".fasta", sequenceOutputDirectory);*/

                    //File ref1SequenceOutput = new File(baseComparisonOutputDirectory,   sequenceInRef1Name + ".fasta");
                    //File ref2SequenceOutput = new File(baseComparisonOutputDirectory, sequenceInRef2Name + ".fasta");

                    //generateFastaForSequence(ReferenceDataSource.of(ref1Path.toPath(), true), sequenceName, new GATKPath(ref1SequenceOutput.toString()));
                    //generateFastaForSequence(ReferenceDataSource.of(ref2Path.toPath(), true), sequenceName, new GATKPath(ref2SequenceOutput.toString()));
                }
            }
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

    public static class FindSNPsOnlyTableWriter extends TableWriter<SNPRecord> {

        public FindSNPsOnlyTableWriter(final Path table, TableColumnCollection columns) throws IOException {
            super(table, columns);
        }

        public FindSNPsOnlyTableWriter(final Writer writer, TableColumnCollection columns) throws IOException {
            super(writer, columns);
        }

        @Override
        protected void composeLine(final SNPRecord record, final DataLine dataLine) {
            dataLine.set("Sequence Name", record.sequence)
                    .set("Position", record.position)
                    .set(record.ref1, record.ref1Allele)
                    .set(record.ref2, record.ref2Allele);
        }
    }
}
