package org.broadinstitute.hellbender.tools.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.SequenceDictionaryValidationArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.file.Path;
import java.util.*;

/**
 * Check a BAM/VCF for compatibility against specified references.
 *
 * <p>This tool generates a table analyzing the compatibility of a sequence file against provided references. The tool works to compare
 * BAM/CRAMs (specified using the -I argument) as well as VCFs (specified using the -V argument) against provided
 * reference(s), specified using the -references-to-compare argument. The table can be directed to a file or standard
 * output using provided command-line arguments.
 *</p>
 *
 * <h3>Input</h3>
 * <p>
 * A BAM/CRAM/VCF and reference(s) for comparison.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A TSV file or output stream displaying the compatibility table.
 * <pre>
 * #Current Reference: reads_data_source_test1_withmd5s_missingchr1.bam
 * Reference	Compatibility
 * hg19mini.fasta	Compatible, the sequence dictionary in reads_data_source_test1_withmd5s_missingchr1.bam is a subset of the hg19mini.fasta reference sequence dictionary. Missing sequence(s): [1]
 * hg19mini_1renamed.fasta	Compatible, the sequence dictionary in reads_data_source_test1_withmd5s_missingchr1.bam is a subset of the hg19mini_1renamed.fasta reference sequence dictionary. Missing sequence(s): [chr1]
 * hg19mini_chr2snp.fasta	Not compatible. Status: [DIFFER_IN_SEQUENCE, DIFFER_IN_SEQUENCES_PRESENT]. Run CompareReferences tool for more information on reference differences.
 * </pre>
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk CheckReferenceCompatibility \
 *   -refcomp reference.fasta \
 *   -I reads.bam \
 *   -O output.table
 * </pre>
 *
 */

@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = ReferenceProgramGroup.class
)

public class CheckReferenceCompatibility extends GATKTool {

    @Argument(fullName = "references-to-compare", shortName = "refcomp", doc = "Reference sequence file(s) to compare.")
    private List<GATKPath> references;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants", optional = true)
    public GATKPath vcfPath;

    /**
     * Output file will be written here.
     *
     * Note: If no output file provided, table will print to standard output.
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "If specified, output reference sequence table in TSV format to this file. Otherwise print it to stdout.", optional = true)
    private GATKPath output;

    private SAMSequenceDictionary dictionaryToCompare;
    private GATKPath dictionaryPath;
    private String dictionaryName;
    private Map<GATKPath, SAMSequenceDictionary> dictionaries = new LinkedHashMap<>();
    private List<CompatibilityRecord> referenceCompatibilities = new ArrayList<>();

    @Override
    protected SequenceDictionaryValidationArgumentCollection getSequenceDictionaryValidationArgumentCollection() {
        return new SequenceDictionaryValidationArgumentCollection.NoValidationCollection();
    }

    @Override
    public void onTraversalStart() {
        initializeSequenceDictionaryForInput();
    }

    private void initializeSequenceDictionaryForInput(){
        // BAMs/CRAMs
        if(hasReads() ^ vcfPath != null){
            if (hasReads()) {
                if (readArguments.getReadPathSpecifiers().size() > 1) {
                    throw new UserException.BadInput("Tool analyzes one BAM at a time.");
                }
                dictionaryToCompare = getHeaderForReads().getSequenceDictionary();
                dictionaryPath = readArguments.getReadPathSpecifiers().get(0);
                dictionaryName = dictionaryPath.toPath().getFileName().toString();
            }
            // VCFs
            else {
                try(final FeatureDataSource<VariantContext> vcfReader = new FeatureDataSource<>(vcfPath.toString())){
                    VCFHeader header = (VCFHeader) vcfReader.getHeader();
                    dictionaryToCompare = header.getSequenceDictionary();
                    dictionaryPath = vcfPath;
                    dictionaryName = dictionaryPath.toPath().getFileName().toString();
                }
            }
            dictionaries.put(dictionaryPath, dictionaryToCompare);
            for(GATKPath path : references){
                dictionaries.put(path, ReferenceDataSource.of(path.toPath()).getSequenceDictionary());
            }
        }
    }

    @Override
    public void traverse() {
        if(dictionaryHasMD5s(dictionaryToCompare)){
            ReferenceSequenceTable table = new ReferenceSequenceTable(dictionaries);
            table.build();

            List<ReferencePair> refPairs = table.compareAgainstKeyReference(dictionaryPath);
            for(ReferencePair pair : refPairs){
                evaluateCompatibility(pair, table);
            }
        }
        else{
            logger.warn("*************************************************************************************************************************");
            logger.warn("* Comparison lacking MD5. All comparisons based on sequence name and length, which could hide mismatching references.  *");
            logger.warn("*************************************************************************************************************************");

            for(Map.Entry<GATKPath, SAMSequenceDictionary> entry : dictionaries.entrySet()){
                if(!entry.getValue().equals(dictionaryToCompare)){
                    evaluateCompatibility(entry.getValue(), entry.getKey());
                }
            }
        }
        writeOutput();
    }

    /**
     * Write compatibility information to output in table format.
     *
     * Note: if no output file specified, displays to standard output.
     *
     */
    private void writeOutput(){
        TableColumnCollection columns = new TableColumnCollection(Arrays.asList("Reference", "Compatibility"));
        try(CheckReferenceCompatibilityTableWriter writer = output == null
                ? new CheckReferenceCompatibilityTableWriter(new OutputStreamWriter(System.out), columns)
                : new CheckReferenceCompatibilityTableWriter(output.toPath(), columns)
        ){
            writer.writeComment(String.format("Current Reference: %s", dictionaryName));
            writer.writeHeaderIfApplies();
            for(CompatibilityRecord record : referenceCompatibilities){
                writer.writeRecord(record);
            }
        }
        catch(IOException exception){
            throw new UserException.CouldNotCreateOutputFile(output, "Failed to write output table.", exception);
        }
    }

    private boolean dictionaryHasMD5s(SAMSequenceDictionary dictionary){
        for (SAMSequenceRecord record : dictionary.getSequences()) {
            String md5FromDict = record.getMd5();
            if (md5FromDict == null) {
                return false;
            }
        }
        return true;
    }

    /** Evaluate the compatibility of a given pair of references when MD5s are present.
     *
     * @param refPair ReferencePair containing the dictionary file and a reference it is being compared against
     * @param table
     */
    private void evaluateCompatibility(ReferencePair refPair, ReferenceSequenceTable table){
        EnumSet<ReferencePair.Status> status = refPair.getAnalysis();
        if(status.contains(ReferencePair.Status.EXACT_MATCH)){
            referenceCompatibilities.add(new CompatibilityRecord(refPair.getRef2(), "Compatible, the sequence dictionaries exactly match"));
        }
        else if(status.contains(ReferencePair.Status.SUBSET) && status.size() == 1){
            referenceCompatibilities.add(new CompatibilityRecord(refPair.getRef2(),String.format("Compatible, the sequence dictionary in %s is a subset of the %s reference sequence dictionary. Missing sequence(s): %s", refPair.getRef1AsString(), refPair.getRef2AsString(), getMissingSequencesIfSubset(dictionaries.get(refPair.getRef2())))));
        }
        else{
            referenceCompatibilities.add(new CompatibilityRecord(refPair.getRef2(), String.format("Not compatible. Status: %s. Run CompareReferences tool for more information on reference differences.", status)));
        }
    }

    /** Evaluate the compatibility of a given pair of references when MD5s are not present.
     *
     * Note: compatibility calls are made based on sequence name and length, alone, and are therefore imprecise.
     *
     * @param referenceDict SAMSequenceDictionary to be compared against the key sequence file
     * @param referenceDictPath path to provided dictionary
     */
    private void evaluateCompatibility(SAMSequenceDictionary referenceDict, GATKPath referenceDictPath){
        String referenceDictName = referenceDictPath.toPath().getFileName().toString();
        SequenceDictionaryUtils.SequenceDictionaryCompatibility compatibilityStatus = SequenceDictionaryUtils.compareDictionaries(referenceDict, dictionaryToCompare, false);
        if(compatibilityStatus.equals(SequenceDictionaryUtils.SequenceDictionaryCompatibility.IDENTICAL)){
            referenceCompatibilities.add(new CompatibilityRecord(referenceDictPath, "All sequence names and lengths match in the sequence dictionaries. Since the MD5s are lacking, we can't confirm there aren't mismatching bases in the references."));
        }
        else if(compatibilityStatus.equals(SequenceDictionaryUtils.SequenceDictionaryCompatibility.SUPERSET)){
            referenceCompatibilities.add(new CompatibilityRecord(referenceDictPath, String.format("All sequence names and lengths present in the sequence dictionaries match, but %s is a subset of %s. Missing sequence(s): %s. Since the MD5s are lacking, we can't confirm there aren't mismatching bases in the references.", dictionaryName, referenceDictName, getMissingSequencesIfSubset(referenceDict))));
        }
        else{
            referenceCompatibilities.add(new CompatibilityRecord(referenceDictPath, String.format("Not compatible. Status: %s. Run CompareReferences tool for more information on reference differences.", compatibilityStatus)));
        }
    }

    private List<String> getMissingSequencesIfSubset(SAMSequenceDictionary referenceDict){
        Set<String> commonSequences = SequenceDictionaryUtils.getCommonContigsByName(dictionaryToCompare, referenceDict);
        List<String> missingSequences = SequenceDictionaryUtils.getContigNamesList(referenceDict);
        missingSequences.removeAll(commonSequences);

        return missingSequences;
    }

    @Override
    public Object onTraversalSuccess() {
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
    }

    /**
     * Minimal class representing a record in the compatibility output table.
     */
    private static class CompatibilityRecord {
        GATKPath ref;
        String compatibilityOutput;

        CompatibilityRecord(GATKPath ref, String compatibility){
            this.ref = ref;
            compatibilityOutput = compatibility;
        }

        public String getRefAsString(){
            return ref.toPath().getFileName().toString();
        }

        public String getCompatibility(){
            return compatibilityOutput;
        }
    }

    /**
     * TableWriter to format and write the table output.
     */
    public static class CheckReferenceCompatibilityTableWriter extends TableWriter<CompatibilityRecord> {

        public CheckReferenceCompatibilityTableWriter(final Path table, TableColumnCollection columns) throws IOException {
            super(table, columns);
        }

        public CheckReferenceCompatibilityTableWriter(final Writer writer, TableColumnCollection columns) throws IOException {
            super(writer, columns);
        }

        @Override
        protected void composeLine(final CompatibilityRecord record, final DataLine dataLine){
            dataLine.set("Reference", record.getRefAsString())
                    .set("Compatibility", record.getCompatibility());
        }
    }
}
