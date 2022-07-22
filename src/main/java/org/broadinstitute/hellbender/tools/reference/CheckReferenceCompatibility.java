package org.broadinstitute.hellbender.tools.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.SequenceDictionaryValidationArgumentCollection;
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

@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = ReferenceProgramGroup.class
)

public class CheckReferenceCompatibility extends GATKTool {

    @Argument(fullName = "references-to-compare", shortName = "refcomp", doc = "Reference sequence file(s) to compare.")
    private List<GATKPath> references;

    /**
     * Output file will be written here.
     *
     * Note: If no output file provided, table will print to standard output.
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "If specified, output reference sequence table in TSV format to this file. Otherwise print it to stdout.", optional = true)
    private GATKPath output;


    private SAMSequenceDictionary dictionaryToCompare;
    private GATKPath dictionaryPath;
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
        if(hasReads()){
            if(readArguments.getReadPathSpecifiers().size() > 1){
                throw new UserException.BadInput("Tool analyzes one BAM at a time.");
            }
            dictionaryToCompare = getHeaderForReads().getSequenceDictionary();
            dictionaryPath = readArguments.getReadPathSpecifiers().get(0);

            dictionaries.put(dictionaryPath, dictionaryToCompare);
            for(GATKPath path : references){
                dictionaries.put(path, ReferenceDataSource.of(path.toPath()).getSequenceDictionary());
            }
        }
        // VCFs
        else{

        }
    }

    @Override
    public void traverse() {
        if(dictionaryHasMD5s(dictionaryToCompare)){
            ReferenceSequenceTable table = new ReferenceSequenceTable(dictionaries);
            table.build();

            List<ReferencePair> refPairs = table.compareAgainstKeyReference(dictionaryPath);
            for(ReferencePair pair : refPairs){
                evaluateCompatibility(pair);
            }
        }
        else{
            logger.warn("****************************");
            logger.warn("* Comparison lacking MD5.  *");
            logger.warn("****************************");

            for(Map.Entry<GATKPath, SAMSequenceDictionary> entry : dictionaries.entrySet()){
                if(!entry.getValue().equals(dictionaryToCompare)){
                    evaluateCompatibility(entry.getValue(), entry.getKey());
                }
            }
        }
        writeOutput();
    }

    public void writeOutput(){
        TableColumnCollection columns = new TableColumnCollection(Arrays.asList("Reference", "Compatibility"));
        try(CheckReferenceCompatibilityTableWriter writer = output == null
                ? new CheckReferenceCompatibilityTableWriter(new OutputStreamWriter(System.out), columns)
                : new CheckReferenceCompatibilityTableWriter(output.toPath(), columns)
        ){
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
        for(SAMSequenceRecord record : dictionary.getSequences()){
            String md5FromDict = record.getMd5();
            if(md5FromDict == null){
                return false;
            }
        }
        return true;
    }

    // for MD5s present
    public void evaluateCompatibility(ReferencePair refPair){
        if(refPair.getAnalysis().contains(ReferencePair.Status.EXACT_MATCH)){
            referenceCompatibilities.add(new CompatibilityRecord(refPair.getRef1(), "Compatible, exact match"));
            //System.out.println(String.format("%s: Compatible, exact match", refPair.getRef1AsString()));
        }
        else if(refPair.getAnalysis().contains(ReferencePair.Status.SUBSET)){
            referenceCompatibilities.add(new CompatibilityRecord(refPair.getRef1(),String.format("Compatible, %s is a subset of %s reference.", refPair.getRef1AsString(), refPair.getRef2AsString())));
            //System.out.println(String.format("%s: Compatible, %s is a subset of %s reference.", refPair.getRef1AsString(), refPair.getRef1AsString(), refPair.getRef2AsString() ));
        }
        else{
            referenceCompatibilities.add(new CompatibilityRecord(refPair.getRef1(),"Not compatible. Run CompareReferences tool for more information on reference differences."));
            //System.out.println(String.format("%s: Not compatible. Run CompareReferences tool for more information on reference differences.", refPair.getRef1AsString()));
        }
    }

    // for no MD5s
    public void evaluateCompatibility(SAMSequenceDictionary dict, GATKPath dictPath){
        String dictName = dictPath.toPath().getFileName().toString();
        SequenceDictionaryUtils.SequenceDictionaryCompatibility compatibilityStatus = SequenceDictionaryUtils.compareDictionaries(dict, dictionaryToCompare, false);
        if(compatibilityStatus.equals(SequenceDictionaryUtils.SequenceDictionaryCompatibility.IDENTICAL)){
            referenceCompatibilities.add(new CompatibilityRecord(dictPath, "All sequence name and lengths match."));
            //System.out.println(String.format("%s: All sequence name and lengths match.", dictName));
        }
        else if(compatibilityStatus.equals(SequenceDictionaryUtils.SequenceDictionaryCompatibility.SUPERSET)){
            referenceCompatibilities.add(new CompatibilityRecord(dictPath, String.format("All present sequence names and lengths match, but %s is a subset of %s.", dictName, dictionaryPath)));
            //System.out.println(String.format("%s: All present sequence names and lengths match, but %s is a subset of %s.", dictName, dictName, dictionaryPath));
        }
        else{
            referenceCompatibilities.add(new CompatibilityRecord(dictPath, "Not compatible. Run CompareReferences tool for more information on reference differences."));
            //System.out.println(String.format("%s: Not compatible. Run CompareReferences tool for more information on reference differences.", dictName));
        }

    }

    @Override
    public Object onTraversalSuccess() {
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
    }

    public static class CompatibilityRecord {
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
            for(int i = 0; i < 2; i++){
                dataLine.set("Reference", record.getRefAsString())
                        .set("Compatibility", record.getCompatibility());
            }
        }
    }



}
