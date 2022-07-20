package org.broadinstitute.hellbender.tools.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.argumentcollections.SequenceDictionaryValidationArgumentCollection;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class CheckReferenceCompatibility extends GATKTool {

    @Argument(fullName = "references-to-compare", shortName = "refcomp", doc = "Reference sequence file(s) to compare.")
    private List<GATKPath> references;

    private SAMSequenceDictionary dictionaryToCompare;
    private GATKPath dictionaryPath;


    @Override
    protected SequenceDictionaryValidationArgumentCollection getSequenceDictionaryValidationArgumentCollection() {
        return new SequenceDictionaryValidationArgumentCollection.NoValidationCollection();
    }

    @Override
    public void onTraversalStart() {
        initializeSequenceDictionaryForInput();
    }

    @Override
    public void traverse() {
        if(dictionaryHasMD5s(dictionaryToCompare)){
            ReferenceSequenceTable table = generateReferenceSequenceTable(dictionaryToCompare);
            List<ReferencePair> refPairs = table.analyzeTable();
        }
        else{
            logger.warn("Comparison lacking MD5.");
            // tap into SequenceDictionaryUtils
        }

        // evaluateStatus();
    }

    @Override
    public Object onTraversalSuccess() {
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
    }

    private void initializeSequenceDictionaryForInput(){
        if(hasReads()){
            if(readArguments.getReadPathSpecifiers().size() > 1){
                throw new UserException.BadInput("Tool analyzes one BAM at a time.");
            }
            dictionaryToCompare = getHeaderForReads().getSequenceDictionary();
            dictionaryPath = readArguments.getReadPathSpecifiers().get(0);
        }
        else{

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

    private ReferenceSequenceTable generateReferenceSequenceTable(SAMSequenceDictionary dictionary){
        Map<GATKPath, SAMSequenceDictionary> dictionaries = new LinkedHashMap<>();
        dictionaries.put(dictionaryPath, dictionary);

        for(GATKPath path : references){
            dictionaries.put(getReferencePath(), ReferenceDataSource.of(path.toPath()).getSequenceDictionary());
        }

        ReferenceSequenceTable table = new ReferenceSequenceTable(dictionaries);

        return table;
    }

}
