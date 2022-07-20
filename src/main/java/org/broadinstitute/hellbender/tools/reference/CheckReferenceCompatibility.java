package org.broadinstitute.hellbender.tools.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.argumentcollections.SequenceDictionaryValidationArgumentCollection;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.List;

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
            // create ref seq table and analyzing
            // ReferenceSequenceTable table = generateReferenceSequenceTable(dictionaryToCompare);

        }
        else{
            // warning: comparison lacking MD5
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
            // dictionaryPath = getHeaderForReads().
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

    /*private ReferenceSequenceTable generateReferenceSequenceTable(SAMSequenceDictionary dictionary){
        ReferenceSequenceTable table = new ReferenceSequenceTable(dictionary, CompareReferences.MD5CalculationMode.USE_DICT);

        return table;
    }*/

}
