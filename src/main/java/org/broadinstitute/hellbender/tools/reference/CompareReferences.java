package org.broadinstitute.hellbender.tools.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.LocalAssembler;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.nio.file.Path;
import java.util.*;

/**
 *
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = ReferenceProgramGroup.class
)
public class CompareReferences extends GATKTool {

    @Argument(fullName = "references-to-compare", shortName = "refcomp", doc = "Reference(s) to compare.")
    private List<GATKPath> references;

    private Map<GATKPath, ReferenceDataSource> referenceSources;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        // add data source for -R reference
        referenceSources = new LinkedHashMap<>();
        referenceSources.put(getReferencePath(), directlyAccessEngineReferenceDataSource());

        // add data sources for remaining references
        for(GATKPath path : references){
            referenceSources.put(path, ReferenceDataSource.of(path.toPath()));
        }
    }

    @Override
    public void traverse() {

      /*  // get dictionaries for both references
        SAMSequenceDictionary dict1 = referenceSources.get(getReferencePath()).getSequenceDictionary();
        SAMSequenceDictionary dict2 = referenceSources.get(references.get(0)).getSequenceDictionary();

        // go through dictionaries and find MD5s for all contigs
        List<String> dict1Md5s = new ArrayList<>();
        List<String> dict2Md5s = new ArrayList<>();*/

        Map<String, List<TableEntry>> table = new LinkedHashMap<>();

        for(Map.Entry<GATKPath, ReferenceDataSource> entry : referenceSources.entrySet()){
            SAMSequenceDictionary dictionary = entry.getValue().getSequenceDictionary();
            for(SAMSequenceRecord record : dictionary.getSequences()){
                System.out.println();
                String name = record.getSequenceName();
                int length = record.getSequenceLength();
                String md5 = record.getMd5();
                String reference = entry.getKey().toPath().getFileName().toString();
                System.out.printf("Name: %s Length: %d md5: %s ref: %s\n", name, length, md5, reference);
                TableEntry newEntry = new TableEntry(name, length, md5, reference);

                // map each MD5 to List of TableEntry objects containing length, md5, and name
                if (!table.containsKey(md5)) {
                    table.put(md5, new ArrayList<TableEntry>());
                }
                table.get(md5).add(newEntry);
            }
        }

        // MD5  Length  Ref1  Ref2
        //System.out.printf("%s\t%s\t%s\t%s\n", "MD5", "Length", , , "");


        // use string format to output as a table
        for(Map.Entry<String, List<TableEntry>> entry : table.entrySet()){
            String currMd5 = entry.getKey();

            for(TableEntry currEntry : entry.getValue()){
                //System.out.printf("%s\t%d\t%s\t");
            }

        }




            // md5
            // length of seq
            // name in ref 1
            // name in ref 2 .. n
            //
           // System.out.printf("%s\t%d\t%s\n", )

    }

    @Override
    public Object onTraversalSuccess() {
        return null;
    }

    @Override
    public void closeTool() {
        for(Map.Entry<GATKPath, ReferenceDataSource> entry : referenceSources.entrySet()){
            if(!entry.getKey().equals(getReferencePath())){
                entry.getValue().close();
            }
        }
    }

    private static class TableEntry {
        private final String sequenceName;
        private final int sequenceLength;
        private final String MD5;
        private final String REFERENCE;

        public TableEntry(String name, int length, String md5, String reference){
            MD5 = md5;
            sequenceLength = length;
            REFERENCE = reference;
            sequenceName = name;
        }

    }


}
