package org.broadinstitute.hellbender.tools.reference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 *
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = ReferenceProgramGroup.class
)
public class CompareReferences extends GATKTool {

    @Argument(fullName = "references-to-compare", shortName = "refcomp", doc = "")
    private List<GATKPath> references;

    private Map<GATKPath, ReferenceDataSource> referenceSources;

    @Override
    public void onTraversalStart() {
        try(ReferenceDataSource testDataSource = ReferenceDataSource.of(references.get(0).toPath())){
            SAMSequenceDictionary dictionary = testDataSource.getSequenceDictionary();
            for(SAMSequenceRecord record : dictionary.getSequences()){
                String name = record.getSequenceName();
                int length = record.getSequenceLength();
                Iterator<Byte> baseIterator = testDataSource.query(new SimpleInterval(name, 1, length));
                while(baseIterator.hasNext()){
                    Byte b = baseIterator.next();

                }
            }
        }
    }

    @Override
    public void traverse() {

    }

    @Override
    public Object onTraversalSuccess() {

    }

    @Override
    public void closeTool() {

    }

}
