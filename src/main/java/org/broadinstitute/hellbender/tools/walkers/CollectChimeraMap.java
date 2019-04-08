package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

@CommandLineProgramProperties(
        summary = "Chimera",
        oneLineSummary = "Chimera",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
public class CollectChimeraMap extends ReadWalker {
    static final int MINIMUM_CHIMERA_FRAGMENT_SIZE = 2000;

    RealMatrix chimeraMap;

    SAMSequenceDictionary dictionary;

    final Map<SAMSequenceRecord, Integer> contigToIndexMap = new HashMap<>();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "")
    private File output;

    @Override
    public void onTraversalStart(){
        dictionary = getBestAvailableSequenceDictionary();
        final AtomicInteger index = new AtomicInteger();
        dictionary.getSequences().stream().filter(s -> isRegularContig(s))
                .forEach(s -> contigToIndexMap.put(s, index.getAndIncrement()));
        chimeraMap = new Array2DRowRealMatrix(contigToIndexMap.size(), contigToIndexMap.size());
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        // TODO: how should I handle unmapped reads?
        if (read.isUnmapped()) { return; }
        if (read.mateIsUnmapped()) { return; }
        if (read.getContig().equals(read.getMateContig()) && read.getFragmentLength() < MINIMUM_CHIMERA_FRAGMENT_SIZE){ return; }
        if (! isRegularContig(dictionary.getSequence(read.getContig()))){ return; }
        if (! isRegularContig(dictionary.getSequence(read.getMateContig()))) { return; }

        final int fromIndex = contigToIndexMap.get(dictionary.getSequence(read.getContig()));
        final int toIndex = contigToIndexMap.get(dictionary.getSequence(read.getMateContig()));

        chimeraMap.addToEntry(fromIndex, toIndex, 1.0);
    }

    @Override
    public Object onTraversalSuccess(){
        try (PrintWriter writer = new PrintWriter(output.getAbsolutePath())){
            for (int i = 0; i < chimeraMap.getRowDimension(); i++){
                writer.println(Arrays.toString(chimeraMap.getRow(i).clone()).replace("[","").replace("]", ""));
            }
        } catch (IOException e){
            throw new UserException("Encountered an exception writing to " + output.getAbsolutePath(), e);
        }

        return("DONE");
    }

    private boolean isRegularContig(final SAMSequenceRecord sequenceRecord){
        final List<String> unwantedKeyWords = Arrays.asList("HLA", "_", "EBV");

        return unwantedKeyWords.stream().noneMatch(w -> sequenceRecord.getSequenceName().contains(w));
    }
}
