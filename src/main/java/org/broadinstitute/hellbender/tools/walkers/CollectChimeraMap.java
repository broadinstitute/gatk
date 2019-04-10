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
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Chimera",
        oneLineSummary = "Chimera",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
public class CollectChimeraMap extends ReadWalker {
    static final int MINIMUM_CHIMERA_FRAGMENT_SIZE = 2000;

    RealMatrix chimeraMap;
    SAMSequenceDictionary dictionary;
    List<SAMSequenceRecord> contigs;
    final Map<SAMSequenceRecord, Integer> contigToIndexMap = new HashMap<>();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "")
    private File output;

    @Argument(fullName = "bam", doc = "bam of chimeric reads")
    public String bam;
    private SAMFileGATKReadWriter bamWriter;

    @Argument(fullName = "normalize", doc = "normalize counts?")
    private boolean normalizeCounts;

    @Override
    public void onTraversalStart(){
        dictionary = getBestAvailableSequenceDictionary();
        final AtomicInteger index = new AtomicInteger();
        contigs = dictionary.getSequences().stream().filter(s -> isRegularContig(s)).collect(Collectors.toList());
        contigs.forEach(s -> contigToIndexMap.put(s, index.getAndIncrement()));
        chimeraMap = new Array2DRowRealMatrix(contigs.size(), contigs.size());
        bamWriter = createSAMWriter(IOUtils.getPath(bam), true);
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (read.isUnmapped()) { return; }
        if (read.mateIsUnmapped()) { return; }
        if (! isRegularContig(dictionary.getSequence(read.getContig()))){ return; }
        if (! isRegularContig(dictionary.getSequence(read.getMateContig()))) { return; }
        if (read.getContig().equals(read.getMateContig()) && read.getFragmentLength() < MINIMUM_CHIMERA_FRAGMENT_SIZE){ return; }

        final int fromIndex = contigToIndexMap.get(dictionary.getSequence(read.getContig()));
        final int toIndex = contigToIndexMap.get(dictionary.getSequence(read.getMateContig()));

        chimeraMap.addToEntry(fromIndex, toIndex, 1.0);
        bamWriter.addRead(read);
    }

    @Override
    public Object onTraversalSuccess(){
        try (PrintWriter writer = new PrintWriter(output.getAbsolutePath())){
            // ### WRITING HEADER ### //
            StringBuilder header = new StringBuilder();
            for (int i = 0; i < contigs.size(); i++){
                final String delimiter = i == contigs.size() - 1 ? "" : "\t";
                header.append(contigs.get(i).getSequenceName() + delimiter);
            }
            writer.println(header.toString());
            // ### END WRITING HEADER ### //

            for (int i = 0; i < chimeraMap.getRowDimension(); i++){
                final long sumContigSizes = contigs.stream().mapToLong(c -> c.getSequenceLength()).sum();
                final StringBuilder row = new StringBuilder();
                for (int j = 0; j < chimeraMap.getColumnDimension(); j++){
                    final String delimiter = j == chimeraMap.getColumnDimension() - 1 ? "" : "\t";
                    if (normalizeCounts){
                        final double fromContigNormalization = (double) contigs.get(i).getSequenceLength() / sumContigSizes;
                        final double toContigNormalization = (double) contigs.get(j).getSequenceLength() / sumContigSizes;
                        final double normalization = fromContigNormalization * toContigNormalization;
                        row.append(chimeraMap.getEntry(i,j)*normalization + delimiter);
                    } else {
                        row.append(chimeraMap.getEntry(i,j) + delimiter);
                    }
                }
                writer.println(row.toString());
            }
        } catch (IOException e){
            throw new UserException("Encountered an exception writing to " + output.getAbsolutePath(), e);
        }

        return("DONE");
    }

    @Override
    public void closeTool() {
        if ( bamWriter != null ) {
            bamWriter.close();
        }
    }

    private boolean isRegularContig(final SAMSequenceRecord sequenceRecord){
        final List<String> unwantedKeyWords = Arrays.asList("HLA", "_", "EBV");

        return unwantedKeyWords.stream().noneMatch(w -> sequenceRecord.getSequenceName().contains(w));
    }
}
