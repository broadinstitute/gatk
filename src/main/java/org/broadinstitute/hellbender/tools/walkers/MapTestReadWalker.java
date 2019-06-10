package org.broadinstitute.hellbender.tools.walkers;


import htsjdk.samtools.util.OverlapDetector;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.CollectReadCounts;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.*;

@CommandLineProgramProperties(
        summary = "TODO",
        oneLineSummary = "TODO",
        programGroup = VariantEvaluationProgramGroup.class
)
public class MapTestReadWalker extends ReadWalker {
    private static final Logger logger = LogManager.getLogger(MapTestReadWalker.class);
    static final String BWA_INDEX_PARAM_FULL_NAME = "bwa-index";
    private static final String MAP_CT_ALPHA_FULL_NAME = "map-ct-alpha";
    private static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 30;

    private Map<SimpleInterval, Integer> preTweakedCounts = new LinkedHashMap<>();
    private Map<SimpleInterval, Integer> postTweakedCounts = new LinkedHashMap<>();
    private OverlapDetector<SimpleInterval> overlapDetector;
    private BwaMemAligner bwaMemAligner;
    private int numTweakedReadsFailMQ = 0;

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output file ")
    public File outputFile;

    @Argument(
            fullName  = MAP_CT_ALPHA_FULL_NAME,
            doc = "alpha")
    public double alpha = 1e-5;

    @Argument(doc = "BWA index image.  Generated from `gatk BwaMemIndexImageCreator -I <fasta_file>`", fullName = BWA_INDEX_PARAM_FULL_NAME)
    private File bwaIndexFile;

    private int readsProcessed = 0;
    private int readsTweaked = 0;

    private final List<GATKRead> readBatch = new ArrayList<>();
    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        readBatch.add(read);
        if (readBatch.size() >= 200000) {
            processReads(readBatch);
            readBatch.clear();
        }
    }

    private void processReads(final List<GATKRead> reads) {
        boolean wasReadTweaked = false;
        final List<byte[]> tweakedReads = new ArrayList<>();
        final Map<byte[], GATKRead> tweakedReadToOriginalRead = new HashMap<>();
        for (final GATKRead read : reads) {
            final byte[] bases = read.getBasesNoCopy();
            final byte[] tweakedbases = new byte[bases.length];
            for (int i = 0; i < bases.length; i++) {
                final double randomSample = Math.random();
                if ((randomSample < alpha) && (bases[i] == BaseUtils.Base.C.base)) {
                    tweakedbases[i] = BaseUtils.Base.T.base;
                    wasReadTweaked = true;
                } else if ((randomSample < alpha) && (bases[i] == BaseUtils.Base.G.base)) {
                    tweakedbases[i] = BaseUtils.Base.A.base;
                    wasReadTweaked = true;
                } else {
                    tweakedbases[i] = bases[i];
                }

            }

            readsProcessed++;

            // Do nothing if the read was not tweaked
            if (wasReadTweaked) {
                tweakedReads.add(tweakedbases);
                tweakedReadToOriginalRead.put(tweakedbases, read);
                readsTweaked++;
            }
        }

        final List<List<BwaMemAlignment>> bwaAlignments = bwaMemAligner.alignSeqs(tweakedReads);

        for (int i = 0; i < bwaAlignments.size(); i++) {
            // For each set of alignments (we know that the first index will be length 1 since we only give one aligned read at a time).
            final List<BwaMemAlignment> bwaAlignmentsForThisRead = bwaAlignments.get(i);

            // Check if the top result for the tweaked read is below MQ30.  If so, score it.
            final BwaMemAlignment firstBwaMemAlignment = bwaAlignmentsForThisRead.get(0);
            if (firstBwaMemAlignment.getMapQual() < 30) {
                numTweakedReadsFailMQ++;
                return;
            }
            final String contig = getBestAvailableSequenceDictionary().getSequences().get(firstBwaMemAlignment.getRefId()).getSequenceName();

            final SimpleInterval postInterval = incrementCountMap(firstBwaMemAlignment, contig, postTweakedCounts);
            final SimpleInterval preInterval = incrementCountMap(tweakedReadToOriginalRead.get(tweakedReads.get(i)), preTweakedCounts);
            if (!preInterval.equals(postInterval)) {
                logger.info("Difference!");
            }
        }
        bwaMemAligner.close();
    }

    private SimpleInterval incrementCountMap(final BwaMemAlignment bwaAlignment, final String contig, final Map<SimpleInterval, Integer> intervalToCountMap) {
        final Set<SimpleInterval> overlaps = overlapDetector.getOverlaps(new SimpleInterval(contig, bwaAlignment.getRefStart(), bwaAlignment.getRefEnd()));
        if (overlaps.iterator().hasNext()) {
            final SimpleInterval simpleInterval = overlaps.iterator().next();
            intervalToCountMap.put(simpleInterval, intervalToCountMap.get(simpleInterval) + 1);
            return simpleInterval;
        } else {
            logger.warn("Could no longer map read to a bin for " + contig + ":" + bwaAlignment.getRefStart() + "-" + bwaAlignment.getRefEnd());
            return null;
        }
    }

    private SimpleInterval incrementCountMap(final GATKRead read, final Map<SimpleInterval, Integer> intervalToCountMap) {
        final Set<SimpleInterval> overlaps = overlapDetector.getOverlaps(new SimpleInterval(read.getAssignedContig(), read.getAssignedStart(), read.getEnd()));
        if (overlaps.iterator().hasNext()) {
            final SimpleInterval simpleInterval = overlaps.iterator().next();
            intervalToCountMap.put(simpleInterval, intervalToCountMap.get(simpleInterval) + 1);
            return simpleInterval;
        } else {
            logger.warn("Could no longer map UNMODIFIED read to a bin for " + read.getAssignedContig() + ":" + read.getStart() + "-" + read.getEnd());
            return null;
        }
    }

    @Override
    public Object onTraversalSuccess() {
        bwaMemAligner.close();
        final String msg1 = "Tweaked " + this.readsTweaked + " reads of " + this.readsProcessed + "\n";
        logger.info(msg1);
        logger.info(numTweakedReadsFailMQ);

        final File postOutputFile = outputFile;
        final File preOutputFile = new File(outputFile.getAbsolutePath() + ".pre.txt");

        writeOutputFile(postOutputFile, postTweakedCounts);
        writeOutputFile(preOutputFile, preTweakedCounts);

        return msg1;
    }

    private void writeOutputFile(final File outputFile, final Map<SimpleInterval, Integer> intervalCountMap) {
        try {
            final FileWriter fileWriter = new FileWriter(outputFile);
            final MapTableWriter writer = new MapTableWriter(fileWriter, new TableColumnCollection(Arrays.asList("CONTIG", "START", "END", "VAL")));
            for (Map.Entry<SimpleInterval, Integer> entry : intervalCountMap.entrySet()) {
                writer.writeRecord(entry);
            }
            writer.close();
        } catch (final IOException ioe) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "Could not create output file.", ioe);
        }
    }

    @Override
    public void onTraversalStart() {
        getTraversalIntervals().forEach(interval -> {
            postTweakedCounts.put(interval, 0);
            preTweakedCounts.put(interval, 0);
        });

        overlapDetector = OverlapDetector.create(getTraversalIntervals());

        bwaMemAligner = new BwaMemAligner(new BwaMemIndex(bwaIndexFile.getAbsolutePath()));
//        bwaMemAligner.setFlagOption(BwaMemAligner.MEM_F_ALL);

    }

    /**
     * Must have this in order to be able to look at the intervals.
     * @return true
     */
    @Override
    public boolean requiresReference(){
        return true;
    }

    @Override
    public boolean requiresIntervals(){
        return true;
    }

    /**
     * Should be same as {@link CollectReadCounts}
     * @return
     */
    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> readFilters = new ArrayList<>(super.getDefaultReadFilters());
        readFilters.add(ReadFilterLibrary.MAPPED);
        readFilters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        readFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        readFilters.add(new MappingQualityReadFilter(DEFAULT_MINIMUM_MAPPING_QUALITY));
        return readFilters;
    }

    private class MapTableWriter extends TableWriter<Map.Entry<SimpleInterval, Integer>> {

        MapTableWriter(final Writer writer, TableColumnCollection tableColumns) throws IOException {
            super(writer, tableColumns);
        }

        @Override
        protected void composeLine(final Map.Entry<SimpleInterval, Integer> record, final DataLine dataLine) {
            // First the Locatable info
            dataLine.set("CONTIG", record.getKey().getContig());
            dataLine.set("START", record.getKey().getStart());
            dataLine.set("END", record.getKey().getEnd());
            dataLine.set("VAL", record.getValue());
        }
    }

}
