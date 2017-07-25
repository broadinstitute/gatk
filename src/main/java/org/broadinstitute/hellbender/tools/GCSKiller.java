package org.broadinstitute.hellbender.tools;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.*;

@CommandLineProgramProperties(
        summary = "GCSKiller",
        oneLineSummary = "GCSKiller",
        programGroup = TestProgramGroup.class
)
public class GCSKiller extends CommandLineProgram {

    @Argument(fullName = "variant", shortName = "V", doc = "")
    private String vcf;

    @Override
    protected Object doWork() {
        try ( final FeatureDataSource<VariantContext> vcfReader = new FeatureDataSource<VariantContext>(vcf, "variants", 0, VariantContext.class, 0, 0) ) {
            final SAMSequenceDictionary vcfDict = ((VCFHeader)vcfReader.getHeader()).getSequenceDictionary();
            final List<SimpleInterval> intervals = IntervalUtils.getAllIntervalsForReference(vcfDict);
            final List<SimpleInterval> queryIntervals = IntervalUtils.cutToShards(intervals, 1000);
            final List<List<SimpleInterval>> groupedQueryIntervals = groupIntervals(queryIntervals, 1000);

            logger.info("Query intervals: " + queryIntervals.size());
            logger.info("Interval groups: " + groupedQueryIntervals.size());

            final ThreadFactory threadFactory = new ThreadFactoryBuilder().setDaemon(true).build();
            final ExecutorService threadPool = Executors.newFixedThreadPool(10, threadFactory);
            final ExecutorCompletionService<Integer> completionService = new ExecutorCompletionService<>(threadPool);
            final List<Future<Integer>> futures = new ArrayList<>();

            int groupCount = 0;
            for ( final List<SimpleInterval> intervalGroup : groupedQueryIntervals ) {
                ++groupCount;
                futures.add(completionService.submit(new QueryWorker(vcf, intervalGroup), groupCount));
            }

            int completedGroups = 0;
            for ( int i = 0; i < groupedQueryIntervals.size(); ++i ) {
                try {
                    final Future<Integer> result = completionService.take();
                    ++completedGroups;
                    logger.info("Worker " + result.get() + " done. " + completedGroups + "/" + groupedQueryIntervals.size() + " total groups complete");
                }
                catch ( InterruptedException e ) {
                    logger.info("Interrupted");
                    e.printStackTrace();
                    Thread.currentThread().interrupt();
                }
                catch ( ExecutionException e ) {
                    logger.info("Error in query worker");
                    e.printStackTrace();
                }
            }

            threadPool.shutdownNow();
        }

        return null;
    }

    private List<List<SimpleInterval>> groupIntervals( final List<SimpleInterval> intervals, final int groupSize ) {
        final List<List<SimpleInterval>> groups = new ArrayList<>();
        final Iterator<SimpleInterval> intervalIter = intervals.iterator();

        while ( intervalIter.hasNext() ) {
            final List<SimpleInterval> currentGroup = new ArrayList<>();
            while ( currentGroup.size() < groupSize && intervalIter.hasNext() ) {
                currentGroup.add(intervalIter.next());
            }
            groups.add(currentGroup);
        }

        return groups;
    }

    private static final class QueryWorker implements Runnable {
        private final String vcf;
        private final List<SimpleInterval> queryIntervals;

        public QueryWorker(final String vcf, List<SimpleInterval> queryIntervals) {
            this.vcf = vcf;
            this.queryIntervals = queryIntervals;
        }

        @Override
        public void run() {
            try ( final FeatureDataSource<VariantContext> workerReader = new FeatureDataSource<>(vcf, "variants", 0, VariantContext.class, 0, 0) ) {
                for ( final SimpleInterval interval : queryIntervals ) {
                    final Iterator<VariantContext> queryIter = workerReader.query(interval);
                    while ( queryIter.hasNext() ) {
                        queryIter.next();
                    }
                }
            }
        }
    }
}
