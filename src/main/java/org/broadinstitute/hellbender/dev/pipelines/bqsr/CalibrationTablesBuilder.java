
package org.broadinstitute.hellbender.dev.pipelines.bqsr;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.dev.tools.walkers.bqsr.BaseRecalibrationArgumentCollection;
import org.broadinstitute.hellbender.dev.tools.walkers.bqsr.BaseRecalibratorWorker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.tools.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

/*
 * Feed me reads and skipIntervals, and I'll get you a RecalibrationTables and a Covariate[].
 *
 * (OK, yes, I need the reference genome too).
 *
 * This class wraps around BaseRecalibratorUprooted, giving it a cozy environment outside of the normal
 * CommandLine framework. It's also completely independent of Dataflow, so it could be tested outside of that framework.
 *
 * Typical use:
 *
 * foo = new CalibrationTablesBuilder(header,reference,args);
 * foo.add(reads, skipIntervals);
 * foo.add(moreReads, skipIntervals);
 * foo.done();
 * return foo.getRecalibrationTables();
 *
 */
public final class CalibrationTablesBuilder {

    private final BaseRecalibratorWorker br;
    private final SAMFileHeader header;
    private final ReferenceDataSource reference;
    private boolean done = false;

    public CalibrationTablesBuilder(final SAMFileHeader header, final String referenceFileName, final BaseRecalibrationArgumentCollection toolArgs) {
        // 1. we'll paste the header onto the new SAMRecord objects, the code needs it.
        this.header = header;
        // 2. the reference we copied
        final File refFile = new File(referenceFileName);
        this.reference = new ReferenceDataSource(refFile);
        // 3. create the class that'll do the actual work
        this.br = BaseRecalibratorWorker.fromArgs(header, toolArgs);
        if (null==br) throw new RuntimeException("invalid tool arguments");
        br.onTraversalStart(refFile);
    }

    public CalibrationTablesBuilder(final SAMFileHeader header, final String referenceFileName, String commandLineArguments) {
        // 1. we'll paste the header onto the new SAMRecord objects, the code needs it.
        this.header = header;
        // 2. the reference we copied
        final File refFile = new File(referenceFileName);
        this.reference = new ReferenceDataSource(refFile);
        // 3. create the class that'll do the actual work
        this.br = BaseRecalibratorWorker.fromCommandLine(header, commandLineArguments, System.out);
        if (null==br) throw new GATKException("invalid tool arguments");
        br.onTraversalStart(refFile);
    }


    /**
     * For each read at this locus, computes the various covariate values and updates the RecalibrationTables as needed.
     * The "skipIntervals" list doesn't need to be complete but it must includes all those that overlap with any read
     * in the "reads" argument.
     *
     * The skipIntervals have to be 1-based, close-ended.
     */
    public void add(final Iterable<GATKRead> reads, final List<SimpleInterval> skipIntervals) {
        if (done) throw new GATKException("Can't call add after done");

        for (final GATKRead read : reads) {
            final SimpleInterval readInterval = read.isUnmapped() ? null : new SimpleInterval(read);
            // TODO: this could probably be sped up by taking advantage of a sorted order.
            final List<SimpleInterval> knownSitesOverlappingReadInterval =
                    skipIntervals.stream().filter(x -> x.overlaps(readInterval))
                            .collect(Collectors.toList());
            br.apply(read, new ReferenceContext(reference, readInterval), knownSitesOverlappingReadInterval);
        }
    }

    /**
     * "Closes" the object, making it immutable. You can't call add anymore, but you can read the tables.
     */
    public void done() {
        done = true;
    }

    /**
     * Returns the recalibrationTables that were computed as a result of the calls to add.
     * Call done() first.
     */
    public RecalibrationTables getRecalibrationTables() {
        if (!done) throw new GATKException("Can't call getRecalibrationTables before done");
        return br.getTables();
    }

    /**
     * Returns the requestedCovariates computed as a result of the calls to add.
     * Call done() first.
     */
    public StandardCovariateList getRequestedCovariates() {
        if (!done) throw new GATKException("Can't call getRequestedCovariates before done");
        return br.getRequestedCovariates();
    }

}
