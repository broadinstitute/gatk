
package org.broadinstitute.hellbender.dev.pipelines.bqsr;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.genomics.gatk.common.GenomicsConverter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.dev.tools.walkers.bqsr.BaseRecalibratorUprooted;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.recalibration.QuantizationInfo;
import org.broadinstitute.hellbender.tools.recalibration.RecalUtils;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.tools.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.tools.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

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
 * return foo.getRecalibrationTables();
 *
 */
public final class CalibrationTablesBuilder {

    private final BaseRecalibratorUprooted br;
    private final SAMFileHeader header;
    private final ReferenceDataSource reference;
    private boolean done = false;

    public CalibrationTablesBuilder(final SAMFileHeader header, final String referenceFileName, final String toolArgs) {
        // 1. we'll paste the header onto the new SAMRecord objects, the code needs it.
        this.header = header;
        // 2. the reference we copied
        final File refFile = new File(referenceFileName);
        this.reference = new ReferenceDataSource(refFile);
        // 3. create the class that'll do the actual work
        this.br = BaseRecalibratorUprooted.fromArgs(header, toolArgs, System.out);
        if (null==br) throw new RuntimeException("invalid tool arguments");
        br.onTraversalStart(refFile);
    }

    /**
     * For each read at this locus, computes the various covariate values and updates the RecalibrationTables as needed.
     * The "skipIntervals" list doesn't need to be complete but it must includes all those that overlap with any read
     * in the "reads" argument.
     *
     * The skipIntervals have to be 1-based, close-ended.
     */
    public void add(final Iterable<Read> reads, final List<SimpleInterval> skipIntervals) {
        if (done) throw new GATKException("Can't call add after done");

        for (final Read r : reads) {
            final SAMRecord sr = GenomicsConverter.makeSAMRecord(r, header);
            final SimpleInterval readInterval = sr.getReadUnmappedFlag() ? null : new SimpleInterval(sr);
            // TODO: this could probably be sped up by taking advantage of a sorted order.
            final List<SimpleInterval> knownSitesOverlappingReadInterval =
                    skipIntervals.stream().filter(x -> x.overlaps(readInterval))
                            .collect(Collectors.toList());
            br.apply(sr, new ReferenceContext(reference, readInterval), knownSitesOverlappingReadInterval);
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
