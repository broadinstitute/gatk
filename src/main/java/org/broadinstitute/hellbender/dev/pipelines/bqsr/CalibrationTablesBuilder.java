
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
 */
public class CalibrationTablesBuilder {

    private final BaseRecalibratorUprooted br;
    private final SAMFileHeader header;
    private final ReferenceDataSource reference;
    private boolean done = false;

    public CalibrationTablesBuilder(SAMFileHeader header, String referenceFileName, String toolArgs) {
        // 1. we'll paste the header onto the new SAMRecord objects, the code needs it.
        this.header = header;
        // 2. the reference we copied
        File refFile = new File(referenceFileName);
        reference = new ReferenceDataSource(refFile);
        // 3. create the class that'll do the actual work
        br = BaseRecalibratorUprooted.fromArgs(header, toolArgs, System.out);
        if (null==br) throw new RuntimeException("invalid tool arguments");
        br.onTraversalStart(refFile);
    }

    /**
     * call this as many times as you want.
     * The skipIntervals have to be 1-based, close-ended.
     */
    public void add(Iterable<Read> reads, List<SimpleInterval> skipIntervals) {
        if (done) throw new GATKException("Can't call add after done");

        for (Read r : reads) {
            SAMRecord sr = GenomicsConverter.makeSAMRecord(r, header);
            final SimpleInterval readInterval = sr.getReadUnmappedFlag() ? null : new SimpleInterval(sr);
            // TODO: this could probably be sped up by taking advantage of a sorted order.
            List<SimpleInterval> knownSitesOverlappingReadInterval =
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
     * Call done() first.
     */
    public RecalibrationTables getRecalibrationTables() {
        if (!done) throw new GATKException("Can't call getRecalibrationTables before done");
        return br.getTables();
    }

    /**
     * Call done() first.
     */
    public Covariate[] getRequestedCovariates() {
        if (!done) throw new GATKException("Can't call getRequestedCovariates before done");
        return br.getRequestedCovariates();
    }

}
