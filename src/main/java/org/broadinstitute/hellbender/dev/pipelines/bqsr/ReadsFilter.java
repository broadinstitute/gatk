package org.broadinstitute.hellbender.dev.pipelines.bqsr;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionTuple;
import com.google.cloud.genomics.gatk.common.GenomicsConverter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;

import java.io.Serializable;

/**
 * A simple filter for reads.
 */
public final class ReadsFilter extends PTransform<PCollection<Read>, PCollection<Read>> implements Serializable {
    private static final long serialVersionUID = 1L;
    
    private final ReadFilter readFilter;
    private final SAMFileHeader header;

    /**
     * ReadsFilter will use the given header to interpret each read, and then let the Read through
     * if the passed filter accepts it.
     */
    public ReadsFilter(final ReadFilter filterToApply, final SAMFileHeader header) {
        if (null==filterToApply) {
            throw new GATKException("Missing argument: filterToApply");
        }
        if (null==header) {
            throw new GATKException("Missing argument: header");
        }
        this.readFilter = filterToApply;
        this.header = header;
    }

    /**
     * Filter out reads we don't want.
     */
    @Override
    public PCollection<Read> apply(PCollection<Read> in) {
        return in.apply(ParDo
                .named(getName())
                .of(new DoFn<Read, Read>() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public void processElement(DoFn<Read,Read>.ProcessContext c) throws Exception {
                        Read r = c.element();
                        // TODO: remove this conversion once we switch to a Read interface
                        final SAMRecord sr = GenomicsConverter.makeSAMRecord(r, header);
                        if (readFilter.test(sr)) {
                            c.output(r);
                        }
                    }
                }));
    }


}
