package org.broadinstitute.hellbender.tools.spark.transforms.bqsr;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * A lightweight wrapper over BaseRecalibrationEngine to make it easier to use from Spark.
 * Takes in reads + contextual data (overlapping reference bases and variants), spits out RecalibrationTables.
 */
public final class BaseRecalibratorEngineSparkWrapper implements Serializable {
    private static final long serialVersionUID = 1L;

    private Broadcast<SAMFileHeader> headerBcast;
    private transient SAMFileHeader header;
    private Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBcast;
    private transient SAMSequenceDictionary referenceSequenceDictionary;

    private RecalibrationArgumentCollection recalArgs;

    private BaseRecalibrationEngine recalibrationEngine;

    // true when we are about to process the first element of a bundle.
    private transient boolean initialized = false;


    // -------------------------------------------------------------
    // public functions (and constructors, those go first)

    /**
     * Takes in reads + contextual data (overlapping reference bases and variants), spits out RecalibrationTables.
     */
    public BaseRecalibratorEngineSparkWrapper(Broadcast<SAMFileHeader> headerBcast, Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBcast, RecalibrationArgumentCollection recalArgs) {
        this.headerBcast = headerBcast;
        this.referenceSequenceDictionaryBcast = referenceSequenceDictionaryBcast;
        this.recalArgs = recalArgs;
        if (this.recalArgs.FORCE_PLATFORM != null) {
            this.recalArgs.DEFAULT_PLATFORM = this.recalArgs.FORCE_PLATFORM;
        }
    }

    // saves to output
    public static void saveTextualReport(String output, SAMFileHeader header, RecalibrationTables rt, RecalibrationArgumentCollection recalArgs) throws IOException {
        OutputStream oStream = BucketUtils.createFile(output);
        QuantizationInfo qi = new QuantizationInfo(rt, recalArgs.QUANTIZING_LEVELS);
        if (recalArgs.FORCE_PLATFORM != null) {
            recalArgs.DEFAULT_PLATFORM = recalArgs.FORCE_PLATFORM;
        }
        StandardCovariateList covariates = new StandardCovariateList(recalArgs, header);
        try ( PrintStream reportStream = new PrintStream(oStream) ) {
            RecalUtils.outputRecalibrationReport(reportStream, recalArgs, qi, rt, covariates);
        }
    }

    public Iterator<RecalibrationTables> apply(Iterator<ContextShard> shards) throws Exception {
        this.header = headerBcast.value();
        this.referenceSequenceDictionary = referenceSequenceDictionaryBcast.value();
        recalibrationEngine = new BaseRecalibrationEngine(recalArgs, header);
        while (shards.hasNext()) {
            ContextShard shard = shards.next();
            for (int i=0; i<shard.reads.size(); i++) {
                final GATKRead read = shard.reads.get(i);
                // Reads are shipped without the header -- put it back in
                ReadUtils.restoreHeaderIfNecessary(read, header);

                final ReadContextData rc = shard.readContext.get(i);
                final Iterable<GATKVariant> variants = rc.getOverlappingVariants();
                final ReferenceBases refBases = rc.getOverlappingReferenceBases();
                final ReferenceDataSource refDS = new ReferenceMemorySource(refBases, referenceSequenceDictionary);
                recalibrationEngine.processRead(read, refDS, variants);
            }
        }
        ArrayList<RecalibrationTables> ret = new ArrayList<>();
        ret.add(recalibrationEngine.getRecalibrationTables());
        return ret.iterator();
    }

}
