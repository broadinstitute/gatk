package org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr;

import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.dataflow.DoFnWLog;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.baq.BAQ;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

/**
 * DoFn for BaseRecalibrator.
 * Takes in reads + contextual data (overlapping reference bases and variants), spits out RecalibrationTables.
 */
public final class BaseRecalibratorFn extends DoFnWLog<KV<GATKRead, ReadContextData>, RecalibrationTables> {
    private static final long serialVersionUID = 1L;

    private PCollectionView<SAMFileHeader> headerView;
    private transient SAMFileHeader header;
    private PCollectionView<SAMSequenceDictionary> referenceSequenceDictionaryView;
    private transient SAMSequenceDictionary referenceSequenceDictionary;

    private RecalibrationArgumentCollection recalArgs;

    private BaseRecalibrationEngine recalibrationEngine;

    // true when we are about to process the first element of a bundle.
    private boolean firstInBundle = true;
    private long nReadsProcessed = 0;

    // -------------------------------------------------------------
    // public functions (and constructors, those go first)

    /**
     * DoFn for BaseRecalibrator.
     * Takes in reads + contextual data (overlapping reference bases and variants), spits out RecalibrationTables.
     */
    public BaseRecalibratorFn(PCollectionView<SAMFileHeader> headerView, PCollectionView<SAMSequenceDictionary> referenceSequenceDictionaryView, RecalibrationArgumentCollection recalArgs) {
        super("BaseRecalibratorFn");
        this.headerView = headerView;
        this.referenceSequenceDictionaryView = referenceSequenceDictionaryView;
        this.recalArgs = recalArgs;
    }

    // This ctor leaves referenceSequenceDictionaryView uninitialized. Caveat emptor.
    private BaseRecalibratorFn(SAMFileHeader header, RecalibrationArgumentCollection recalArgs) {
        super("BaseRecalibratorFn");
        this.header = header;
        this.recalArgs = recalArgs;
    }

    // saves to BRAC.BRAC.RAC.RECAL_TABLE_FILE.getName()
    public static void SaveTextualReport(File output, SAMFileHeader header, RecalibrationTables rt, RecalibrationArgumentCollection recalArgs) throws IOException {
        QuantizationInfo qi = new QuantizationInfo(rt, recalArgs.QUANTIZING_LEVELS);
        BaseRecalibratorFn worker = new BaseRecalibratorFn(header, recalArgs);
        worker.onTraversalStart();
        worker.saveReport(output, rt, qi, worker.recalibrationEngine.getCovariates());
    }

    @Override
    public void startBundle(Context c) throws Exception {
        super.startBundle(c);
        // can't set things up here because we don't have our side input yet.
        // So instead we're doing it in processElement.
        nReadsProcessed = 0;
        // resetting firstInBundle in case we have multiple bundles via the same instance.
        firstInBundle = true;
    }

    @Override
    public void processElement(ProcessContext c) throws Exception {

        // can't set things up in startBundle because it doesn't have our side input yet.
        // So instead we're doing it here in processElement.
        if (firstInBundle) {
            init(c);
            onTraversalStart();
            firstInBundle = false;
        }
        GATKRead r = c.element().getKey();
        ReadContextData rc = c.element().getValue();
        Iterable<Variant> variants = rc.getOverlappingVariants();
        final ReferenceBases refBases = rc.getOverlappingReferenceBases();
        ReferenceDataSource refDS = new ReferenceMemorySource(refBases, referenceSequenceDictionary);

        recalibrationEngine.processRead(r, refDS, variants);
        nReadsProcessed++;
    }

    @Override
    public void finishBundle(Context c) throws Exception {
        bunny.stepEnd("finishBundle");
        if (null!=recalibrationEngine) {
            final RecalibrationTables recalibrationTables = recalibrationEngine.getRecalibrationTables();
            if (null==recalibrationTables) {
                throw new GATKException.ShouldNeverReachHereException("null recalibrationTables");
            }
            c.output(recalibrationTables);
        } else {
            // this can happen, but only if we processed zero read in this bundle
            if (nReadsProcessed!=0) {
                throw new GATKException.ShouldNeverReachHereException("null==recalibrationEngine yet nReadsProcessed!=0");
            }
        }
        bunny.stepEnd("Processed "+nReadsProcessed+" reads.");
        super.finishBundle(c);
    }

    // -------------------------------------------------------------
    // helper functions (non-public)

    private void init(ProcessContext c) {
        header = c.sideInput(headerView);
        referenceSequenceDictionary = c.sideInput(referenceSequenceDictionaryView);
    }

    private void onTraversalStart() {
        if (recalArgs.FORCE_PLATFORM != null)
            recalArgs.DEFAULT_PLATFORM = recalArgs.FORCE_PLATFORM;

        recalibrationEngine = new BaseRecalibrationEngine(recalArgs, header);
    }

    private void saveReport( File output, RecalibrationTables rt, QuantizationInfo quantizationInfo, StandardCovariateList finalRequestedCovariates ) throws FileNotFoundException {
        try ( PrintStream reportStream = new PrintStream(output) ) {
            RecalUtils.outputRecalibrationReport(reportStream, recalArgs, quantizationInfo, rt, finalRequestedCovariates, recalArgs.SORT_BY_ALL_COLUMNS);
        }
    }
}
