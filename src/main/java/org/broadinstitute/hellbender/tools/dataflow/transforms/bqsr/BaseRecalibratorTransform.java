package org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr;

import com.google.cloud.dataflow.sdk.coders.Coder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.dataflow.DoFnWLog;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Base Quality Score Recalibration, phase 1
 * This contains the Dataflow-specific bits. Most of the work is done by BaseRecalibratorFn.
 * The Hellbender-specific bits are in BaseRecalibratorDataflow2.
 */
public final class BaseRecalibratorTransform extends PTransform<PCollection<KV<GATKRead, ReadContextData>>, PCollection<BaseRecalOutput>> {
    private static final long serialVersionUID = 1L;
    private PCollectionView<SAMFileHeader> headerView;
    private PCollectionView<SAMSequenceDictionary> refDictionary;
    RecalibrationArgumentCollection recalArgs;

    public BaseRecalibratorTransform(final PCollectionView<SAMFileHeader> headerView, PCollectionView<SAMSequenceDictionary> refDictionary, RecalibrationArgumentCollection recalArgs) {
        this.headerView = headerView;
        this.refDictionary = refDictionary;
        this.recalArgs = recalArgs;
    }

    @Override
    public PCollection<BaseRecalOutput> apply( PCollection<KV<GATKRead, ReadContextData>> input ) {
        PCollection<RecalibrationTables> oneStatPerWorker =
            input.apply(ParDo.named("BaseRecalibrator")
            .withSideInputs(headerView, refDictionary)
            .of(new BaseRecalibratorFn(headerView, refDictionary, recalArgs)));
        PCollection<RecalibrationTables> aggregate = aggregateStatistics(oneStatPerWorker);
        return addQuantizationInfo(headerView, recalArgs, aggregate);
    }

    @Override
    public Coder<BaseRecalOutput> getDefaultOutputCoder() {
        return SerializableCoder.of(BaseRecalOutput.class);
    }

    /**
     * Merge the statistics from each block. The resulting "collection" contains a single element, with the answer.
     */
    private static PCollection<RecalibrationTables> aggregateStatistics(final PCollection<RecalibrationTables> tables) {
        return tables
            // aggregate
            .apply(Combine.globally(new RecalibrationTablesMerger()))
                // call finalize on the result
            .apply(ParDo
                .named("finalizeRecalTables")
                .of(new DoFnWLog<RecalibrationTables, RecalibrationTables>("finalizeRecalTables") {
                    private static final long serialVersionUID = 1L;

                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        RecalibrationTables tables = c.element();
                        if (null==tables) {
                            // the merger may return null when there are no inputs at all. In that case we don't want to
                            // crash (though it's really an edge case).
                            log.warn("No recalibration tables!");
                        } else {
                            // normal case: recalibrate
                            BaseRecalibrationEngine.finalizeRecalibrationTables(tables);
                        }
                        c.output(tables);
                    }
                }));
    }

    /**
    * addQuantizationInfo takes the computed RecalibrationTable and adds the QuantizationInfo and RequestedCovariates objects.
    * We call this triplet "BaseRecalOutput". It contains everything we need from phase 1 to continue onto phase 2 of BQSR.
    */
    private static PCollection<BaseRecalOutput> addQuantizationInfo(PCollectionView<SAMFileHeader> headerView, RecalibrationArgumentCollection recalArgs, PCollection<RecalibrationTables> recal) {
        return recal.apply(ParDo
            .named("addQuantizationInfo")
            .withSideInputs(headerView)
            .of(new DoFnWLog<RecalibrationTables, BaseRecalOutput>("addQuantizationInfo") {
                private static final long serialVersionUID = 1L;

                @Override
                public void processElement(ProcessContext c) throws IOException {
                    RecalibrationTables rt = c.element();
                    SAMFileHeader header = c.sideInput(headerView);
                    //BaseRecalOutput ret = new BaseRecalOutput(rt, baseRecalibratorWorker.getQuantizationInfo(rt), baseRecalibratorWorker.getRequestedCovariates());
                    // Saving and loading back the report actually changes it. So we have to do it.
                    // TODO(issue#799): Figure out what it changes, and just do that instead of doing the whole rigamarole.
                    File temp = IOUtils.createTempFile("temp-recalibrationtable-", ".tmp");
                    try {
                        BaseRecalibratorFn.SaveTextualReport(temp, header, rt, recalArgs);
                        BaseRecalOutput ret = new BaseRecalOutput(temp);
                        c.output(ret);
                    } catch (FileNotFoundException e) {
                        throw new GATKException("can't find my own temporary file", e);
                    } catch (IOException e) {
                        throw new GATKException("unable to save temporary report to " + temp.getPath(), e);
                    }
                }
            }));
    }

}



