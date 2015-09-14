package org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;

/**
 * Functions to load a BaseRecalOutput.
 * Either from GCS on the worker, or from local on the client.
 */
public final class BaseRecalOutputSource implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Recalibration report on GCS/HDFS -> PCollection of a single BaseRecalOutput.
     * The loading is done at the worker.
     *
     * @param pipeline the pipeline, with authentication information.
     * @param GCSFileName the path to the recalibration report. Must start with "gs://"
     */
    static public PCollection<BaseRecalOutput> of(final Pipeline pipeline, String GCSFileName) {
        return pipeline.apply("calibration report name", Create.of(GCSFileName))
                .apply(ParDo.of(new DoFn<String, BaseRecalOutput>() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public void processElement(ProcessContext c) {
                        final String fname = c.element();
                        File dest = IOUtils.createTempFile("temp-BaseRecal-", ".tmp");
                        try {
                            BucketUtils.copyFile(fname, c.getPipelineOptions(), dest.getPath());
                        } catch (IOException x) {
                            throw new GATKException("Unable to download recalibration table from '" + fname + "'.", x);
                        }
                        c.output(new BaseRecalOutput(dest));
                    }

                }).named("ingest calibration report"));
    }

    /**
     * Recalibration report -> PCollection of a single BaseRecalOutput.
     *
     * If the recalibration report is on remote storage, then the loading will be done at the worker.
     * Otherwise, it'll be done as this method is called (presumably that's on the client).
     *
     * @param pipeline the pipeline, with authentication information.
     * @param path the Recalibration report
     */
    public static PCollection<BaseRecalOutput> loadFileOrRemote(final Pipeline pipeline, String path) {
        if (BucketUtils.isRemoteStorageUrl(path)) {
            return BaseRecalOutputSource.of(pipeline, path);
        } else{
            final BaseRecalOutput recalInfo = new BaseRecalOutput(new File(path));
            return pipeline.apply("recal_file ingest", Create.of(recalInfo));
        }
    }

}
