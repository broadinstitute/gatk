package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public final class PSScoreArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "URI for the taxonomic scores output",
            fullName = "scoresOutputPath")
    public String scoresPath;

    @Argument(doc = "URI to the reference taxonomy database build using PathSeqBuildReferenceTaxonomy",
            fullName = "taxonomicDatabasePath")
    public String taxonomyDatabasePath;

    @Argument(doc = "Alignment coverage score threshold, as a fraction of the read length (between 0 and 1)",
            fullName = "scoreMinCoverage",
            minValue = 0.0,
            maxValue = 1.0,
            optional = true)
    public double minCoverage = 0.90;

    @Argument(doc = "Alignment identity score threshold, as a fraction of the read length (between 0 and 1)",
            fullName = "scoreMinIdentity",
            minValue = 0.0,
            maxValue = 1.0,
            optional = true)
    public double minIdentity = 0.90;

    @Argument(doc = "Write accessions found in the reads header but not the taxonomy database to this file",
            fullName = "scoreWarningsFile",
            optional = true)
    public String headerWarningFile = null;

    @Argument(doc = "Estimated reads per Spark partition",
            fullName = "scoreReadsPerPartition",
            minValue = 1,
            optional = true)
    public int readsPerPartition = 200000;

}