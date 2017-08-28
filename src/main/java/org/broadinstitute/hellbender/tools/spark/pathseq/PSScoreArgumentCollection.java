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

    @Argument(doc = "Alignment identity score threshold, as a fraction of the read length (between 0 and 1).",
            fullName = "scoreMinIdentity",
            minValue = 0.0,
            maxValue = 1.0,
            optional = true)
    public double minIdentity = 0.90;

    @Argument(doc = "Identity margin, as a fraction of the best hit (between 0 and 1). Alignments will be counted when " +
            "their number of matches is at least (1 - identityMargin) x 100% of the best hit.",
            fullName = "identityMargin",
            minValue = 0.0,
            maxValue = 1.0,
            optional = true)
    public double identityMargin = 0.02;

    @Argument(doc = "Divide abundance scores by each taxon's reference genome length (in millions)",
            fullName = "divideByGenomeLength",
            optional = true)
    public boolean divideByGenomeLength = false;

    @Argument(doc = "If true, normalized abundance scores will be reported as a percentage within their kingdom.",
            fullName = "normalizeByKingdom",
            optional = true)
    public boolean normalizeByKingdom = true;

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