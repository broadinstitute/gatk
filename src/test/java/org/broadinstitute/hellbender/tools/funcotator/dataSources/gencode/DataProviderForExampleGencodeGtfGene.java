package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import htsjdk.tribble.annotation.Strand;
import org.broadinstitute.hellbender.utils.codecs.gencode.*;

import java.util.Collections;

/**
 * A class to hold a method that provides a valid {@link GencodeGtfGeneFeature} object for testing.
 * Created by jonn on 10/3/17.
 */
public class DataProviderForExampleGencodeGtfGene {

    public static GencodeGtfGeneFeature createGencodeGtfGeneFeature() {

        // Create the Features as they exist in the test file:
        GencodeGtfFeatureBaseData data;

        data = new GencodeGtfFeatureBaseData(1, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.GENE,
                1, 3000, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", null, GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", null, null, null, -1, null, GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED, null, null);
        final GencodeGtfGeneFeature gene = (GencodeGtfGeneFeature)GencodeGtfFeature.create(data);

        // ======================

        data = new GencodeGtfFeatureBaseData(2, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.TRANSCRIPT,
                1, 1000, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", -1, null, GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfTranscriptFeature transcript1 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(3, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                1, 200, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", 1, "TEST_EXON1", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon1 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(4, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.START_CODON,
                99, 101, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", 1, "TEST_EXON1", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfStartCodonFeature startCodon1 = (GencodeGtfStartCodonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(5, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                201, 400, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", 2, "TEST_EXON2", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon2 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(6, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.CDS,
                201, 400, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", 2, "TEST_EXON2", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfCDSFeature cds1 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(7, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                401, 600, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", 3, "TEST_EXON3", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon3 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(8, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.CDS,
                401, 600, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", 3, "TEST_EXON3", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfCDSFeature cds2 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(9, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                601, 800, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", 4, "TEST_EXON4", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon4 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(10, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.CDS,
                601, 800, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", 4, "TEST_EXON4", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfCDSFeature cds3 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(data);


        data = new GencodeGtfFeatureBaseData(11, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                801, 1000, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", 5, "TEST_EXON5", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon5 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(12, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.STOP_CODON,
                900, 902, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", 5, "TEST_EXON5", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfStopCodonFeature stopCodon1 = (GencodeGtfStopCodonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(13, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.UTR,
                1, 98, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", 1, "TEST_EXON1", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfUTRFeature utr1 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(14, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.UTR,
                903, 1000, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT1", 5, "TEST_EXON5", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfUTRFeature utr2 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(data);

        // ======================

        data = new GencodeGtfFeatureBaseData(15, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.TRANSCRIPT,
                1001, 2000, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", -1, null, GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfTranscriptFeature transcript2 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(16, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                1001, 1200, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 1, "TEST_EXON1", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon6 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(17, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.START_CODON,
                1099, 1101, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 1, "TEST_EXON1", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfStartCodonFeature startCodon2 = (GencodeGtfStartCodonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(18, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.CDS,
                1099, 1200, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 1, "TEST_EXON1", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfCDSFeature cds4 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(19, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                1201, 1400, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 2, "TEST_EXON2", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon7 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(20, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.CDS,
                1201, 1400, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 2, "TEST_EXON2", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfCDSFeature cds5 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(21, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                1401, 1600, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 3, "TEST_EXON3", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon8 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(22, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.CDS,
                1401, 1600, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 3, "TEST_EXON3", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfCDSFeature cds6 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(23, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                1601, 1800, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 4, "TEST_EXON4", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon9 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(24, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.CDS,
                1601, 1800, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 4, "TEST_EXON4", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfCDSFeature cds7 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(data);


        data = new GencodeGtfFeatureBaseData(25, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                1801, 2000, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 5, "TEST_EXON5", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon10 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(26, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.CDS,
                1801, 1899, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 5, "TEST_EXON5", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfCDSFeature cds8 = (GencodeGtfCDSFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(27, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.STOP_CODON,
                1900, 1902, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 5, "TEST_EXON5", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfStopCodonFeature stopCodon2 = (GencodeGtfStopCodonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(28, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.UTR,
                1001, 1098, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 1, "TEST_EXON1", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfUTRFeature utr3 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(29, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.UTR,
                1903, 2000, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT2", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT2", 5, "TEST_EXON5", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfUTRFeature utr4 = (GencodeGtfUTRFeature) GencodeGtfFeature.create(data);
        
        // ======================

        data = new GencodeGtfFeatureBaseData(30, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.TRANSCRIPT,
                2001, 3000, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT3", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT3", -1, null, GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfTranscriptFeature transcript3 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(31, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                2001, 2200, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT3", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT3", 1, "TEST_EXON1", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon11 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(32, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                2201, 2400, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT3", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT3", 2, "TEST_EXON2", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon12 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(33, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                2401, 2600, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT3", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT3", 3, "TEST_EXON3", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon13 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(34, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                2601, 2800, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT3", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT3", 4, "TEST_EXON4", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon14 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);

        data = new GencodeGtfFeatureBaseData(35, "chr1", GencodeGtfFeature.AnnotationSource.ENSEMBL, GencodeGtfFeature.FeatureType.EXON,
                2801, 3000, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT3", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING,
                null, "TEST_GENE", GencodeGtfFeature.GeneTranscriptType.PROTEIN_CODING, null, "TEST_TRANSCRIPT3", 5, "TEST_EXON5", GencodeGtfFeature.LocusLevel.AUTOMATICALLY_ANNOTATED,
                Collections.emptyList(),
                null
        );
        final GencodeGtfExonFeature exon15 = (GencodeGtfExonFeature) GencodeGtfFeature.create(data);
        
        // ======================
        // Here we'll set the version number of each feature.
        // Normally this is set in the decode method, so we need to do it manually.

        gene.setUcscGenomeVersion("hg38");

        transcript1.setUcscGenomeVersion("hg38");
        exon1.setUcscGenomeVersion("hg38");
        exon2.setUcscGenomeVersion("hg38");
        exon3.setUcscGenomeVersion("hg38");
        exon4.setUcscGenomeVersion("hg38");
        exon5.setUcscGenomeVersion("hg38");
        cds1.setUcscGenomeVersion("hg38");
        cds2.setUcscGenomeVersion("hg38");
        cds3.setUcscGenomeVersion("hg38");
        startCodon1.setUcscGenomeVersion("hg38");
        stopCodon1.setUcscGenomeVersion("hg38");
        utr1.setUcscGenomeVersion("hg38");
        utr2.setUcscGenomeVersion("hg38");

        transcript1.setUcscGenomeVersion("hg38");
        exon6.setUcscGenomeVersion("hg38");
        exon7.setUcscGenomeVersion("hg38");
        exon8.setUcscGenomeVersion("hg38");
        exon9.setUcscGenomeVersion("hg38");
        exon10.setUcscGenomeVersion("hg38");
        cds4.setUcscGenomeVersion("hg38");
        cds5.setUcscGenomeVersion("hg38");
        cds6.setUcscGenomeVersion("hg38");
        cds7.setUcscGenomeVersion("hg38");
        cds8.setUcscGenomeVersion("hg38");
        startCodon2.setUcscGenomeVersion("hg38");
        stopCodon2.setUcscGenomeVersion("hg38");
        utr3.setUcscGenomeVersion("hg38");
        utr4.setUcscGenomeVersion("hg38");

        transcript3.setUcscGenomeVersion("hg38");
        exon11.setUcscGenomeVersion("hg38");
        exon12.setUcscGenomeVersion("hg38");
        exon13.setUcscGenomeVersion("hg38");
        exon14.setUcscGenomeVersion("hg38");
        exon15.setUcscGenomeVersion("hg38");

        // ======================

        // Aggregate the Features as they should be:
        exon1.setStartCodon(startCodon1);
        exon2.setCds(cds1);
        exon3.setCds(cds2);
        exon4.setCds(cds3);
        exon5.setStopCodon(stopCodon1);

        transcript1.addExon(exon1);
        transcript1.addExon(exon2);
        transcript1.addExon(exon3);
        transcript1.addExon(exon4);
        transcript1.addExon(exon5);
        transcript1.addUtr(utr1);
        transcript1.addUtr(utr2);

        gene.addTranscript(transcript1);

        // ======================

        exon6.setStartCodon(startCodon2);
        exon6.setCds(cds4);
        
        exon7.setCds(cds5);
        exon8.setCds(cds6);
        exon9.setCds(cds7);

        exon10.setStopCodon(stopCodon2);
        exon10.setCds(cds8);

        transcript2.addExon(exon6);
        transcript2.addExon(exon7);
        transcript2.addExon(exon8);
        transcript2.addExon(exon9);
        transcript2.addExon(exon10);
        transcript2.addUtr(utr3);
        transcript2.addUtr(utr4);
        gene.addTranscript(transcript2);

        // ======================

        transcript3.addExon(exon11);
        transcript3.addExon(exon12);
        transcript3.addExon(exon13);
        transcript3.addExon(exon14);
        transcript3.addExon(exon15);
        gene.addTranscript(transcript3);

        // ======================

        return gene;
    }


}
