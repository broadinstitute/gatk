package org.broadinstitute.hellbender.tools.walkers.conversion;

import htsjdk.samtools.util.Interval;
import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfFeature;

import java.util.HashMap;

public class GtfToBed extends FeatureWalker <GencodeGtfFeature>{

    @Argument(shortName = "G", fullName = "gencode_gtf_file", doc = "Gencode GTF file")
    public GATKPath inputFile;

    HashMap<String, GtfInfo> idToInfo = new HashMap<>();

    @Override
    protected boolean isAcceptableFeatureType(Class<? extends Feature> featureType) {
        return featureType.isAssignableFrom(GencodeGtfFeature.class);
    }

    //Apply runs PER line
    @Override
    public void apply(GencodeGtfFeature feature, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        int geneStart = 0;
        int geneEnd = 0;
        String featureType = String.valueOf(feature.getFeatureType());
        if (featureType.equals("gene")) {
            geneStart = feature.getStart();
            geneEnd = feature.getEnd();
            Interval interval = new Interval(feature.getContig(), geneStart, geneEnd);
            GtfInfo gtfInfo = new GtfInfo(GtfInfo.Type.GENE, feature.getGeneName(), interval);
            idToInfo.put(feature.getGeneId(), gtfInfo);
        } else if (featureType.equals("transcript")) {
            Interval interval = new Interval(feature.getContig(), feature.getStart(), feature.getEnd());
            GtfInfo gtfInfo = new GtfInfo(GtfInfo.Type.TRANSCRIPT, feature.getGeneName(), interval);
            idToInfo.put(feature.getTranscriptId(), gtfInfo);
            if(feature.getStart() < idToInfo.get(feature.getGeneId()).getStart()){
                geneStart = feature.getStart();
                Interval geneInterval = new Interval(feature.getContig(), geneStart, geneEnd); //interval class is immutable, so I need to create a new interval
                GtfInfo gtfGeneInfo = new GtfInfo(GtfInfo.Type.GENE, feature.getGeneName(), geneInterval);
                idToInfo.put(feature.getGeneId(), gtfGeneInfo);
            }
            if(feature.getEnd() > idToInfo.get(feature.getGeneId()).getEnd()){
                geneEnd = feature.getEnd();
                Interval geneInterval = new Interval(feature.getContig(), geneStart, geneEnd); //interval class is immutable, so I need to create a new interval
                GtfInfo gtfGeneInfo = new GtfInfo(GtfInfo.Type.GENE, feature.getGeneName(), geneInterval);
                idToInfo.put(feature.getGeneId(), gtfGeneInfo);
            }
        }
    }

    @Override
    public GATKPath getDrivingFeaturePath() {
        return inputFile;
    }
}

