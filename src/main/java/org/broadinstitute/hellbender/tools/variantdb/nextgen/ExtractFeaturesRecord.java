package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.samtools.util.Locatable;
import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;

import java.util.Objects;

public class ExtractFeaturesRecord implements Locatable {

    private final long location;
    private final String contig;
    private final int start;
    private final int end;
    private final String ref;
    private final String allele;
    private final Double rawQual;
    private final Double refAD;                     // nullable, if null -> 0
    private final Float asMQRankSum;                // nullable
    private final String asMQRankSumFreqTable;      // nullable
    private final Float asReadPosRankSum;           // nullable
    private final String asReadPosRankSumFreqTable; // nullable
    private final Double rawMQ;
    private final Double rawAD;
    private final Double rawADGT1;
    private final int sbRefPlus;
    private final int sbRefMinus;
    private final int sbAltPlus;
    private final int sbAltMinus;
    private final int numHetSamples;
    private final int numHomvarSamples;
    private final int distinctAlleles;
    private final int hqGenotypeSamples;
    private final Double qualApprox;

    // FEATURE_EXTRACT_FIELDS = Arrays.asList(
    //      LOCATION_FIELD_NAME,
    //      REF_ALLELE_FIELD_NAME,
    //      "allele",
    //      RAW_QUAL,
    //      "ref_ad",
    //      AS_MQRankSum,
    //      "AS_MQRankSum_ft",
    //      AS_ReadPosRankSum,
    //      "AS_ReadPosRankSum_ft",
    //      RAW_MQ,
    //      RAW_AD,
    //      "RAW_AD_GT_1",
    //      "SB_REF_PLUS",
    //      "SB_REF_MINUS",
    //      "SB_ALT_PLUS",
    //      "SB_ALT_MINUS",
    //      "num_het_samples",
    //      "num_homvar_samples");


    public ExtractFeaturesRecord(GenericRecord genericRecord) {
        this.location = Long.parseLong(genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME).toString());
        this.contig = SchemaUtils.decodeContig(location);
        this.start = SchemaUtils.decodePosition(location);
        this.end = start;

        this.ref = Objects.toString(genericRecord.get("ref"), null);
        this.allele = Objects.toString(genericRecord.get("allele"), null);
        this.rawQual = Double.valueOf(genericRecord.get(SchemaUtils.RAW_QUAL).toString());
        this.rawMQ = Double.valueOf(genericRecord.get(SchemaUtils.RAW_MQ).toString());
        this.rawAD = Double.valueOf(genericRecord.get(SchemaUtils.RAW_AD).toString());
        this.rawADGT1 = Double.valueOf(genericRecord.get("RAW_AD_GT_1").toString());
        this.qualApprox = Double.valueOf(genericRecord.get("sum_qualapprox").toString());

        this.sbRefPlus = Double.valueOf(genericRecord.get("SB_REF_PLUS").toString()).intValue();
        this.sbRefMinus = Double.valueOf(genericRecord.get("SB_REF_MINUS").toString()).intValue();
        this.sbAltPlus = Double.valueOf(genericRecord.get("SB_ALT_PLUS").toString()).intValue();
        this.sbAltMinus = Double.valueOf(genericRecord.get("SB_ALT_MINUS").toString()).intValue();
        this.numHetSamples = Integer.parseInt(genericRecord.get("num_het_samples").toString());
        this.numHomvarSamples = Integer.parseInt(genericRecord.get("num_homvar_samples").toString());
        
        this.distinctAlleles = Integer.parseInt(genericRecord.get("distinct_alleles").toString());
        this.hqGenotypeSamples = Integer.parseInt(genericRecord.get("hq_genotype_samples").toString());

        // nullable fields - TODO double check that these are the only nullable fields
        Object asMQRankSumNullable = genericRecord.get(SchemaUtils.AS_MQRankSum);
        this.asMQRankSum = ( asMQRankSumNullable == null ) ? null : Float.valueOf(Objects.toString(asMQRankSumNullable));
        this.asMQRankSumFreqTable = Objects.toString(genericRecord.get(SchemaUtils.AS_MQRankSum + "_ft"));

        Object asReadPosRankSumNullable = genericRecord.get(SchemaUtils.AS_ReadPosRankSum);
        this.asReadPosRankSum = ( asReadPosRankSumNullable == null ) ? null : Float.valueOf(Objects.toString(asReadPosRankSumNullable));
        this.asReadPosRankSumFreqTable = Objects.toString(genericRecord.get(SchemaUtils.AS_ReadPosRankSum + "_ft"));

        // if ref_ad is not defined, set it to zero
        Object refADNullable = genericRecord.get("ref_ad");
        this.refAD = ( refADNullable == null ) ? 0 : Double.valueOf(Objects.toString(refADNullable));

    }

    @Override
    public String getContig() { return this.contig; }

    @Override
    public int getStart() { return this.start; }

    @Override
    public int getEnd() { return this.end; }

    public long getLocation() { return this.location; }

    public String getRef() { return this.ref; }

    public String getAllele() { return this.allele; }

    public Double getRawQual() { return this.rawQual; }

    public Double getRefAD() { return this.refAD; }

    public Float getAsMQRankSum() { return this.asMQRankSum; }

    public String getAsMQRankSumFreqTable() { return this.asMQRankSumFreqTable; }

    public Float getAsReadPosRankSum() { return this.asReadPosRankSum; }

    public String getAsReadPosRankSumFreqTable() { return this.asReadPosRankSumFreqTable; }

    public Double getRawMQ() { return this.rawMQ; }

    public Double getRawAD() { return this.rawAD; }

    public Double getRawADGT1() { return this.rawADGT1; }

    public Double getQualApprox() { return this.qualApprox; }

    public int getSbRefPlus() { return this.sbRefPlus; }

    public int getSbRefMinus() { return this.sbRefMinus; }

    public int getSbAltPlus() { return this.sbAltPlus; }

    public int getSbAltMinus() { return this.sbAltMinus; }

    public int getNumHetSamples() { return this.numHetSamples; }

    public int getNumHomvarSamples() { return this.numHomvarSamples; }

    public int getDistinctAlleles() { return this.distinctAlleles; }

    public int getHqGenotypeSamples() { return this.hqGenotypeSamples; }
    
}
