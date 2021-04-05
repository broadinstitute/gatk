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
    private final String refAD;
    private final Float asMQRankSum;
//    private final String asMQRankSumFt;
    private final Float asReadPosRankSum;
//    private final String asReadPosRankSumFt;
    private final Double rawMQ;
    private final Double rawAD;
    private final Double rawADGT1;
    private final int sbRefPlus;
    private final int sbRefMinus;
    private final int sbAltPlus;
    private final int sbAltMinus;
    private final int numHetSamples;
    private final int numHomvarSamples;

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

        // TODO check which if any of these are nullable
        this.ref = Objects.toString(genericRecord.get("ref"), null);
        this.allele = Objects.toString(genericRecord.get("allele"), null);
        this.rawQual = Double.valueOf(genericRecord.get(SchemaUtils.RAW_QUAL).toString());
        this.refAD = Objects.toString(genericRecord.get("ref_ad"), null);
        this.asMQRankSum = Float.parseFloat(Objects.toString(genericRecord.get(SchemaUtils.AS_MQRankSum), null));
//        this.asMQRankSumFt = "";
        this.asReadPosRankSum = Float.parseFloat(Objects.toString(genericRecord.get(SchemaUtils.AS_ReadPosRankSum), null));
//        this.asReadPosRankSumFt = "";
        this.rawMQ = Double.valueOf(genericRecord.get(SchemaUtils.RAW_MQ).toString());
        this.rawAD = Double.valueOf(genericRecord.get(SchemaUtils.RAW_AD).toString());
        this.rawADGT1 = Double.valueOf(genericRecord.get("RAW_AD_GT_1").toString());
        this.sbRefPlus = Double.valueOf(genericRecord.get("SB_REF_PLUS").toString()).intValue();
        this.sbRefMinus = Double.valueOf(genericRecord.get("SB_REF_MINUS").toString()).intValue();
        this.sbAltPlus = Double.valueOf(genericRecord.get("SB_ALT_PLUS").toString()).intValue();
        this.sbAltMinus = Double.valueOf(genericRecord.get("SB_ALT_MINUS").toString()).intValue();
        this.numHetSamples = Integer.parseInt(genericRecord.get("num_het_samples").toString());
        this.numHomvarSamples = Integer.parseInt(genericRecord.get("num_homvar_samples").toString());
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

    public String getRefAD() { return this.refAD; }

    public Float getAsMQRankSum() { return this.asMQRankSum; }

    public Float getAsReadPosRankSum() { return this.asReadPosRankSum; }

    public Double getRawMQ() { return this.rawMQ; }

    public Double getRawAD() { return this.rawAD; }

    public Double getRawADGT1() { return this.rawADGT1; }

    public int getSbRefPlus() { return this.sbRefPlus; }

    public int getSbRefMinus() { return this.sbRefMinus; }

    public int getSbAltPlus() { return this.sbAltPlus; }

    public int getSbAltMinus() { return this.sbAltMinus; }

    public int getNumHetSamples() { return this.numHetSamples; }

    public int getNumHomvarSamples() { return this.numHomvarSamples; }
}
