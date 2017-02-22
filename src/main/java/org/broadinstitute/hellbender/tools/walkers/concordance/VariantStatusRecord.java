package org.broadinstitute.hellbender.tools.walkers.concordance;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by Takuto Sato on 2/8/17.
 */
public class VariantStatusRecord {
    private static final String CHROMOSOME_COLUMN_NAME = "CHROM";
    private static final String START_POSITION_COLUMN_NAME = "START";
    private static final String END_POSITION_COLUMN_NAME = "END";
    private static final String REF_ALLELE_COLUMN_NAME = "REF";
    private static final String ALT_ALLELE_COLUMN_NAME = "ALT";
    private static final String VARIANT_TYPE_COLUMN_NAME = "TYPE";
    private static final String ALLELE_FRACTION_COLUMN_NAME = "TUMOR_ALT_AF";
    private static final String TRUTH_STATUS_COLUMN_NAME = "TRUTH_STATUS";
    private static final String[] VARIANT_TABLE_COLUMN_HEADERS = {CHROMOSOME_COLUMN_NAME, START_POSITION_COLUMN_NAME, END_POSITION_COLUMN_NAME,
            REF_ALLELE_COLUMN_NAME, ALT_ALLELE_COLUMN_NAME, VARIANT_TYPE_COLUMN_NAME, ALLELE_FRACTION_COLUMN_NAME, TRUTH_STATUS_COLUMN_NAME};

    Status truthStatus;
    VariantContext variantContext;
    double altAlleleFraction;

    // Possible statuses
    public enum Status {
        TRUE_POSITIVE, FALSE_POSITIVE, FALSE_NEGATIVE;
    }

    public VariantStatusRecord(final VariantContext variantContext, final Status truthStatus, final double altAllelefraction){
        this.truthStatus = truthStatus;
        this.variantContext = variantContext;
        this.altAlleleFraction = altAllelefraction;
    }

    public Status getTruthStatus() { return truthStatus; }

    public VariantContext getVariantContext() { return variantContext; }

    public double getAltAlleleFraction(){ return altAlleleFraction; }

    public static class Writer extends TableWriter<VariantStatusRecord> {
        private Writer(final File output) throws IOException {
            super(output, new TableColumnCollection(VariantStatusRecord.VARIANT_TABLE_COLUMN_HEADERS));
        }

        @Override
        protected void composeLine(final VariantStatusRecord record, final DataLine dataLine) {
            // Avoid printing brackets [ ] when we convert an array to a string
            StringBuilder altAlleleBuilder = new StringBuilder();
            for (Allele allele : record.getVariantContext().getAlternateAlleles()){
                altAlleleBuilder.append(allele.toString());
            }

            // And remove the asterisk (*) next to the ref allele
            String refAllele = record.getVariantContext().getReference().toString();
            refAllele = refAllele.substring(0, refAllele.length() - 1);

            dataLine.set(VariantStatusRecord.CHROMOSOME_COLUMN_NAME, record.getVariantContext().getContig())
                    .set(VariantStatusRecord.START_POSITION_COLUMN_NAME, record.getVariantContext().getStart())
                    .set(VariantStatusRecord.END_POSITION_COLUMN_NAME, record.getVariantContext().getEnd())
                    .set(VariantStatusRecord.REF_ALLELE_COLUMN_NAME, refAllele)
                    .set(VariantStatusRecord.ALT_ALLELE_COLUMN_NAME, altAlleleBuilder.toString())
                    .set(VariantStatusRecord.VARIANT_TYPE_COLUMN_NAME, record.getVariantContext().getType().toString())
                    .set(VariantStatusRecord.ALLELE_FRACTION_COLUMN_NAME, record.getAltAlleleFraction()) // TODO: must be able to retrieve the allele fraction from tumor
                    .set(VariantStatusRecord.TRUTH_STATUS_COLUMN_NAME, record.getTruthStatus().toString());
        }
    }

    public static Writer getWriter(final File outputTable){
        try {
            Writer writer = new Writer(outputTable);
            return writer;
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception opening %s.", outputTable), e);
        }
    }

    public static class Reader extends TableReader<VariantStatusRecord> {
        public Reader(final File output) throws IOException {
            super(output);
        }

        @Override
        protected VariantStatusRecord createRecord(final DataLine dataLine) {
            final String chr = dataLine.get(CHROMOSOME_COLUMN_NAME);
            final long start = Long.parseLong(dataLine.get(START_POSITION_COLUMN_NAME));
            final long end = Long.parseLong(dataLine.get(END_POSITION_COLUMN_NAME));
            final Allele refAllele = Allele.create(dataLine.get(REF_ALLELE_COLUMN_NAME), true);
            final List<String> alts = Arrays.asList(dataLine.get(ALT_ALLELE_COLUMN_NAME).split(","));
            final List<Allele> altAlleles = alts.stream().map(str -> Allele.create(str, false)).collect(Collectors.toList());
            final List<Allele> alleles = new ArrayList<>();
            alleles.add(refAllele);
            alleles.addAll(altAlleles);
            final VariantContext.Type type = VariantContext.Type.valueOf(dataLine.get(VARIANT_TYPE_COLUMN_NAME));
            final double alleleFraction = Double.parseDouble(dataLine.get(ALLELE_FRACTION_COLUMN_NAME));
            final Status status = Status.valueOf(dataLine.get(TRUTH_STATUS_COLUMN_NAME));

            VariantContextBuilder builder = new VariantContextBuilder();
            VariantContext variant = builder.chr(chr).start(start).stop(end).alleles(alleles).make();

            return new VariantStatusRecord(variant, status, alleleFraction);
        }

    }

    // TODO: write a reader and add tests
    // TODO: handle * next to the ref allele, brackets around the alt alleles, etc.


}