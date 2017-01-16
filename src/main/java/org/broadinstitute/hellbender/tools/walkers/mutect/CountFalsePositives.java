package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;


@CommandLineProgramProperties(
        summary = "Count PASS (false positive) variants in a vcf file for Mutect2 NA12878 normal-normal evaluation",
        oneLineSummary = "Count PASS variants",
        programGroup = VariantProgramGroup.class

)
/**
 * Created by tsato on 12/28/16.
 */
public class CountFalsePositives extends CommandLineProgram {
    @Argument(doc = "a line-delimited list of mutect2 vcf files", fullName= "", shortName = "V", optional = false)
    protected File vcfList;

    @Argument(doc = "output file", fullName= "output", shortName = "O", optional = false)
    protected File output;

    // TODO: eventually use tumor and normal sample names instead of the file name. To do so we must extract them from the vcf, which I don't know how.
    @Override
    public Object doWork() {
        Utils.regularReadableUserFile(vcfList);
        try ( FalsePositiveTableWriter writer = new FalsePositiveTableWriter(output) ) {
            List<String> vcfFilePaths = Files.readAllLines(Paths.get(vcfList.toString()));
            List<FalsePositiveRecord> falsePositiveRecords = vcfFilePaths.stream().map(str -> new FalsePositiveRecord(new File(str))).collect(Collectors.toList());
            writer.writeAllRecords(falsePositiveRecords);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading %s, %s, or files listed in them", vcfList.toString(), output));
        }

        return "SUCCESS";
    }

    private final String ID_COLUMN_NAME = "filename";
    private final String SNP_COLUMN_NAME = "snp";
    private final String INDEL_COLUMN_NAME = "indel";

    private class FalsePositiveTableWriter extends TableWriter<FalsePositiveRecord> {
        private FalsePositiveTableWriter(final File output) throws IOException {
            super(output, new TableColumnCollection(ID_COLUMN_NAME, SNP_COLUMN_NAME, INDEL_COLUMN_NAME));
        }

        @Override
        protected void composeLine(final FalsePositiveRecord falsePositiveRecord, final DataLine dataLine) {
            dataLine.set(ID_COLUMN_NAME, falsePositiveRecord.getId())
                    .set(SNP_COLUMN_NAME, falsePositiveRecord.getSnpFalsePositives())
                    .set(INDEL_COLUMN_NAME, falsePositiveRecord.getIndelFalsePositives());
        }
    }

    private class FalsePositiveRecord {
        private String id;
        private long snpFalsePositives;
        private long indelFalsePositives;

        private FalsePositiveRecord(final File mutect2Vcf) {
            Utils.nonNull(mutect2Vcf, "the vcf file cannot be null");

            id = FilenameUtils.removeExtension(mutect2Vcf.getName());
            snpFalsePositives = StreamSupport.stream(new FeatureDataSource<VariantContext>(mutect2Vcf).spliterator(), false)
                    .filter(vc -> vc.getFilters().isEmpty())
                    .filter(vc -> vc.isSNP())
                    .count();
            indelFalsePositives = StreamSupport.stream(new FeatureDataSource<VariantContext>(mutect2Vcf).spliterator(), false)
                    .filter(vc -> vc.getFilters().isEmpty())
                    .filter(vc -> !vc.isSNP())
                    .count();
        }

        public String getId(){
            return id;
        }

        public long getSnpFalsePositives(){
            return snpFalsePositives;
        }

        public long getIndelFalsePositives(){
            return indelFalsePositives;
        }

    }
}
