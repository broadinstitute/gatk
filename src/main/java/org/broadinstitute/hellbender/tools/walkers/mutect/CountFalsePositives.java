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
            List<File> vcfFiles = vcfFilePaths.stream().map(File::new).collect(Collectors.toList());
            writer.writeAllRecords(vcfFiles);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading (or reading a file in) %s or %s", vcfList.toString(), output));
        }

        return "SUCCESS";
    }

    /***
     *
     * @param m2vcf a mutect2 vcf file
     * @return a pair of SNP (left) and INDEL (right) false positives
     */
    private static Pair<Long, Long> countFalsePositives(File m2vcf) {
        final long snpFalsePositives = StreamSupport.stream(new FeatureDataSource<VariantContext>(m2vcf).spliterator(), false)
                .filter(vc -> vc.getFilters().isEmpty())
                .filter(vc -> vc.isSNP())
                .count();
        final long indelFalsePositives = StreamSupport.stream(new FeatureDataSource<VariantContext>(m2vcf).spliterator(), false)
                .filter(vc -> vc.getFilters().isEmpty())
                .filter(vc -> !vc.isSNP())
                .count();
        return Pair.of(snpFalsePositives, indelFalsePositives);
    }

    private class FalsePositiveTableWriter extends TableWriter<File> {
        private FalsePositiveTableWriter(final File output) throws IOException {
            super(output, new TableColumnCollection("filename", "snp", "indel"));
        }

        @Override
        protected void composeLine(final File vcf, final DataLine dataLine) {
            Utils.nonNull(vcf, "the vcf file cannot be null");
            Utils.nonNull(dataLine, "dataLine cannot be null");

            Pair<Long, Long> falsePositives = countFalsePositives(vcf);

            dataLine.set("filename", FilenameUtils.removeExtension(vcf.getName()))
                    .set("snp", falsePositives.getLeft())
                    .set("indel",falsePositives.getRight());
        }
    }
}
