package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Stream;
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

    @Override
    public Object doWork() {
        Utils.regularReadableUserFile(vcfList);
        String header = String.format("%s\t%s\t%s\n", "sample", "fps", "snp-or-indel");
        StringBuilder sb = new StringBuilder(header);

        try (
                Stream<String> vcfFiles = Files.lines(Paths.get(vcfList.toString()));
                PrintWriter outputWriter = new PrintWriter(new FileWriter(output))
        ) {
            sb = vcfFiles.reduce(sb,
                    (acc, f) -> acc.append(countFalsePositivesAndBuildString(new File(f))),
                    StringBuilder::append);
            outputWriter.print(sb.toString());
        } catch (IOException e){
            // TODO: catching the exception here and Utils.regularReadableUserFile() are redundant
            throw new UserException(String.format("Encountered an IO exception while reading (or reading a file in) %s or %s", vcfList.toString(), output));
        }

        return "SUCCESS";
    }

    private static String countFalsePositivesAndBuildString(File m2vcf) {
        final String sampleName = FilenameUtils.removeExtension(m2vcf.getName());

        final long snpFalsePositives = StreamSupport.stream(new FeatureDataSource<VariantContext>(m2vcf).spliterator(), false)
                .filter(vc -> vc.getFilters().isEmpty())
                .filter(vc -> vc.isSNP())
                .count();
        final long indelFalsePositives = StreamSupport.stream(new FeatureDataSource<VariantContext>(m2vcf).spliterator(), false)
                .filter(vc -> vc.getFilters().isEmpty())
                .filter(vc -> !vc.isSNP())
                .count();

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%s\t%d\t%s\n", sampleName, snpFalsePositives, "snp"));
        sb.append(String.format("%s\t%d\t%s\n", sampleName, indelFalsePositives, "indel"));
        return sb.toString();
    }
}
