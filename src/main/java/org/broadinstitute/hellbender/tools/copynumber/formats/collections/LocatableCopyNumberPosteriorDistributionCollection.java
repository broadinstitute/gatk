package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LocatableCopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class LocatableCopyNumberPosteriorDistributionCollection {

    private final String sample;
    private final List<LocatableCopyNumberPosteriorDistribution> records;

    public LocatableCopyNumberPosteriorDistributionCollection(final File file) {
        final String path = file.getAbsolutePath();
        if (path.toLowerCase().endsWith(FileExtensions.VCF.toLowerCase()) || path.toLowerCase().endsWith(FileExtensions.COMPRESSED_VCF.toLowerCase())) {
            final VCFFileReader reader = new VCFFileReader(file, false);
            final List<String> samples = reader.getFileHeader().getSampleNamesInOrder();
            Utils.validate(samples.size() > 0, "No samples found in vcf: " + file.getAbsolutePath());
            this.sample = samples.get(0);

            this.records = Utils.stream(reader.iterator())
                    .map(record -> parseRecord(record, sample)).collect(Collectors.toList());
        } else {
            throw new UserException.BadInput("CNV intervals file must be " + FileExtensions.VCF + " or " + FileExtensions.COMPRESSED_VCF);
        }
    }

    private static LocatableCopyNumberPosteriorDistribution parseRecord(final VariantContext variantContext,
                                                                        final String sample) {
        final SimpleInterval interval = new SimpleInterval(variantContext.getContig(), variantContext.getStart(), variantContext.getEnd());
        final GenotypesContext genotypesContext = variantContext.getGenotypes();
        if (genotypesContext.isEmpty()) {
            throw new UserException.BadInput("No genotypes found in variant context " + variantContext.getID());
        }
        final Genotype genotype = genotypesContext.get(sample);
        if (!genotype.hasExtendedAttribute("CNLP")) {
            throw new UserException.BadInput("Copy number genotype not found in variant context " + variantContext.getID());
        }
        final String[] phredScaledLikelihoodStrings = ((String) genotype.getExtendedAttribute("CNLP")).split(",");

        //Posteriors reported as integer phred-scaled likelihoods and need to be renormalized
        final double[] approximatePosteriors = new double[phredScaledLikelihoodStrings.length];
        double total = 0;
        for (int i = 0; i < phredScaledLikelihoodStrings.length; i++) {
            final double logLikelihood = -Double.valueOf(phredScaledLikelihoodStrings[i]) / (10.0 * Math.log(10.0));
            approximatePosteriors[i] = Math.max(Math.exp(logLikelihood), Double.MIN_NORMAL);
            total += approximatePosteriors[i];
        }

        final Map<IntegerCopyNumberState,Double> copyNumberLogPosteriors = new HashMap<>(SVUtils.hashMapCapacity(phredScaledLikelihoodStrings.length));
        for (int i = 0; i < phredScaledLikelihoodStrings.length; i++) {
            copyNumberLogPosteriors.put(new IntegerCopyNumberState(i), Math.min(Math.log(approximatePosteriors[i] / total), Double.MIN_NORMAL));
        }
        return new LocatableCopyNumberPosteriorDistribution(copyNumberLogPosteriors, sample, interval);
    }

    public List<LocatableCopyNumberPosteriorDistribution> getRecords() {
        return records;
    }

    public String getSampleName() {
        return sample;
    }

    public OverlapDetector<LocatableCopyNumberPosteriorDistribution> createNewOverlapDetector() {
        return OverlapDetector.create(records);
    }
}