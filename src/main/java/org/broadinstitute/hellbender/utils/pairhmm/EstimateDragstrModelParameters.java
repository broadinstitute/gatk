package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.samtools.*;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.LineReader;
import org.apache.commons.collections4.EnumerationUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import static org.broadinstitute.hellbender.utils.pairhmm.DragstrConstants.*;

import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.*;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

@CommandLineProgramProperties(
        programGroup = ReferenceProgramGroup.class,
        summary = "Determine the presence of STR in a reference sequence",
        oneLineSummary = "Determines the presence of STR in a reference sequence"
)
public class EstimateDragstrModelParameters extends GATKTool {

    @ArgumentCollection
    private final DragstrCasesSamplerArgumentCollection dragstrCasesSamplerArgumentCollection = new DragstrCasesSamplerArgumentCollection();
    private Pattern ZIP_ENTRY_NAME_REGEXP = Pattern.compile("^(\\d+)/(\\d+).bin$");

    @ArgumentCollection
    public DragstrModelEstimatorArgumentCollection dragstrModelEstimatorArgumentCollection = new DragstrModelEstimatorArgumentCollection();

    @Argument(doc = "location of the .zip file containing the locations to sample from. This zip file is to be generated using SampleSitesForDRAGstrModel tool ",
              fullName = SAMPLING_LOCI_ARGUMENT_FULL_NAME)
    private String samplingLoci = null;

    @Argument(fullName= MAX_PERIOD_ARGUMENT_FULL_NAME,
            doc="maximum STR period sampled", optional = true, minValue = 1, maxValue = 10)
    private int maxPeriod = DEFAULT_MAX_PERIOD;

    @Argument(fullName= MAX_REPEATS_ARGUMENT_FULL_NAME,
            doc="maximum STR repeat sampled", optional = true, minValue = 1, maxValue = 20)
    private int maxRepeat = DEFAULT_MAX_REPEATS;


    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              doc = "name of the output parameter file")
    private String outputPath;

    @Override
    public boolean requiresReads() {
        return true;
    }



    @Override
    public void traverse() {
        final DragstrModelEstimator estimator = new DragstrModelEstimator(dragstrModelEstimatorArgumentCollection);
        final DragstrCasesSampler sampler = new DragstrCasesSampler(dragstrCasesSamplerArgumentCollection, directlyAccessEngineReferenceDataSource(), directlyAccessEngineReadsDataSource());
        final DragstrModelEstimator.Estimate estimate = estimator.createEstimate(maxPeriod, maxRepeat);
        try (final ZipFile zipIn = stageZipInput()) {
            checkSameSampleDictionary(zipIn);
            final List<ZipEntry> zipEntries = EnumerationUtils.toList(zipIn.entries());
            for (int period = 1; period <= maxPeriod; period++) {
                final DragstrModelEstimator.PeriodCases periodCases = estimator.createPeriodCases(period, maxRepeat, 10000);
                final String periodEntryPrefix = "" + period + "/";
                for (final ZipEntry entry : zipEntries) {
                    final String entryName = entry.getName();
                    if (entryName.startsWith(periodEntryPrefix)) {
                        final Matcher matcher = ZIP_ENTRY_NAME_REGEXP.matcher(entryName);
                        if (matcher.matches()) {
                            final int repeat = Integer.parseInt(matcher.group(2));
                            try (final BinaryTableReader<DragstrLocus> locusReader = DragstrLocus.binaryReader(zipIn.getInputStream(entry))) {
                                final DragstrModelEstimator.RepeatCases repeatCases = periodCases.getRepeatCases(repeat);
                                final List<DragstrLocus> loci = locusReader.readAll();
                                sampler.sample(repeatCases, loci);
                            } catch (final IOException ex) {
                                throw new UserException.CouldNotReadInputFile(samplingLoci, "could not read entry " + entryName + " from " + samplingLoci, ex);
                            }
                        }
                    }
                }
                logger.info("Estimating gap penalty and het prior for period " + period);
                final long start = System.currentTimeMillis();
                estimator.estimate(estimate, periodCases);
                logger.info("Estimate done in " + ((System.currentTimeMillis() - start) / 1000) + " seconds");
                logger.debug("GP: " + Arrays.toString(estimate.gp[period - 1]));
                logger.debug("Pr: " + Arrays.toString(estimate.ph_het_variant[period - 1]));
            }
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(samplingLoci, "could not open zip file " + samplingLoci, ex);
        }
        final DragstrParams params = DragstrParams.fromEstimate(estimate);
        params.print(outputPath);
    }

    private void checkSameSampleDictionary(final ZipFile zipIn) {
        final ZipEntry zipEntry = zipIn.getEntry("reference.dict");
        if (zipEntry != null) {
            try (final LineReader reader = new BufferedLineReader(zipIn.getInputStream(zipEntry))) {
                final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
                final SAMSequenceDictionary zipDictionary = codec.decode(reader, zipIn.getName() + "/reference.dict").getSequenceDictionary();
                final SAMSequenceDictionary currentDictionary = getBestAvailableSequenceDictionary();
                if (!areThemCompatibleDictionaries(zipDictionary, currentDictionary)) {
                   throw new UserException("the input reference's dictionary does not match the loci zip dictionary");
                };
                logger.info("Loci zip file reference dictionary check was successful");
            } catch (final IOException ex) {
                throw new GATKException("could not open the reference dictionary inside the input loci zip ", ex);
            }
        } else {
            logger.warn("Missing reference dicitonary input sites zip file, we will proceed without checkin that" +
                    " that the reference provided is the same that was used to create this loci zip");
        }
    }

    private boolean areThemCompatibleDictionaries(final SAMSequenceDictionary a, final SAMSequenceDictionary b) {
        if (a == b) {
            return true;
        } else {
            final List<SAMSequenceRecord> aSeqs = a.getSequences();
            final List<SAMSequenceRecord> bSeqs = b.getSequences();
            if (aSeqs.size() != bSeqs.size()) {
                return false;
            } else {
                for (int i = 0; i < aSeqs.size(); i++) {
                    if (!aSeqs.get(i).isSameSequence(bSeqs.get(i))) {
                        return false;
                    }
                }
                return true;
            }
        }
    }

    private ZipFile stageZipInput() throws IOException {
        final File staged = BucketUtils.stageFile(samplingLoci, null, true, true);
        return new ZipFile(staged) {
                public void close() throws IOException {
                    super.close();
                    if (!staged.delete()) {
                        throw new GATKException("could not delete temporary copy of zip input at " + staged);
                    }
                }
        };
    }
}
