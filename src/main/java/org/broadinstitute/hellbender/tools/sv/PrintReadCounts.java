package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.codecs.copynumber.SimpleCountCodec;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

/**
 * <p>Prints (and optionally subsets) an rd (DepthEvidence) file or a counts file
 * as one or more (for multi-sample DepthEvidence files) counts files for CNV determination.</p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         A locus-sorted DepthEvidence file (name ends with ".rd.txt", ".rd.txt.gz", or "rd.bci"),
 *         or a counts file (name ends with ".counts.tsv", or ".counts.tsv.gz").
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         A counts file (or files) for use by one of the CNV-calling tools.
 *     </li>
 * </ul>
 * <p>The output files are named "{output-prefix}{sample-name}.counts.tsv"</p>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk PrintCNVCounts \
 *       -F input.rd.txt.gz \
 *       -L chr1
 * </pre>
 *
 * @author Ted Sharpe &lt;tsharpe@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Prints count files for CNV determination.",
        oneLineSummary = "Prints count files for CNV determination.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public class PrintReadCounts extends FeatureWalker<Feature> {
    public static final String INPUT_ARGNAME = "input-counts";
    public static final String OUTPUT_PREFIX_ARGNAME = "output-prefix";
    public static final String OUTPUT_FILES_ARGNAME = "output-file-list";

    @Argument(
            doc = "Input file of counts",
            fullName = INPUT_ARGNAME,
            shortName = StandardArgumentDefinitions.FEATURE_SHORT_NAME
    )
    private GATKPath inputPath;

    @Argument(
            doc = "Output file path prefix. Paths have the form \"{output-prefix}{sample-name}.counts.tsv\". " +
                    "Default is the current working directory.",
            fullName = OUTPUT_PREFIX_ARGNAME,
            optional = true
    )
    private String outputPrefix = "";

    @Argument(
            doc = "Two-column list of outputs in the form <sample-name>\\t<output file path>",
            fullName = OUTPUT_FILES_ARGNAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = true
    )
    private String outputFilesFilename;

    private BufferedWriter[] writers;
    private List<String> sampleNames;

    @Override
    public GATKPath getDrivingFeaturePath() {
        return inputPath;
    }

    @Override
    protected boolean isAcceptableFeatureType( final Class<? extends Feature> featureType ) {
        return featureType.equals(DepthEvidence.class) || featureType.equals(SimpleCount.class);
    }

    @Override
    public void onTraversalStart() {
        final Object header = getDrivingFeaturesHeader();
        if ( header instanceof SampleLocatableMetadata ) {
            initializeWriter((SampleLocatableMetadata)header);
        } else if ( header instanceof SVFeaturesHeader ) {
            initializeWriters((SVFeaturesHeader)header);
        } else {
            throw new UserException("Input file has no header.");
        }
    }

    @Override
    public void apply( final Feature feature,
                       final ReadsContext readsContext,
                       final ReferenceContext referenceContext,
                       final FeatureContext featureContext ) {
        if ( feature instanceof DepthEvidence ) {
            apply((DepthEvidence)feature);
        } else {
            apply((SimpleCount)feature);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        final int nSamples = sampleNames.size();
        for ( int idx = 0; idx != nSamples; ++idx ) {
            try {
                writers[idx].close();
            } catch ( final IOException ioe ) {
                throw new UserException("Failed to close output file for sample " +  sampleNames.get(idx), ioe);
            }
        }
        return null;
    }

    private void initializeWriter( final SampleLocatableMetadata metadata ) {
        final String sampleName = metadata.getSampleName();
        sampleNames = Collections.singletonList(sampleName);
        writers = new BufferedWriter[1];
        try ( final BufferedWriter outputFilesWriter =
                      outputFilesFilename == null ? null :
                              new BufferedWriter(new FileWriter(outputFilesFilename)) ) {
            writers[0] = createWriter(sampleName, metadata.toHeader(), outputFilesWriter);
        } catch ( final IOException ioe ) {
            throw new UserException("Can't open list of output files " + outputFilesFilename, ioe);
        }
    }

    private void initializeWriters( final SVFeaturesHeader featuresHeader ) {
        sampleNames = featuresHeader.getSampleNames();
        SAMSequenceDictionary dictionary = featuresHeader.getDictionary();
        if ( dictionary == null ) {
            dictionary = getMasterSequenceDictionary();
            if ( dictionary == null ) {
                throw new UserException("No dictionary available.  Supply one with --sequence-dictonary.");
            }
        }
        final int nSampleNames = sampleNames.size();
        writers = new BufferedWriter[nSampleNames];
        try ( final BufferedWriter outputFilesWriter =
                      outputFilesFilename == null ? null :
                              new BufferedWriter(new FileWriter(outputFilesFilename)) ) {
            for ( int idx = 0; idx != nSampleNames; ++idx ) {
                final SAMFileHeader samFileHeader = new SAMFileHeader(dictionary);
                final SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord("GATKCopyNumber");
                final String sampleName = sampleNames.get(idx);
                readGroupRecord.setSample(sampleName);
                samFileHeader.addReadGroup(readGroupRecord);
                writers[idx] = createWriter(sampleName, samFileHeader, outputFilesWriter);
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't open list of output files " + outputFilesFilename, ioe);
        }
    }

    private BufferedWriter createWriter( final String sampleName,
                                         final SAMFileHeader header,
                                         final BufferedWriter outputFilesWriter ) {
        final String filename = outputPrefix + sampleName + ".counts.tsv";
        if ( outputFilesWriter != null ) {
            try {
                outputFilesWriter.write(sampleName + "\t" + filename);
                outputFilesWriter.newLine();
            } catch ( final IOException ioe ) {
                throw new UserException("Can't write to list of output files " + outputFilesFilename, ioe);
            }
        }

        final BufferedWriter writer;
        try {
            writer = new BufferedWriter(new FileWriter(filename));
        } catch ( final IOException ioe ) {
            throw new UserException("Can't open " + filename + " for output.", ioe);
        }
        try {
            new SAMTextHeaderCodec().encode(writer, header, true);
        } catch ( final RuntimeIOException ioe ) {
            throw new UserException("Can't write header to " + filename, ioe);
        }
        try {
            writer.write("CONTIG\tSTART\tEND\tCOUNT");
            writer.newLine();
        } catch ( final IOException ioe ) {
            throw new UserException("Can't writer column header to " + filename, ioe);
        }
        return writer;
    }

    private void apply( final SimpleCount simpleCount ) {
        final BufferedWriter writer = writers[0];
        try {
            writer.write(SimpleCountCodec.encode(simpleCount));
            writer.newLine();
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write record to output.", ioe);
        }
    }

    private void apply( final DepthEvidence depthEvidence ) {
        // don't use DepthEvidence codec to encode--we want 1-based coordinates
        final String intervalFields =
                depthEvidence.getContig() + "\t" + depthEvidence.getStart() + "\t" + depthEvidence.getEnd() + "\t";
        final int[] counts = depthEvidence.getCounts();
        for ( int idx = 0; idx != writers.length; ++idx ) {
            try {
                writers[idx].write(intervalFields + counts[idx]);
                writers[idx].newLine();
            } catch ( final IOException ioe ) {
                throw new UserException("Can't write to output file for sample " + sampleNames.get(idx), ioe);
            }
        }
    }
}
