package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import org.apache.parquet.filter2.predicate.Operators;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.codecs.gtf.Gff3Codec;
import org.broadinstitute.hellbender.utils.codecs.gtf.GffCodec;
import org.broadinstitute.hellbender.utils.codecs.gtf.GtfFeature;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.*;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "SortGff feature file",
        oneLineSummary = "SortGffFeatureFile",
        programGroup = OtherProgramGroup.class
)
public class SortGff extends FeatureWalker<GtfFeature> {

    @Argument(shortName = "F", fullName = "feature_file", doc = "Feature file (eg., VCF or BED file)")
    public File featuresFile;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file")
    private File outputFile = null;


    SortingCollection<GtfFeature> sorter = SortingCollection.newInstance(
            GtfFeature.class,
            new Gff3SortingCollectionCodec(),
            (GtfFeature f1, GtfFeature f2) -> f1.getContig().compareTo(f2.getContig()) == 0 ? f1.getStart() - f2.getStart() : f1.getContig().compareTo(f2.getContig()),
            500000
            );

    @Override
    protected boolean isAcceptableFeatureType(Class<? extends Feature> featureType) {
        return featureType.isAssignableFrom(GtfFeature.class);
    }

    @Override
    public void apply(final GtfFeature feature, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext )
    {
        sorter.add(feature);
    }

    @Override
    public File getDrivingFeatureFile() {
        return featuresFile;
    }

    @Override
    public Object onTraversalSuccess() {
        final Gff3SortingCollectionCodec outCodec = new Gff3SortingCollectionCodec();
        try (final PrintStream outputStream = new PrintStream(outputFile)) {
            outputStream.println("##gff-version 3");
            outCodec.setOutputStream(outputStream);
            for(final GtfFeature feature: sorter) {
                outCodec.encode(feature);
            }
        } catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(outputFile, e);
        }

        return null;
    }


//    class Gff3SortingCollectionCodec implements SortingCollection.Codec<GtfFeature> {
//        ObjectOutputStream os;
//        ObjectInputStream is;
//
//        @Override
//        public SortingCollection.Codec<GtfFeature> clone() {
//            return new Gff3SortingCollectionCodec();
//        }
//
//        @Override
//        public void setOutputStream(final OutputStream os) {
//            try {
//                this.os = new ObjectOutputStream(os);
//            } catch (final IOException ex) {
//                throw new GATKException("cannot open output stream in SortGff");
//            }
//        }
//
//        @Override
//        public void setInputStream(final InputStream is) {
//            try {
//                this.is = new ObjectInputStream(is);
//            } catch (final IOException ex) {
//                throw new GATKException("cannot open input stream in SortGff");
//            }
//        }
//
//        @Override
//        public GtfFeature decode() {
//            try {
//                final SerializableGtfFeature serializableGtfFeature = (SerializableGtfFeature)is.readObject();
//                return new GtfFeature(serializableGtfFeature.contig,
//                        serializableGtfFeature.source,
//                        serializableGtfFeature.type,
//                        serializableGtfFeature.start,
//                        serializableGtfFeature.end,
//                        serializableGtfFeature.strand,
//                        serializableGtfFeature.phase,
//                        serializableGtfFeature.attributes);
//            } catch (final Exception ex) {
//                throw new GATKException("error reading from input stream in SortGff");
//            }
//        }
//
//        @Override
//        public void encode(final GtfFeature val) {
//            try {
//                final SerializableGtfFeature serializableGtfFeature = new SerializableGtfFeature(
//                        val.getContig(),
//                        val.getSource(),
//                        val.getType(),
//                        val.getStart(),
//                        val.getEnd(),
//                        val.getStrand(),
//                        val.getPhase(),
//                        val.getAttributes());
//                os.writeObject(serializableGtfFeature);
//            } catch (final Exception ex) {
//                throw new GATKException("error writing to output stream in SortGff");
//            }
//        }
//
//        public class SerializableGtfFeature implements Serializable {
//            public final String contig;
//            public final String source;
//            public final String type;
//            public final int start;
//            public final int end;
//            public final Strand strand;
//            public final int phase;
//            public final Map<String, String> attributes;
//
//            SerializableGtfFeature(final String contig, final String source, final String type, final int start, final int end, final Strand strand, final int phase, final Map<String, String> attributes) {
//                this.contig = contig;
//                this.source = source;
//                this.type = type;
//                this.start = start;
//                this.end = end;
//                this.phase = phase;
//                this.strand = strand;
//                this.attributes = attributes;
//            }
//        }
//    }

    class Gff3SortingCollectionCodec implements SortingCollection.Codec<GtfFeature> {
        PrintStream writer;
        LineIterator lineIteratorIn;
        Gff3Codec gff3Codec;


        Gff3SortingCollectionCodec() {
            gff3Codec = new Gff3Codec();
        }

        @Override
        public SortingCollection.Codec<GtfFeature> clone() {
            return new Gff3SortingCollectionCodec();
        }

        @Override
        public void setOutputStream(final OutputStream os) {
            writer = new PrintStream(os);
        }

        @Override
        public void setInputStream(final InputStream is) {
            lineIteratorIn = new LineIteratorImpl(new SynchronousLineReader(is));
        }

        @Override
        public GtfFeature decode() {
            if (!lineIteratorIn.hasNext()) {
                return null;
            }
            return gff3Codec.decode(lineIteratorIn);
        }

        @Override
        public void encode(final GtfFeature val) {
            final String lineNoAttributes = String.join(Gff3Codec.FIELD_DELIMITER,
                                                        new String[] {
                                                                val.getContig(),
                                                                val.getSource(),
                                                                val.getType(),
                                                                Integer.toString(val.getStart()),
                                                                Integer.toString(val.getEnd()),
                                                                ".",
                                                                val.getStrand().toString(),
                                                                val.getPhase()<0 ? "." : Integer.toString(val.getPhase())
                                                            }
                                                        );
            final List<String> attributesStrings = val.getAttributes().entrySet().stream().map(e -> String.join("=", new String[] {e.getKey(), e.getValue()})).collect(Collectors.toList());
            final String attributesString = String.join(";", attributesStrings);

            final String lineString = lineNoAttributes + Gff3Codec.FIELD_DELIMITER + attributesString;
            writer.println(lineString);
        }


    }
}
