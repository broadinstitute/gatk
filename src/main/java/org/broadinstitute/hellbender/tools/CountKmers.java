package org.broadinstitute.hellbender.tools;

import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.WorkflowProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMap;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.*;
import java.nio.file.Path;
import java.util.Map;

/**
 * Count the number of kmers in a FASTA-specified set of kmers found in the reads of a SAM/BAM/CRAM file.
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A single BAM file</li>
 *     <li>A FASTA file to parse for the set of kmers that interest us</li>
 * </ul>
 *
 * <h3>Example</h3>
 *
 * <pre>
 *   gatk CountKmers \
 *     -I input_reads.bam
 *     --interesting-kmers kmers.fa
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Count the number of kmers in a FASTA-specified set of kmers found in " +
                "the reads of a SAM/BAM/CRAM file",
        oneLineSummary = "Count interesting kmers in reads",
        programGroup = CoverageAnalysisProgramGroup.class
)
@DocumentedFeature
@WorkflowProperties
public final class CountKmers extends ReadWalker {
    @Argument(doc="FASTA file of interesting kmers", fullName="fasta-file")
    private GATKPath fastaPath;

    @Argument(doc="output file of kmer counts", fullName="output-file",
                shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME )
    private GATKPath outputPath;

    private static final int DEFAULT_KMER_SIZE = 31;
    @Argument(doc="kmer size", fullName="kmer-size", optional=true)
    private int kSize = DEFAULT_KMER_SIZE;

    private KmerCounter kmerCounter;

    @Override public void onTraversalStart() {
        kmerCounter = new KmerCounter(kSize);
        kmerCounter.parseFASTAintoKmers(fastaPath, logger);
    }

    @Override public void apply( final GATKRead read,
                                 final ReferenceContext referenceContext,
                                 final FeatureContext featureContext ) {
        kmerCounter.countKmers(read.getBasesNoCopy());
    }

    @Override public Object onTraversalSuccess() {
        kmerCounter.write(outputPath);
        return null;
    }

    /** A number representing a fixed-length sequence of nucleotides */
    public static class Kmer {
        private final long kVal;

        public Kmer( final long kVal ) {
            this.kVal = kVal;
        }

        public Kmer( final Kmer kmer ) {
            this.kVal = kmer.kVal;
        }

        public boolean equals( final Kmer that ) {
            return this.kVal == that.kVal;
        }

        @Override public boolean equals( final Object obj ) {
            if ( this == obj ) return true;
            if ( !(obj instanceof Kmer) ) return false;
            return equals((Kmer)obj);
        }

        @Override public int hashCode() { return (int)SVUtils.fnvLong64(kVal); }

        public String toString( final int kSizeArg ) {
            final StringBuilder sb = new StringBuilder(kSizeArg);

            long val = kVal;
            int kSize = kSizeArg;
            while ( --kSize >= 0 ) {
                sb.append("ACGT".charAt((int)(val & 3)));
                val >>= 2;
            }

            // we built the string in least-significant to most-significant bit order.  reverse it.
            return sb.reverse().toString();
        }
    }

    /** A kmer with a immutable sequence name, and a mutable count */
    public static class KmerCount extends Kmer implements Map.Entry<Kmer, Integer> {
        private final String seqName;
        private int count;

        public KmerCount( final Kmer kmer, final String seqName ) {
            super(kmer);
            this.seqName = seqName;
            count = 0;
        }

        @Override public Kmer getKey() { return new Kmer(this); }
        @Override public Integer getValue() { return count; }
        @Override public Integer setValue( final Integer count ) {
            final int oldCount = this.count;
            this.count = count;
            return oldCount;
        }

        public String getSeqName() { return seqName; }
        public int getCount() { return count; }
        public void incrementCount() { count += 1; }
    }

    /** Builder that creates a kmer when enough calls have been accumulated */
    public static class KmerBuilder {
        private final int kSize;
        private final int kRCShift;
        private final long kMask;
        private final long kCanonicalBit;
        private int kCount;
        private long kVal;
        private long kValRC;

        public KmerBuilder( final int kSize ) {
            this.kSize = kSize;
            this.kRCShift = 2 * kSize - 2;
            this.kMask = (1L << (2 * kSize)) - 1L;
            this.kCanonicalBit = 1L << kSize;
        }

        public Kmer addCall( final char nextChar ) {
            kVal <<= 2;
            kValRC >>= 2;
            switch ( nextChar ) {
                case 'a': case 'A':
                    kValRC |= (3L << kRCShift);
                    break;
                case 'c': case 'C':
                    kVal |= 1L;
                    kValRC |= (2L << kRCShift);
                    break;
                case 'g': case 'G':
                    kVal |= 2L;
                    kValRC |= (1L << kRCShift);
                    break;
                case 't': case 'T': case 'u': case 'U':
                    kVal |= 3L;
                    break;
                default:
                    kCount = -1; // not a valid call, reset the counter
            }
            if ( ++kCount < kSize ) {
                return null;
            }
            return new Kmer(((kVal & kCanonicalBit) == 0 ? kVal : kValRC) & kMask);
        }
    }

    public static class KmerCounter {
        private final int kSize;
        private final HopscotchMap<Kmer, Integer, KmerCount> kmerMap;
        private static final int EOF = -1;
        private static final int DEFAULT_MAP_SIZE = 1000000;

        public KmerCounter( final int kSize ) {
            if ( kSize < 1 || kSize > 31 || (kSize & 1) == 0 ) {
                throw new UserException("kmer-size must be between 1 and 31 and odd");
            }
            this.kSize = kSize;
            this.kmerMap = new HopscotchMap<>(DEFAULT_MAP_SIZE);
        }

        public void parseFASTAintoKmers( final GATKPath fastaPath, final Logger logger ) {
            final Path path = fastaPath.toPath();
            try ( final Reader reader = IOUtils.makeReaderMaybeGzipped(path) ) {
                int nextChar = reader.read();
                if ( nextChar == EOF ) {
                    throw new GATKException("fasta-file " + path + " is empty");
                }
                if ( nextChar != '>' ) {
                    throw new GATKException("fasta-file " + path + " does not start with a comment");
                }
                final StringBuilder seqNameBuilder = new StringBuilder();
                while ( nextChar != EOF ) {

                    // We've read a '>'.  Parse the comment line.
                    while ( (nextChar = reader.read()) != '\n' ) {
                        if ( nextChar == '\r' ) {
                            continue;
                        }
                        if ( nextChar == EOF ) {
                            throw new GATKException("fasta-file " + path + " ends with a comment");
                        }
                        seqNameBuilder.append((char)nextChar);
                    }
                    final String seqName = seqNameBuilder.toString();

                    // We've read the comment line.  Now read the base calls.
                    final KmerBuilder kmerBuilder = new KmerBuilder(kSize);
                    while ( (nextChar = reader.read()) != '>' && nextChar != EOF ) {
                        if ( Character.isWhitespace(nextChar) ) {
                            continue;
                        }
                        final Kmer kmer = kmerBuilder.addCall((char)nextChar);
                        if ( kmer != null ) {
                            if ( !kmerMap.add(new KmerCount(kmer, seqName)) ) {
                                final KmerCount kmerCount = kmerMap.find(kmer);
                                logger.warn("The kmer " + kmer.toString(kSize) +
                                        " appears first in sequence " + kmerCount.getSeqName() +
                                        ".  Its appearance in sequence " + seqName +
                                        " will be ignored.");
                            }
                        }
                    }
                    seqNameBuilder.setLength(0);
                }
            } catch ( final IOException ioe ) {
                throw new GATKException("unable to read fasta-file " + path, ioe);
            }
        }

        public void countKmers( final byte[] calls ) {
            final KmerBuilder kmerBuilder = new KmerBuilder(kSize);
            for ( int idx = 0; idx != calls.length; ++idx ) {
                final Kmer kmer = kmerBuilder.addCall((char)(calls[idx] & 0x00FF));
                if ( kmer != null ) {
                    final KmerCount kmerCount = kmerMap.find(kmer);
                    if ( kmerCount != null ) {
                        kmerCount.incrementCount();
                    }
                }
            }
        }

        public void write( final GATKPath outputPath ) {
            try ( final BufferedWriter out = IOUtils.makeWriterMaybeGzipped(outputPath) ) {
                for ( final KmerCount kmerCount : kmerMap ) {
                    out.append(kmerCount.toString(kSize))
                            .append('\t')
                            .append(kmerCount.getSeqName())
                            .append('\t')
                            .append(Integer.toString(kmerCount.getCount()));
                    out.newLine();
                }
            } catch ( final IOException ioe ) {
                throw new UserException("can't write output to " + outputPath.toPath(), ioe);
            }
        }
    }
}
