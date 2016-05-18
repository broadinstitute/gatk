package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
* Mock-up active region data used in testing.
*
* @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
*/
public final class AssemblyRegionTestDataSet {

    private String reference;
    private String[] haplotypeCigars;
    private List<String> haplotypeStrings;
    private String[] readCigars;
    private byte[] bq;
    private byte[] dq;
    private byte[] iq;
    private int kmerSize;
    private List<Haplotype> haplotypeList;
    private List<GATKRead> readList;
    private AssemblyResultSet assemblyResultSet;
    private String stringRepresentation;
    private List<List<Civar.ElementOffset>> readEventOffsetList;
    private GenomeLocParser genomeLocParser;

    /** Create a new active region data test set */
    public AssemblyRegionTestDataSet(final int kmerSize, final String reference, final String[] haplotypes,
                                   final String[] readCigars, final byte[] bq, final byte[] dq, final byte[] iq) {
        this.reference = reference;
        this.haplotypeCigars = haplotypes;
        this.readCigars = readCigars;
        this.bq = bq;
        this.dq = dq;
        this.iq = iq;
        this.kmerSize = kmerSize;
        this.genomeLocParser = new GenomeLocParser(ArtificialReadUtils.createArtificialSamHeader(1, 1, reference.length()).getSequenceDictionary());
    }

    public String getReference() {
        return reference;
    }

    public String toString() {
        if (stringRepresentation == null)
            return super.toString();
        else return stringRepresentation;
    }

    public AssemblyResultSet assemblyResultSet() {
        if (assemblyResultSet == null) {
            final ReadThreadingGraph rtg = new ReadThreadingGraph(kmerSize);
            rtg.addSequence("anonymous", this.getReference().getBytes(), true);
            for (final String haplotype : this.haplotypesStrings()) {
                rtg.addSequence("anonymous", haplotype.getBytes(), false);
            }
            rtg.buildGraphIfNecessary();
            if (rtg.hasCycles())
                throw new RuntimeException("there is cycles in the reference with kmer size " + kmerSize + ". Don't use this size for the benchmark or change the reference");

            List<Haplotype> haplotypeList = this.haplotypeList();

            assemblyResultSet = new AssemblyResultSet();
            final AssemblyResult ar = new AssemblyResult((haplotypeList.size() > 1 ?
                    AssemblyResult.Status.ASSEMBLED_SOME_VARIATION : AssemblyResult.Status.JUST_ASSEMBLED_REFERENCE),rtg.toSequenceGraph(), rtg);

            for (final Haplotype h : haplotypeList)
                assemblyResultSet.add(h, ar);
        }
        return assemblyResultSet;
    }

    public List<String> haplotypesStrings() {
        if (haplotypeStrings != null) {
            return haplotypeStrings;
        }
        final List<String> result = new ArrayList<>(haplotypeCigars.length);
        String reference = this.reference;
        for (final String cigar : haplotypeCigars) {
            if (cigar.matches("^Civar:.*$")) {
                stringRepresentation = cigar.substring(6);
                result.addAll(expandAllCombinations(cigar.substring(6),reference));
            } else if (cigar.matches("^.*\\d+.*$")) {
                result.add(applyCigar(reference, cigar,0,true));
            } else {
                result.add(cigar);
            }
        }
        haplotypeStrings = result;
        return result;
    }

    private List<String> expandAllCombinations(final String cigarString, final String reference) {
        final Civar civar = Civar.fromCharSequence(cigarString);
        final List<Civar> unrolledCivars = civar.optionalizeAll().unroll();
        List<String> result = new ArrayList<>(unrolledCivars.size());
        for (final Civar c : unrolledCivars) {
            result.add(c.applyTo(reference));
        }
        return result;
    }

    private List<Haplotype> expandAllHaplotypeCombinations(final String civarString, final String reference) {
        final Civar civar = Civar.fromCharSequence(civarString);
        final List<Civar> unrolledCivars = civar.optionalizeAll().unroll();
        List<Haplotype> result = new ArrayList<>(unrolledCivars.size());
        for (final Civar c : unrolledCivars) {
            final String baseString = c.applyTo(reference);
            final Haplotype haplotype = new Haplotype(baseString.getBytes(),baseString.equals(reference));
            haplotype.setGenomeLocation(genomeLocParser.createGenomeLoc("1",1,reference.length()));
            try {
            haplotype.setCigar(c.toCigar(reference.length()));
            } catch (final RuntimeException ex) {
                c.applyTo(reference);
                c.toCigar(reference.length());
                throw new RuntimeException("" + c + " " + ex.getMessage(),ex);
            }
            result.add(haplotype);
        }
        return result;
    }


    public List<Haplotype> haplotypeList() {
        if (haplotypeList == null) {

          final List<Haplotype> result = new ArrayList<>(haplotypeCigars.length);
          final String reference = this.reference;
          for (final String cigar : haplotypeCigars) {
              if (cigar.matches("^Civar:.*$")) {
                  stringRepresentation = cigar.substring(6);
                  result.addAll(expandAllHaplotypeCombinations(cigar.substring(6), reference));
              } else if (cigar.matches("^.*\\d+.*$")) {
                  result.add(cigarToHaplotype(reference, cigar, 0, true));
              } else {
                  final Haplotype h = new Haplotype(cigar.getBytes());
                  h.setGenomeLocation(genomeLocParser.createGenomeLoc("chr1",1,reference.length()));
                  result.add(h);
              }
          }
          haplotypeList = result;
        }
        return haplotypeList;
    }


    SAMSequenceDictionary artificialSAMSequenceDictionary() {
        return new SAMSequenceDictionary(Collections.singletonList(new SAMSequenceRecord("00", reference.length())));
    }

    SAMFileHeader artificialSAMFileHeader() {
        return ArtificialReadUtils.createArtificialSamHeader(artificialSAMSequenceDictionary());
    }

    public List<GATKRead> readList() {
        if (readList == null) {
            final SAMFileHeader header = artificialSAMFileHeader();
            readList = new ArrayList<>(readCigars.length);
            final List<String> haplotypes = haplotypesStrings();
            int count = 0;
            for (final String descr : readCigars) {
                String sequence;
                if (descr.matches("^\\d+:\\d+:.+$")) {
                    final String[] parts = descr.split(":");
                    int allele = Integer.parseInt(parts[0]);
                    int offset = Integer.parseInt(parts[1]);
                    final String cigar = parts[2];
                    final String base = allele == 0 ? reference : haplotypes.get(allele - 1);
                    sequence = applyCigar(base, cigar, offset, false);
                    final GATKRead samRecord = ArtificialReadUtils.createArtificialRead(header, "read_" + count, 0, 1, sequence.getBytes(), Arrays.copyOf(bq, sequence.length()));
                    readList.add(samRecord);
                } else if (descr.matches("^\\*:\\d+:\\d+$")) {
                    int readCount = Integer.parseInt(descr.split(":")[1]);
                    int readLength = Integer.parseInt(descr.split(":")[2]);
                    readList.addAll(generateSamRecords(haplotypes, readCount, readLength, header, count));
                } else {
                    sequence = descr;
                    final GATKRead samRecord = ArtificialReadUtils.createArtificialRead(header, "read_" + count, 0, 1, sequence.getBytes(), Arrays.copyOf(bq, sequence.length()));
                    readList.add(samRecord);
                }
                count = readList.size();
            }
        }
        return readList;
    }

    public List<List<Civar.ElementOffset>> readEventOffsetList() {
        if (haplotypeCigars.length != 1 || !haplotypeCigars[0].startsWith("Civar:"))
            throw new UnsupportedOperationException();
        if (readEventOffsetList == null) {
            final Civar civar = Civar.fromCharSequence(haplotypeCigars[0].substring(6));
            final List<Civar> unrolledCivars = civar.optionalizeAll().unroll();

            readEventOffsetList = new ArrayList<>(readCigars.length);
            int count = 0;
            for (final String descr : readCigars) {
                if (descr.matches("^\\d+:\\d+:.+$")) {
                    throw new UnsupportedOperationException();
                } else if (descr.matches("^\\*:\\d+:\\d+$")) {
                    int readCount = Integer.parseInt(descr.split(":")[1]);
                    int readLength = Integer.parseInt(descr.split(":")[2]);
                    readEventOffsetList.addAll(generateElementOffsetRecords(haplotypesStrings(), unrolledCivars, readCount, readLength, count));
                } else {
                    throw new UnsupportedOperationException();
                }
                count = readEventOffsetList.size();
            }
            readEventOffsetList = Collections.unmodifiableList(readEventOffsetList);
        }
        return readEventOffsetList;
    }

    public List<String> readStrings() {
        final List<String> result = new ArrayList<>(readCigars.length);
        final List<String> haplotypes = haplotypesStrings();
        for (final String descr : readCigars) {
            String sequence;
            if (descr.matches("^\\d+:\\d+:.+$")) {
                final String[] parts = descr.split(":");
                int allele = Integer.parseInt(parts[0]);
                int offset = Integer.parseInt(parts[1]);
                final String cigar = parts[2];
                final String base = allele == 0 ? reference : haplotypes.get(allele - 1);
                sequence = applyCigar(base, cigar, offset, false);
                result.add(sequence);
            } else if (descr.matches("\\*:^\\d+:\\d+")) {
                int readCount = Integer.parseInt(descr.split(":")[1]);
                int readLength = Integer.parseInt(descr.split(":")[2]);
                result.addAll(generateReads(haplotypes, readCount, readLength));
            } else {
                sequence = descr;
                result.add(sequence);
            }
        }
        return result;
    }

    private List<String> generateReads(final List<String> haplotypes, final int readCount, final int readLength) {
        final List<String> result = new ArrayList<>(readCount);
        for (int i = 0; i < readCount; i++) {
            int hi = i % haplotypes.size();
            final String h = haplotypes.get(hi);
            int offset = i % h.length() - readLength;
            result.add(h.substring(offset,offset + readLength));
        }
        return result;
    }

    private List<GATKRead> generateSamRecords(final List<String> haplotypes, final int readCount, final int readLength, final SAMFileHeader header, final int idStart) {
        int id = idStart;
        final List<GATKRead> result = new ArrayList<>(readCount);
        for (int i = 0; i < readCount; i++) {
            int hi = i % haplotypes.size();
            final String h = haplotypes.get(hi);
            int offset = h.length() <= readLength ? 0 : i % (h.length() - readLength);
            int to = Math.min(h.length(), offset + readLength);
            byte[] bases = h.substring(offset,to).getBytes();
            byte[] quals = Arrays.copyOf(bq, to - offset);
            final GATKRead samRecord = ArtificialReadUtils.createArtificialRead(header, "read_" + id++, 0, offset + 1, bases, quals);
            result.add(samRecord);
        }
        return result;
    }


    private List<List<Civar.ElementOffset>> generateElementOffsetRecords(final List<String> haplotypes, final List<Civar> unrolledCivars, final int readCount, final int readLength, final int count) {

        final List<List<Civar.ElementOffset>> result = new ArrayList<>(readCount);
        for (int i = 0; i < readCount; i++) {
            int hi = i % unrolledCivars.size();
            final Civar c = unrolledCivars.get(hi);
            final String h = haplotypes.get(hi);
            int offset = h.length() <= readLength ? 0 : i % (h.length() - readLength);
            int to = Math.min(h.length(), offset + readLength);
            result.add(c.eventOffsets(reference,offset,to));
        }
        return result;
    }

    private static final Pattern cigarPattern = Pattern.compile("(\\d+)([=A-Z])");


    private Haplotype cigarToHaplotype(final String reference, final String cigar, final int offset, final boolean global) {
        final String sequence = applyCigar(reference,cigar,offset,global);
        final Haplotype haplotype = new Haplotype(sequence.getBytes(),reference.equals(sequence));
        haplotype.setGenomeLocation(genomeLocParser.createGenomeLoc("chr1",1,reference.length()));
        haplotype.setCigar(Civar.fromCharSequence(cigar).toCigar(reference.length()));
        return haplotype;
    }

    private String applyCigar(final String reference, final String cigar, final int offset, final boolean global) {
        final Matcher pm = cigarPattern.matcher(cigar);
        StringBuffer sb = new StringBuffer();
        int index = offset;
        while (pm.find()) {
            int length = Integer.parseInt(pm.group(1));
            char operator = pm.group(2).charAt(0);
            switch (operator) {
                case '=' :
                    try {
                      sb.append(reference.substring(index, index + length));
                    } catch (Exception e) {
                      throw new RuntimeException(" " + index + " " + (index + length) + " " + reference.length() + " " + cigar,e);
                    }
                    index += length; break;
                case 'D' :
                    index += length; break;
                case 'I' :
                    String insert = cigar.substring(pm.end(),pm.end() + length).toUpperCase();
                    sb.append(insert); break;
                case 'V' :
                    sb.append(transversionV(reference.charAt(index))); index++; break;
                case 'W' :
                        sb.append(transversionW(reference.charAt(index))); index++; break;
                case 'T' :
                    sb.append(transition(reference.charAt(index))); index++; break;
                default:
                    throw new UnsupportedOperationException("cigar operator " + operator + " not supported.");
            }
        }
        if (global && index != reference.length()) {
            throw new RuntimeException(" haplotype cigar does not explain reference length (" + index + " != " + reference.length() + ") on cigar " + cigar);
        } else if (index > reference.length()) {
            throw new RuntimeException(" index beyond end ");
        }
        return sb.toString();
    }

    protected int kmerSize() {
        return kmerSize;
    }

    private char transversionV(final char c) {
        switch (Character.toUpperCase(c)) {
            case 'A': return 'C';
            case 'G': return 'T';
            case 'C': return 'A';
            case 'T': return 'G';
            default:
                return c;
        }

    }

    private char transversionW(final char c) {
        switch (Character.toUpperCase(c)) {
            case 'A': return 'T';
            case 'G': return 'C';
            case 'T': return 'A';
            case 'C': return 'G';
            default:
                return c;
        }

    }

    private char transition(final char c) {
        switch (Character.toUpperCase(c)) {
            case 'A': return 'G';
            case 'G': return 'A';
            case 'T': return 'C';
            case 'C': return 'T';
            default:
                return c;
        }

    }
}
