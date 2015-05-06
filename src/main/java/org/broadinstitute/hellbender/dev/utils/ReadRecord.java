package org.broadinstitute.hellbender.dev.utils;

import com.google.api.services.genomics.model.CigarUnit;
import com.google.api.services.genomics.model.LinearAlignment;
import com.google.api.services.genomics.model.Position;
import com.google.api.services.genomics.model.Read;
import com.google.common.collect.Iterables;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMTagUtil;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * A wrapper around Read that pretends to be a SAMRecord.
 * The advantage of doing it this way is that we can report an error if we attempt to access a field that
 * we don't actually have.
 *
 * Same rules as SAMRecord: don't modify the return values, even if the typesystem would let you.
 */
public class ReadRecord extends SAMRecord {

  private enum Tag {

    AM(Type.INTEGER),
    AS(Type.INTEGER),
    BC(Type.STRING),
    BQ(Type.STRING),
    CC(Type.STRING),
    CM(Type.INTEGER),
    CO(Type.STRING),
    CP(Type.INTEGER),
    CQ(Type.STRING),
    CS(Type.STRING),
    CT(Type.STRING),
    E2(Type.STRING),
    FI(Type.INTEGER),
    FS(Type.STRING),
    FZ(Type.SHORT_ARRAY),
    H0(Type.INTEGER),
    H1(Type.INTEGER),
    H2(Type.INTEGER),
    HI(Type.INTEGER),
    IH(Type.INTEGER),
    LB(Type.STRING),
    MC(Type.STRING),
    MD(Type.STRING),
    MQ(Type.INTEGER),
    NH(Type.INTEGER),
    NM(Type.INTEGER),
    OC(Type.STRING),
    OP(Type.INTEGER),
    OQ(Type.STRING),
    PG(Type.STRING),
    PQ(Type.INTEGER),
    PT(Type.STRING),
    PU(Type.STRING),
    Q2(Type.STRING),
    QT(Type.STRING),
    R2(Type.STRING),
    RG(Type.STRING),
    RT(Type.STRING),
    SA(Type.STRING),
    SM(Type.INTEGER),
    TC(Type.INTEGER),
    U2(Type.STRING),
    UQ(Type.INTEGER);

    private enum Type {

      INTEGER {
        @Override Object convert(List<String> list) {
          return Integer.valueOf(Integer.parseInt(Iterables.getOnlyElement(list)));
        }
      },

      SHORT_ARRAY {
        @Override Object convert(List<String> list) {
          int size = list.size();
          short[] array = new short[size];
          int i = 0;
          for (String string : list) {
            array[i++] = Short.parseShort(string);
          }
          return array;
        }
      },

      STRING {
        @Override Object convert(List<String> list) {
          return Iterables.getOnlyElement(list);
        }
      };

      abstract Object convert(List<String> list);
    }

    private final Type type;

    private Tag(Type type) {
      this.type = type;
    }

    Object convert(List<String> list) {
      return type.convert(list);
    }
  }

  private static byte[] integerListToByteArray(List<Integer> list) {
        int size = list.size();
        byte[] array = new byte[size];
        int i = 0;
        for (Integer j : list) {
          array[i++] = j.byteValue();
        }
        return array;
      }

  private static HashMap<String,CigarOperator> cigarEls = new HashMap<String,CigarOperator>()
    {{
      put("ALIGNMENT_MATCH", htsjdk.samtools.CigarOperator.M);
      put("CLIP_HARD", htsjdk.samtools.CigarOperator.H);
      put("CLIP_SOFT", htsjdk.samtools.CigarOperator.S);
      put("DELETE", htsjdk.samtools.CigarOperator.D);
      put("INSERT", htsjdk.samtools.CigarOperator.I);
      put("PAD", htsjdk.samtools.CigarOperator.P);
      put("SEQUENCE_MATCH", null);
      put("SEQUENCE_MISMATCH", null);
      put("SKIP", htsjdk.samtools.CigarOperator.N);
    }};

  static public htsjdk.samtools.CigarElement genToSam(CigarUnit cu) {
    htsjdk.samtools.CigarOperator op = cigarEls.get(cu.getOperation());
    return new htsjdk.samtools.CigarElement(cu.getOperationLength().intValue(), op);
  }

  static public CigarUnit samToGen(htsjdk.samtools.CigarElement sam) {
    return new CigarUnit()
        .setOperationLength(new Long(sam.getLength()))
        .setOperation(sam.getOperator().toString());
  }

  public ReadRecord(Read read, SAMFileHeader header) {
    super(header);

    // look at the read and fill our own fields.
    String readGroupId = read.getReadGroupId();
    if (null != readGroupId) {
      this.setAttribute(SAMTag.RG.name(), readGroupId);
    }
    String fragmentName = read.getFragmentName();
    if (null != fragmentName) {
        this.setReadName(fragmentName);
      }
    Boolean properPlacement = read.getProperPlacement();
    if (null != properPlacement) {
      this.setProperPairFlag(properPlacement);
      }
    Boolean duplicateFragment = read.getDuplicateFragment();
    if (null != duplicateFragment) {
      this.setDuplicateReadFlag(duplicateFragment);
      }
    Integer fragmentLength = read.getFragmentLength();
    if (null != fragmentLength) {
      this.setInferredInsertSize(fragmentLength);
      }
    Integer numberReads = read.getNumberReads();
    if (null != numberReads && 1 < numberReads) {
      this.setReadPairedFlag(true);
      Integer readNumber = read.getReadNumber();
        if (null != readNumber) {
            switch (readNumber.intValue()) {
              case 0:
                this.setFirstOfPairFlag(true);
                  break;
            case 1:
              this.setSecondOfPairFlag(true);
              break;
               default:
                   throw new IllegalStateException();
              }
          }
      }
    Boolean failedVendorQualityChecks = read.getFailedVendorQualityChecks();
    if (null != failedVendorQualityChecks) {
        this.setReadFailsVendorQualityCheckFlag(failedVendorQualityChecks);
      }
    LinearAlignment alignment = read.getAlignment();
    if (null != alignment) {
      Position position = alignment.getPosition();
      if (null != position) {
        String referenceName = position.getReferenceName();
        if (null != referenceName) {
          this.setReferenceName(referenceName);
        }
        Long alignmentStart = position.getPosition();
        if (null != alignmentStart) {
          // convert from 0-based to 1-based
          this.setAlignmentStart(alignmentStart.intValue()+1);
        }
        Boolean reverseStrand = position.getReverseStrand();
        if (null != reverseStrand) {
          this.setReadNegativeStrandFlag(reverseStrand);
        }
      }
      Integer mappingQuality = alignment.getMappingQuality();
      if (null != mappingQuality) {
        this.setMappingQuality(mappingQuality);
      }
      List<CigarUnit> cigar = alignment.getCigar();
      if (null != cigar) {
        this.setCigar(
            new htsjdk.samtools.Cigar(read.getAlignment().getCigar().stream().map(ReadRecord::genToSam).collect(Collectors.toList())));
      }
    }
    Boolean secondaryAlignment = read.getSecondaryAlignment();
    if (null != secondaryAlignment) {
      this.setNotPrimaryAlignmentFlag(secondaryAlignment);
    }
    Boolean supplementaryAlignment = read.getSupplementaryAlignment();
    if (null != supplementaryAlignment) {
      this.setSupplementaryAlignmentFlag(supplementaryAlignment);
    }
    String alignedSequence = read.getAlignedSequence();
    if (null != alignedSequence) {
      this.setReadString(alignedSequence);
    }
    List<Integer> alignedQuality = read.getAlignedQuality();
    if (null != alignedQuality) {
      this.setBaseQualities(integerListToByteArray(alignedQuality));
    }
    if (read.getNumberReads()>1) {
      // 0x1: template having multiple segments in sequencing
      this.setReadPairedFlag(true);
      // 0x80:  the last segment in the template
      this.setSecondOfPairFlag(read.getReadNumber() + 1 == read.getNumberReads());
      // 0x8: next segment in the template unmapped
      this.setMateUnmappedFlag(read.getNextMatePosition()==null);
    }
    Position nextMatePosition = read.getNextMatePosition();
    if (null != nextMatePosition) {
      String referenceName = nextMatePosition.getReferenceName();
      if (null != referenceName) {
        this.setMateReferenceName(referenceName);
      }
      Long position = nextMatePosition.getPosition();
      if (null != position) {
        // convert from 0-based to 1-based
        this.setMateAlignmentStart(position.intValue()+1);
      }
      Boolean reverseStrand = nextMatePosition.getReverseStrand();
      if (null != reverseStrand) {
        this.setMateNegativeStrandFlag(reverseStrand);
      }
    }
    Map<String, List<String>> info = read.getInfo();
    if (null != info) {
      for (Map.Entry<String, List<String>> entry : info.entrySet()) {
        String key = entry.getKey();
        List<String> value = entry.getValue();
        try {
          Tag tag = Tag.valueOf(key);
          this.setAttribute(
              key, null == tag ? Iterables.getOnlyElement(value) : tag.convert(value));
        } catch (IllegalArgumentException x) {
          // X,Y, or Z tag. Not sure what to do.
        }
      }
    }
    // we leave "validationstringency" at the default value since there is no equivalent in Read.
  }

  public htsjdk.samtools.SAMReadGroupRecord getReadGroup() {
    final String rgId = (String)getAttribute(SAMTagUtil.getSingleton().RG);
    if (rgId == null) {
      return null;
    }
    return getHeader().getReadGroup(rgId);
  }

}