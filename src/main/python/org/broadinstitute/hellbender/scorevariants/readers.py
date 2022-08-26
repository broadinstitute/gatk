from pysam import VariantRecord
from pysam import VariantFile, AlignmentFile
from Bio import SeqIO
from scorevariants.encoders import (
    AnnotationEncoder,
    ReferenceEncoder,
    VariantTypeEncoder,
    ReadTensorEncoder,
)
import torch
from typing import Dict, List, Optional
from functools import partial


def parse_reference_dict(reference_file: str):
    """Parses a FASTA File into a dictionary

    Args:
        reference_file: str or File-like object in fasta format

    """
    return SeqIO.to_dict(SeqIO.parse(reference_file, format="fasta"))


def determine_polymorphic_type(variant: VariantRecord):
    """
    Motivated by gatk htsjdk:
    https://github.com/samtools/htsjdk/blob/b5af659e6c71addab4baeed297aadf4468d94bd0/src/main/java/htsjdk/variant/variantcontext/VariantContext.java#L1468

    """
    variant_type = None
    for allele in variant.alts or []:
        if allele == variant.ref:
            continue

        biallelic_type = _type_of_biallelic_variant(variant.ref, allele)

        if variant_type is None:
            variant_type = biallelic_type
        elif biallelic_type != variant_type:
            return "MIXED"

        if variant_type is None:
            return "NO_VARIATION"

    return variant_type


def _type_of_biallelic_variant(ref, allele):
    """
    Determine the variant type of a single allele.
    """
    if len(ref) == len(allele):
        if len(allele) == 1:
            return "SNP"
        else:
            return "MNP"
    return "INDEL"

_torch_dtypes = {
    "type": torch.long,
    "reference": torch.float32,
    "best_practices": torch.float32,
    "read_tensor": torch.float32,
}

    
class _BaseReader:
    def __init__(self, variant_file, annotation_encoder, type_encoder, reference_dict):
        self.variant_file = variant_file
        if annotation_encoder is None:
            annotation_encoder = AnnotationEncoder()
        self.annotation_encoder = annotation_encoder
        if type_encoder is None:
            type_encoder = VariantTypeEncoder()
        self.type_encoder = type_encoder
        self.reference = reference_dict
        self.determine_type = determine_polymorphic_type
        
    
    def get_variants(self) -> List[VariantRecord]:
        """Gets all variants in the variant file

        Returns:
            list of VariantRecord containing all variants
            from the vcf

        """
        self.variant_file.reset()
        variants = [v for v in self.variant_file]
        self.variant_file.reset()
        return variants
    
    
    def __iter__(self):
        """
        Iterates over the variants in the file.

        """
        for variant in self.variant_file:
            yield self.variant_representation(variant)
            
        
    def variant_representation(
            self, variant: VariantRecord, alt_sep: str = None
            ) -> Dict:
        """Gets the representation and metadata for a variant

        Args:
            variant:
                A pysam.VariantRecord containing the variant to encode
            alt_sep:
                Optional str for separating the alts. default = ', '
        Returns:
            A dict containing the following keys:
            * chrom : str for the chromosome of the variant
            * pos : int for the position of the variant (1-based)
            * ref : str for the reference sequence of the variant
            * alt : str for the alt alleles
            * type : torch.long for the variant type (SNP/INDEL)
            * best_practices: torch.float32 for the annotation vector
            * reference: torch.float32 for the reference tensor

        """
        if alt_sep is None:
            alt_sep = ", "
        chrom = variant.chrom
        # convert from pysam 0-based index to 1-based index
        pos = variant.pos
        ref = variant.ref
        alt = alt_sep.join(variant.alts or [])
        type_ = self.type_encoder(self.determine_type(variant))
        annotations = self.annotation_encoder(variant)
        representations = [
            ("chrom", chrom),
            ("pos", pos),
            ("ref", ref),
            ("alt", alt),
            ("type", type_),
            ("best_practices", annotations),
        ]
        specific_representations = self._get_specific_representations(variant)
        representations.extend(specific_representations)

        out_values = dict()
        for key, value in representations:
            if key in _torch_dtypes:
                value = torch.tensor(value, dtype=_torch_dtypes[key])
            out_values[key] = value

        return out_values
    

class ReferenceTensorReader(_BaseReader):
    def __init__(
        self,
        variant_file: VariantFile,
        reference_dict: Dict[str, str],
        reference_encoder=None,
        annotation_encoder=None,
        variant_type_encoder=None,
    ):
        """Reads variants from a vcf and and encodes
        the reference sequence as a tensor, as well as their metadata.

        Args:
            variant_file:
                pysam VariantFile containing variants to encode
            reference_dict:
                dict containing the sequence (value) for
                each contig (key)
            reference_encoder:
                Optional Encoder to use encode the reference sequence
            annotation_encoder:
                Optional Encoder used to encode the variant annotations
            variant_type_encoder:
                Optional Encoder used to encode the variant type
                (e.g., SNP v. INDEL)

        Attributes:
            determine_type:
                function used to determine the type of the variant.
                returns a string.

        """
        super().__init__(variant_file, annotation_encoder, variant_type_encoder, reference_dict)

        if reference_encoder is None:
            reference_encoder = ReferenceEncoder()
        self.reference_encoder = reference_encoder


    @classmethod
    def from_files(cls, vcf: str, fasta_file: str, **kwargs):
        """Creates a ReferenceTensorReader from paths to files

        Args:
            vcf: str or file-like object to vcf file
            fasta_file : str of file-like object to fasta file
            **kwargs : arugments for ReferenceTensorReader constructor

        Returns:
            An instnace of ReferenceTensorReader

        """
        variant_file = VariantFile(vcf)
        reference_dict = parse_reference_dict(fasta_file)
        return ReferenceTensorReader(variant_file, reference_dict, **kwargs)

    def _get_specific_representations(self, variant):
        reference_tensor = self.reference_encoder(self.reference, variant)
        return [("reference", reference_tensor) ]
    
    
class TensorReader(_BaseReader):
    def __init__(
        self,
        variant_file: VariantFile,
        reference_dict: Dict[str, str],
        read_encoder=None,
        reference_encoder=None,
        annotation_encoder=None,
        variant_type_encoder=None,
    ):
        """Reads variants from a vcf and and encodes
        the reference sequence as a tensor, as well as their metadata.

        Args:
            variant_file:
                pysam VariantFile containing variants to encode
            reference_dict:
                dict containing the sequence (value) for
                each contig (key)
            read_encoder:
                Optional Encoder to use to encode the reads to tensor
            reference_encoder:
                Optional Encoder to use encode the reference sequence
            annotation_encoder:
                Optional Encoder used to encode the variant annotations
            variant_type_encoder:
                Optional Encoder used to encode the variant type
                (e.g., SNP v. INDEL)

        Attributes:
            determine_type:
                function used to determine the type of the variant.
                returns a string.

        """
        super().__init__(variant_file, annotation_encoder, variant_type_encoder, reference_dict)
        self.additional_encoders = {
            'reference': partial(reference_encoder, self.reference) if callable(reference_encoder) else None,
            'read_tensor': read_encoder,
        }

    @classmethod
    def from_files(cls, vcf: str, fasta_file: str, alignment_file: Optional[str] = None, mode='rb', force_include_reference=False, **kwargs):
        """Creates a ReferenceTensorReader from paths to files

        Args:
            vcf: str or file-like object to vcf file
            fasta_file : str of file-like object to fasta file
            alignemnt_file: str or file-like object to an alignment file
            mode: str. mode for opening alignment file. default = 'rb'
            force_include_reference: bool. if alignment_file is supplied,
                still read the reference tensor
            **kwargs : arugments for ReferenceTensorReader constructor

        Returns:
            An instnace of TensorReader. Read tensors are given in alignment file
            is specified. Otherwise, Reference tensors are given.

        """
        variant_file = VariantFile(vcf)
        reference_dict = parse_reference_dict(fasta_file)
        
        def _add_ref_encoder(kw):
            ref_encoder = ReferenceEncoder()
            kw['reference_encoder'] = ref_encoder
        
        cls_kwargs = dict()
        if alignment_file is None:
            _add_ref_encoder(cls_kwargs)
        else:
            alignment_object = AlignmentFile(alignment_file, mode)
            read_encoder = ReadTensorEncoder(alignment_object, reference_dict)
            cls_kwargs['read_encoder'] = read_encoder
            if force_include_reference:
                _add_ref_encoder(cls_kwargs)

        return cls(
                variant_file,
                reference_dict,
                **cls_kwargs,
                **kwargs,
            )

    def _get_specific_representations(self, variant):
        return [(k, encode(variant)) for k, encode in self.additional_encoders.items() if encode is not None]
    
    
