import logging
import vcf
from typing import List, Tuple

_logger = logging.getLogger(__name__)


# TODO: for now I'm going to do the lazy thing and just traverse the VCF each time for each sample
def read_sample_segments_and_calls(intervals_vcf: str,
                                   clustered_vcf: str,
                                   sample_name: str,
                                   contig: str) -> List[Tuple[int, int, int]]:
    """
    Get the segmentation "path" to use for calculating qualities based on the VCF with clustered breakpoints
    :param intervals_vcf:
    :param clustered_vcf:
    :param sample_name:
    :param contig:
    :return: {copy number, start index, stop index (inclusive)}
    """
    intervals = vcf.Reader(filename=intervals_vcf)
    intervals2 = vcf.Reader(filename=intervals_vcf)
    segments = vcf.Reader(filename=clustered_vcf)

    path: List[Tuple[int, int, int]] = []
    segment_start_index = 0
    segment_end_index = 0

    # A record corresponds to [CHROM,POS,REF,ALT]
    try:
        interval_start_iter = iter(intervals.fetch(contig))
        interval_end_iter = iter(intervals2.fetch(contig))
    except ValueError:
        print('ERROR: could not fetch intervals')
        raise
    else:
        start_interval = next(interval_start_iter)
        end_interval = next(interval_end_iter)
        intervals_copy_number = try_getting_format_attribute(end_interval, sample_name, 'CN')

    try:
        segments_iter = iter(segments.fetch(contig))
    except ValueError:
        print('WARN: no segments found on contig {0}'.format(contig))
        return path
    else:
        segments_rec = next(segments_iter)
        segment_copy_number = try_getting_format_attribute(segments_rec, sample_name, 'CN')

    # we assume segments are sorted by start, but may be overlapping
    while segments_rec is not None and start_interval is not None:
        # make sure interval start matches
        while start_interval is not None and start_interval.POS < segments_rec.POS:
            try:
                start_interval = next(interval_start_iter)
                segment_start_index += 1
                end_interval = next(interval_end_iter)
                segment_end_index += 1
            except StopIteration:
                print('ERROR: ran out of intervals with unmatched segments remaining')
                raise
        # once start matches, move the interval end
        while end_interval is not None and try_getting_info_attribute(segments_rec, 'END') > \
                try_getting_info_attribute(end_interval, 'END'):
            try:
                end_interval = next(interval_end_iter)
                segment_end_index += 1
                intervals_copy_number = try_getting_format_attribute(end_interval, sample_name, 'CN')
            except StopIteration:
                print('WARN: ran out of intervals with segment end unmatched')
                end_interval = None

        # add the segment
        if segment_end_index < segment_start_index:
            print('Sample {0} contains segment at {1}:{2} with end index greater than start index'.format(sample_name, contig, segments_rec.POS))
        path.append((segment_copy_number, segment_start_index, segment_end_index))
        # do this the dumb way because each reader gets the same iterator
        segment_end_index = 0
        interval_end_iter = iter(intervals2.fetch(contig))
        end_interval = next(interval_end_iter)

        # get the next segment
        try:
            segments_rec = next(segments_iter)
            segment_copy_number = try_getting_format_attribute(segments_rec, sample_name, 'CN')
        except StopIteration:
            segments_rec = None
            segments_iter = None

    return path


def try_getting_info_attribute(record,
                               attribute: str) -> int:
    try:
        value = record.INFO[attribute]
    except AttributeError:
        print('No {} field for record at position:{}'.format(attribute, record.POS))
    else:
        return value


def try_getting_format_attribute(record,
                                 sample_name: str,
                                 attribute: str) -> int:
    try:
        value = record.genotype(sample_name)[attribute]
    except AttributeError:
        print('No {} field for {} intervals at position:{}'.format(attribute, sample_name, record.POS))
    else:
        return value
