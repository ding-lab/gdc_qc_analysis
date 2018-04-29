import argparse
import gzip
from pathlib import Path
from maf_utils import MC3MAF


def main(args):
    maf_reader = MC3MAF(Path(args.maf_pth))
    with gzip.open(args.out_pth, 'wt') as f:
        for record in maf_reader:
            chrom, start, end = \
                record.chromosome, int(record.start_position), int(record.end_position)
            # convert to 0-based half open interval
            start -= 1
            print(chrom, start, end, sep="\t", file=f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract genomic location from given MAF to gzip'd BED format."
    )
    parser.add_argument('maf_pth', help="Path to the MAF file")
    parser.add_argument('out_pth', help="Path to the output gzip'd BED file")
    args = parser.parse_args()
    main(args)
