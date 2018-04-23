import argparse
from collections import namedtuple
import gzip


def maf_reader(pth):
    with open(pth) as f:
        header = next(f)[:-1].split('\t')
        # Generate the record based on the header
        Record = namedtuple('Record', header)
        for line in f:
            yield Record._make(line[:-1].split('\t'))


def main(args):
    with gzip.open(args.out_pth,'wt') as f:
        for record in maf_reader(args.maf_pth):
            chrom, start, end = \
                record.Chromosome, int(record.Start_Position), int(record.End_Position)
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
