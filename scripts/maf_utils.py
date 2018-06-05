from collections import namedtuple
from pathlib import Path
import gzip


class MAF:
    """
    General purpose of MAF reader.

    It assumes the file begins wtih some comments lines (starting with `#`),
    then the column header, and the actual variant records in TSV format.

    Arguments:
        pth (pathlib.Path or str): Path object to the MAF file.
    """
    def __init__(self, pth):
        pth = Path(pth)
        self.pth = pth
        if pth.suffix == '.gz':
            self._file = gzip.open(str(pth), 'rt')
        else:
            self._file = open(str(pth))

        # A reader wrapping the underlying file object
        # which also returns the line number
        self._reader = enumerate(self._file, 1)

        # Header comments appear before column
        self.header_comments = []
        self.raw_columns = self.read_header()

        # Set up columns
        self.columns = self.make_columns(self.raw_columns)
        self._record_cls = self.make_record_class()

    def read_header(self):
        """
        Read the header comments and return the parsed the column header.
        """
        line_no, line = next(self._reader)
        while line.startswith('#'):
            self.header_comments.append(line.rstrip('\n'))
            line_no, line = next(self._reader)

        # Treat the first noncomment line as columns
        return line.rstrip('\n').split('\t')

    def make_columns(self, raw_columns):
        """Define the columns a variant record should store."""
        return [c.lower() for c in raw_columns]

    def make_record_class(self):
        """Define the record class."""
        return namedtuple(f'MAFRecord', self.columns)

    def make_record(self, vals):
        """Given the MAF record values from a row, construct the record"""
        return self._record_cls._make(vals)

    def __iter__(self):
        return self

    def __next__(self):
        line_no, line = next(self._reader)
        cols = line.rstrip('\n').split('\t')
        return self.make_record(cols)


class MC3MAF(MAF):
    """
    MC3 MAF reader.

    It renames the VEP strand column to strand_vep, and renames the chromosome
    name to include 'chr' prefix.
    """
    def make_columns(self, raw_columns):
        # Strand is duplicated
        renamed_cols = []
        for c in raw_columns:
            if c == 'STRAND':
                renamed_cols.append('strand_vep')
            else:
                renamed_cols.append(c.lower())
        # Add line number as a column
        renamed_cols.append('raw_file_line_number')
        return renamed_cols

    def __next__(self):
        line_no, line = next(self._reader)
        cols = line.rstrip('\n').split('\t')
        r = self._record_cls(*cols, line_no)
        # Rename chromosome
        r = r._replace(chromosome=f'chr{r.chromosome}')
        return r


class GDCMAF(MAF):
    """
    GDC MAF reader

    Based on the file name, it will determine the cancer type and variant caller
    of the MAF and passes them as additional columns.
    """
    def __init__(self, pth):
        super().__init__(pth)
        _, cancer_type, caller, *__ = self.pth.name.split('.')
        self.cancer_type = cancer_type
        self.caller = caller

    def make_columns(self, raw_columns):
        columns = super().make_columns(raw_columns)
        # Add cancer type and caller information
        columns.extend(['cancer_type', 'caller', 'raw_file_line_number'])
        return columns

    def __next__(self):
        line_no, line = next(self._reader)
        cols = line.rstrip('\n').split('\t')
        return self._record_cls(*cols, self.cancer_type, self.caller, line_no)
