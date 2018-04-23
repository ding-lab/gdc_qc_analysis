from collections import namedtuple
import gzip


class MAF:
    def __init__(self, pth):
        self.pth = pth
        if pth.suffix == '.gz':
            self._file = gzip.open(pth, 'rt')
        else:
            self._file = open(pth)

        # Header is the comments before column
        self.header = []
        raw_columns = self.read_header()

        # Set up columns
        self.columns = self.make_columns(raw_columns)
        self._record = self.make_record_class()

    def read_header(self):
        line = next(self._file)
        while line.startswith('#'):
            self.header.append(line.rstrip())
            line = next(self._file)

        # Treat the first noncomment line as columns
        return line.rstrip().split('\t')

    def make_columns(self, raw_columns):
        return [c.lower() for c in raw_columns]

    def make_record_class(self):
        return namedtuple(f'MAFRecord', self.columns)

    def make_record(self, vals):
        """Given the MAF record values from a row, construct the record"""
        return self._record._make(vals)

    def __iter__(self):
        return self

    def __next__(self):
        cols = next(self._file).rstrip().split('\t')
        return self.make_record(cols)


class MC3MAF(MAF):
    def make_columns(self, raw_columns):
        # Strand is duplicated
        renamed_cols = []
        for c in raw_columns:
            if c == 'STRAND':
                renamed_cols.append('strand_vep')
            else:
                renamed_cols.append(c.lower())
        return renamed_cols


class GDCMAF(MAF):
    """GDCMAF will detect"""
    def __init__(self, pth):
        super().__init__(pth)
        _, cancer_type, caller, *__ = self.pth.name.split('.')
        self.cancer_type = cancer_type
        self.caller = caller

    def make_columns(self, raw_columns):
        columns = super().make_columns(raw_columns)
        # Add cancer type and caller information
        columns.extend(['cancer_type', 'caller'])
        return columns

    def __next__(self):
        cols = next(self._file).rstrip().split('\t')
        return self._record(*cols, self.cancer_type, self.caller)
