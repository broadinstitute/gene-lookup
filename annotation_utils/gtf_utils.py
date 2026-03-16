import gzip

VALID_FEATURE_TYPES = {"gene", "transcript", "CDS", "UTR", "exon", "start_codon", "stop_codon"}


def parse_gtf(gtf_path, feature_type=None):
    """Parse a gtf file and return a generator of records

    Args:
        gtf_path (str): path of the gtf file
        feature_type (str): if not None, only keep features of this type. Allowed values are:
            'gene', 'transcript', 'CDS', 'UTR', 'exon', 'start_codon', 'stop_codon'
    """

    if feature_type and feature_type not in VALID_FEATURE_TYPES:
        raise ValueError(f"Invalid feature_type: {feature_type}. Expecting one of: {VALID_FEATURE_TYPES}")

    fopen = gzip.open if gtf_path.endswith("gz") else open
    with fopen(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                print(f"WARNING: unable to parse line: {line}")
                continue

            if feature_type is not None and fields[2] != feature_type:
                continue

            record = {
                "chrom": fields[0],
                "source": fields[1],
                "feature": fields[2],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
            }

            info = fields[8]
            record.update({k: v.strip(';" ') for k, v in [x.split(" ") for x in info.split("; ")]})

            yield record
