"""Encode confirmed reference calls in a merged SV VCF before PLINK import."""

import re
import sys


MISSING_DIPLOID_GT = re.compile(rb"([/|]?)\.([/|])\.")
LEADING_MISSING_GT = re.compile(rb"(^|\t)([/|]?)\.([/|])\.(?=[:\t]|$)")


def normalize_record(line):
    """Return the VCF line and the number of fully missing diploid GTs changed."""
    if line.startswith(b"#"):
        return line, 0

    body = line.rstrip(b"\r\n")
    columns = body.split(b"\t", 9)
    if len(columns) != 10:
        return line, 0

    formats = columns[8].split(b":")
    if b"GT" not in formats:
        return line, 0

    gt_index = formats.index(b"GT")
    if gt_index == 0:
        columns[9], changed = LEADING_MISSING_GT.subn(rb"\g<1>\g<2>0\g<3>0", columns[9])
    else:
        samples = columns[9].split(b"\t")
        changed = 0
        for index, sample in enumerate(samples):
            fields = sample.split(b":")
            if gt_index >= len(fields):
                continue
            gt = MISSING_DIPLOID_GT.fullmatch(fields[gt_index])
            if gt:
                fields[gt_index] = gt[1] + b"0" + gt[2] + b"0"
                samples[index] = b":".join(fields)
                changed += 1
        columns[9] = b"\t".join(samples)

    return (b"\t".join(columns) + line[len(body):], changed) if changed else (line, 0)


def main():
    changed = 0
    for line in sys.stdin.buffer:
        record, count = normalize_record(line)
        sys.stdout.buffer.write(record)
        changed += count
    print(f"Encoded {changed:,} fully missing diploid GTs as reference.", file=sys.stderr)


if __name__ == "__main__":
    main()
