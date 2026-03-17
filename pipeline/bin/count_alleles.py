#!/usr/bin/env python3
"""Count reference and alternative alleles at variant sites.

Usage: count_alleles.py <sorted.bam> <variants.vcf.gz> <reference.fa> <output.csv> [--threads N]

Replaces GATK3 ASEReadCounter with a pysam-based implementation.
For each variant site in the VCF, counts reads supporting the reference
and alternative alleles from the BAM pileup.

Output columns match GATK ASEReadCounter format:
  contig, position, variantID, refAllele, altAllele,
  refCount, altCount, totalCount, lowMAPQDepth, lowBaseQDepth,
  rawDepth, otherBases, improperPairs
"""

import argparse
from concurrent.futures import ProcessPoolExecutor

import pysam


def _process_chunk(chunk_args):
    """Process a chunk of VCF records in a worker process."""
    bam_path, ref_path, records, min_mapq, min_baseq = chunk_args
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref = pysam.FastaFile(ref_path)
    results = []

    for chrom, pos, variant_id, ref_allele, alt_allele in records:
        ref_count = 0
        alt_count = 0
        other_count = 0
        low_mapq = 0
        low_baseq = 0
        raw_depth = 0

        # pysam pileup uses 0-based coordinates
        for pileup_col in bam.pileup(
            chrom,
            pos - 1,
            pos,
            truncate=True,
            stepper="all",
            min_base_quality=0,
            min_mapping_quality=0,
        ):
            if pileup_col.reference_pos != pos - 1:
                continue

            for read in pileup_col.pileups:
                if read.is_del or read.is_refskip:
                    continue

                raw_depth += 1
                alignment = read.alignment

                if alignment.mapping_quality < min_mapq:
                    low_mapq += 1
                    continue

                base_qual = alignment.query_qualities[read.query_position]
                if base_qual < min_baseq:
                    low_baseq += 1
                    continue

                base = alignment.query_sequence[read.query_position].upper()
                if base == ref_allele.upper():
                    ref_count += 1
                elif base == alt_allele.upper():
                    alt_count += 1
                else:
                    other_count += 1

        total_count = ref_count + alt_count

        results.append(
            f"{chrom}\t{pos}\t{variant_id}\t{ref_allele}\t{alt_allele}\t"
            f"{ref_count}\t{alt_count}\t{total_count}\t{low_mapq}\t{low_baseq}\t"
            f"{raw_depth}\t{other_count}\t0"
        )

    bam.close()
    ref.close()
    return results


def count_alleles(bam_path, vcf_path, ref_path, out_path, min_mapq=10, min_baseq=20, threads=1):
    vcf = pysam.VariantFile(vcf_path)

    # Collect all VCF records
    records = []
    for rec in vcf:
        chrom = rec.chrom
        pos = rec.pos  # 1-based in VCF
        ref_allele = rec.alleles[0]
        alt_alleles = [a for a in rec.alleles[1:] if a != "."]
        if not alt_alleles:
            continue
        alt_allele = alt_alleles[0]

        # Only handle SNPs
        if len(ref_allele) != 1 or len(alt_allele) != 1:
            continue

        variant_id = rec.id if rec.id else "."
        records.append((chrom, pos, variant_id, ref_allele, alt_allele))

    vcf.close()

    if threads <= 1 or len(records) == 0:
        # Single-threaded: no subprocess overhead
        all_results = _process_chunk((bam_path, ref_path, records, min_mapq, min_baseq))
    else:
        # Split records into chunks for parallel processing
        chunk_size = max(1, len(records) // threads)
        chunks = []
        for i in range(0, len(records), chunk_size):
            chunk = records[i : i + chunk_size]
            chunks.append((bam_path, ref_path, chunk, min_mapq, min_baseq))

        all_results = []
        with ProcessPoolExecutor(max_workers=threads) as executor:
            for chunk_results in executor.map(_process_chunk, chunks):
                all_results.extend(chunk_results)

    with open(out_path, "w") as out:
        out.write(
            "contig\tposition\tvariantID\trefAllele\taltAllele\t"
            "refCount\taltCount\ttotalCount\tlowMAPQDepth\tlowBaseQDepth\t"
            "rawDepth\totherBases\timproperPairs\n"
        )
        for line in all_results:
            out.write(line + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count alleles at variant sites")
    parser.add_argument("bam", help="Sorted BAM file")
    parser.add_argument("vcf", help="VCF file with variant sites")
    parser.add_argument("ref", help="Reference FASTA")
    parser.add_argument("output", help="Output CSV path")
    parser.add_argument("--threads", type=int, default=1, help="Number of parallel workers")
    args = parser.parse_args()

    count_alleles(args.bam, args.vcf, args.ref, args.output, threads=args.threads)
