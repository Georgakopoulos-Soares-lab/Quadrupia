if __name__ == "__main__":

    from pathlib import Path
    import sys
    import argparse
    import csv
    from coverage_extractor import extract_id

    extractions = Path(sys.argv[1]).resolve()
    suffix = sys.argv[2]
    assemblies = sys.argv[3]

    extraction_files = {extract_id(file): file for file in extractions.glob(f"*{suffix}")}
    with open("design.csv", mode="w") as g, open(assemblies, mode="r", encoding="utf-8") as f:
        writer = csv.DictWriter(g, delimiter="\t", fieldnames=["accession_id", "gff", "extraction"])
        writer.writeheader()
        for line in f:
            line = line.strip()
            if line.count("\t") > 0:
                line = line.split("\t")[0]
            accession_id = extract_id(line)
            # accession_id = Path(line).name.split('.gff')[0]
            print(accession_id)
            gff_file = Path(line.replace(".fna", ".gff"))
            if gff_file.is_file() and accession_id in extraction_files:
                writer.writerow({"accession_id": accession_id,
                                 "gff": gff_file,
                                 "extraction": extraction_files[accession_id]
                                })
