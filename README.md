# HelixVault

HelixVault is a DNA-based steganography and data integrity tool designed as a production-ready Python CLI for encoding and recovering secret payloads within FASTA-hosted genomic sequences.

The system combines:
- DNA bit encoding using the required A/C/G/T mapping,
- biological constraint enforcement for GC content and homopolymer control,
- Hamming-style parity protection for mutation resilience,
- FASTA steganographic concealment in non-coding or junk DNA regions,
- CLI commands for encode, corrupt, and extract workflows.

## Repository structure

- `dna_encoder.py` — core DNA bitstream encoding, decoding, and biological validation.
- `error_correction.py` — Hamming-based DNA parity correction layer.
- `steg_engine.py` — host-sequence processing, junk-region detection, and payload hiding.
- `simulator.py` — biological mutation simulator.
- `main.py` — polished command-line interface.
- `tests/test_helixvault.py` — end-to-end verification of the encode/hide/corrupt/extract cycle.

## Quick start

Encode a file into a host FASTA and output a steganographic FASTA alongside a retrieval key:

```bash
python main.py encode-hide secret.txt host.fasta stego.fasta key.json --seed 42
```

Corrupt the steganographic FASTA to simulate biological mutations:

```bash
python main.py corrupt stego.fasta corrupted.fasta --mutation-rate 0.02 --seed 7
```

Recover the secret file from the corrupted FASTA:

```bash
python main.py extract-repair corrupted.fasta key.json restored_secret.txt
```

## Notes

The system intentionally uses only the Python standard library to keep the implementation portable and easy to audit.
