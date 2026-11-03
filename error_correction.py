from __future__ import annotations

from random import Random
from typing import Iterable, List


class HelixVaultEncodingError(ValueError):
    """Raised when DNA encoding or decoding fails."""


class BiologicalConstraintError(HelixVaultEncodingError):
    """Raised when the DNA payload violates biological constraints."""


class DNAEncoder:
    """Convert binary data to DNA bases and enforce biological constraints.

    Base mapping:
    - 00 -> A
    - 01 -> C
    - 10 -> G
    - 11 -> T
    """

    BASE_TO_BITS = {"A": 0b00, "C": 0b01, "G": 0b10, "T": 0b11}
    BITS_TO_BASE = {0b00: "A", 0b01: "C", 0b10: "G", 0b11: "T"}
    ALTERNATE_BASES = {"A": ("C", "G", "T"), "C": ("A", "G", "T"), "G": ("A", "C", "T"), "T": ("A", "C", "G")}

    @classmethod
    def _validate_dna(cls, sequence: str) -> str:
        normalized = sequence.upper()
        invalid = set(normalized) - set(cls.BASE_TO_BITS)
        if invalid:
            raise HelixVaultEncodingError(
                f"Sequence contains invalid DNA bases: {sorted(invalid)}"
            )
        return normalized

    @staticmethod
    def _byte_to_bits(data: bytes) -> List[int]:
        bits: List[int] = []
        for byte in data:
            for shift in range(7, -1, -1):
                bits.append((byte >> shift) & 1)
        return bits

    @staticmethod
    def _bits_to_byte(bits: Iterable[int]) -> int:
        value = 0
        for bit in bits:
            value = (value << 1) | int(bit)
        return value

    @classmethod
    def bits_to_dna(cls, bits: Iterable[int]) -> str:
        sequence: List[str] = []
        bit_buffer: List[int] = []
        for bit in bits:
            bit_buffer.append(int(bit))
            if len(bit_buffer) == 2:
                pair = (bit_buffer[0] << 1) | bit_buffer[1]
                sequence.append(cls.BITS_TO_BASE[pair])
                bit_buffer = []
        if bit_buffer:
            raise HelixVaultEncodingError(
                "Bit stream length is odd and cannot be mapped to DNA bases."
            )
        return "".join(sequence)

    @classmethod
    def dna_to_bits(cls, sequence: str) -> List[int]:
        normalized = cls._validate_dna(sequence)
        bits: List[int] = []
        for base in normalized:
            pair = cls.BASE_TO_BITS[base]
            bits.append((pair >> 1) & 1)
            bits.append(pair & 1)
        return bits

    @staticmethod
    def _gc_content(sequence: str) -> float:
        if not sequence:
            return 0.0
        gc_count = sum(1 for base in sequence if base in {"G", "C"})
        return gc_count / len(sequence)

    @staticmethod
    def _has_homopolymer_violation(sequence: str) -> bool:
        if len(sequence) < 4:
            return False
        current = sequence[0]
        run_length = 1
        for base in sequence[1:]:
            if base == current:
                run_length += 1
                if run_length > 3:
                    return True
            else:
                current = base
                run_length = 1
        return False

    @staticmethod
    def _run_length(sequence: List[str], index: int) -> int:
        base = sequence[index]
        left = index
        while left > 0 and sequence[left - 1] == base:
            left -= 1
        right = index
        while right + 1 < len(sequence) and sequence[right + 1] == base:
            right += 1
        return right - left + 1

    @classmethod
    def _base_index(cls, base: str) -> int:
        return cls.BASE_TO_BITS[base.upper()]

    @classmethod
    def _index_to_base(cls, value: int) -> str:
        return cls.BITS_TO_BASE[value & 0b11]

    @classmethod
    def _scramble_sequence(cls, sequence: str, seed: int) -> str:
        rng = Random(seed)
        out: List[str] = []
        for index, base in enumerate(sequence):
            base_index = cls._base_index(base)
            mask = rng.getrandbits(2)
            new_index = base_index ^ mask
            out.append(cls._index_to_base(new_index))
        return "".join(out)

    @classmethod
    def _repair_biological_constraints(cls, sequence: str, seed: int) -> str:
        candidate = list(sequence)
        for attempt in range(1, 2000):
            gc_ratio = cls._gc_content("".join(candidate))
            if 0.40 <= gc_ratio <= 0.60 and not cls._has_homopolymer_violation(
                "".join(candidate)
            ):
                return "".join(candidate)

            rng = Random(seed + attempt)
            for index, base in enumerate(candidate):
                run = cls._run_length(candidate, index)
                if run > 3:
                    options = ["A", "C", "G", "T"]
                    options = [item for item in options if item != base]
                    candidate[index] = rng.choice(options)
                    continue

                gc_ratio = cls._gc_content("".join(candidate))
                if gc_ratio < 0.40 and base in {"A", "T"}:
                    candidate[index] = "G" if rng.random() < 0.5 else "C"
                elif gc_ratio > 0.60 and base in {"G", "C"}:
                    candidate[index] = "A" if rng.random() < 0.5 else "T"

            # Deterministic XOR-like scrambling from the seed keeps the process reproducible.
            candidate = list(cls._scramble_sequence("".join(candidate), seed + attempt * 7))

        raise BiologicalConstraintError(
            "Unable to satisfy the GC-content and homopolymer constraints on the DNA payload."
        )

    @classmethod
    def enforce_constraints(cls, sequence: str, seed: int = 0) -> str:
        sanitized = cls._validate_dna(sequence)
        if not sanitized:
            return ""
        gc_ratio = cls._gc_content(sanitized)
        if 0.40 <= gc_ratio <= 0.60 and not cls._has_homopolymer_violation(sanitized):
            return sanitized
        return cls._repair_biological_constraints(sanitized, seed=seed)

    @classmethod
    def bytes_to_dna(cls, payload: bytes, seed: int = 0) -> str:
        bits = cls._byte_to_bits(payload)
        dna = cls.bits_to_dna(bits)
        if not dna:
            return ""
        return cls.enforce_constraints(dna, seed=seed)

    @classmethod
    def dna_to_bytes(cls, sequence: str) -> bytes:
        bits = cls.dna_to_bits(sequence)
        if not bits:
            return b""
        if len(bits) % 8 != 0:
            bits = bits[: -(len(bits) % 8)]
        byte_values: List[int] = []
        for index in range(0, len(bits), 8):
            byte_values.append(cls._bits_to_byte(bits[index : index + 8]))
        return bytes(byte_values)


def bits_to_dna(bits: Iterable[int]) -> str:
    """Public helper to convert bits into DNA."""
    return DNAEncoder.bits_to_dna(bits)


def dna_to_bits(sequence: str) -> List[int]:
    """Public helper to convert DNA back into bits."""
    return DNAEncoder.dna_to_bits(sequence)


def bytes_to_dna(payload: bytes, seed: int = 0) -> str:
    """Public helper to encode raw bytes into DNA."""
    return DNAEncoder.bytes_to_dna(payload, seed=seed)


def dna_to_bytes(sequence: str) -> bytes:
    """Public helper to decode DNA back into bytes."""
    return DNAEncoder.dna_to_bytes(sequence)


__all__ = [
    "DNAEncoder",
    "HelixVaultEncodingError",
    "BiologicalConstraintError",
    "bytes_to_dna",
    "dna_to_bits",
    "dna_to_bytes",
    "bits_to_dna",
]


if __name__ == "__main__":
    payload = b"HelixVault"
    encoded = DNAEncoder.bytes_to_dna(payload, seed=10)
    recovered = DNAEncoder.dna_to_bytes(encoded)
    print(encoded)
    print(recovered)
    assert recovered == payload
