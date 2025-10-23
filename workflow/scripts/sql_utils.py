#!/usr/bin/env python3
"""Utility helpers for accessing the pcDBG SQLite databases."""

from __future__ import annotations

import sqlite3
from functools import lru_cache
from itertools import islice
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple

from loguru import logger

SQLITE_MAX_VARIABLES = 999


def _chunked(iterable: Sequence[str], size: int) -> Iterator[List[str]]:
    """Yield ``size``-sized lists from *iterable* while preserving order."""
    it = iter(iterable)
    while True:
        window = list(islice(it, size))
        if not window:
            return
        yield window


def _normalise_ids(unitig_ids: Iterable[str]) -> List[str]:
    """Return a list of unitig identifiers with duplicates removed."""
    seen = set()
    ordered_ids: List[str] = []
    for unitig_id in unitig_ids:
        if unitig_id in seen:
            continue
        seen.add(unitig_id)
        ordered_ids.append(unitig_id)
    return ordered_ids


def lookup_unitigs_batch(
    conn: sqlite3.Connection,
    unitig_ids: Iterable[str],
    output_neighbors: bool = True,
) -> Dict[str, Dict[str, object]]:
    """
    Retrieve metadata for *unitig_ids* from ``unitigs`` (and optionally
    ``adjacency``) tables.

    Parameters
    ----------
    conn:
        SQLite connection.
    unitig_ids:
        Iterable of unitig identifiers to materialise.
    output_neighbors:
        When ``True`` include adjacency lists for both directions.

    Returns
    -------
    dict
        Mapping of unitig id to ``{"length": int, "colors": str,
        "neighbors": {"+": [(neighbor_id, neighbor_dir), ...], "-": [...]}}``.
    """

    ids = _normalise_ids(unitig_ids)
    if not ids:
        return {}

    cursor = conn.cursor()
    unitig_info: Dict[str, Dict[str, object]] = {
        unitig_id: {
            "length": None,
            "colors": None,
            "neighbors": {"+": [], "-": []},
        }
        for unitig_id in ids
    }

    for chunk in _chunked(ids, SQLITE_MAX_VARIABLES):
        placeholders = ", ".join("?" for _ in chunk)
        cursor.execute(
            f"SELECT unitig_id, length, colors FROM unitigs "
            f"WHERE unitig_id IN ({placeholders})",
            chunk,
        )
        for unitig_id, length, colors in cursor.fetchall():
            unitig_info[unitig_id]["length"] = length
            unitig_info[unitig_id]["colors"] = colors

    missing = [uid for uid, meta in unitig_info.items() if meta["length"] is None]
    if missing:
        preview = ", ".join(missing[:5])
        suffix = "..." if len(missing) > 5 else ""
        logger.warning(
            "Requested {} unitig ids missing from database ({}{})",
            len(missing),
            preview,
            suffix,
        )
        for uid in missing:
            unitig_info.pop(uid, None)

    if output_neighbors and unitig_info:
        ids_with_data = list(unitig_info.keys())
        for chunk in _chunked(ids_with_data, SQLITE_MAX_VARIABLES):
            placeholders = ", ".join("?" for _ in chunk)
            cursor.execute(
                f"""
                SELECT unitig_id, unitig_direction, neighbor_id, neighbor_direction
                  FROM adjacency
                 WHERE unitig_id IN ({placeholders})
                """,
                chunk,
            )
            for unitig_id, unitig_dir, neighbor_id, neighbor_dir in cursor.fetchall():
                meta = unitig_info.get(unitig_id)
                if meta is None:
                    continue
                neighbours = meta["neighbors"]
                if unitig_dir not in neighbours:
                    neighbours[unitig_dir] = []
                neighbours[unitig_dir].append((neighbor_id, neighbor_dir))

    return unitig_info


@lru_cache(maxsize=1024)
def lookup_unitig_cached(
    conn: sqlite3.Connection,
    unitig_id: str,
    output_neighbors: bool = True,
) -> Dict[str, object]:
    """Cache-friendly wrapper around :func:`lookup_unitigs_batch`."""

    result = lookup_unitigs_batch(conn, [unitig_id], output_neighbors)
    if not result:
        raise KeyError(f"Unitig {unitig_id} not present in the SQLite database")
    return result[unitig_id]


def lookup_unitig_sequences(
    conn: sqlite3.Connection, unitig_ids: Iterable[str]
) -> Dict[str, str]:
    """Fetch the DNA sequence for an iterable of unitig identifiers."""

    ids = _normalise_ids(unitig_ids)
    if not ids:
        return {}

    cursor = conn.cursor()
    sequences: Dict[str, str] = {}
    for chunk in _chunked(ids, SQLITE_MAX_VARIABLES):
        placeholders = ", ".join("?" for _ in chunk)
        cursor.execute(
            f"SELECT unitig_id, sequence FROM unitigs "
            f"WHERE unitig_id IN ({placeholders})",
            chunk,
        )
        for unitig_id, sequence in cursor.fetchall():
            if sequence is None:
                continue
            sequences[unitig_id] = sequence

    missing = set(ids) - sequences.keys()
    if missing:
        preview = ", ".join(sorted(missing)[:5])
        suffix = "..." if len(missing) > 5 else ""
        logger.warning(
            "Sequences for {} unitigs are missing in the SQLite database ({}{})",
            len(missing),
            preview,
            suffix,
        )
    return sequences


def lookup_unitig_sequence(conn: sqlite3.Connection, unitig_id: str) -> str | None:
    """Retrieve the sequence of a single unitig (or ``None`` if absent)."""

    sequences = lookup_unitig_sequences(conn, [unitig_id])
    return sequences.get(unitig_id)
