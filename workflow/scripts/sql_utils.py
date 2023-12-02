#!/usr/bin/env python3

import sqlite3
from functools import lru_cache, partial
from loguru import logger


def lookup_unitigs_batch(conn, unitig_ids, output_neighbors=True):
    unitig_ids_tuple = tuple(unitig_ids)
    unitig_query = (
        "SELECT unitig_id, length, colors FROM unitigs WHERE unitig_id IN ({})"
        .format(", ".join(["?"] * len(unitig_ids))))
    neighbor_query = "SELECT unitig_id, unitig_direction, neighbor_id, neighbor_direction FROM adjacency WHERE unitig_id IN ({})".format(
        ", ".join(["?"] * len(unitig_ids)))

    cur = conn.cursor()

    # Fetch unitig data
    cur.execute(unitig_query, unitig_ids_tuple)
    unitigs_data = cur.fetchall()
    unitig_info = {
        unitig_id: {
            "length": length,
            "colors": colors,
            "neighbors": {
                "+": [],
                "-": []
            }
        }
        for unitig_id, length, colors in unitigs_data
    }

    # Fetch neighbor data
    if output_neighbors:
        cur.execute(neighbor_query, unitig_ids_tuple)
        neighbors_data = cur.fetchall()
        for unitig_id, unitig_dir, neighbor_id, neighbor_dir in neighbors_data:
            if unitig_id in unitig_info:
                direction = ("+" if unitig_dir == "+" else "-"
                             )  # Adjust based on how you store direction
                unitig_info[unitig_id]["neighbors"][direction].append(
                    (neighbor_id, neighbor_dir))

    return unitig_info


@lru_cache(maxsize=1024)
def lookup_unitig_cached(conn, unitig_id, output_neighbors=True):
    return lookup_unitigs_batch(conn, [unitig_id], output_neighbors)[unitig_id]


def lookup_unitig_sequence(conn, unitig_id):
    """
    Retrieves the sequence of a unitig given its unitig_id.
    """
    cur = conn.cursor()
    cur.execute("SELECT sequence FROM unitigs WHERE unitig_id = ?",
                (unitig_id, ))
    result = cur.fetchone()
    if result:
        return result[0]  # Return the sequence
    else:
        return None  # No sequence found for the given unitig_id
