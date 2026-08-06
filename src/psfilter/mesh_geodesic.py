"""Memory-conscious local geodesic distances on a mesh graph.

The graph is built from the edges of an openCARP ``.elem`` file. Distances
are evaluated with bounded, target-aware Dijkstra searches rather than an
all-pairs distance matrix. This is suitable for local phase-singularity
matching on large meshes because only vertices reachable within the supplied
radius are explored.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
import heapq
from pathlib import Path
import time
from typing import Iterable

import numpy as np
from scipy.sparse import coo_matrix


PathLike = str | Path


_ELEMENT_EDGE_PATTERNS: dict[str, tuple[tuple[int, int], ...]] = {
    "Ln": ((0, 1),),
    "Tr": ((0, 1), (1, 2), (2, 0)),
    "Qd": ((0, 1), (1, 2), (2, 3), (3, 0)),
    "Tt": (
        (0, 1), (0, 2), (0, 3),
        (1, 2), (1, 3), (2, 3),
    ),
    "Py": (
        (0, 1), (1, 2), (2, 3), (3, 0),
        (0, 4), (1, 4), (2, 4), (3, 4),
    ),
    "Pr": (
        (0, 1), (1, 2), (2, 0),
        (3, 4), (4, 5), (5, 3),
        (0, 3), (1, 4), (2, 5),
    ),
    "Hx": (
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ),
}

_ELEMENT_NODE_COUNTS = {
    "Ln": 2,
    "Tr": 3,
    "Qd": 4,
    "Tt": 4,
    "Py": 5,
    "Pr": 6,
    "Hx": 8,
}


@dataclass(frozen=True)
class MeshGraphInfo:
    """Static information about a mesh geodesic graph."""

    n_vertices: int
    n_directed_edges: int
    n_undirected_edges: int
    element_count: int
    build_seconds: float
    compact_memory_bytes: int


class MeshGeodesicGraph:
    """CSR mesh graph with bounded local Dijkstra queries.

    Parameters are normally created through :meth:`from_elem_file`.
    The class keeps a bounded LRU cache of exact finite pair distances.
    Failed radius-limited searches are not cached because a later query may
    use a larger radius.
    """

    def __init__(
        self,
        *,
        coordinates: np.ndarray,
        indptr: np.ndarray,
        indices: np.ndarray,
        weights: np.ndarray,
        element_count: int = 0,
        build_seconds: float = 0.0,
        cache_size: int = 200_000,
    ) -> None:
        coordinates = np.asarray(coordinates, dtype=float)
        indptr = np.asarray(indptr)
        indices = np.asarray(indices)
        weights = np.asarray(weights)

        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError(
                "coordinates must have shape (n_vertices, 3)."
            )
        if indptr.ndim != 1 or indptr.size != coordinates.shape[0] + 1:
            raise ValueError("Invalid CSR indptr array.")
        if indices.ndim != 1 or weights.ndim != 1:
            raise ValueError("CSR indices and weights must be one-dimensional.")
        if indices.size != weights.size:
            raise ValueError("CSR indices and weights must have equal length.")
        if np.any(indices < 0) or np.any(indices >= coordinates.shape[0]):
            raise ValueError("CSR graph contains an invalid vertex index.")
        if np.any(weights < 0) or not np.all(np.isfinite(weights)):
            raise ValueError("Graph weights must be non-negative and finite.")

        cache_size = int(cache_size)
        if cache_size < 0:
            raise ValueError("cache_size must be non-negative.")

        self.coordinates = coordinates
        self.indptr = indptr
        self.indices = indices
        self.weights = weights
        self.cache_size = cache_size
        self._cache: OrderedDict[tuple[int, int], float] = OrderedDict()

        # Reusable arrays avoid allocating one million-entry arrays for every
        # local query. Only touched entries are initialized for each query.
        self._distance_workspace = np.empty(coordinates.shape[0], dtype=float)
        self._visit_mark = np.zeros(coordinates.shape[0], dtype=np.uint32)
        self._query_mark = np.uint32(0)

        self.distance_requests = 0
        self.dijkstra_runs = 0
        self.cache_hits = 0
        self.cache_misses = 0
        self.visited_vertices = 0
        self.relaxed_edges = 0

        compact_memory_bytes = int(
            self.indptr.nbytes
            + self.indices.nbytes
            + self.weights.nbytes
            + self._distance_workspace.nbytes
            + self._visit_mark.nbytes
        )
        self.info = MeshGraphInfo(
            n_vertices=int(coordinates.shape[0]),
            n_directed_edges=int(indices.size),
            n_undirected_edges=int(indices.size // 2),
            element_count=int(element_count),
            build_seconds=float(build_seconds),
            compact_memory_bytes=compact_memory_bytes,
        )

    @classmethod
    def from_elem_file(
        cls,
        coordinates: np.ndarray,
        elem_path: PathLike,
        *,
        cache_size: int = 200_000,
        edge_chunk_size: int = 500_000,
        weight_dtype: np.dtype | type = np.float32,
        weight_row_block: int = 100_000,
    ) -> "MeshGeodesicGraph":
        """Build an undirected weighted mesh graph from an openCARP file.

        Duplicate edges shared by adjacent elements are collapsed by the CSR
        conversion. Edge weights are Euclidean edge lengths; shortest paths
        along these edges approximate the surface geodesic distance.
        """
        start = time.perf_counter()
        coordinates = np.asarray(coordinates, dtype=float)

        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError(
                "coordinates must have shape (n_vertices, 3)."
            )
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("coordinates contain NaN or infinite values.")

        n_vertices = int(coordinates.shape[0])
        elem_path = Path(elem_path)
        if not elem_path.is_file():
            raise FileNotFoundError(f"Mesh element file not found: {elem_path}")

        edge_chunk_size = int(edge_chunk_size)
        if edge_chunk_size < 1:
            raise ValueError("edge_chunk_size must be at least 1.")

        source_chunks: list[np.ndarray] = []
        target_chunks: list[np.ndarray] = []
        source_buffer: list[int] = []
        target_buffer: list[int] = []

        def flush_buffers() -> None:
            if not source_buffer:
                return
            source_chunks.append(np.asarray(source_buffer, dtype=np.int32))
            target_chunks.append(np.asarray(target_buffer, dtype=np.int32))
            source_buffer.clear()
            target_buffer.clear()

        with elem_path.open("r", encoding="utf-8") as handle:
            first_line = handle.readline().strip()
            try:
                expected_elements = int(first_line)
            except ValueError as exc:
                raise ValueError(
                    f"First line of {elem_path} must be the element count."
                ) from exc

            element_count = 0
            for line_number, raw_line in enumerate(handle, start=2):
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split()
                element_type = parts[0]
                node_count = _ELEMENT_NODE_COUNTS.get(element_type)
                edge_pattern = _ELEMENT_EDGE_PATTERNS.get(element_type)

                if node_count is None or edge_pattern is None:
                    raise ValueError(
                        f"Unsupported element type {element_type!r} at "
                        f"{elem_path}:{line_number}."
                    )
                if len(parts) < node_count + 1:
                    raise ValueError(
                        f"Malformed {element_type} element at "
                        f"{elem_path}:{line_number}."
                    )

                try:
                    nodes = [int(value) for value in parts[1:1 + node_count]]
                except ValueError as exc:
                    raise ValueError(
                        f"Non-integer node index at {elem_path}:{line_number}."
                    ) from exc

                if min(nodes) < 0 or max(nodes) >= n_vertices:
                    raise ValueError(
                        f"Element at {elem_path}:{line_number} contains a "
                        f"vertex outside 0 to {n_vertices - 1}."
                    )

                for local_a, local_b in edge_pattern:
                    node_a = nodes[local_a]
                    node_b = nodes[local_b]
                    if node_a == node_b:
                        continue
                    source_buffer.append(node_a)
                    target_buffer.append(node_b)

                element_count += 1
                if len(source_buffer) >= edge_chunk_size:
                    flush_buffers()

        flush_buffers()

        if element_count != expected_elements:
            raise ValueError(
                f"Expected {expected_elements} elements in {elem_path}, "
                f"but parsed {element_count}."
            )
        if not source_chunks:
            raise ValueError(f"No mesh edges were found in {elem_path}.")

        source = np.concatenate(source_chunks)
        target = np.concatenate(target_chunks)
        directed_source = np.concatenate((source, target))
        directed_target = np.concatenate((target, source))

        # Values are placeholders. CSR conversion collapses duplicate edge
        # positions; exact edge lengths are computed afterward.
        adjacency = coo_matrix(
            (
                np.ones(directed_source.size, dtype=np.uint8),
                (directed_source, directed_target),
            ),
            shape=(n_vertices, n_vertices),
        ).tocsr()
        adjacency.setdiag(0)
        adjacency.eliminate_zeros()
        adjacency.sort_indices()

        indptr = adjacency.indptr.copy()
        indices = adjacency.indices.copy()
        weights = np.empty(indices.size, dtype=weight_dtype)

        weight_row_block = max(1, int(weight_row_block))
        for row_start in range(0, n_vertices, weight_row_block):
            row_end = min(n_vertices, row_start + weight_row_block)
            edge_start = int(indptr[row_start])
            edge_end = int(indptr[row_end])
            if edge_end == edge_start:
                continue

            counts = np.diff(indptr[row_start:row_end + 1])
            rows = np.repeat(
                np.arange(row_start, row_end, dtype=indices.dtype),
                counts,
            )
            cols = indices[edge_start:edge_end]
            delta = coordinates[rows] - coordinates[cols]
            block_weights = np.sqrt(
                np.einsum("ij,ij->i", delta, delta)
            )
            weights[edge_start:edge_end] = block_weights.astype(
                weight_dtype,
                copy=False,
            )

        build_seconds = time.perf_counter() - start
        return cls(
            coordinates=coordinates,
            indptr=indptr,
            indices=indices,
            weights=weights,
            element_count=element_count,
            build_seconds=build_seconds,
            cache_size=cache_size,
        )

    def _pair_key(self, source: int, target: int) -> tuple[int, int]:
        source = int(source)
        target = int(target)
        return (source, target) if source <= target else (target, source)

    def _cache_get(self, source: int, target: int) -> float | None:
        if self.cache_size == 0:
            return None
        key = self._pair_key(source, target)
        value = self._cache.get(key)
        if value is None:
            return None
        self._cache.move_to_end(key)
        self.cache_hits += 1
        return float(value)

    def _cache_put(self, source: int, target: int, distance: float) -> None:
        if self.cache_size == 0 or not np.isfinite(distance):
            return
        key = self._pair_key(source, target)
        self._cache[key] = float(distance)
        self._cache.move_to_end(key)
        while len(self._cache) > self.cache_size:
            self._cache.popitem(last=False)

    def _next_query_mark(self) -> np.uint32:
        next_mark = int(self._query_mark) + 1
        if next_mark >= np.iinfo(np.uint32).max:
            self._visit_mark.fill(0)
            next_mark = 1
        self._query_mark = np.uint32(next_mark)
        return self._query_mark

    def distances_to_targets(
        self,
        source: int,
        targets: Iterable[int] | np.ndarray,
        *,
        max_distance: float,
    ) -> np.ndarray:
        """Return exact geodesic distances up to ``max_distance``.

        Targets that are disconnected or farther than the limit receive
        ``np.inf``. The search stops as soon as every unresolved target has
        been reached or the heap exceeds the radius.
        """
        source = int(source)
        target_array = np.asarray(list(targets), dtype=int).reshape(-1)
        max_distance = float(max_distance)

        if source < 0 or source >= self.info.n_vertices:
            raise ValueError(f"Invalid source vertex: {source}")
        if np.any(target_array < 0) or np.any(
            target_array >= self.info.n_vertices
        ):
            raise ValueError("A target vertex is outside the mesh.")
        if not np.isfinite(max_distance) or max_distance < 0:
            raise ValueError("max_distance must be non-negative and finite.")

        self.distance_requests += int(target_array.size)
        result = np.full(target_array.size, np.inf, dtype=float)
        if target_array.size == 0:
            return result

        unresolved_positions: dict[int, list[int]] = {}
        for position, target in enumerate(target_array):
            target = int(target)
            if target == source:
                result[position] = 0.0
                continue

            cached = self._cache_get(source, target)
            if cached is not None:
                if cached <= max_distance:
                    result[position] = cached
                continue

            self.cache_misses += 1
            unresolved_positions.setdefault(target, []).append(position)

        if not unresolved_positions:
            return result

        self.dijkstra_runs += 1
        query_mark = self._next_query_mark()
        distances = self._distance_workspace
        marks = self._visit_mark
        distances[source] = 0.0
        marks[source] = query_mark
        heap: list[tuple[float, int]] = [(0.0, source)]
        remaining = set(unresolved_positions)

        while heap and remaining:
            current_distance, vertex = heapq.heappop(heap)
            if current_distance > max_distance:
                break
            if marks[vertex] != query_mark:
                continue
            if current_distance != distances[vertex]:
                continue

            self.visited_vertices += 1

            if vertex in remaining:
                for position in unresolved_positions[vertex]:
                    result[position] = current_distance
                self._cache_put(source, vertex, current_distance)
                remaining.remove(vertex)
                if not remaining:
                    break

            edge_start = int(self.indptr[vertex])
            edge_end = int(self.indptr[vertex + 1])
            self.relaxed_edges += edge_end - edge_start

            for edge_index in range(edge_start, edge_end):
                neighbour = int(self.indices[edge_index])
                candidate = current_distance + float(self.weights[edge_index])
                if candidate > max_distance:
                    continue
                if (
                    marks[neighbour] != query_mark
                    or candidate < distances[neighbour]
                ):
                    distances[neighbour] = candidate
                    marks[neighbour] = query_mark
                    heapq.heappush(heap, (candidate, neighbour))

        return result

    def distance(
        self,
        source: int,
        target: int,
        *,
        max_distance: float,
    ) -> float:
        """Return one bounded geodesic distance or ``np.inf``."""
        return float(
            self.distances_to_targets(
                source,
                [target],
                max_distance=max_distance,
            )[0]
        )

    def clear_cache(self) -> None:
        """Clear cached finite pair distances."""
        self._cache.clear()

    def statistics(self) -> dict[str, int | float]:
        """Return build and query diagnostics."""
        return {
            "n_vertices": self.info.n_vertices,
            "n_directed_edges": self.info.n_directed_edges,
            "n_undirected_edges": self.info.n_undirected_edges,
            "element_count": self.info.element_count,
            "graph_build_seconds": self.info.build_seconds,
            "compact_memory_bytes": self.info.compact_memory_bytes,
            "cache_entries": len(self._cache),
            "cache_size_limit": self.cache_size,
            "distance_requests": self.distance_requests,
            "dijkstra_runs": self.dijkstra_runs,
            "cache_hits": self.cache_hits,
            "cache_misses": self.cache_misses,
            "visited_vertices": self.visited_vertices,
            "relaxed_edges": self.relaxed_edges,
        }
