"""Common tooling abstractions for agent nodes."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, Tuple, Any


@dataclass(frozen=True)
class Node:
    """Lightweight wrapper describing an executable tool node."""

    name: str
    entry_point: Callable[..., Dict[str, Any]]
    description: str
    dependencies: Tuple[str, ...] = field(default_factory=tuple)
    aliases: Tuple[str, ...] = field(default_factory=tuple)
    tags: Tuple[str, ...] = field(default_factory=tuple)
    run_in_second_loop: bool = False

    def __call__(self, *args, **kwargs) -> Dict[str, Any]:
        return self.entry_point(*args, **kwargs)

    def iter_bindings(self) -> Iterable[Tuple[str, Callable[..., Dict[str, Any]]]]:
        """Yield `(key, callable)` pairs for the registry map."""
        yield self.name, self.entry_point
        for alias in self.aliases:
            yield alias, self.entry_point

    def iter_catalog_entries(self) -> Iterable[Tuple[str, str]]:
        """Yield `(key, description)` for the catalog lookup."""
        yield self.name, self.description
        for alias in self.aliases:
            yield alias, self.description

