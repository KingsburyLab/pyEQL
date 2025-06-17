"""Reference data source parsers."""

from typing import Any

from pyEQL.benchmark import _BenchmarkEntry


def parse_crc(d: dict[str, Any]) -> _BenchmarkEntry:
    """Parse data from CRC."""
    # read solution from "Mol. form" key

    # mean activity coefficient from <i>gamma</i>(0.1 m)-like keys

    # electrical conductivity from <i>kappa</i>(0.5%)-like keys

    # molar electrical conductivity from <i>Λ</i>/S cm<sup>2</sup> mol<sup>-1</sup><br/>0.0001 M-like keys


def parse_idst(d: dict[str, Any]) -> _BenchmarkEntry:
    # molal concentration, molal activity/osmotic coefficient
    pass
