"""Portable, offline HTML and deterministic SVG for the interaction atlas."""
from __future__ import annotations

import html
import json
import math
from importlib.resources import files

from .catalog import SECTORS, new_ledger, validate_catalog, validate_ledger


def vertex_svg(vertex, catalog):
    """Draw an ordered, all-incoming vertex; labels retain UFO leg indices.

    This is a vertex glyph, not a generated scattering diagram. Fermion arrows
    follow fermion flow: toward the vertex for particles, away for antiparticles.
    """
    esc = html.escape
    particles = [catalog["particles"][p] for p in vertex["particles"]]
    title = f"{vertex['id']}: " + ", ".join(p["name"] for p in particles)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 420 280" role="img" aria-label="{esc(title)}">',
             '<rect width="420" height="280" rx="14" fill="#f6f8fc"/>']
    cx, cy = 210, 140
    for index, p in enumerate(particles):
        angle = -math.pi / 2 + 2 * math.pi * index / len(particles)
        dx, dy = math.cos(angle), math.sin(angle)
        ex, ey = cx + 94 * dx, cy + 94 * dy
        color = "#a34718" if abs(p["color"]) != 1 else "#215f91"
        if p["spin"] == 3:
            pts = []
            # Gluon coils have an additional longitudinal oscillation.
            for step in range(181):
                t = step / 180
                transverse = 4 * math.sin(14 * math.pi * t)
                along = 94 * t + (3 * math.sin(28 * math.pi * t) if p["pdg_code"] == 21 else 0)
                pts.append(f"{cx+along*dx-transverse*dy:.2f},{cy+along*dy+transverse*dx:.2f}")
            parts.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color}" stroke-width="2"/>')
        else:
            dash = ' stroke-dasharray="7 5"' if p["spin"] == 1 else (
                ' stroke-dasharray="2 4"' if p["spin"] == -1 else "")
            parts.append(f'<path d="M {cx} {cy} L {ex:.2f} {ey:.2f}" stroke="{color}" stroke-width="2"{dash}/>')
            if p["spin"] == 2:
                direction = -1 if p["pdg_code"] > 0 else 1
                tx, ty = cx + 52 * dx, cy + 52 * dy
                tip = (tx + direction * 6 * dx, ty + direction * 6 * dy)
                back = (tx - direction * 6 * dx, ty - direction * 6 * dy)
                points = [tip, (back[0] - 4 * dy, back[1] + 4 * dx),
                          (back[0] + 4 * dy, back[1] - 4 * dx)]
                parts.append(f'<polygon points="{" ".join(f"{x:.2f},{y:.2f}" for x,y in points)}" fill="{color}"/>')
        lx, ly = cx + 119 * dx, cy + 119 * dy
        label = esc(f"{index+1}: {p['name']}")
        parts.append(f'<text x="{lx:.2f}" y="{ly:.2f}" text-anchor="middle" dominant-baseline="middle" font-family="sans-serif" font-size="14" fill="#172a3b">{label}</text>')
    parts.append('<circle cx="210" cy="140" r="5" fill="#172a3b"/></svg>')
    return "".join(parts)


def render_html(catalog, ledger=None):
    """Render a self-contained searchable atlas, with no CDN or local server."""
    validate_catalog(catalog)
    ledger = new_ledger(catalog) if ledger is None else ledger
    validate_ledger(ledger, catalog)
    payload = {
        "catalog": catalog, "ledger": ledger, "sectors": SECTORS,
        "svgs": {v["id"]: vertex_svg(v, catalog) for v in catalog["vertices"]},
        "guide": json.loads(files(__package__).joinpath("_data/guide.json").read_text(encoding="utf-8")),
    }
    # A JSON string must not be able to close its enclosing script element.
    encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False).replace("<", "\\u003c").replace("&", "\\u0026")
    template = files(__package__).joinpath("_data/atlas.html").read_text(encoding="utf-8")
    return template.replace("__ATLAS_DATA__", encoded)
