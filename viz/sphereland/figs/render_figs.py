"""Render pictures of Sphereland's instruments from the page itself (written 2026-11-16 for the README of mtft 1.1.0).

    python viz/sphereland/figs/render_figs.py --fonts <node_modules> --out <folder> [--only 22 10 18 24] [--scale 2]

Needs Playwright for Python and Google Chrome.  The page asks fonts.googleapis.com for its three typefaces
(Bricolage Grotesque, IBM Plex Sans, IBM Plex Mono); this script answers those requests from local copies, the
fontsource packages from the npm registry, so that a machine without access to Google still draws the page in its own
type:

    npm install @fontsource/ibm-plex-sans@5.3.0 @fontsource/ibm-plex-mono@5.3.0 @fontsource-variable/bricolage-grotesque@5.3.0

Nothing else on the page is changed: each picture is a screenshot of one instrument as the page draws it when it opens.
The page's console errors are reported, and a run with any error exits with status 1.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

PAGE = Path(__file__).resolve().parents[1] / "index.html"
INSTRUMENTS = {
    "1": "viewPanel", "2": "tablePanel", "3": "holesPanel", "4": "shellPanel", "5": "ladderPanel", "6": "motionPanel",
    "7": "responsePanel", "8": "knotPanel", "9": "slicePanel", "10": "modPanel", "11": "qubitPanel", "12": "gapPanel",
    "13": "kakPanel", "14": "clockPanel", "15": "besselPanel", "16": "refPanel", "17": "casimirPanel",
    "18": "gravityPanel", "19": "kroneckerPanel", "20": "poissonPanel", "21": "aetherPanel", "22": "flipPanel",
    "23": "stripPanel", "24": "curveGapPanel", "25": "relPanel",
}
FONT_CSS = [("@fontsource/ibm-plex-sans", f"{w}.css") for w in (400, 500, 600)] + \
           [("@fontsource/ibm-plex-mono", f"{w}.css") for w in (400, 500)] + \
           [("@fontsource-variable/bricolage-grotesque", "opsz.css")]
LOCAL = "https://fonts.gstatic.com/local/"
# the pictures in the README: file name -> (instrument, figure), each drawn above its caption
README_PICTURES = {
    "i22_two_fillings.png": ("22", 1),
    "i22_which_filling_wins.png": ("22", 3),
    "i10_modular_knot.png": ("10", 1),
    "i24_modular_spectrum.png": ("24", 1),
    "i18_resonance_knot.png": ("18", 3),
}


def font_css(modules: Path) -> str:
    """The fontsource CSS for the page's weights, every subset, with its files served under LOCAL."""
    out = []
    for package, name in FONT_CSS:
        css = (modules / package / name).read_text(encoding="utf-8")
        css = css.replace("'Bricolage Grotesque Variable'", "'Bricolage Grotesque'")
        out.append(re.sub(r"url\(\./files/([^)]+)\)", lambda m, p=package: f"url({LOCAL}{p}/files/{m.group(1)})", css))
    return "\n".join(out)


def main(argv=None) -> int:
    from playwright.sync_api import sync_playwright

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fonts", required=True, type=Path, help="a node_modules folder holding the three fontsource packages")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--only", nargs="*", default=None, help="instrument numbers (default: all 25)")
    ap.add_argument("--scale", type=float, default=2.0)
    ap.add_argument("--width", type=int, default=1280)
    ap.add_argument("--chrome", default="/usr/bin/google-chrome")
    ap.add_argument("--figures", action="store_true", help="also photograph each <figure> of the chosen instruments")
    ap.add_argument("--no-captions", action="store_true", help="photograph a figure above its <figcaption>")
    ap.add_argument("--readme", action="store_true", help="only the README's pictures, under their names, with PROVENANCE.txt")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    css = font_css(args.fonts)
    errors, fonts_served, record = [], [], {}

    def serve_css(route):
        route.fulfill(status=200, content_type="text/css", body=css, headers={"Access-Control-Allow-Origin": "*"})

    def serve_font(route):
        rel = route.request.url[len(LOCAL):]
        path = (args.fonts / rel).resolve()
        if args.fonts.resolve() not in path.parents or not path.is_file():
            route.fulfill(status=404, body="")
            return
        fonts_served.append(rel)
        route.fulfill(status=200, content_type="font/woff2", body=path.read_bytes(), headers={"Access-Control-Allow-Origin": "*"})

    with sync_playwright() as p:
        browser = p.chromium.launch(executable_path=args.chrome, headless=True)
        context = browser.new_context(viewport={"width": args.width, "height": 1000}, device_scale_factor=args.scale,
                                      color_scheme="light", reduced_motion="no-preference")
        page = context.new_page()
        page.on("console", lambda m: errors.append(f"console.{m.type}: {m.text}") if m.type == "error" else None)
        page.on("pageerror", lambda e: errors.append(f"pageerror: {e}"))
        page.route("https://fonts.googleapis.com/**", serve_css)
        page.route(LOCAL + "**", serve_font)
        page.goto(PAGE.as_uri(), wait_until="load")
        page.evaluate("document.fonts.ready.then(() => true)")
        page.wait_for_timeout(1500)
        if args.readme:
            for name, (number, k) in README_PICTURES.items():
                fig = page.locator(f"#{INSTRUMENTS[number]} figure").nth(k - 1)
                fig.scroll_into_view_if_needed()
                page.wait_for_timeout(600)
                fbox = fig.bounding_box()
                top = fig.locator("figcaption").first.bounding_box()["y"]
                clip = {"x": fbox["x"], "y": fbox["y"], "width": fbox["width"], "height": top - fbox["y"]}
                page.screenshot(path=str(args.out / name), clip=clip)
                record[name] = {"instrument": number, "figure": k, "panel": INSTRUMENTS[number],
                                "css_px": [round(clip["width"]), round(clip["height"])],
                                "sha256": hashlib.sha256((args.out / name).read_bytes()).hexdigest()}
        chosen = [] if args.readme else (args.only or list(INSTRUMENTS))
        for number in chosen:
            element = page.locator(f"#{INSTRUMENTS[number]}")
            element.scroll_into_view_if_needed()
            page.wait_for_timeout(600)                     # let an instrument that draws on view finish its first frame
            target = args.out / f"sphereland_{int(number):02d}.png"
            element.screenshot(path=str(target))
            box = element.bounding_box()
            record[number] = {"panel": INSTRUMENTS[number], "file": target.name, "css_px": [round(box["width"]), round(box["height"])],
                              "sha256": hashlib.sha256(target.read_bytes()).hexdigest()}
            if args.figures:
                figures = page.locator(f"#{INSTRUMENTS[number]} figure")
                for k in range(figures.count()):
                    fig = figures.nth(k)
                    fig.scroll_into_view_if_needed()
                    page.wait_for_timeout(300)
                    ftarget = args.out / f"sphereland_{int(number):02d}_fig{k + 1}.png"
                    fbox = fig.bounding_box()
                    caption = fig.locator("figcaption")
                    if args.no_captions and caption.count():
                        top = caption.first.bounding_box()["y"]
                        clip = {"x": fbox["x"], "y": fbox["y"], "width": fbox["width"], "height": max(1.0, top - fbox["y"])}
                        page.screenshot(path=str(ftarget), clip=clip)
                        fbox = clip
                    else:
                        fig.screenshot(path=str(ftarget))
                    record[f"{number}.{k + 1}"] = {"panel": INSTRUMENTS[number], "figure": k + 1, "file": ftarget.name,
                                                   "css_px": [round(fbox["width"]), round(fbox["height"])],
                                                   "sha256": hashlib.sha256(ftarget.read_bytes()).hexdigest()}
        families = page.evaluate("Array.from(document.fonts).filter(f => f.status === 'loaded').map(f => f.family + ' ' + f.weight)")
        version = browser.version
        browser.close()
    summary = {"page": PAGE.name, "page_sha256": hashlib.sha256(PAGE.read_bytes()).hexdigest(), "chrome": version,
               "scale": args.scale, "viewport_width": args.width, "fonts_loaded": sorted(set(families)),
               "font_files_served": len(set(fonts_served)), "errors": errors, "pictures": record}
    print(json.dumps(summary, indent=1, ensure_ascii=False))
    if args.readme and not errors:
        lines = [f"Pictures of Sphereland's instruments for the README, rendered from {PAGE.name} (SHA-256 {summary['page_sha256']})",
                 f"by headless Google Chrome {version} through Playwright, at {args.scale:g}x, viewport {args.width} px wide, light scheme;",
                 "fonts: " + ", ".join(summary["fonts_loaded"]) + " (fontsource 5.3.0 from npm, served in place of Google Fonts).",
                 "Each picture is one <figure> of the page as it opens, above its caption. Made by",
                 "    python viz/sphereland/figs/render_figs.py --fonts <node_modules> --out viz/sphereland/figs --readme", ""]
        lines += [f"{name}  instrument {r['instrument']}, figure {r['figure']} (#{r['panel']}), {r['css_px'][0]} x {r['css_px'][1]} css px, SHA-256 {r['sha256']}"
                  for name, r in record.items()]
        (args.out / "PROVENANCE.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
