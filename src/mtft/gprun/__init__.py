"""mtft.gprun — a local job runner for long PARI/GP computations.

    py -m mtft.gprun

Opens a browser page on 127.0.0.1 with a drop zone, a Run/Stop toggle, a
live log tail, and a download button.  Pure standard library: no Flask, no
pip installs, nothing leaves the machine.

Why this exists.  The RT-1B holdout run is a multi-hour job, and the whole
arc kept tripping over one failure: what was computed and what the staged
artifact reproduces drifted apart.  So every run here is stamped
automatically:

  * the script is SHA-256'd **before** launch and the hash is written into
    the output header, so a result can never be separated from the exact
    bytes that produced it;
  * a copy of the script is frozen inside the run directory;
  * gp's version, the start and end times, the wall clock and the exit
    code all land in the same header;
  * output is flushed to disk continuously, so an interrupted run leaves a
    readable partial log rather than nothing.

Run directories are `mtft_runs/<timestamp>_<name>/` containing
`script.gp`, `output.txt` and `meta.json`.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import threading
import time
import webbrowser
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

__all__ = ["Job", "serve", "find_gp", "main"]

RUNS = Path("mtft_runs")


def find_gp() -> str | None:
    """Locate the gp executable, honouring MTFT_GP if set."""
    env = os.environ.get("MTFT_GP")
    if env and Path(env).exists():
        return env
    for name in ("gp", "gp.exe"):
        found = shutil.which(name)
        if found:
            return found
    for guess in (r"C:\Program Files\PARI\gp.exe",
                  r"C:\Program Files (x86)\PARI\gp.exe",
                  "/usr/bin/gp", "/usr/local/bin/gp", "/opt/homebrew/bin/gp"):
        if Path(guess).exists():
            return guess
    return None


def _gp_version(gp: str) -> str:
    try:
        out = subprocess.run([gp, "--version"], capture_output=True, text=True,
                             timeout=20)
        return (out.stdout + out.stderr).strip().splitlines()[0]
    except Exception:
        return "unknown"


class Job:
    """One GP run: frozen script, streamed output, stamped provenance."""

    def __init__(self, script_text: str, name: str, gp: str):
        self.name = "".join(c for c in name if c.isalnum() or c in "._-") or "script"
        self.gp = gp
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        self.dir = RUNS / f"{stamp}_{Path(self.name).stem}"
        self.dir.mkdir(parents=True, exist_ok=True)
        self.script = self.dir / "script.gp"
        self.script.write_text(script_text, encoding="utf-8")
        self.sha256 = hashlib.sha256(script_text.encode("utf-8")).hexdigest()
        self.out = self.dir / "output.txt"
        self.meta = self.dir / "meta.json"
        self.proc = None
        self.started = None
        self.finished = None
        self.returncode = None
        self._cancel = False
        self._lock = threading.Lock()

    # ---------------------------------------------------------------- run
    def start(self):
        self.started = time.time()
        header = [
            "=" * 72,
            "mtft.gprun",
            f"script        : {self.name}",
            f"sha256        : {self.sha256}",
            f"gp            : {_gp_version(self.gp)}",
            f"started (UTC) : {datetime.now(timezone.utc).isoformat()}",
            "",
            "Verify the hash above against your freeze manifest before",
            "reading any result. A result separated from its script hash is",
            "not reproducible.",
            "=" * 72,
            "",
        ]
        self.out.write_text("\n".join(header), encoding="utf-8")
        self._write_meta("running")
        threading.Thread(target=self._run, daemon=True).start()

    def _run(self):
        try:
            with open(self.out, "a", encoding="utf-8", errors="replace") as fh:
                self.proc = subprocess.Popen(
                    [self.gp, "-q", str(self.script)],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, errors="replace", bufsize=1,
                    cwd=str(self.script.parent))
                # Stop may have been pressed before Popen returned; honour it,
                # otherwise an early Stop silently does nothing and the job
                # runs on with no way to reach it.
                if self._cancel:
                    self.proc.terminate()
                for line in self.proc.stdout:
                    fh.write(line)
                    fh.flush()          # partial runs must stay readable
                self.proc.wait()
                self.returncode = self.proc.returncode
        except Exception as exc:                       # pragma: no cover
            with open(self.out, "a", encoding="utf-8") as fh:
                fh.write(f"\n*** runner error: {exc}\n")
            self.returncode = -1
        self.finished = time.time()
        with open(self.out, "a", encoding="utf-8") as fh:
            fh.write("\n" + "=" * 72 + "\n")
            fh.write(f"finished (UTC): {datetime.now(timezone.utc).isoformat()}\n")
            fh.write(f"wall clock    : {self.elapsed_str()}\n")
            fh.write(f"exit code     : {self.returncode}\n")
            fh.write(f"sha256        : {self.sha256}\n")
            if self.returncode not in (0, None):
                fh.write("NOTE: nonzero exit. Treat the output above as a "
                         "PARTIAL run, not a result.\n")
            fh.write("=" * 72 + "\n")
        self._write_meta("finished")

    def stop(self):
        self._cancel = True
        if self.proc and self.proc.poll() is None:
            self.proc.terminate()
            try:
                self.proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.proc.kill()

    # -------------------------------------------------------------- state
    @property
    def running(self):
        return self.proc is not None and self.proc.poll() is None

    def elapsed(self):
        if not self.started:
            return 0.0
        return (self.finished or time.time()) - self.started

    def elapsed_str(self):
        s = int(self.elapsed())
        return f"{s // 3600:d}h {(s % 3600) // 60:02d}m {s % 60:02d}s"

    def tail(self, offset=0):
        try:
            data = self.out.read_text(encoding="utf-8", errors="replace")
        except FileNotFoundError:
            return "", 0
        return data[offset:], len(data)

    def _write_meta(self, state):
        self.meta.write_text(json.dumps({
            "name": self.name, "sha256": self.sha256, "state": state,
            "started": self.started, "finished": self.finished,
            "returncode": self.returncode, "elapsed_seconds": self.elapsed(),
        }, indent=2), encoding="utf-8")


PAGE = """<!doctype html><meta charset=utf-8>
<title>mtft.gprun</title>
<style>
 :root{--bg:#12100e;--fg:#e8e3db;--dim:#8a8177;--acc:#c9782e;--ok:#5f9e5f;--bad:#b4534b}
 *{box-sizing:border-box} body{margin:0;background:var(--bg);color:var(--fg);
  font:14px/1.55 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
 header{padding:18px 24px;border-bottom:1px solid #2a2622;display:flex;
  align-items:baseline;gap:16px}
 h1{font-size:15px;margin:0;letter-spacing:.14em;text-transform:uppercase}
 .sub{color:var(--dim);font-size:12px}
 main{padding:24px;max-width:1100px}
 #drop{border:1px dashed #3a352f;padding:34px;text-align:center;color:var(--dim);
  border-radius:3px;transition:.15s;cursor:pointer}
 #drop.hot{border-color:var(--acc);color:var(--fg);background:#191612}
 .row{display:flex;gap:12px;align-items:center;margin:18px 0;flex-wrap:wrap}
 button{font:inherit;padding:9px 20px;border:1px solid #3a352f;background:#1b1815;
  color:var(--fg);border-radius:3px;cursor:pointer;letter-spacing:.06em}
 button:hover:not(:disabled){border-color:var(--acc)}
 button:disabled{opacity:.35;cursor:default}
 button.go{background:var(--acc);border-color:var(--acc);color:#140f09;font-weight:600}
 button.stop{background:var(--bad);border-color:var(--bad);color:#fff}
 pre{background:#0d0b09;border:1px solid #241f1b;padding:16px;border-radius:3px;
  height:60vh;overflow:auto;white-space:pre-wrap;font-size:12.5px;margin:0}
 .k{color:var(--dim)} .v{color:var(--fg)} .hash{font-size:11px;color:var(--dim);
  word-break:break-all}
 .pill{padding:3px 10px;border-radius:2px;font-size:11px;letter-spacing:.1em;
  text-transform:uppercase}
 .idle{background:#241f1b;color:var(--dim)} .run{background:var(--acc);color:#140f09}
 .done{background:var(--ok);color:#0b140b} .fail{background:var(--bad);color:#fff}
</style>
<header><h1>mtft.gprun</h1><span class=sub id=gpv></span></header>
<main>
 <div id=drop>drop a <b>.gp</b> file here &mdash; or click to choose</div>
 <input type=file id=file accept=".gp,.txt" hidden>
 <div class=row>
  <button class=go id=go disabled>Run</button>
  <button class=stop id=halt disabled>Stop</button>
  <button id=dl disabled>Download output.txt</button>
  <span class="pill idle" id=state>idle</span>
  <span class=k>elapsed</span><span class=v id=el>0h 00m 00s</span>
 </div>
 <div class=row><span class=k>script</span><span class=v id=nm>&mdash;</span></div>
 <div class=row><span class=k>sha256</span><span class="hash" id=sh>&mdash;</span></div>
 <pre id=log>Waiting for a script.

Every run is hashed before launch and the hash is written into the output
header. Check it against your freeze manifest before trusting a result.</pre>
</main>
<script>
let text=null,name=null,off=0,timer=null;
const $=i=>document.getElementById(i);
const drop=$('drop');
drop.onclick=()=>$('file').click();
drop.ondragover=e=>{e.preventDefault();drop.classList.add('hot')};
drop.ondragleave=()=>drop.classList.remove('hot');
drop.ondrop=e=>{e.preventDefault();drop.classList.remove('hot');take(e.dataTransfer.files[0])};
$('file').onchange=e=>take(e.target.files[0]);
function take(f){if(!f)return;const r=new FileReader();
 r.onload=()=>{text=r.result;name=f.name;$('nm').textContent=name;
  $('go').disabled=false;drop.textContent=name+' \\u2014 ready';
  $('log').textContent='Loaded '+name+' ('+text.length+' bytes). Press Run.';};
 r.readAsText(f);}
$('go').onclick=async()=>{
 $('go').disabled=true;off=0;$('log').textContent='';
 const r=await fetch('/api/run',{method:'POST',headers:{'X-Name':name},body:text});
 const j=await r.json();
 if(j.error){$('log').textContent='ERROR: '+j.error;$('go').disabled=false;return}
 $('sh').textContent=j.sha256;$('halt').disabled=false;poll();
 timer=setInterval(poll,1200);};
$('halt').onclick=()=>fetch('/api/stop',{method:'POST'});
$('dl').onclick=()=>location='/api/download';
async function poll(){
 const r=await fetch('/api/status?offset='+off);const j=await r.json();
 if(j.chunk){$('log').textContent+=j.chunk;off=j.offset;$('log').scrollTop=1e9;}
 $('el').textContent=j.elapsed;
 const s=$('state');
 if(j.running){s.className='pill run';s.textContent='running';}
 else if(j.started){const bad=j.returncode!==0;
  s.className='pill '+(bad?'fail':'done');
  s.textContent=bad?'partial / failed':'complete';
  $('halt').disabled=true;$('go').disabled=false;$('dl').disabled=false;
  if(timer){clearInterval(timer);timer=null;}}
}
fetch('/api/gp').then(r=>r.json()).then(j=>{
 $('gpv').textContent=j.gp?j.version+'  \\u00b7  '+j.gp:'gp NOT FOUND \\u2014 set MTFT_GP';
 if(!j.gp)$('gpv').style.color='#b4534b';});
</script>"""


class _Handler(BaseHTTPRequestHandler):
    job: Job | None = None
    gp: str | None = None

    def log_message(self, *a):
        pass

    def _send(self, code, body, ctype="application/json"):
        raw = body if isinstance(body, bytes) else body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)

    def do_GET(self):
        u = urlparse(self.path)
        if u.path == "/":
            return self._send(200, PAGE, "text/html; charset=utf-8")
        if u.path == "/api/gp":
            return self._send(200, json.dumps(
                {"gp": _Handler.gp,
                 "version": _gp_version(_Handler.gp) if _Handler.gp else ""}))
        if u.path == "/api/status":
            j = _Handler.job
            if not j:
                return self._send(200, json.dumps(
                    {"running": False, "started": False, "elapsed": "0h 00m 00s"}))
            off = int(parse_qs(u.query).get("offset", ["0"])[0])
            chunk, total = j.tail(off)
            return self._send(200, json.dumps(
                {"running": j.running, "started": True, "chunk": chunk,
                 "offset": total, "elapsed": j.elapsed_str(),
                 "returncode": j.returncode, "dir": str(j.dir)}))
        if u.path == "/api/download":
            j = _Handler.job
            if not j:
                return self._send(404, b"no run", "text/plain")
            return self._send(200, j.out.read_bytes(),
                              "text/plain; charset=utf-8")
        self._send(404, b"not found", "text/plain")

    def do_POST(self):
        u = urlparse(self.path)
        if u.path == "/api/run":
            if _Handler.job and _Handler.job.running:
                return self._send(200, json.dumps({"error": "a run is already active"}))
            if not _Handler.gp:
                return self._send(200, json.dumps(
                    {"error": "gp not found; install PARI/GP or set MTFT_GP"}))
            n = int(self.headers.get("Content-Length", 0))
            text = self.rfile.read(n).decode("utf-8", "replace")
            name = self.headers.get("X-Name", "script.gp")
            job = Job(text, name, _Handler.gp)
            _Handler.job = job
            job.start()
            return self._send(200, json.dumps(
                {"sha256": job.sha256, "dir": str(job.dir)}))
        if u.path == "/api/stop":
            if _Handler.job:
                _Handler.job.stop()
            return self._send(200, json.dumps({"stopped": True}))
        self._send(404, b"not found", "text/plain")


def serve(port=8731, open_browser=True):
    """Start the local runner UI."""
    _Handler.gp = find_gp()
    srv = ThreadingHTTPServer(("127.0.0.1", port), _Handler)
    url = f"http://127.0.0.1:{port}/"
    print(f"mtft.gprun  ->  {url}")
    print(f"gp: {_Handler.gp or 'NOT FOUND - install PARI/GP or set MTFT_GP'}")
    print(f"runs land in: {RUNS.resolve()}")
    print("Ctrl-C to stop the server (a running job is terminated with it).")
    if open_browser:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        if _Handler.job:
            _Handler.job.stop()
        print("\nstopped.")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    port = 8731
    if "--port" in argv:
        port = int(argv[argv.index("--port") + 1])
    serve(port=port, open_browser="--no-browser" not in argv)


if __name__ == "__main__":
    main()
