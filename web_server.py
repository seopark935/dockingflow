"""Serve the DockingFlow GUI over HTTP, for running it on a headless docking server.

`python3 gui.py --web` calls `serve()`, which exposes the same `PipelineAPI`
the desktop window uses:

  - `GET /` serves `gui_assets/index.html` with `window.DF_WEB = true`
    injected, which switches the page's `pywebview.api.<method>(...)` calls
    to HTTP and its file dialogs to an in-page picker.
  - `POST /api/<method>` calls `PipelineAPI.<method>(*args)` with the JSON
    array body as arguments and returns the result as JSON.

Security: the API can start runs and delete directories as the user who
started the server, so by default it binds to 127.0.0.1 only (reach it from
a laptop through an SSH tunnel), and every API call must carry a random
per-process token, handed out in the startup URL. The token matters on a
shared server, where other local users can reach 127.0.0.1 too.

Standard library only, so nothing needs installing on the server.
"""
from __future__ import annotations

import hmac
import json
import secrets
import socket
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

# PipelineAPI methods the browser must not call: desktop-only plumbing.
NOT_EXPOSED = {"set_window", "browse_file", "browse_folder"}

TOKEN_HEADER = "X-DockingFlow-Token"
MAX_BODY_BYTES = 1_000_000


def exposed_methods(api: Any) -> set[str]:
    return {
        name
        for name in dir(api)
        if not name.startswith("_") and name not in NOT_EXPOSED and callable(getattr(api, name))
    }


def make_handler(api: Any, html_path: Path, token: str) -> type[BaseHTTPRequestHandler]:
    methods = exposed_methods(api)

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, format: str, *args: Any) -> None:
            pass  # the page polls every ~700ms; per-request logging would flood the log

        def _send(self, code: int, body: bytes, content_type: str) -> None:
            self.send_response(code)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def _send_json(self, code: int, obj: Any) -> None:
            self._send(code, json.dumps(obj).encode(), "application/json")

        def do_GET(self) -> None:
            if self.path.split("?", 1)[0] != "/":
                self._send(404, b"not found", "text/plain")
                return
            html = html_path.read_text(encoding="utf-8")
            html = html.replace("<head>", "<head>\n<script>window.DF_WEB = true;</script>", 1)
            self._send(200, html.encode(), "text/html; charset=utf-8")

        def do_POST(self) -> None:
            if not hmac.compare_digest(self.headers.get(TOKEN_HEADER, ""), token):
                self._send_json(403, {"ok": False, "message": "Missing or wrong access token. Open the URL printed at startup."})
                return
            if not self.path.startswith("/api/"):
                self._send_json(404, {"ok": False, "message": "not found"})
                return
            name = self.path[len("/api/"):]
            if name not in methods:
                self._send_json(404, {"ok": False, "message": f"unknown API method: {name}"})
                return

            length = int(self.headers.get("Content-Length") or 0)
            if length > MAX_BODY_BYTES:
                self._send_json(413, {"ok": False, "message": "request too large"})
                return
            try:
                args = json.loads(self.rfile.read(length) or b"[]")
                if not isinstance(args, list):
                    raise ValueError("body must be a JSON array of arguments")
                result = getattr(api, name)(*args)
            except Exception as exc:  # PipelineAPI methods catch their own errors; this is bad input
                self._send_json(400, {"ok": False, "message": f"{name}: {exc}"})
                return
            self._send_json(200, result)

    return Handler


def serve(api: Any, html_path: Path, host: str = "127.0.0.1", port: int = 8765) -> None:
    """Serve the GUI until interrupted. Prints the access URL and SSH-tunnel command."""
    token = secrets.token_urlsafe(16)
    try:
        server = ThreadingHTTPServer((host, port), make_handler(api, html_path, token))
    except OSError as exc:
        sys.exit(f"Can't listen on {host}:{port} ({exc}). Is the GUI already running? Try another --port.")
    server.daemon_threads = True

    hostname = socket.gethostname()
    print(f"DockingFlow web GUI running on {hostname}, port {port}.", flush=True)
    if host in ("127.0.0.1", "localhost"):
        print("From your laptop, open an SSH tunnel (leave it running):", flush=True)
        print(f"  ssh -N -L {port}:localhost:{port} <you>@{hostname}", flush=True)
        print("then open in your browser:", flush=True)
        print(f"  DOCKINGFLOW_URL=http://localhost:{port}/?token={token}", flush=True)
    else:
        print(f"WARNING: listening on {host}, not just localhost; anyone with the URL can control runs.", flush=True)
        print(f"  DOCKINGFLOW_URL=http://{hostname}:{port}/?token={token}", flush=True)
    print("Runs keep going if the browser or tunnel closes; stopping this process stops them.", flush=True)

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
