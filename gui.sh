#!/bin/bash
# DockingFlow GUI launcher for the docking server (e.g. over MobaXterm).
#
#   bash gui.sh           open the GUI (starts the GUI server first if needed)
#   bash gui.sh browser   same, but in a browser window instead of the Tk one
#   bash gui.sh status    is the GUI server running? print its URL
#   bash gui.sh stop      stop the GUI server (also stops a run in progress;
#                         starting a run again later resumes where it left off)
#   bash gui.sh setup     optional: install a browser-engine window (pywebview +
#                         Qt, into .venv) -- only needed if the server has
#                         neither tkinter nor a web browser
#
# How it works: the GUI server (`gui.py --web`) runs detached in the
# background and owns any docking run, so runs keep going when you close
# the window or MobaXterm. The window is just a view onto it, shown on your
# own screen through MobaXterm's built-in X server (X11 forwarding). It's a
# Tkinter window (gui_tk.py) when Python has tkinter -- much faster over X11
# than a browser, which has to send every frame as pixels.

cd "$(dirname "$0")" || exit 1
PIDFILE=.gui_server.pid
LOG=gui_server.log
PORT="${DOCKINGFLOW_PORT:-8765}"
PYTHON="${PYTHON:-python3}"
BROWSER_PROFILE="$HOME/.dockingflow-browser"

running() { [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; }
url() { grep -m1 'DOCKINGFLOW_URL=' "$LOG" 2>/dev/null | sed 's/.*DOCKINGFLOW_URL=//'; }

start_server() {
    running && return 0
    nohup "$PYTHON" -u gui.py --web --port "$PORT" > "$LOG" 2>&1 &
    echo $! > "$PIDFILE"
    for _ in $(seq 20); do
        [ -n "$(url)" ] && break
        sleep 0.5
    done
    if ! running || [ -z "$(url)" ]; then
        echo "The GUI server failed to start. Its log ($LOG):"
        cat "$LOG"
        rm -f "$PIDFILE"
        exit 1
    fi
    echo "Started the DockingFlow GUI server (pid $(cat "$PIDFILE"), log: $LOG)."
}

no_window_help() {
    cat <<EOF

Couldn't open a window: $1

The GUI server is running, so you can still reach it from your own browser:
  1. In MobaXterm, open Tools > Network > MobaSSHTunnel (port forwarding) > New SSH tunnel:
       Local port forwarding, "My computer" port $PORT,
       SSH server = this server (your usual login),
       Remote server = localhost, port $PORT.
     Save, then click Start (the play button).
  2. Open this URL in your browser:
       $(url)
EOF
}

open_window() {
    local u
    u="$(url)"
    if [ -z "$DISPLAY" ]; then
        no_window_help "X11 forwarding is off (DISPLAY isn't set). In MobaXterm, edit the session:
  SSH > Advanced SSH settings > tick 'X11-Forwarding', and check the 'X server' button
  (top right) is on. Then reconnect and run 'bash gui.sh' again."
        return 1
    fi

    echo "Opening the GUI window on your screen (it may take a few seconds over the network)..."

    # 1. The Tkinter window: fast over X11, and part of Python itself.
    #    ('bash gui.sh browser' skips it, to use a browser-based window instead.)
    if [ "$WINDOW" != browser ]; then
        local py
        for py in "$PYTHON" /usr/bin/python3; do
            if "$py" -c "import tkinter" 2>/dev/null; then
                "$py" gui_tk.py "$u" 2>>"$LOG" && return
                echo "(the Tk window failed; see $LOG -- trying a browser-based window instead)"
                break
            fi
        done
        [ -n "$py" ] && ! "$py" -c "import tkinter" 2>/dev/null && echo \
            "(Python's tkinter isn't installed here -- ask your admin for the 'python3-tk' package" \
            "for the fast window. Falling back to a slower browser-based window.)"
    fi

    # 2. A browser-engine window, if 'bash gui.sh setup' was run.
    local vpy=".venv/bin/python3"
    if [ -x "$vpy" ] && "$vpy" -c "import webview" 2>/dev/null; then
        "$vpy" gui.py --viewer "$u" 2>>"$LOG" && return
        echo "(the desktop window failed; see $LOG -- trying a browser instead)"
    fi

    # 3. A browser installed on the server, as an app-style window.
    local b
    for b in chromium chromium-browser google-chrome google-chrome-stable; do
        if command -v "$b" >/dev/null 2>&1; then
            "$b" --app="$u" --user-data-dir="$BROWSER_PROFILE" --no-first-run --disable-gpu \
                --disable-dev-shm-usage 2>>"$LOG" >/dev/null && return
        fi
    done
    if command -v firefox >/dev/null 2>&1; then
        mkdir -p "$BROWSER_PROFILE-firefox"
        firefox --no-remote --profile "$BROWSER_PROFILE-firefox" "$u" 2>>"$LOG" >/dev/null && return
    fi

    no_window_help "no web browser (chromium/chrome/firefox) is installed on this server.
Tip: 'bash gui.sh setup' installs a small desktop window instead (no admin rights needed)."
    return 1
}

WINDOW=""
case "${1:-open}" in
    stop)
        if running; then kill "$(cat "$PIDFILE")"; echo "Stopped the GUI server."; else echo "Not running."; fi
        rm -f "$PIDFILE"
        ;;
    status)
        if running; then
            echo "GUI server running (pid $(cat "$PIDFILE")): $(url)"
        else
            echo "GUI server not running."
        fi
        ;;
    setup)
        echo "Installing pywebview + Qt into .venv (roughly 200 MB download)..."
        [ -x .venv/bin/python3 ] || "$PYTHON" -m venv .venv || exit 1
        .venv/bin/python3 -m pip install --upgrade pip >/dev/null
        .venv/bin/python3 -m pip install "pywebview[qt]" || exit 1
        echo "Done. Run 'bash gui.sh' to open the GUI."
        ;;
    open|browser)
        WINDOW="$1"
        start_server
        if open_window; then
            echo
            echo "Window closed. The GUI server keeps running (and so does any docking run)."
        fi
        echo "Reopen with 'bash gui.sh'; stop everything with 'bash gui.sh stop'."
        ;;
    *)
        echo "Usage: bash gui.sh [open|browser|status|stop|setup]"
        exit 1
        ;;
esac
