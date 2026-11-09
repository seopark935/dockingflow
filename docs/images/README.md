# Screenshots for the manual

`docs/MANUAL.md` refers to the images below. They are not in the repository:
screenshots have to be taken from a real session, and they go stale whenever
the window changes, so whoever updates the manual retakes them.

## How to capture

The DockingFlow window is drawn on **your own PC** (MobaXterm forwards it over
X11), so capture it there:

- **Windows**: press `Win` + `Shift` + `S`, drag around the window, then paste
  into Paint and save as PNG.
- **macOS**: `Cmd` + `Shift` + `4`, then press space and click the window.
- **On the server itself** (if it has a desktop and ImageMagick):
  `import -window root ~/shot.png`

Save each file with the exact name below, in this folder, then view the manual
in any Markdown viewer to check it reads well.

## Wanted images

| File | What it should show |
|---|---|
| `01-overview.png` | The whole window just after `bash gui.sh`, before anything is configured. |
| `02-running.png` | A run in progress: tranche table part-done, log scrolling, a few hits listed. |
| `03-ligands-zinc.png` | Tab 1 with the source set to ZINC22 and a preset chosen, showing the live count. |
| `04-ligands-chembl.png` | Tab 1 with the source set to ChEMBL, showing the weight/logP boxes, the fetch dropdown and the Check Open Babel button. |
| `05-targets.png` | Tab 2 with a real receptor loaded and a box filled in by "Center on known site". |
| `06-resources.png` | Tab 4 showing the detected server and the plan line (`mpirun -np …`). |
| `07-details.png` | The Details… window for one tranche, showing the summary counts. |

Crop to the window, keep the full width, and avoid including anything
confidential (paths with personal directories are fine to blur).
