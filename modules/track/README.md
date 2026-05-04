# track

This module generates both UCSC track lists and `igv.js` HTML pages from coverage tracks.

## Outputs

- `ucsc_track.txt` for UCSC Genome Browser style track loading.
- `igv_track.html` for interactive viewing with `igv.js`.

## Inputs

- BigWig files for standard coverage tracks.
- Plus and minus bedGraph files for strand-specific iCLIP tracks.

## Notes

- The HTML output references tracks relative to the output page, so keep the HTML next to the track files or serve the directory over HTTP.
- When loading locally, a browser may block `file://` access to track data; a simple local web server avoids that issue.
